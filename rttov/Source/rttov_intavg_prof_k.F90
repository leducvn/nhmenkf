SUBROUTINE rttov_intavg_prof_k( &
            & opts,      &
            & chanprof,  &
            & kni,       &
            & kno,       &
            & ProfIn,    &
            & ProfIn_k,  &
            & ProfOut,   &
            & ProfOut_k)
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!-----------------------------------------------------------------------
! usage of arguments
!              set of profile elements
!              -----------------------
! nprofiles     number of profiles
! kni           number of source levs
! kno           number of target levs
! ProfIn        source for interpolator
! ProfIn_k      source for interpolator
! ProfOut       target for interpolator
! ProfOut_k     target for interpolator
!
!
!--------------------------------------------------------------------------
! application within RTTOV (for NWP SAF)
!   - profiles on USER levels -> prediction of gas/atm opdeps on COEF levels
!   - gas/atm opdeps on COEF levels -> transmittances on USER levels for RTE
!---------------------------------------------------------------------------
!
!     History:
!
!     Version   Date      Comment
!     -------   ----      -------
!
!     1         12/2006   main code intavg restructured Peter Rayer
!
!     2         10/2007   Optimisation Deborah Salmond and Pascal Brunel
!
!     3         11/2007   layeravg now deals better with output levels near
!                         and beyond the upper and lower boundaries of the
!                         input levels, so there is now no need for a special
!                         case here - A. Geer
!
!     4         01/2008   lgradp option for TL/AD for user pressure levels,
!                         Niels Bormann
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
  USE rttov_types, ONLY : rttov_options, rttov_chanprof, profile_type
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_options ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno                 ! number of levels
  TYPE(profile_type)  , INTENT(IN)    :: ProfIn   (:)             ! atmospheric profiles
  TYPE(profile_type)  , INTENT(INOUT) :: ProfIn_k (size(chanprof))
  TYPE(profile_type)  , INTENT(INOUT) :: ProfOut  (size(ProfIn)  )! atmospheric profiles
  TYPE(profile_type)  , INTENT(INOUT) :: ProfOut_k(size(chanprof))
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_k.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, kn, istart(kno, size(ProfIn)              ), iend (kno, size(ProfIn))
  INTEGER(KIND=jpim) :: prof
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles                                                                                                  ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                                                                  ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=jprb) :: pvlev  (kni, size(ProfIn))
  REAL   (KIND=jprb) :: pvlev_k(kni)
  REAL   (KIND=jprb) :: ppo    (kno, size(ProfIn))
  REAL   (KIND=jprb) ::      &
    & zlnpi  (kni, size(ProfIn)), zlnpi_k(kni), zpz   (kni, kno         , size(ProfIn)), zpz_k(kni, kno         )
  REAL   (KIND=jprb) :: zlnpo  (kno, size(ProfIn)), zlnpo_k(kno)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_K', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = size(ProfIn)
  nchannels         = size(chanprof)
  IF (opts%lgradp) THEN
    llevels_different = .TRUE.
    zpz_k             = 0._jprb
  ELSE
    DO prof = 1, nprofiles - 1
      DO jo = 1, kni
        IF (ProfIn(prof)%p(jo) /= ProfIn(prof + 1)%p(jo)) THEN
          llevels_different = .TRUE.
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  DO prof = 1, nprofiles
    IF (llevels_different .OR. prof == 1) THEN
! source levels
      pvlev(1:kni, prof) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni, prof) = log(pvlev(1:kni, prof))
! target levels
      ppo(1:kno, prof)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno, prof) = log(ppo(1:kno, prof))
!--------------------------------------------------------------------------
!1. profile interpolation user -> coef
!--------------------------------------------------------------------------
! replaced by profile assignments
!1.2 - interpolation procedure for profile elements (user -> coef)
!------------------------------------------------------------------
! for rttov same set of weights for all elements, so one call suffices
      CALL rttov_layeravg( &
            & zlnpo(1, prof),  &
            & zlnpi(1, prof),  &
            & kno,             &
            & kni,             &
            & zpz(1, 1, prof), &
            & istart(1, prof), &
            & iend(1, prof))
    ENDIF
  ENDDO
  DO kn = 1, nchannels
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF
    IF (opts%lgradp) THEN
      DO jo = 1, kno
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = 0.0_JPRB
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
          & ProfOut_k(kn)%t(jo) * ProfIn(prof)%t(istart(jo, prof):iend(jo, prof))
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
          & ProfOut_k(kn)%q(jo) * ProfIn(prof)%q(istart(jo, prof):iend(jo, prof))
        IF (opts%ozone_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%o3(jo) * ProfIn(prof)%o3(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%clw_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%clw(jo) * ProfIn(prof)%clw(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%co2_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%co2(jo) * ProfIn(prof)%co2(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%n2o_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%n2o(jo) * ProfIn(prof)%n2o(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%co_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%co(jo) * ProfIn(prof)%co(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%ch4_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            & ProfOut_k(kn)%ch4(jo) * ProfIn(prof)%ch4(istart(jo, prof):iend(jo, prof))
        ENDIF
      ENDDO
    ENDIF
    DO jo = 1, kno
      ProfIn_k(kn)%t(istart(jo, prof):iend(jo, prof)) = ProfIn_k(kn)%t(istart(jo, prof):iend(jo, prof)) +      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%t(jo)
      ProfIn_k(kn)%q(istart(jo, prof):iend(jo, prof)) = ProfIn_k(kn)%q(istart(jo, prof):iend(jo, prof)) +      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%q(jo)
      IF (opts%ozone_data) ProfIn_k(kn)%o3(istart(jo, prof):iend(jo, prof))  =      &
        & ProfIn_k(kn)%o3(istart(jo, prof):iend(jo, prof))     &
   &  +                       &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%o3(jo)
      IF (opts%clw_data  ) ProfIn_k(kn)%clw(istart(jo, prof):iend(jo, prof)) =      &
        & ProfIn_k(kn)%clw(istart(jo, prof):iend(jo, prof))     &
   &  +                      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%clw(jo)
      IF (opts%co2_data  ) ProfIn_k(kn)%co2(istart(jo, prof):iend(jo, prof)) =      &
        & ProfIn_k(kn)%co2(istart(jo, prof):iend(jo, prof))     &
   &  +                      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%co2(jo)
      IF (opts%n2o_data  ) ProfIn_k(kn)%n2o(istart(jo, prof):iend(jo, prof)) =      &
        & ProfIn_k(kn)%n2o(istart(jo, prof):iend(jo, prof))     &
   &  +                      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%n2o(jo)
      IF (opts%co_data   ) ProfIn_k(kn)%co(istart(jo, prof):iend(jo, prof))  =      &
        & ProfIn_k(kn)%co(istart(jo, prof):iend(jo, prof))     &
   &  +                       &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%co(jo)
      IF (opts%ch4_data  ) ProfIn_k(kn)%ch4(istart(jo, prof):iend(jo, prof)) =      &
        & ProfIn_k(kn)%ch4(istart(jo, prof):iend(jo, prof))     &
   &  +                      &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%ch4(jo)
    ENDDO
! target levels
    IF (opts%lgradp) THEN
      zlnpi_k = 0.0_JPRB
      zlnpo_k = 0.0_JPRB
      CALL rttov_layeravg_k( &
            & zlnpo(1, prof),  &
            & zlnpo_k,         &
            & zlnpi(1, prof),  &
            & zlnpi_k,         &
            & kno,             &
            & kni,             &
            & zpz(1, 1, prof), &
            & zpz_k,           &
            & istart(1, prof), &
            & iend(1, prof))
      pvlev_k(1:kni)        = 1 / pvlev(1:kni, prof) * zlnpi_k(1:kni)
      ProfIn_k(kn)%p(1:kni) = ProfIn_k(kn)%p(1:kni) + pvlev_k(1:kni)
    ENDIF
  ENDDO
! chan loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_K', 1_jpim, ZHOOK_HANDLE)
!------------------------------------------------------------
END SUBROUTINE rttov_intavg_prof_k
