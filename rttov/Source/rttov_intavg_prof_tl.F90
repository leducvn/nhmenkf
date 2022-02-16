SUBROUTINE rttov_intavg_prof_tl( &
            & opts,       &
            & kni,        &
            & kno,        &
            & ProfIn,     &
            & ProfIn_tl,  &
            & ProfOut,    &
            & ProfOut_tl)
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
! ProfIn_tl     source for interpolator
! ProfOut       target for interpolator
! ProfOut_tl    target for interpolator
!
!--------------------------------------------------------------------------
! application within RTTOV (for NWP SAF)
!   - profiles on USER levels -> prediction of gas/atm opdeps on COEF levels
!   - gas/atm opdeps on COEF levels -> transmittances on USER levels for RTE
!---------------------------------------------------------------------------
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
  USE rttov_types, ONLY : rttov_options, profile_type
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_options), INTENT(IN)    :: opts
  INTEGER(KIND=jpim) , INTENT(IN)    :: kni , kno               ! number of levels
  TYPE(profile_type) , INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(profile_type) , INTENT(IN)    :: ProfIn_tl (size(ProfIn))
  TYPE(profile_type) , INTENT(IN)    :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut_tl(size(ProfIn))
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_tl.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, istart  (kno), iend(kno     )
  INTEGER(KIND=jpim) :: prof
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles                                                     ! Number of profiles
  REAL   (KIND=jprb) :: pvlev   (kni)
  REAL   (KIND=jprb) :: pvlev_tl(kni)
  REAL   (KIND=jprb) :: ppo     (kno)
  REAL   (KIND=jprb) :: zlnpi   (kni), zlnpi_tl(kni), zpz (kni, kno), zpz_tl(kni, kno)
  REAL   (KIND=jprb) :: zlnpo   (kno), zlnpo_tl(kno)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_TL', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = size(ProfIn)
  IF (opts%lgradp) THEN
    llevels_different = .TRUE.
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
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = log(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = log(ppo(1:kno))
!--------------------------------------------------------------------------
!1. profile interpolation user -> coef
!--------------------------------------------------------------------------
! replaced by profile assignments
!1.2 - interpolation procedure for profile elements (user -> coef)
!------------------------------------------------------------------
! for rttov same set of weights for all elements, so one call suffices
      IF (opts%lgradp) THEN
        pvlev_tl(1:kni) = ProfIn_tl(prof)%p(1:kni)
        zlnpi_tl(1:kni) = 1. / pvlev(1:kni) * pvlev_tl(1:kni)
        zlnpo_tl(1:kno) = 0.0_JPRB
        CALL rttov_layeravg_tl( &
              & zlnpo,    &
              & zlnpo_tl, &
              & zlnpi,    &
              & zlnpi_tl, &
              & kno,      &
              & kni,      &
              & zpz,      &
              & zpz_tl,   &
              & istart,   &
              & iend)
      ELSE
        CALL rttov_layeravg( &
              & zlnpo,  &
              & zlnpi,  &
              & kno,    &
              & kni,    &
              & zpz,    &
              & istart, &
              & iend)
      ENDIF
    ENDIF
    DO jo = 1, kno
! TL
      ProfOut_tl(prof)%t(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%t(istart(jo):iend(jo))))
      ProfOut_tl(prof)%q(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%q(istart(jo):iend(jo))))
      IF (opts%ozone_data) THEN
        ProfOut_tl(prof)%o3(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%o3(istart(jo):iend(jo))))
      ENDIF
      IF (opts%clw_data) THEN
        ProfOut_tl(prof)%clw(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%clw(istart(jo):iend(jo))))
      ENDIF
      IF (opts%co2_data) THEN
        ProfOut_tl(prof)%co2(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%co2(istart(jo):iend(jo))))
      ENDIF
      IF (opts%n2o_data) THEN
        ProfOut_tl(prof)%n2o(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%n2o(istart(jo):iend(jo))))
      ENDIF
      IF (opts%co_data) THEN
        ProfOut_tl(prof)%co(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%co(istart(jo):iend(jo))))
      ENDIF
      IF (opts%ch4_data) THEN
        ProfOut_tl(prof)%ch4(jo) = sum(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%ch4(istart(jo):iend(jo))))
      ENDIF
      IF (opts%lgradp) THEN
        ProfOut_tl(prof)%t(jo) =      &
          & ProfOut_tl(prof)%t(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%t(istart(jo):iend(jo))))
        ProfOut_tl(prof)%q(jo) =      &
          & ProfOut_tl(prof)%q(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%q(istart(jo):iend(jo))))
        IF (opts%ozone_data) THEN
          ProfOut_tl(prof)%o3(jo) =      &
            & ProfOut_tl(prof)%o3(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%o3(istart(jo):iend(jo))))
        ENDIF
        IF (opts%clw_data) THEN
          ProfOut_tl(prof)%clw(jo) =      &
            & ProfOut_tl(prof)%clw(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%clw(istart(jo):iend(jo))))
        ENDIF
        IF (opts%co2_data) THEN
          ProfOut_tl(prof)%co2(jo) =      &
            & ProfOut_tl(prof)%co2(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%co2(istart(jo):iend(jo))))
        ENDIF
        IF (opts%n2o_data) THEN
          ProfOut_tl(prof)%n2o(jo) =      &
            & ProfOut_tl(prof)%n2o(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%n2o(istart(jo):iend(jo))))
        ENDIF
        IF (opts%co_data) THEN
          ProfOut_tl(prof)%co(jo) =      &
            & ProfOut_tl(prof)%co(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%co(istart(jo):iend(jo))))
        ENDIF
        IF (opts%ch4_data) THEN
          ProfOut_tl(prof)%ch4(jo) =      &
            & ProfOut_tl(prof)%ch4(jo) + sum(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%ch4(istart(jo):iend(jo))))
        ENDIF
      ENDIF
    ENDDO
! target levels
  ENDDO
! profiles
!------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_TL', 1_jpim, ZHOOK_HANDLE)
!
END SUBROUTINE rttov_intavg_prof_tl
