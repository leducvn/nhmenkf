SUBROUTINE rttov_intavg_prof_ad( &
            & opts,       &
            & kni,        &
            & kno,        &
            & ProfIn,     &
            & ProfIn_ad,  &
            & ProfOut,    &
            & ProfOut_ad)
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
! ProfIn_ad     source for interpolator
! ProfOut       target for interpolator
! ProfOut_ad    target for interpolator
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
  TYPE(profile_type) , INTENT(INOUT) :: ProfIn_ad (size(ProfIn))
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut_ad(size(ProfIn))
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_ad.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, istart  (kno), iend(kno     )
  INTEGER(KIND=jpim) :: prof
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles                                                     ! Number of profiles
  REAL   (KIND=jprb) :: pvlev   (kni)
  REAL   (KIND=jprb) :: pvlev_ad(kni)
  REAL   (KIND=jprb) :: ppo     (kno)
  REAL   (KIND=jprb) :: zlnpi   (kni), zlnpi_ad(kni), zpz (kni, kno), zpz_ad(kni, kno)
  REAL   (KIND=jprb) :: zlnpo   (kno), zlnpo_ad(kno)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_AD', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = size(ProfIn)
  IF (opts%lgradp) THEN
    llevels_different = .TRUE.
    zpz_ad            = 0._jprb
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
      CALL rttov_layeravg( &
            & zlnpo,  &
            & zlnpi,  &
            & kno,    &
            & kni,    &
            & zpz,    &
            & istart, &
            & iend)
    ENDIF
    IF (opts%lgradp) THEN
      DO jo = 1, kno
        zpz_ad(istart(jo):iend(jo), jo) = 0.0_JPRB
        zpz_ad(istart(jo):iend(jo), jo) =      &
          & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%t(jo) * ProfIn(prof)%t(istart(jo):iend(jo))
        zpz_ad(istart(jo):iend(jo), jo) =      &
          & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%q(jo) * ProfIn(prof)%q(istart(jo):iend(jo))
        IF (opts%ozone_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%o3(jo) * ProfIn(prof)%o3(istart(jo):iend(jo))
        ENDIF
        IF (opts%clw_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%clw(jo) * ProfIn(prof)%clw(istart(jo):iend(jo))
        ENDIF
        IF (opts%co2_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%co2(jo) * ProfIn(prof)%co2(istart(jo):iend(jo))
        ENDIF
        IF (opts%n2o_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%n2o(jo) * ProfIn(prof)%n2o(istart(jo):iend(jo))
        ENDIF
        IF (opts%co_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%co(jo) * ProfIn(prof)%co(istart(jo):iend(jo))
        ENDIF
        IF (opts%ch4_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
            & zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%ch4(jo) * ProfIn(prof)%ch4(istart(jo):iend(jo))
        ENDIF
      ENDDO
    ENDIF
    DO jo = 1, kno
      ProfIn_ad(prof)%t(istart(jo):iend(jo)) =      &
        & ProfIn_ad(prof)%t(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%t(jo)
      ProfIn_ad(prof)%q(istart(jo):iend(jo)) =      &
        & ProfIn_ad(prof)%q(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%q(jo)
      IF (opts%ozone_data) ProfIn_ad(prof)%o3(istart(jo):iend(jo))  =      &
   & ProfIn_ad(prof)%o3(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%o3(jo)
      IF (opts%clw_data  ) ProfIn_ad(prof)%clw(istart(jo):iend(jo)) =      &
   & ProfIn_ad(prof)%clw(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%clw(jo)
      IF (opts%co2_data  ) ProfIn_ad(prof)%co2(istart(jo):iend(jo)) =      &
   & ProfIn_ad(prof)%co2(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%co2(jo)
      IF (opts%n2o_data  ) ProfIn_ad(prof)%n2o(istart(jo):iend(jo)) =      &
   & ProfIn_ad(prof)%n2o(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%n2o(jo)
      IF (opts%co_data   ) ProfIn_ad(prof)%co(istart(jo):iend(jo))  =      &
   & ProfIn_ad(prof)%co(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%co(jo)
      IF (opts%ch4_data  ) ProfIn_ad(prof)%ch4(istart(jo):iend(jo)) =      &
   & ProfIn_ad(prof)%ch4(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%ch4(jo)
    ENDDO
! target levels
    IF (opts%lgradp) THEN
      zlnpi_ad = 0.0_JPRB
      zlnpo_ad = 0.0_JPRB
      CALL rttov_layeravg_ad( &
            & zlnpo,    &
            & zlnpo_ad, &
            & zlnpi,    &
            & zlnpi_ad, &
            & kno,      &
            & kni,      &
            & zpz,      &
            & zpz_ad,   &
            & istart,   &
            & iend)
      pvlev_ad(1:kni)          = 1 / pvlev(1:kni) * zlnpi_ad(1:kni)
      ProfIn_ad(prof)%p(1:kni) = ProfIn_ad(prof)%p(1:kni) + pvlev_ad(1:kni)
    ENDIF
  ENDDO
! profiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_AD', 1_jpim, ZHOOK_HANDLE)
!------------------------------------------------------------
END SUBROUTINE rttov_intavg_prof_ad
