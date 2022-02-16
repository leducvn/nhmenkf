!
SUBROUTINE rttov_intavg_chan_ad( &
            & lgradp,      &
            & sun,         &
            & kni,         &
            & kno,         &
            & chanprof,    &
            & ProfIn,      &
            & ProfOut,     &
            & ProfOut_ad,  &
            & OpdepIn,     &
            & OpdepIn_ad,  &
            & OpdepOut,    &
            & OpdepOut_ad)
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
!              set of optical depths
!              ---------------------------
! nchannels,     number of channels this profile
! sun           solar switch by channel
! kni           number of source levs
! kno           number of target levs
! nprofiles     number of profiles
! lprofiles     channels -> profiles
! ProfIn        for target p-levs & addsolar
! ProfOut       for source p-levs
! ProfOut_ad    for source p-levs
! OpdepIn       direct for source
! OpdepOut      direct for target
! OpdepIn_ad    target for interpolator
! OpdepOut_ad   source for interpolator
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
!     3         11/2007   Extrapolate optical depths near surface so as to
!                         replicate behaviour of rttov_transmit.F90 when
!                         RTTOV levels are used throughout. - Alan Geer
!                         Niels Bormann
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
!
! imported type definitions:
  USE rttov_types, ONLY : rttov_chanprof, profile_type, opdp_path_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: lgradp
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: sun(size(chanprof))     ! switch for solar beam
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno                ! number of levels
  TYPE(profile_type  ), INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut   (size(profin))! atmospheric profiles
  TYPE(profile_type  ), INTENT(INOUT) :: ProfOut_ad(size(profin))! atmospheric profiles
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn                 ! optical depths
  TYPE(opdp_path_type), INTENT(INOUT) :: OpdepIn_ad
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut                ! optical depths
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut_ad
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_ad.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, kn, istart  (kno, size(ProfIn)), iend(kno, size(ProfIn))
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: nprofiles                                                                                    ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                                                    ! Number of radiances computed (channels used * profiles)
  LOGICAL(KIND=jplm) :: llevels_different
  REAL(KIND=jprb)    :: pvlev (kni)
  REAL(KIND=jprb)    :: zlnpi (kni)                   , zpz     (kni, kno, size(ProfIn))
  REAL(KIND=jprb)    :: zlnpo (kno)
  REAL(KIND=jprb)    :: zdp   (kno)
  REAL(KIND=jprb)    :: ppo   (kno)                   , ppo_ad  (kno)
  REAL(KIND=jprb)    :: zdp_ad(kno)                   , zlnpo_ad(kno)                   , zlnpi_ad(kni)
  REAL(KIND=jprb)    :: zpz_ad(kni, kno, size(ProfIn))
  REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
! initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_AD', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = size(ProfIn)
  nchannels         = size(chanprof)
  IF (lgradp) THEN
    llevels_different = .TRUE.
  ELSE
    DO prof = 1, nprofiles - 1
      DO jo = 1, kno
        IF (ProfOut(prof)%p(jo) /= ProfOut(prof + 1)%p(jo)) THEN
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
! interpolation
!--------------------
      CALL rttov_layeravg( &
            & zlnpo,           &
            & zlnpi,           &
            & kno,             &
            & kni,             &
            & zpz(1, 1, prof), &
            & istart(1, prof), &
            & iend(1, prof))
! The layer averaging tends towards a constant value when target levels are
! below the lowest RTTOV level. We need instead to extrapolate optical depths
! so as replicate behaviour of rttov_transmit.F90 when RTTOV levels are used
! throughout.
      WHERE (ppo(:) > pvlev(kni))
        zdp(:)                = (ppo(:) - pvlev(kni)) / (pvlev(kni) - pvlev(kni - 1))
        zpz(kni - 1, :, prof) =  - zdp(:)
        zpz(kni, :, prof)     = zdp(:) + 1.0_jprb
        istart(:, prof)       = kni - 1
        iend(:, prof)         = kni
      ENDWHERE
    ENDIF
  ENDDO
  IF (lgradp) THEN
    zpz_ad = 0.0_JPRB
    DO kn = 1, nchannels! loop over channels for profile passed in this time
      prof = chanprof(kn)%prof
      DO jo = 1, kno! o/p levels
! for viewing path of radiation generated by atmosphere
! nb for solar beam
! weights for any channel of sun_level same as for any channel of atm_level
! nb there could be a problem with variable incidence angles from raytracing
        zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) = zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) +      &
          & OpdepIn%atm_level(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_ad%atm_level(jo, kn)
        IF (sun(kn)) zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) =      &
   & zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) +              &
   & OpdepIn%sun_level(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_ad%sun_level(jo, kn)
      ENDDO
! target levels
    ENDDO
! channels
  ENDIF
  DO kn = 1, nchannels! loop over channels for profile passed in this time
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF
    DO jo = 1, kno
! for fwd model weights are given by zpz (ignore zps and zpzps)
      OpdepIn_ad%atm_level(istart(jo, prof):iend(jo, prof), kn) =      &
        & OpdepIn_ad%atm_level(istart(jo, prof):iend(jo, prof), kn) +  &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_ad%atm_level(jo, kn)
      IF (sun(kn)) OpdepIn_ad%sun_level(istart(jo, prof):iend(jo, prof), kn) =      &
   & OpdepIn_ad%sun_level(istart(jo, prof):iend(jo, prof), kn) +              &
   & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_ad%sun_level(jo, kn)
    ENDDO
! target levels
  ENDDO
! channels
  IF (lgradp) THEN
    DO prof = 1, nprofiles
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = log(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = log(ppo(1:kno))
      ppo_ad(:)    = 0.0_JPRB
      WHERE (ppo(:) > pvlev(kni))
        zdp_ad(:)                = zpz_ad(kni, :, prof)
        zpz_ad(kni, :, prof)     = 0.0_JPRB
        zdp_ad(:)                = zdp_ad(:) - zpz_ad(kni - 1, :, prof)
        zpz_ad(kni - 1, :, prof) = 0.0_JPRB
        ppo_ad(:)                = zdp_ad(:) / (pvlev(kni) - pvlev(kni - 1))
      ENDWHERE
      zlnpi_ad = 0.0_JPRB
      zlnpo_ad = 0.0_JPRB
      CALL rttov_layeravg_ad( &
            & zlnpo,              &
            & zlnpo_ad,           &
            & zlnpi,              &
            & zlnpi_ad,           &
            & kno,                &
            & kni,                &
            & zpz(1, 1, prof),    &
            & zpz_ad(1, 1, prof), &
            & istart(1, prof),    &
            & iend(1, prof))
      ppo_ad(1:kno)             = ppo_ad(1:kno) + 1._JPRB / ppo(1:kno) * zlnpo_ad(1:kno)
      ProfOut_ad(prof)%p(1:kno) = ProfOut_ad(prof)%p(1:kno) + ppo_ad(1:kno)
    ENDDO
  ENDIF
!
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_AD', 1_jpim, ZHOOK_HANDLE)
!
END SUBROUTINE rttov_intavg_chan_ad
