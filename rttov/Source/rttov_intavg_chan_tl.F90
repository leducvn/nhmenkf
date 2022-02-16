!
SUBROUTINE rttov_intavg_chan_tl( &
            & lgradp,      &
            & sun,         &
            & kni,         &
            & kno,         &
            & chanprof,    &
            & ProfIn,      &
            & ProfOut,     &
            & ProfOut_tl,  &
            & OpdepIn,     &
            & OpdepIn_tl,  &
            & OpdepOut,    &
            & OpdepOut_tl)
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
!                set of optical depths
!                ---------------------------
! nchannels,     number of channels this profile
! sun            solar switch by channel
! kni            number of source levs
! kno            number of target levs
! nprofiles      number of profiles
! lprofiles      channels -> profiles
! ProfIn         for source p-levs & addsolar
! ProfOut        for target p-levs
! ProfOut_tl     for target p-levs
! OpdepIn        direct for source
! OpdepIn_tl     source for interpolator
! OpdepOut       direct for target
! OpdepOut_tl    target for interpolator
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
!
!     4         01/2008   lgradp option for TL/AD for user pressure levels,
!                         Niels Bormann
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
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
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut_tl(size(ProfIn))! TL of atmospheric profiles
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn                 ! optical depths
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn_tl
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut                ! optical depths
  TYPE(opdp_path_type), INTENT(INOUT) :: OpdepOut_tl
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_tl.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, kn, istart(kno, size(ProfIn)), iend(kno, size(ProfIn))
  INTEGER(KIND=jpim) :: prof
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles                                                                                    ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                                                    ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=jprb) :: pvlev   (kni)
  REAL   (KIND=jprb) :: zlnpi   (kni)                   , zpz     (kni, kno, size(ProfIn))
  REAL   (KIND=jprb) :: zpz_tl  (kni, kno, size(ProfIn))
  REAL   (KIND=jprb) :: ppo     (kno)
  REAL   (KIND=jprb) :: zlnpi_tl(kni)                   , zlnpo_tl(kno)                   , ppo_tl(kno)
  REAL   (KIND=jprb) :: zlnpo   (kno)
  REAL   (KIND=jprb) :: zdp     (kno)                   , zdp_tl  (kno)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_TL', 0_jpim, ZHOOK_HANDLE)
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
!-------------------------------------------
!2. optical depth interpolation coef -> user
!-------------------------------------------
!2.2 - interpolation procedure for optical depths (coef -> user)
!     nb same Subroutine layeravg called by dir/tl/ad/k
!---------------------------------------------------------------
  DO prof = 1, nprofiles
    IF (llevels_different .OR. prof == 1) THEN
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = log(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = log(ppo(1:kno))
      IF (lgradp) THEN
        ppo_tl(1:kno)   = ProfOut_tl(prof)%p(1:kno)
        zlnpo_tl(1:kno) = 1._JPRB / ppo(1:kno) * ppo_tl(1:kno)
        zlnpi_tl        = 0.0_JPRB
        CALL rttov_layeravg_tl( &
              & zlnpo,              &
              & zlnpo_tl,           &
              & zlnpi,              &
              & zlnpi_tl,           &
              & kno,                &
              & kni,                &
              & zpz(1, 1, prof),    &
              & zpz_tl(1, 1, prof), &
              & istart(1, prof),    &
              & iend(1, prof))
      ELSE
        CALL rttov_layeravg( &
              & zlnpo,           &
              & zlnpi,           &
              & kno,             &
              & kni,             &
              & zpz(1, 1, prof), &
              & istart(1, prof), &
              & iend(1, prof))
      ENDIF
! The layer averaging tends towards a constant value when target levels are
! below the lowest RTTOV level. We need instead to extrapolate optical depths
! so as replicate behaviour of rttov_transmit.F90 when RTTOV levels are used
! throughout.
      IF (lgradp) THEN
        WHERE (ppo(:) > pvlev(kni))
          zdp_tl(:)                = ppo_tl(:) / (pvlev(kni) - pvlev(kni - 1))
          zpz_tl(kni - 1, :, prof) =  - zdp_tl(:)
          zpz_tl(kni, :, prof)     = zdp_tl(:)
        ENDWHERE
      ENDIF
      WHERE (ppo(:) > pvlev(kni))
        zdp(:)                = (ppo(:) - pvlev(kni)) / (pvlev(kni) - pvlev(kni - 1))
        zpz(kni - 1, :, prof) =  - zdp(:)
        zpz(kni, :, prof)     = zdp(:) + 1.0_jprb
        istart(:, prof)       = kni - 1
        iend(:, prof)         = kni
      ENDWHERE
    ENDIF
  ENDDO
  DO kn = 1, nchannels! loop over channels for profile passed in this time
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF
    DO jo = 1, kno! o/p levels
! for viewing path of radiation generated by atmosphere
! nb for solar beam
! weights for any channel of sun_level same as for any channel of atm_level
! nb there could be a problem with variable incidence angles from raytracing
      OpdepOut_tl%atm_level(jo, kn) =      &
        & sum(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn_tl%atm_level(istart(jo, prof):iend(jo, prof), kn))
      IF (sun(kn)) OpdepOut_tl%sun_level(jo, kn) =      &
   & sum(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn_tl%sun_level(istart(jo, prof):iend(jo, prof), kn))
    ENDDO
! target levels
  ENDDO
! channels
  IF (lgradp) THEN
    DO kn = 1, nchannels! loop over channels for profile passed in this time
      prof = chanprof(kn)%prof
      DO jo = 1, kno! o/p levels
! for viewing path of radiation generated by atmosphere
! nb for solar beam
! weights for any channel of sun_level same as for any channel of atm_level
! nb there could be a problem with variable incidence angles from raytracing
        OpdepOut_tl%atm_level(jo, kn) = OpdepOut_tl%atm_level(jo, kn) +      &
          & sum(zpz_tl(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn%atm_level(istart(jo, prof):iend(jo, prof), kn))
        IF (sun(kn)) OpdepOut_tl%sun_level(jo, kn) = OpdepOut_tl%sun_level(jo, kn) +      &
   & sum(zpz_tl(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn%sun_level(istart(jo, prof):iend(jo, prof), kn))
      ENDDO
! target levels
    ENDDO
! channels
  ENDIF
!------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_TL', 1_jpim, ZHOOK_HANDLE)
!
END SUBROUTINE rttov_intavg_chan_tl
