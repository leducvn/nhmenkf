!
SUBROUTINE rttov_intavg_chan_k( &
            & lgradp,     &
            & sun,        &
            & kni,        &
            & kno,        &
            & chanprof,   &
            & ProfIn,     &
            & ProfOut,    &
            & ProfOut_k,  &
            & OpdepIn,    &
            & OpdepIn_k,  &
            & OpdepOut,   &
            & OpdepOut_k)
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
! OpdepIn        direct for source
! OpdepIn_k      source for interpolator
! OpdepOut       direct for target
! OpdepOut_k     target for interpolator
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
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: sun(size(chanprof))      ! switch for solar beam
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno                 ! number of levels
  TYPE(profile_type  ), INTENT(IN)    :: ProfIn   (:)             ! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut  (size(ProfIn)  )! atmospheric profiles
  TYPE(profile_type  ), INTENT(INOUT) :: ProfOut_k(size(chanprof))! atmospheric profiles
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn                  ! optical depths
  TYPE(opdp_path_type), INTENT(INOUT) :: OpdepIn_k
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut                 ! optical depths
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut_k
!INTF_END
#include "rttov_layeravg.h"
#include "rttov_layeravg_k.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, kn, istart (kno, size(ProfIn)), iend(kno, size(ProfIn))
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: nprofiles                                                              ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                              ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=jprb) :: pvlev(kni)
  REAL   (KIND=jprb) :: zlnpi(kni), zpz    (kni, kno, size(ProfIn))
  REAL   (KIND=jprb) :: zlnpo(kno)
  REAL   (KIND=jprb) :: zdp  (kno)
  REAL   (KIND=jprb) :: ppo  (kno), ppo_k  (kno)
  REAL   (KIND=jprb) :: zdp_k(kno), zlnpo_k(kno)                   , zlnpi_k(kni)
  REAL(KIND=jprb), ALLOCATABLE :: zpz_k(:, :, :, :)
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
! initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_K', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = size(ProfIn)
  nchannels         = size(chanprof)
  IF (lgradp) THEN
    llevels_different = .TRUE.
    ALLOCATE (zpz_k(kni, kno, nchannels, nprofiles))
  ELSE
    ALLOCATE (zpz_k(1, 1, 1, 1))
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
!
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
    zpz_k = 0.0_JPRB
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
        zpz_k(istart(jo, prof):iend(jo, prof), jo, kn, prof) = zpz_k(istart(jo, prof):iend(jo, prof), jo, kn, prof) +      &
          & OpdepIn%atm_level(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_k%atm_level(jo, kn)
        IF (sun(kn)) zpz_k(istart(jo, prof):iend(jo, prof), jo, kn, prof) =      &
   & zpz_k(istart(jo, prof):iend(jo, prof), jo, kn, prof) +              &
   & OpdepIn%sun_level(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_k%sun_level(jo, kn)
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
      OpdepIn_k%atm_level(istart(jo, prof):iend(jo, prof), kn) =      &
        & OpdepIn_k%atm_level(istart(jo, prof):iend(jo, prof), kn) +  &
        & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_k%atm_level(jo, kn)
      IF (sun(kn)) OpdepIn_k%sun_level(istart(jo, prof):iend(jo, prof), kn) =      &
   & OpdepIn_k%sun_level(istart(jo, prof):iend(jo, prof), kn) +              &
   & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_k%sun_level(jo, kn)
    ENDDO
! target levels
  ENDDO
! channels
  IF (lgradp) THEN
    DO kn = 1, nchannels! loop over channels for profile passed in this time
      IF (llevels_different) THEN
        prof = chanprof(kn)%prof
      ELSE
        prof = 1
      ENDIF
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = log(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = log(ppo(1:kno))
      ppo_k(:)     = 0.0_JPRB
      WHERE (ppo(:) > pvlev(kni))
        zdp_k(:)                    = zpz_k(kni, :, kn, prof)
        zpz_k(kni, :, kn, prof)     = 0.0_JPRB
        zdp_k(:)                    = zdp_k(:) - zpz_k(kni - 1, :, kn, prof)
        zpz_k(kni - 1, :, kn, prof) = 0.0_JPRB
        ppo_k(:)                    = zdp_k(:) / (pvlev(kni) - pvlev(kni - 1))
      ENDWHERE
      zlnpi_k = 0.0_JPRB
      zlnpo_k = 0.0_JPRB
      CALL rttov_layeravg_k( &
            & zlnpo,                 &
            & zlnpo_k,               &
            & zlnpi,                 &
            & zlnpi_k,               &
            & kno,                   &
            & kni,                   &
            & zpz(1, 1, prof),       &
            & zpz_k(1, 1, kn, prof), &
            & istart(1, prof),       &
            & iend(1, prof))
      ppo_k(1:kno)           = ppo_k(1:kno) + 1._JPRB / ppo(1:kno) * zlnpo_k(1:kno)
      ProfOut_k(kn)%p(1:kno) = ProfOut_k(kn)%p(1:kno) + ppo_k(1:kno)
    ENDDO
  ENDIF
  DEALLOCATE (zpz_k)
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_K', 1_jpim, ZHOOK_HANDLE)
!
!
END SUBROUTINE rttov_intavg_chan_k
