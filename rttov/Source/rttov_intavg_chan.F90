!
SUBROUTINE rttov_intavg_chan( &
            & lgradp,   &
            & sun,      &
            & kni,      &
            & kno,      &
            & chanprof, &
            & ProfIn,   &
            & ProfOut,  &
            & OpdepIn,  &
            & OpdepOut)
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
!               set of optical depths
!               ---------------------------
! nchannels      number of channels
! sun            solar switch by channel
! kni            number of source levs
! kno            number of target levs
! nprofiles      number of profiles
! lprofiles      channel -> profiles
! ProfIn         need source p-levels
! ProfOut        need target p-levels & addsolar
! OpdepIn        source for interpolator
! OpdepOut       target for interpolator
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
!     1         12/2006   main code rttov_intavg restructured Peter Rayer
!
!     2         10/2007   Optimisation Deborah Salmond and Pascal Brunel
!
!     3         11/2007   Extrapolate optical depths near surface so as to
!                         approximate behaviour of rttov_transmit.F90 when
!                         RTTOV levels are used throughout. - Alan Geer
!     4         01/2008   Cleaning, Niels Bormann
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!-----------------------------------------------------------------------
! imported type definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & rttov_chanprof, &
       & profile_type,   &
       & opdp_path_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: lgradp
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: sun(size(chanprof))! switch for solar beam
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno           ! number of levels
  TYPE(profile_type  ), INTENT(IN)    :: ProfIn (:)         ! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut(:)         ! atmospheric profiles
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn            ! optical depths
  TYPE(opdp_path_type), INTENT(INOUT) :: OpdepOut           ! optical depths
!INTF_END
#include "rttov_layeravg.h"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, kn, istart(kno, size(ProfIn)), iend(kno, size(ProfIn))
  INTEGER(KIND=jpim) :: prof
  LOGICAL(KIND=jplm) :: llevels_different
  REAL   (KIND=jprb) :: pvlev(kni)
  REAL   (KIND=jprb) :: ppo  (kno)
  REAL   (KIND=jprb) :: zlnpi(kni), zpz(kni, kno, size(ProfIn))
  REAL   (KIND=jprb) :: zlnpo(kno)
  REAL   (KIND=jprb) :: zdp  (kno)
  INTEGER(KIND=jpim) :: nprofiles                                                 ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                 ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN', 0_jpim, ZHOOK_HANDLE)
  nprofiles         = size(ProfIn)
  nchannels         = size(chanprof)
  llevels_different = .FALSE.
  DO prof = 1, nprofiles - 1
    DO jo = 1, kno
      IF (ProfOut(prof)%p(jo) /= ProfOut(prof + 1)%p(jo)) THEN
        llevels_different = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDDO
!2.2 - interpolation procedure for optical depths (coef -> user)
!---------------------------------------------------------------
!     nb same Subroutine rttov_layeravg called by dir/tl/ad/k
  DO prof = 1, nprofiles
    IF (llevels_different .OR. prof == 1) THEN
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = log(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = log(ppo(1:kno))
! for rttov same set of weights for any channel, so one call suffices
! for direct run weights are given by zpz
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
  DO kn = 1, nchannels! loop over channels for profile passed in this time
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF
    DO jo = 1, kno! o/p levels
! for viewing path of radiation generated by atmosphere
      OpdepOut%atm_level(jo, kn) =      &
        & sum(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn%atm_level(istart(jo, prof):iend(jo, prof), kn))
! for solar beam
! weights for any channel of sun_level same as for any channel of atm_level
! nb there could be a problem with variable incidence angles from raytracing
      IF (sun(kn)) OpdepOut%sun_level(jo, kn) =      &
   & sum(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepIn%sun_level(istart(jo, prof):iend(jo, prof), kn))
    ENDDO
! target levels
  ENDDO
! channels
!------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN', 1_jpim, ZHOOK_HANDLE)
!
END SUBROUTINE rttov_intavg_chan
