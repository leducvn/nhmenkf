!    Compute principal component scores
SUBROUTINE rttov_pcscores_k( &
            & opts,        &
            & chanprof,    &
            & chanprof_pc, &
            & pccomp,      &
            & pccomp_k,    &
            & coef_pccomp, &
            & clear_k_pc)
!     Description:
!     To compute principal component scores
!
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
!     Method:
!     See: User's manual and scientific report for PC_RTTOV
!          (Available from EUMETSAT)
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           02/12/2009   Oroginal version: Marco Matricardi. ECMWF.
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_pccomp,      &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:)     ! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)     ! Channel indices
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_k
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  REAL(KIND=jprb)        , INTENT(INOUT) :: clear_k_pc(:, :, :)
!INTF_END
  INTEGER(KIND=jpim) :: j, k, chan, chan_pc, prof
  INTEGER(KIND=jpim) :: nprofiles                ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim) :: npcscores_p
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   CALCULATE PRINCIPAL COMPONENT SCORES
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 0_jpim, ZHOOK_HANDLE)
  npcscores_p = size(clear_k_pc, 1)
  nprofiles   = size(clear_k_pc, 3)
  nchannels   = size(chanprof)
  DO prof = 1, nprofiles
    DO j = npcscores_p, 1,  - 1
      chan_pc = chanprof_pc(j)%chan
      DO k = nchannels / nprofiles, 1,  - 1
        chan = chanprof(k)%chan
        clear_k_pc(j, k, prof) = clear_k_pc(j, k, prof) +                                                                   &
          & pccomp_k%pcscores(j + (prof - 1) * npcscores_p) * coef_pccomp%pcreg(opts%ipcreg)%coefficients(k, chan_pc) /  &
          & coef_pccomp%noise(chan)
      ENDDO
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores_k
