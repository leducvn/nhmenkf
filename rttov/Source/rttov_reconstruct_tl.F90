!    Reconstruc radiances from principal component scores
SUBROUTINE rttov_reconstruct_tl( &
            & chanprof_in, &
            & chanprof_pc, &
            & pccomp,      &
            & pccomp_tl,   &
            & coef_pccomp)
!     Description:
!     Reconstructs radiances from principal component scores
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
!     1           02/12/2009   Orginal version: Marco Matricardi. ECMWF.
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_tl
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
  INTEGER(KIND=jpim) :: i, j, k, prof, chan, chan_pc, prof_pc
  REAL   (KIND=jprb) :: rad
  INTEGER(KIND=jpim) :: nchannels_rec                                    ! Number of output radiances
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   CALCULATE PRINCIPAL COMPONENT SCORES
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels_rec = size(chanprof_in)
  npcscores     = size(chanprof_pc)
  nprofiles = chanprof_in(nchannels_rec)%prof
  k = 1
  DO prof = 1, nprofiles
    DO i = 1, nchannels_rec/nprofiles
      chan = chanprof_in(i)%chan
      rad = 0
      DO j = 1 + (prof-1) * npcscores/nprofiles, prof * npcscores/nprofiles
        chan_pc = chanprof_pc(j)%chan
        rad = rad + coef_pccomp%eigenvectors(chan, chan_pc) * pccomp_tl%pcscores(j)
      ENDDO
      pccomp_tl%clear_pccomp(k) = rad * coef_pccomp%noise_in(chan)
      k = k+1
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_reconstruct_tl
