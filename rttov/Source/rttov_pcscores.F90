!    Compute principal component scores
SUBROUTINE rttov_pcscores( &
            & opts,         &
            & chanprof,     &
            & chanprof_pc,  &
            & pccomp,       &
            & coef_pccomp,  &
            & radiancedata)
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
       & rttov_coef_pccomp, &
       & radiance_Type,     &
       & profile_Type
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:)! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)! Channel indices
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(radiance_Type    ), INTENT(INOUT) :: radiancedata
!INTF_END
  INTEGER(KIND=jpim) :: i, j, k, chan, chan_pc, prof, prof_pc
  INTEGER(KIND=jpim) :: nchannels                                 ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   CALCULATE PRINCIPAL COMPONENT SCORES
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  npcscores = size(chanprof_pc)
  nprofiles = chanprof(nchannels)%prof
  
  pccomp%pcscores = 0._jprb
  k = 1
  DO prof = 1, nprofiles
    DO i = 1, nchannels/nprofiles
      chan = chanprof(i)%chan
      DO j = 1 + (prof-1) * npcscores/nprofiles, prof * npcscores/nprofiles
        chan_pc = chanprof_pc(j)%chan
        pccomp%pcscores(j) = pccomp%pcscores(j) +      &
          & coef_pccomp%pcreg(opts%ipcreg)%coefficients(i, chan_pc) * radiancedata%clear(k) / coef_pccomp%noise(chan)
      ENDDO
      k = k + 1
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores
