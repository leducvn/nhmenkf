! $Id: ropp_qc_bgqc.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* QC/ropp_qc_bgqc *
!
! NAME
!    ropp_qc_bgqc - Background quality control.
!
! SYNOPSIS
!    call ropp_qc_bgqc(obs, config, diag)
! 
! DESCRIPTION
!    This subroutine performs a background quality control on bending
!    angle or refractivity observation data.
!
! INPUTS
!    obs         Observation vector
!    config      Configuration options
!    diag        Diagnostic output structure
!
! OUTPUT
!    diag        Diagnostic structure with updated variables
!
! NOTES
!    ropp_qc_bgqc requires that some fields in the diag structure have
!    been properly filled, e.g. by a previous call to ropp_qc_OmB().
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_qc_OmB
!    ropp_qc_pge
!
! REFERENCES
!
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_bgqc_1dbangle(obs, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV
  USE ropp_fm_types
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DBangle)              :: obs
  TYPE(VarConfig)                :: config
  TYPE(VarDiag)                  :: diag

  INTEGER                        :: n_good
  INTEGER                        :: n_bad = 0
  INTEGER, DIMENSION(:), POINTER :: idx => null()

  CHARACTER(len =  10)           :: istr
  CHARACTER(len =  10)           :: rstr
  CHARACTER(len = 256)           :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_bgqc')

! 1.3 Find bad data points
! ------------------------

  n_good =  COUNT(obs%bangle > ropp_MDTV .AND. obs%weights > 0.0_wp)
  idx    => WHERE(obs%bangle > ropp_MDTV .AND. obs%weights > 0.0_wp .AND. &
              ABS(diag%OmB) > MIN(config%bgqc_reject_factor * diag%OmB_sigma, &
              config%bgqc_reject_factor), n_bad)

! 1.4 Set weights of rejected data points to 0
! --------------------------------------------

  IF (config % bgqc_apply) THEN
     IF (n_bad > 0) THEN
        obs % weights(idx) = 0.0_wp
        WRITE(istr, '(i10)')   n_bad                      
        istr = ADJUSTL(istr)
        WRITE(rstr, '(f10.1)') config%bgqc_reject_factor  
        rstr = ADJUSTL(rstr)
        CALL message(msg_info, &
            "Background quality control removes " // TRIM(istr) //         &
            " bending angle data points from the \n" //                    &
            "   observations as their deviation from the background " //   &
            " exceeds " // TRIM(rstr) // "\n" //                           & 
            "   times the expected (1-sigma) error.\n")
     ELSE
        CALL message(msg_info, &
             "Background quality control lets all bending angle values pass.\n")
     ENDIF
  ENDIF

! 1.5 Set diagnostic variables
! ----------------------------

  diag % n_data        = n_good
  diag % n_bgqc_reject = n_bad

  IF (config % bgqc_apply) THEN
     IF (REAL(n_bad)/REAL(n_good) < 0.01 * config%bgqc_reject_max_percent) THEN
        diag % ok = .TRUE.
     ELSE
        diag % ok = .FALSE.
        WRITE(rstr, '(f10.1)') 100. * REAL(n_bad)/REAL(n_good)  
        rstr = ADJUSTL(rstr)
        CALL message(msg_error, &
             "Background quality control has rejected " // TRIM(rstr) //    &
             "% of the data points of this profile.\n" //                   &
             "   The entire profile is flagged as ""not ok"".\n")
     ENDIF
  ENDIF

! 1.6 Clean up
! ------------

  DEALLOCATE(idx)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_bgqc_1dbangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_bgqc_1drefrac(obs, config, diag)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV
  USE ropp_fm_types
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DRefrac)              :: obs
  TYPE(VarConfig)                :: config
  TYPE(VarDiag)                  :: diag

  INTEGER                        :: n_good
  INTEGER                        :: n_bad = 0
  INTEGER, DIMENSION(:), POINTER :: idx => null()

  CHARACTER(len =  10)           :: istr
  CHARACTER(len =  10)           :: rstr
  CHARACTER(len = 256)           :: routine

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_bgqc')

! 2.3 Find bad data points
! ------------------------

  n_good =  COUNT(obs%refrac > ropp_MDTV .AND. obs%weights > 0.0_wp)
  idx    => WHERE(obs%refrac > ropp_MDTV .AND. obs%weights > 0.0_wp .AND. &
              ABS(diag%OmB) > MIN(config%bgqc_reject_factor * diag%OmB_sigma, &
              config%bgqc_reject_factor), n_bad)

! 2.4 Set weights of rejected data points to 0
! --------------------------------------------

  IF (config % bgqc_apply) THEN
     IF (n_bad > 0) THEN
        obs % weights(idx) = 0.0_wp
        WRITE(istr, '(i10)')   n_bad                      
        istr = ADJUSTL(istr)
        WRITE(rstr, '(f10.1)') config%bgqc_reject_factor  
        rstr = ADJUSTL(rstr)
        CALL message(msg_info, &
             "Background quality control removes " // TRIM(istr) //         &
             " refractivity data points from the\n" //                      &
             "   observations as their deviation from the background " //   &
             " exceeds " // TRIM(rstr) // "\n" //                           &
             "   times the expected (1-sigma) error.\n")
     ELSE
        CALL message(msg_info, &
             "Background quality control lets all refractivity values pass.\n")
     ENDIF
  ENDIF

! 2.5 Set diagnostic variables
! ----------------------------

  diag % n_data        = n_good
  diag % n_bgqc_reject = n_bad

  IF (config % bgqc_apply) THEN
     IF (REAL(n_bad)/REAL(n_good) < 0.01 * config%bgqc_reject_max_percent) THEN
        diag % ok = .TRUE.
     ELSE
        diag % ok = .FALSE. 
        WRITE(rstr, '(f10.1)') 100. * REAL(n_bad)/REAL(n_good)  
        rstr = ADJUSTL(rstr)
        CALL message(msg_error, &
             "Background quality control has rejected " // TRIM(rstr)  //  &
             "% of the data points of this profile.\n" //                  &
             "   The entire profile is flagged as ""not ok"".\n")
     ENDIF
  ENDIF

! 2.6 Clean up
! ------------

  DEALLOCATE(idx)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_bgqc_1drefrac
