! $Id: ropp_qc_ion.f90 2280 2009-10-21 16:32:29Z frhl $

!****s* QC/ropp_qc_ion_bangle *
!
! NAME
!    ropp_qc_ion_bangle - Check the (noisy) L1 and L2 bending angle profiles
!                         if we're modelling them directly using -direct_ion.
!
! SYNOPSIS
!    CALL ropp_qc_ion_bangle(obs, config, diag)
!
! DESCRIPTION
!    This subroutine rejects elements of the L1 and L2 profiles that are 
!    too wiggly, as measured by the second derivative.
!
! INPUTS
!    obs         Obs1dBangle structure, containing {L1 and L2}
!    config      Configuration options (not currently used)
!
! OUTPUT
!    obs         Updated Obs1dBangle structure
!    diag        Diagnostics structure (not currently used)
!
! NOTES
!    Only applied if direct modelling of the ionopshere is being made,
!    ie, if -direct_ion is active.
!    Config and diag are included for possible future use.
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_qc_bgqc
!    ropp_qc_ion_bg
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

SUBROUTINE ropp_qc_ion_bangle(obs, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_1dvar_types
! USE ropp_qc, not_this => ropp_qc_ion_bangle

  IMPLICIT NONE

  TYPE(Obs1dBangle), INTENT(inout)    :: obs
  TYPE(VarConfig), INTENT(in)         :: config
  TYPE(VarDiag), INTENT(inout)        :: diag

  INTEGER                             :: j1, j2

  INTEGER                             :: n_steep1, n_steep2
  REAL(wp)                            :: eps, slope1, slope2
  REAL(wp), PARAMETER                 :: slope_thresh=1.0e-3_wp

  CHARACTER(len=256)                  :: routine
  CHARACTER(len=5)                    :: str_const

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_ion_bangle')

! 1.3 First derivative check
! --------------------------

  IF (obs%n_L1 < 3) THEN

    CALL message(msg_warn, &
      "Not enough points to apply first derivative check ... bypassing.")

  ELSE
    
    eps = TINY(1.0_wp)
    
    n_steep1 = 0  ;  n_steep2 = 0

! 1.3.1 bangle_L1 checks

    DO j1=2,obs%n_L1-1

      IF ((ALL(obs%bangle(j1-1:j1+1) > ropp_MDTV)) .AND. &
          (ALL(obs%impact(j1-1:j1+1) > ropp_MDTV))) THEN
        slope1 = 2.0_wp*(obs%bangle(j1-1) - obs%bangle(j1)) / &
                        (obs%bangle(j1-1) + obs%bangle(j1) + eps) / &
                        (obs%impact(j1-1) - obs%impact(j1) + eps)
        slope2 = 2.0_wp*(obs%bangle(j1+1) - obs%bangle(j1)) / &
                        (obs%bangle(j1+1) + obs%bangle(j1) + eps) / &
                        (obs%impact(j1+1) - obs%impact(j1) + eps)
        IF ((ABS(slope1) > slope_thresh) .AND. &
            (ABS(slope2) > slope_thresh)) THEN
          n_steep1 = n_steep1 + 1
          obs%weights(j1) = 0.0_wp
        END IF
      END IF

    END DO

    j1 = 1
    slope2 = 2.0_wp*(obs%bangle(j1+1) - obs%bangle(j1)) / &
                    (obs%bangle(j1+1) + obs%bangle(j1) + eps) / &
                    (obs%impact(j1+1) - obs%impact(j1) + eps)
    IF ((ABS(slope2) > slope_thresh)) THEN
      n_steep1 = n_steep1 + 1
      obs%weights(j1) = 0.0_wp
    END IF

    j1 = obs%n_L1
    slope1 = 2.0_wp*(obs%bangle(j1-1) - obs%bangle(j1)) / &
                    (obs%bangle(j1-1) + obs%bangle(j1) + eps) / &
                    (obs%impact(j1-1) - obs%impact(j1) + eps)
    IF ((ABS(slope1) > slope_thresh)) THEN
      n_steep1 = n_steep1 + 1
      obs%weights(j1) = 0.0_wp
    END IF

    IF (n_steep1 > 0) THEN
      WRITE (str_const, '(I5)') n_steep1
      CALL message(msg_info, &
        "First derivative check removes " // str_const // &
        " points from L1 bending angle profile.\n")
    END IF

! 1.3.2 bangle_L2 checks

    DO j2=obs%n_L1+2,2*obs%n_L1-1

      IF ((ALL(obs%bangle(j2-1:j2+1) > ropp_MDTV)) .AND. &
          (ALL(obs%impact(j2-1:j2+1) > ropp_MDTV))) THEN
        slope1 = 2.0_wp*(obs%bangle(j2-1) - obs%bangle(j2)) / &
                        (obs%bangle(j2-1) + obs%bangle(j2) + eps) / &
                        (obs%impact(j2-1) - obs%impact(j2) + eps)
        slope2 = 2.0_wp*(obs%bangle(j2+1) - obs%bangle(j2)) / &
                        (obs%bangle(j2+1) + obs%bangle(j2) + eps) / &
                        (obs%impact(j2+1) - obs%impact(j2) + eps)
        IF ((ABS(slope1) > slope_thresh) .AND. &
            (ABS(slope2) > slope_thresh)) THEN
          n_steep2 = n_steep2 + 1
          obs%weights(j2) = 0.0_wp
        END IF
      END IF

    END DO

    j2 = obs%n_L1 + 1
    slope2 = 2.0_wp*(obs%bangle(j2+1) - obs%bangle(j2)) / &
                    (obs%bangle(j2+1) + obs%bangle(j2) + eps) / &
                    (obs%impact(j2+1) - obs%impact(j2) + eps)
    IF ((ABS(slope2) > slope_thresh)) THEN
      n_steep2 = n_steep2 + 1
      obs%weights(j2) = 0.0_wp
    END IF

    j2 = 2*obs%n_L1
    slope1 = 2.0_wp*(obs%bangle(j2-1) - obs%bangle(j2)) / &
                    (obs%bangle(j2-1) + obs%bangle(j2) + eps) / &
                    (obs%impact(j2-1) - obs%impact(j2) + eps)
    IF ((ABS(slope1) > slope_thresh)) THEN
      n_steep2 = n_steep2 + 1
      obs%weights(j2) = 0.0_wp
    END IF

    IF (n_steep2 > 0) THEN
      WRITE (str_const, '(I5)') n_steep2
      CALL message(msg_info, &
        "First derivative check removes " // str_const // &
        " points from L2 bending angle profile.\n")
    END IF

  END IF

! 1.4 Return to sender
! --------------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_ion_bangle
