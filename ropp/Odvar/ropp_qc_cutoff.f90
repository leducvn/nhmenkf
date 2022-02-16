! $Id: ropp_qc_cutoff.f90 2280 2009-10-21 16:32:29Z frhl $

!****s* QC/ropp_qc_cutoff *
!
! NAME
!    ropp_qc_cutoff - Down-weight data outside required observation
!                     height range [min_1dvar_height to max_1dvar_height]
!
! SYNOPSIS
!    call ropp_qc_cutoff(obs, config)
!
! DESCRIPTION
!    This subroutine down-weights all observations outside the height range
!    interval specified by configuration parameters min_1dvar_height and
!    max_1dvar_height
!
! INPUTS
!    obs         Observation vector
!    config      Configuration options
!
! OUTPUT
!    obs         Updated observation vector
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!
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

SUBROUTINE ropp_qc_cutoff_bangle(obs, config)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_cutoff_bangle
  IMPLICIT NONE

  TYPE(Obs1DBangle)                     :: obs
  TYPE(VarConfig)                       :: config
  CHARACTER(len =  256)                 :: routine
  CHARACTER(len=5)                      :: istr, rstr
  INTEGER                               :: ncut

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_cutoff')

! 1.3 Observation data height range check
! ---------------------------------------

! 1.3.1 Valid bending angle values

  ncut = COUNT((obs%impact - obs%r_curve) < config%min_1dvar_height*1000.0_wp)
  IF (ncut > 0) THEN
    WRITE (istr, '(I5)') ncut
    WRITE (rstr, '(F5.1)') config%min_1dvar_height
    CALL message(msg_info, TRIM(ADJUSTL(istr)) // " points rejected for " // &
      "being below minimum impact height of " // &
      TRIM(ADJUSTL(rstr)) // " km.\n")
    WHERE((obs%impact - obs%r_curve) < config%min_1dvar_height*1000.0_wp)
      obs%weights = 0.0_wp
    END WHERE
  END IF

  ncut = COUNT((obs%impact - obs%r_curve) > config%max_1dvar_height*1000.0_wp)
  IF (ncut > 0) THEN
    WRITE (istr, '(I5)') ncut
    WRITE (rstr, '(F5.1)') config%max_1dvar_height
    CALL message(msg_info, TRIM(ADJUSTL(istr)) // " points rejected for " // &
      "being above maximum impact height of " // &
      TRIM(ADJUSTL(rstr)) // " km.\n")
    WHERE((obs%impact - obs%r_curve) > config%max_1dvar_height*1000.0_wp)
      obs%weights = 0.0_wp
    END WHERE
  END IF

! 1.3.2 Reset routine name

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_cutoff_bangle

!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_cutoff_refrac(obs, config)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_cutoff_refrac

  IMPLICIT NONE

  TYPE(Obs1DRefrac)                     :: obs
  TYPE(VarConfig)                       :: config
  CHARACTER(len =  256)                 :: routine
  CHARACTER(len=5)                      :: istr, rstr
  INTEGER                               :: ncut

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_cutoff')

! 2.3 Observation data height range check
! ---------------------------------------

! 2.3.1 Valid refractivity values

  ncut = COUNT(obs%geop < config%min_1dvar_height*1000.0_wp)
  IF (ncut > 0) THEN
    WRITE (istr, '(I5)') ncut
    WRITE (rstr, '(F5.1)') config%min_1dvar_height
    CALL message(msg_info, TRIM(ADJUSTL(istr)) // " points rejected for " // &
      "being below minimum geopotential height of " // &
      TRIM(ADJUSTL(rstr)) // " km.\n")
    WHERE(obs%geop < config%min_1dvar_height*1000.0_wp)
      obs%weights = 0.0_wp
    END WHERE
  END IF

  ncut = COUNT(obs%geop > config%max_1dvar_height*1000.0_wp)
  IF (ncut > 0) THEN
    WRITE (istr, '(I5)') ncut
    WRITE (rstr, '(F5.1)') config%max_1dvar_height
    CALL message(msg_info, TRIM(ADJUSTL(istr)) // " points rejected for " // &
      "being above maximum geopotential height of " // &
      TRIM(ADJUSTL(rstr)) // " km.\n")
    WHERE(obs%geop > config%max_1dvar_height*1000.0_wp)
      obs%weights = 0.0_wp
    END WHERE
  END IF

! 2.3.2 Reset routine name

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_cutoff_refrac
