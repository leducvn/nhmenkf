! $Id: ropp_qc_pge.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* QC/ropp_qc_pge *
!
! NAME
!    ropp_qc_pge - Probability of gross error.
!
! SYNOPSIS
!    call rop_qc_pge(obs, config, diag)
! 
! DESCRIPTION
!    This subroutine calculates the Probability of Gross Error for bending
!    angle or refractivity observation data, and optionally uses the result
!    for quality control.
!
! INPUTS
!    obs          Observation vector
!    config       Configuration options
!    diag         Diagnostic structure
!
! OUTPUT
!    diag         Diagnostic structure with updated PGE variable
!
! NOTES
!    ropp_qc_pge requires that some fields in the diag structure have
!    been properly filled, e.g. by a previous call to ropp_qc_OmB().
!
! SEE ALSO
!    ropp_qc_OmB
!    ropp_qc_bgqc
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

SUBROUTINE ropp_qc_pge_1dbangle(obs, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV
  USE ropp_fm_types
  USE ropp_fm_constants
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DBangle)              :: obs
  TYPE(VarConfig)                :: config
  TYPE(VarDiag)                  :: diag

  CHARACTER(len = 256)           :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_pge')

! 1.3 Calculate gamma (if required)
! ---------------------------------

  IF (config % pge_d > ropp_MDTV) THEN
     diag % pge_gamma = SQRT(2.0_wp*pi) * config % pge_fg /      &
                        (2.0_wp * (1.0_wp - config % pge_fg) * config % pge_d)
  ELSE
     diag % pge_gamma = config % pge_fg
  ENDIF

! 1.4 Calculate PGE
! -----------------

  ALLOCATE(diag % pge(SIZE(obs % bangle)))
  ALLOCATE(diag % pge_weights(SIZE(obs % bangle)))

  diag % pge = diag % pge_gamma /     &
               (diag % pge_gamma + EXP(- 0.5_wp * (diag%OmB/diag%OmB_sigma)**2))
  diag % pge_weights = 1.0_wp - diag % pge

! 1.5 Set PGE weights
! -------------------

  IF (config % pge_apply) THEN
     WHERE (obs % weights > 0)
        obs % weights = diag % pge_weights
     END WHERE
     CALL message(msg_info, &
          "Observational weights were modified with weights " //   &
          "based on the Probability of Gross Error.\n")
  ENDIF

! 1.6 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_pge_1dbangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_pge_1drefrac(obs, config, diag)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV
  USE ropp_fm_types
  USE ropp_fm_constants
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DRefrac)              :: obs
  TYPE(VarConfig)                :: config
  TYPE(VarDiag)                  :: diag

  CHARACTER(len = 256)           :: routine

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_pge')

! 2.3 Calculate gamma (if required)
! ---------------------------------

  IF (config % pge_d > ropp_MDTV) THEN
     diag % pge_gamma = SQRT(2.0_wp * pi) * config % pge_fg /    &
                        (2.0_wp * (1.0_wp - config % pge_fg) * config % pge_d)
  ELSE
     diag % pge_gamma = config % pge_fg
  ENDIF

! 2.4 Calculate PGE
! -----------------

  ALLOCATE(diag % pge(SIZE(obs % refrac)))
  ALLOCATE(diag % pge_weights(SIZE(obs % refrac)))

  diag % pge = diag % pge_gamma /    &
               (diag % pge_gamma + EXP(- 0.5_wp * (diag%OmB/diag%OmB_sigma)**2))
  diag % pge_weights = 1.0_wp - diag % pge

! 2.5 Set PGE weights
! -------------------

  IF (config % pge_apply) THEN
     WHERE (obs % weights > 0)
        obs % weights = diag % pge_weights
     END WHERE
     CALL message(msg_info, &
          "Observational weights were modified with weights " //  &
          "based on the Probability of Gross Error.\n")
  ENDIF

! 2.6 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_pge_1drefrac
