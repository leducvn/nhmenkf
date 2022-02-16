! $Id: ropp_qc_omb.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* QC/ropp_qc_OmB *
!
! NAME
!    ropp_qc_OmB - Calculate O minus B.
!
! SYNOPSIS
!    call ropp_qc_OmB(obs, bg, diag)
!
! DESCRIPTION
!    This subroutine calculates the O minus B (i.e., the difference between
!    the actual observations and the background forward modelled into
!    observation space) along with the expected uncertainty of this difference.
!    The results are placed in the diag structure.
!
! INPUTS
!    obs      Observation vector
!    bg       Background vector
!
! OUTPUT
!    diag     Diagnostic structure updated with O - B field
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

SUBROUTINE ropp_qc_OmB_bangle(obs, bg, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE messages
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_OmB_bangle

  IMPLICIT NONE

  TYPE(Obs1DBangle)                     :: obs
  TYPE(State1DFM)                       :: bg
  TYPE(VarConfig)                       :: config
  TYPE(VarDiag)                         :: diag

  INTEGER                               :: i, m, n

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: B
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: R
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: OmB_covar

  CHARACTER(len =  256)                 :: routine

  LOGICAL :: dummy
  dummy = config%compress !! fix nag 'unused dummy variable'

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_OmB')

! 1.3 Allocate memory
! -------------------

  n = SIZE(bg % state)
  m = SIZE(obs % bangle)

  IF (ASSOCIATED(diag % OmB))       DEALLOCATE(diag % OmB)
  IF (ASSOCIATED(diag % OmB_sigma)) DEALLOCATE(diag % OmB_sigma)
  IF (ASSOCIATED(diag % B_sigma))   DEALLOCATE(diag % B_sigma)

  ALLOCATE(diag % OmB(m), diag % OmB_sigma(m), diag % B_sigma(m))

  ALLOCATE(B(n, n), K(m, n), R(m, m), OmB_covar(m, m))

! 1.4 Forward model the background and its errors
! -----------------------------------------------

  diag % bg_bangle = obs
  CALL ropp_fm_bangle_1d(bg, diag % bg_bangle)
  CALL ropp_fm_bangle_1d_grad(bg, diag % bg_bangle, K)

! 1.5 Calculate forward modelled background error
! -----------------------------------------------

  B = bg%cov

  CALL matrix_toast(B, K, R)

  diag % bg_bangle % cov = R

! 1.6 Calculate O - B
! -------------------

  WHERE (obs % bangle > ropp_MDTV)
     diag % OmB = obs % bangle - diag % bg_bangle % bangle
  ELSEWHERE
     diag % OmB = ropp_MDFV
  END WHERE

! 1.7 Calculate expected O - B sigma
! ----------------------------------

  OmB_covar = obs % cov + diag % bg_bangle % cov

  DO i = 1, m
    diag % B_sigma(i) = SQRT(R(i,i))
    diag % OmB_sigma(i) = SQRT(OmB_covar(i,i))
  END DO

! 1.8 Clean up
! ------------

  DEALLOCATE(K, B, R, OmB_covar)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_OmB_bangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_OmB_refrac(obs, bg, config, diag)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE messages
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_OmB_refrac

  IMPLICIT NONE

  TYPE(Obs1DRefrac)                     :: obs
  TYPE(State1DFM)                       :: bg
  TYPE(VarConfig)                       :: config
  TYPE(VarDiag)                         :: diag

  INTEGER                               :: i, m, n

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: B
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: R
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: OmB_covar

  CHARACTER(len =  256)                 :: routine

  LOGICAL :: dummy
  dummy = config%compress !! fix nag 'unused dummy variable'

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_OmB')

! 2.3 Allocate memory
! -------------------

  n = SIZE(bg % state)
  m = SIZE(obs % refrac)

  IF (ASSOCIATED(diag % OmB))       DEALLOCATE(diag % OmB)
  IF (ASSOCIATED(diag % OmB_sigma)) DEALLOCATE(diag % OmB_sigma)
  IF (ASSOCIATED(diag % B_sigma)) DEALLOCATE(diag % B_sigma)

  ALLOCATE(diag % OmB(m), diag % OmB_sigma(m), diag % B_sigma(m))

  ALLOCATE(B(n, n), K(m, n), R(m, m), OmB_covar(m, m))

! 2.4 Forward model the background and its errors
! -----------------------------------------------

  diag % bg_refrac = obs
  IF (bg%new_ref_op) THEN
    CALL ropp_fm_refrac_1d_new(bg, diag % bg_refrac)
  ELSE
    CALL ropp_fm_refrac_1d(bg, diag % bg_refrac)
  END IF
  CALL ropp_fm_refrac_1d_grad(bg, diag % bg_refrac, K)

! 2.5 Calculate forward modelled background error
! -----------------------------------------------

  B = bg%cov

  CALL matrix_toast(B, K, R)

  diag % bg_refrac % cov = R

! 2.6 Calculate O - B
! -------------------

  WHERE (obs % refrac > ropp_MDTV)
     diag % OmB = obs % refrac - diag % bg_refrac % refrac
  ELSEWHERE
     diag % OmB = ropp_MDFV
  END WHERE

! 2.7 Calculate expected O - B sigma
! ---------------------------------

  OmB_covar = obs % cov + diag % bg_refrac % cov

  DO i = 1, m
    diag % B_sigma(i) = SQRT(R(i,i))
    diag % OmB_sigma(i) = SQRT(OmB_covar(i,i))
  END DO

! 2.8 Clean up
! ------------

  DEALLOCATE(K, B, R, OmB_covar)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_OmB_refrac
