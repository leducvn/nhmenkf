! $Id: ropp_1dvar_diagnostics.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Diagnostics/ropp_1dvar_diagnostics *
!
! NAME
!    ropp_1dvar_diagnostics - 1DVar postprocessing and diagnostics.
!
! SYNOPSIS
!    call rop_1dvar_diagnostics(obs, state, config, diag)
!
! DESCRIPTION
!    This subroutine performs some postprocessing and diagnostic
!    calculations for a 1DVar. In particular,
!
!       - O - A differences for the solution
!       - error estimates for the solution
!
!    are calculated and stored into the respective members of the
!    diagnostic structure 'diag' as well as into the state vector's
!    error covariance matrix.
!
! INPUTS
!    obs
!    bg
!    config
!
! OUTPUT
!    bg
!    diag
!
! NOTES
!    The calculation of the expected errors in O - A is currently open.
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

SUBROUTINE ropp_1dvar_diag_1dbangle(obs, state, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV, ropp_MDFV
  USE matrix
  USE ropp_fm
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DBangle), INTENT(inout)      :: obs
  TYPE(State1DFM),   INTENT(inout)      :: state
  TYPE(VarConfig),   INTENT(in)         :: config
  TYPE(VarDiag),     INTENT(inout)      :: diag

  INTEGER                               :: m, n, i
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Rm1, Bm1, K, Kt, S, R, B
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KtRm1K
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: OmA_covar

  CHARACTER(len =  256) :: outstr
  CHARACTER(len =  256) :: routine

! LOGICAL :: dummy ! Commented at 21 July, 2016
! dummy = config%compress !! fix nag 'unused dummy variable' ! Commented at 21 July, 2016

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_diagnostics')

! 1.2 Calculate forward modelled solution, its gradient and O - A
! ---------------------------------------------------------------

  n = SIZE(state % state)
  m = SIZE(obs % bangle)

  ALLOCATE(diag % OmA(m), diag % OmA_sigma(m))
  ALLOCATE(Rm1(m, m), Bm1(n, n), K(m, n), Kt(n, m), KtRm1K(n, n), S(n,n))
  ALLOCATE(R(m,m), B(n,n), OmA_covar(m, m))

  diag % res_bangle = obs

  CALL ropp_fm_bangle_1d(state, diag % res_bangle)

  CALL ropp_fm_bangle_1d_grad(state, diag % res_bangle, K)

  WHERE (obs % bangle > ropp_MDTV)
     diag % OmA = (obs % bangle - diag % res_bangle % bangle)
  ELSEWHERE
     diag % OmA = ropp_MDFV
  END WHERE

! 1.3 Calculate error covariance of the solution
! ----------------------------------------------

  Bm1 = matrix_invert(state % cov)
  Rm1 = matrix_invert(obs % cov)
  Kt  = TRANSPOSE(K)
  CALL matrix_toast(Rm1, Kt, KtRm1K)

  S = Bm1 + KtRm1K
  state % cov = matrix_invert(S)

! 1.4 Calculate expected O-A sigma
! --------------------------------

  B = state%cov
  CALL matrix_toast(B, K, R)
  diag%res_bangle%cov = R
  OmA_covar = obs%cov + diag%res_bangle%cov
  DO i=1,m
    diag%Oma_sigma(i) = SQRT(OmA_covar(i,i))
  ENDDO

! 1.5 Calculate observation contribution to cost function
! -------------------------------------------------------

  ALLOCATE(diag%J_obs(m))
  diag%J_obs = 0.5_wp * (diag%OmA*obs%weights) *    &
                 matrix_solve(obs % cov, (diag%OmA*obs%weights))

! 1.6 Quality control
! -------------------

  IF ( diag%n_iter >= 25 ) THEN
    diag%ok = .FALSE.
    WRITE(outstr,'(A,I2,A,F6.3)') "Warning: Too many iterations, niter = ", &
       diag%n_iter, " 2J/m = ", diag%J_scaled
    CALL message(msg_diag, TRIM(outstr))
  ENDIF

  IF ( diag%J_scaled >= 5 ) THEN
    diag%ok = .FALSE.
    WRITE(outstr,'(A,I2,A,F6.3)') "Warning: 2J/m greater than 5 , niter = ", &
       diag%n_iter, " 2J/m = ", diag%J_scaled
    CALL message(msg_diag, TRIM(outstr))
  ENDIF

! 1.7 Clean up
! ------------

  DEALLOCATE(Rm1, Bm1, K, Kt, KtRm1K, S, R, B, OmA_covar)
  CALL message_set_routine(routine)


END SUBROUTINE ropp_1dvar_diag_1dbangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_diag_1drefrac(obs, state, config, diag)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE arrays
  USE messages
  USE ropp_utils, ONLY: ropp_MDTV, ropp_MDFV
  USE matrix
  USE ropp_fm
  USE ropp_1dvar_types

  IMPLICIT NONE

  TYPE(Obs1DRefrac), INTENT(inout)      :: obs
  TYPE(State1DFM),   INTENT(inout)      :: state
  TYPE(VarConfig),   INTENT(in)         :: config
  TYPE(VarDiag),     INTENT(inout)      :: diag

  INTEGER                               :: m, n, i
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Rm1, Bm1, K, Kt, S, R, B
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KtRm1K
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: OmA_covar
  CHARACTER(len =  256) :: outstr
  CHARACTER(len =  256) :: routine

! LOGICAL :: dummy ! Commented at 21 July, 2016
! dummy = config%compress !! fix nag 'unused dummy variable' ! Commented at 21 July, 2016

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_diagnostics')

! 2.2 Calculate forward modelled solution, it's gradient and O - A
! ----------------------------------------------------------------

  n = SIZE(state % state)
  m = SIZE(obs % refrac)

  ALLOCATE(diag % OmA(m), diag % OmA_sigma(m))
  ALLOCATE(Rm1(m, m), Bm1(n, n), K(m, n), Kt(n, m), KtRm1K(n, n), S(n,n))
  ALLOCATE(R(m,m), B(n,n), OmA_covar(m, m))

  diag % res_refrac = obs

  IF (state%new_ref_op) THEN
    CALL ropp_fm_refrac_1d_new(state, diag % res_refrac)
  ELSE
    CALL ropp_fm_refrac_1d(state, diag % res_refrac)
  END IF
  CALL ropp_fm_refrac_1d_grad(state, diag % res_refrac, K)

  WHERE (obs % refrac > ropp_MDTV)
     diag % OmA = (obs % refrac - diag % res_refrac % refrac)
  ELSEWHERE
     diag % OmA = ropp_MDFV
  END WHERE

! 2.2 Calculate error covariance of the solution
! ----------------------------------------------

  Bm1 = matrix_invert(state % cov)
  Rm1 = matrix_invert(obs % cov)
  Kt  = TRANSPOSE(K)
  CALL matrix_toast(Rm1, Kt, KtRm1K)

  S = Bm1 + KtRm1K
  state % cov = matrix_invert(S)

! 2.3 Calculate expected O-A sigma
! --------------------------------

  B = state % cov
  CALL matrix_toast(B, K, R)
  diag%res_refrac%cov = R

  OmA_covar = obs%cov + diag%res_refrac%cov
  DO i=1,m
    diag%OmA_sigma(i) = SQRT(OmA_covar(i,i))
  ENDDO

! 2.4 Calculate observation contribution to cost function
! -------------------------------------------------------

  ALLOCATE(diag%J_obs(m))
  diag%J_obs = 0.5_wp * (diag%OmA*obs%weights) *      &
                  matrix_solve(obs % cov, (diag%OmA*obs%weights))

! 2.5 Quality control
! -------------------

  IF ( diag%n_iter >= 25 ) THEN
    diag%ok = .FALSE.
    WRITE(outstr,'(A,I2,A,F6.3)') "Warning: Too many iterations, niter = ", &
       diag%n_iter, " 2J/m = ", diag%J_scaled
    CALL message(msg_diag, TRIM(outstr))
  ENDIF

  IF ( diag%J_scaled >= 5 ) THEN
    diag%ok = .FALSE.
    WRITE(outstr,'(A,I2,A,F6.3)') "Warning: 2J/m greater than 5 , niter = ", &
       diag%n_iter, " 2J/m = ", diag%J_scaled
    CALL message(msg_diag, TRIM(outstr))
  ENDIF

! 2.6 Clean up
! ------------

  DEALLOCATE(Rm1, Bm1, K, Kt, KtRm1K, S, OmA_covar)
  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_diag_1drefrac
