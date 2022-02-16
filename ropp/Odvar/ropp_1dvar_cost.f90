! $Id: ropp_1dvar_cost.f90 4452 2015-01-29 14:42:02Z idculv $
                                    
!****f* 1DVar/ropp_1dvar_cost
!
! NAME
!    ropp_1dvar_cost - Evaluate a cost function for a one-dimensional 
!                      variational retrieval of radio occultation data.
!
! SYNOPSIS
!    call ropp_1dvar_cost(yo, bg, control, precon, J, control_ad, config, indic)
!         
! 
! DESCRIPTION
!    This subroutine evaluates a quadratic cost function for a
!    variational data assimilation procedure. 
!
!    More specifically, this routine calculates a cost function
!
!           1  /         |  -1 |         \
!       J = - < y - H(x) | O   | y - H(x) > + 
!           2  \         |     |         /
!
!                       1  /       |  -1 |       \
!                       - < x - x  | B   | x - x  >
!                       2  \     b |     |      b/
!
!    where the background state x_b is given by the state vector
!    state.
!
!    On exit, control_ad contains the gradient of J with respect to the 
!    variables in control.
!
! INPUTS
!    yo                 Observation data structure.
!    bg                 Background data structure.
!    control            Control variable structure.
!    precon             Preconditioning matrix.
!    J                  Cost function value.
!    config             Configuration structure.
!    indic              State variable used for the communication with the
!                         calling program and the minimiser. The following
!                         values are recognised as input:
!                           1       Initialise arrays if additional 
!                                     convergence checks are to be applied.
!                           other:  Ignored
!                         The following values are output if additional
!                         convergence checks are to be applied:
!                           0       Iteration converged
!                           4       Iteration not converged
!
! OUTPUT
!    J 
!    control_ad
!
! NOTES
!    If called with indic < 0 (which should never happen during the 
!    minimization), variables for the convergence checks are reset. The 
!    calculations for cost function and its gradient are nevertheless 
!    performed. This allows the simultanious initialization
!    of the cost function minimization and the convergence checks with a single
!    call to ropp_1dvar_cost.
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
! 1. Bending angle cost function
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_cost_bangle(yo, bg, control, precon, J, control_ad,    &
                                  config, indic)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE messages
  USE arrays
  USE ropp_fm
! USE ropp_1dvar, not_this => ropp_1dvar_cost_bangle
  USE ropp_1dvar
  USE matrix
  USE ropp_fm_types
  
  IMPLICIT NONE

  TYPE(Obs1dBangle),      INTENT(inout)     :: yo
  TYPE(State1dFM),        INTENT(inout)     :: bg
  TYPE(State1dFM),        INTENT(in)        :: control
  TYPE(matrix_sq),        INTENT(in)        :: precon
  REAL(wp),               INTENT(out)       :: J
  REAL(wp), DIMENSION(:), INTENT(out)       :: control_ad
  TYPE(VarConfig),        INTENT(in)        :: config
  INTEGER,                INTENT(inout)     :: indic

  TYPE(State1dFM)                           :: x
  TYPE(State1dFM)                           :: x_ad
  TYPE(Obs1dBangle)                         :: y
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: delta_x
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: delta_y
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: y_ad
  REAL(wp)                                  :: J_ad
  REAL(wp),                            SAVE :: J_last
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: delta_J
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: state_last
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: state_sigma
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: delta_state
  INTEGER,                             SAVE :: i_pointer
  INTEGER,                             SAVE :: n_iter
  INTEGER                                   :: i, m

  CHARACTER(len =   4)                      :: it_str
  CHARACTER(len =  15)                      :: ch_str, co_str
  CHARACTER(len = 256)                      :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_cost')

! 1.3 Useful variables
! --------------------

  m = SIZE(yo % bangle)

  ALLOCATE(delta_x(SIZE(bg%state)))
  ALLOCATE(delta_y(SIZE(yo%bangle)))
  ALLOCATE(y_ad(SIZE(yo%bangle)))

! 1.4 Initialise rolling buffers and pointer for convergence checks
! -----------------------------------------------------------------

  IF (config % conv_check_apply .AND. indic == 1) THEN

     CALL message(msg_info, &
          'Using absolute decrease of cost function and relative change \n '// &
          'of state vector as additional convergence criteria.\n')

     i_pointer = 0
     n_iter    = 0
     J_last    = 0.0_wp

     IF (ALLOCATED(delta_J)) DEALLOCATE(delta_J)
     ALLOCATE(delta_J(config%conv_check_n_previous))
     delta_J(:) = 0.0_wp

     IF (ALLOCATED(state_last)) DEALLOCATE(state_last)
     ALLOCATE(state_last(SIZE(bg%state)))
     state_last(:) = 0.0_wp
     
     IF (ALLOCATED(state_sigma)) DEALLOCATE(state_sigma)
     ALLOCATE(state_sigma(SIZE(bg%state)))
     DO i = 1, SIZE(state_sigma)
        ! Direct readout from matrix_pp type
        state_sigma(i) = SQRT(bg%cov%d(i + i*(i-1)/2)) 
     ENDDO

     IF (ALLOCATED(delta_state)) DEALLOCATE(delta_state)
     ALLOCATE(delta_state(config%conv_check_n_previous))
     delta_state(:) = 0.0_wp

  ENDIF

  indic = 4

! 1.5 Reset adjoint variables
! ---------------------------

  ALLOCATE(x_ad%state(SIZE(control%state)))
  ALLOCATE(x_ad%pres(control%n_lev))
  ALLOCATE(x_ad%temp(control%n_lev))
  ALLOCATE(x_ad%shum(control%n_lev))
  ALLOCATE(x_ad%geop(control%n_lev))

  n_iter     = n_iter + 1

  control_ad = 0.0_wp
  x_ad%state = 0.0_wp
  x_ad%pres  = 0.0_wp
  x_ad%temp  = 0.0_wp
  x_ad%shum  = 0.0_wp
  x_ad%geop  = 0.0_wp
  x_ad%ne_max = 0.0_wp
  x_ad%h_peak = 0.0_wp
  x_ad%h_width = 0.0_wp
  y_ad       = 0.0_wp
  J_ad       = 1.0_wp

! 1.6 Forward model and cost function
! -----------------------------------

! 1.6.1 Preconditioning

  x = control

  IF (config % use_precond) THEN
     CALL ropp_control2state(precon, control % state, x % state)
  ENDIF
  
! 1.6.2 Calculate pseudo observations

  ALLOCATE(y % bangle(m))

  y % g_sfc        = yo % g_sfc
  y % r_earth      = yo % r_earth
  y % r_curve      = yo % r_curve
  y % undulation   = yo % undulation

  CALL copy_alloc(yo % impact,  y % impact)
  CALL copy_alloc(yo % weights, y % weights)
  
! 1.6.3 Forward model

  IF(ASSOCIATED(x%ak))THEN
     CALL ropp_fm_state2state_ecmwf(x)
  ELSE
     CALL ropp_fm_state2state_meto(x)
  ENDIF

  CALL ropp_fm_bangle_1d(x, y)

! 1.6.4 Compute cost function

  delta_x = (x%state  - bg%state)
  delta_y = (y%bangle - yo%bangle) * y%weights

  J = 0.5_wp * DOT_PRODUCT(delta_y, matrix_solve(yo%cov, delta_y)) &
    + 0.5_wp * DOT_PRODUCT(delta_x, matrix_solve(bg%cov, delta_x))

! 1.7 Adjoint code for gradient
! -----------------------------

  x_ad%state = x_ad%state + J_ad*matrix_solve(bg%cov, delta_x)
  y_ad = y_ad + J_ad*matrix_solve(yo%cov, delta_y)
  J_ad = 0.0_wp

  y_ad = y_ad * y%weights

  CALL ropp_fm_bangle_1d_ad(x, x_ad, y, y_ad)

  IF(ASSOCIATED(x%ak))THEN
     CALL ropp_fm_state2state_ecmwf_ad(x, x_ad)
  ELSE
     CALL ropp_fm_state2state_meto_ad(x, x_ad)
  ENDIF

  IF (config % use_precond) THEN
     CALL ropp_control2state_ad(precon, control_ad, x_ad%state)
  ELSE
     control_ad = control_ad + x_ad%state
     x_ad%state = 0.0_wp
  ENDIF

! 1.8 Keep track of current state
! -------------------------------

  IF (config % conv_check_apply .AND. indic > 0) THEN

     i_pointer = MOD(i_pointer + 1, config % conv_check_n_previous)
     IF (i_pointer == 0) i_pointer = config % conv_check_n_previous

     ! Check that cost function is decreasing
     IF(n_iter > config % conv_check_n_previous .AND.    &
        (J - J_last) > maxval(delta_J)) THEN  
       CALL message(msg_info,    &
           'Convergence assumed to be achieved as the cost function is ' //  &
           'increasing \n ')
       J = J_last
       x%state = state_last
       n_iter = n_iter - 1
       indic = 0
       RETURN
     ENDIF

     delta_J(i_pointer)     = J_last - J                                     
     J_last = J
     delta_state(i_pointer) = MAXVAL(ABS(state_last - x%state)/state_sigma)  
     state_last = x%state

     indic = 4

  ENDIF

! 1.9 Check for convergence
! -------------------------

  IF (config % conv_check_apply .AND. indic == 4) THEN

     IF (config%minropp%impres == 0) THEN
        WRITE(it_str, '(i4)') n_iter
        WRITE(ch_str, '(g15.5)') J                         
        ch_str = ADJUSTL(ch_str)
        IF (n_iter > 1) THEN
           WRITE(co_str, '(g15.5)') delta_state(i_pointer)  
           co_str = ADJUSTL(co_str)
        ELSE
           co_str = ' -             '
        ENDIF

        CALL message(msg_cont, &
             '   n_iter = ' // it_str // '   J = ' // ch_str //    &
             '   max(relative change in state) = ' // co_str)
     ENDIF
 
     IF (MAXVAL(delta_state) < config%conv_check_max_delta_state .AND.      &
         n_iter > config % conv_check_n_previous) THEN
        WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_state  
        ch_str = ADJUSTL(ch_str)
        WRITE(it_str, '(i2)')    config%conv_check_n_previous       
        it_str = ADJUSTL(it_str)
        CALL message(msg_cont, '')
        CALL message(msg_info,    &
           'Convergence assumed to be achieved as the state vector did ' // &
           'not change by more\n   ' // 'than ' // TRIM(ch_str) // ' ' //   &
           'relative to the assumed background errors for the last ' //     &
           TRIM(it_str) // ' iterations.\n')
        indic = 0
     ELSE IF (MAXVAL(ABS(delta_J)) < config%conv_check_max_delta_J .AND.    &
        n_iter > config % conv_check_n_previous) THEN
        WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_J  
        ch_str = ADJUSTL(ch_str)
        WRITE(it_str, '(i2)')    config%conv_check_n_previous   
        it_str = ADJUSTL(it_str)
        CALL message(msg_cont, '')
        CALL message(msg_info,    &
           'Convergence assumed to be achieved as the cost function did ' // &
           'not change by more\n   ' // 'than ' // TRIM(ch_str) //           &
           ' for the last ' // TRIM(it_str) // ' iterations.\n')
        indic = 0
     ELSE
        indic = 4
     ENDIF 

  ENDIF

! 1.10 Clean up
! -------------
  
  CALL ropp_fm_free(y)
  CALL ropp_fm_free(x_ad)

  IF (config % conv_check_apply .AND. indic == 0) THEN
     DEALLOCATE(delta_state)
     DEALLOCATE(state_sigma)
     DEALLOCATE(state_last)
     DEALLOCATE(delta_J)
  ENDIF

  DEALLOCATE(delta_x)
  DEALLOCATE(delta_y)
  DEALLOCATE(y_ad)
  
  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_cost_bangle


!-------------------------------------------------------------------------------
! 2. Refractivity cost function
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_cost_refrac(yo, bg, control, precon, J, control_ad,   &
                                  config, indic)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE messages
  USE arrays
  USE ropp_fm
! USE ropp_1dvar, not_this => ropp_1dvar_cost_refrac
  USE ropp_1dvar
  USE matrix
  USE ropp_fm_types
  
  IMPLICIT NONE

  TYPE(Obs1dRefrac),      INTENT(inout)     :: yo
  TYPE(State1dFM),        INTENT(inout)     :: bg
  TYPE(State1dFM),        INTENT(in)        :: control
  TYPE(matrix_sq),        INTENT(in)        :: precon
  REAL(wp),               INTENT(out)       :: J
  REAL(wp), DIMENSION(:), INTENT(out)       :: control_ad
  TYPE(VarConfig),        INTENT(in)        :: config
  INTEGER,                INTENT(inout)     :: indic

  TYPE(State1dFM)                           :: x
  TYPE(State1dFM)                           :: x_ad
  TYPE(Obs1dRefrac)                         :: y
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: delta_x
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: delta_y
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: y_ad
  REAL(wp)                                  :: J_ad
  REAL(wp),                            SAVE :: J_last
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: delta_J
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: state_last
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: state_sigma
  REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: delta_state
  INTEGER,                             SAVE :: i_pointer
  INTEGER,                             SAVE :: n_iter
  INTEGER                                   :: i, m

  CHARACTER(len =   4)                      :: it_str
  CHARACTER(len =  15)                      :: ch_str, co_str
  CHARACTER(len = 256)                      :: routine

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_cost')

! 2.3 Useful variables
! --------------------

  m = SIZE(yo % refrac)

  ALLOCATE(delta_x(SIZE(bg%state)))
  ALLOCATE(delta_y(SIZE(yo%refrac)))
  ALLOCATE(y_ad(SIZE(yo%refrac)))

! 2.4 Initialise rolling buffers and pointer for convergence checks
! -----------------------------------------------------------------

  IF (config % conv_check_apply .AND. indic == 1) THEN

     CALL message(msg_info, &
          'Using absolute decrease of cost function and relative change \n' // &
          'of state vector as additional convergence criterium.\n')

     i_pointer = 0
     n_iter    = 0
     J_last    = 0.0_wp

     IF (ALLOCATED(delta_J)) DEALLOCATE(delta_J)
     ALLOCATE(delta_J(config%conv_check_n_previous))
     delta_J(:) = 0.0_wp
     
     IF (ALLOCATED(state_last)) DEALLOCATE(state_last)
     ALLOCATE(state_last(SIZE(bg%state)))
     state_last(:) = 0.0_wp
     
     IF (ALLOCATED(state_sigma)) DEALLOCATE(state_sigma)
     ALLOCATE(state_sigma(SIZE(bg%state)))
     DO i = 1, SIZE(state_sigma)
        ! Direct readout from matrix_pp type
        state_sigma(i) = SQRT(bg%cov%d(i + i*(i-1)/2)) 
     ENDDO
     
     IF (ALLOCATED(delta_state)) DEALLOCATE(delta_state)
     ALLOCATE(delta_state(config%conv_check_n_previous))
     delta_state(:) = 0.0_wp

  ENDIF

  indic = 4

! 2.5 Reset adjoint variables
! ---------------------------
  
  ALLOCATE(x_ad%state(SIZE(control%state)))
  ALLOCATE(x_ad%pres(control%n_lev))
  ALLOCATE(x_ad%temp(control%n_lev))
  ALLOCATE(x_ad%shum(control%n_lev))
  ALLOCATE(x_ad%geop(control%n_lev))

  n_iter     = n_iter + 1
  control_ad = 0.0_wp
  x_ad%state = 0.0_wp
  x_ad%pres  = 0.0_wp
  x_ad%temp  = 0.0_wp
  x_ad%shum  = 0.0_wp
  x_ad%geop  = 0.0_wp
  y_ad       = 0.0_wp
  J_ad       = 1.0_wp

! 2.6 Forward model and cost function
! -----------------------------------

! 2.6.1 Preconditioning

  x = control

  IF (config % use_precond) THEN
     CALL ropp_control2state(precon, control % state, x % state)
  ENDIF

! 2.6.2 Calculate pseudo observations
 
  ALLOCATE(y % refrac(m))
  CALL copy_alloc(yo % geop,    y % geop)
  CALL copy_alloc(yo % weights, y % weights)
  
! 2.6.3 Forward model
  
  IF(ASSOCIATED(x%ak))THEN
     CALL ropp_fm_state2state_ecmwf(x)
  ELSE
     CALL ropp_fm_state2state_meto(x)
  ENDIF

  IF (x%new_ref_op) THEN
    CALL ropp_fm_refrac_1d_new(x, y)
  ELSE
    CALL ropp_fm_refrac_1d(x, y)
  END IF

! 2.6.4 Compute cost function

  delta_x = (x%state  - bg%state)
  delta_y = (y%refrac - yo%refrac) * y%weights

  J = 0.5_wp * DOT_PRODUCT(delta_y, matrix_solve(yo%cov, delta_y)) &
       + 0.5_wp * DOT_PRODUCT(delta_x, matrix_solve(bg%cov, delta_x))

! 2.7 Adjoint code for gradient
! -----------------------------
  
  x_ad%state = x_ad%state + J_ad*matrix_solve(bg%cov, delta_x)
  y_ad = y_ad + J_ad*matrix_solve(yo%cov, delta_y)
  J_ad = 0.0_wp
  
  y_ad = y_ad * y%weights

  IF (x%new_ref_op) THEN
    CALL ropp_fm_refrac_1d_new_ad(x, x_ad, y, y_ad)
  ELSE
    CALL ropp_fm_refrac_1d_ad(x, x_ad, y, y_ad)
  END IF

  IF(ASSOCIATED(x%ak))THEN
     CALL ropp_fm_state2state_ecmwf_ad(x, x_ad)
  ELSE
     CALL ropp_fm_state2state_meto_ad(x, x_ad)
  ENDIF

  IF (config % use_precond) THEN
     CALL ropp_control2state_ad(precon, control_ad, x_ad%state)
  ELSE
     control_ad = control_ad + x_ad%state
     x_ad%state = 0.0_wp
  ENDIF

! 2.8 Keep track of current state
! -------------------------------

  IF (config % conv_check_apply .AND. indic > 0) THEN

     i_pointer = MOD(i_pointer + 1, config % conv_check_n_previous)
     IF (i_pointer == 0) i_pointer = config % conv_check_n_previous
     
     ! Check that cost function is decreasing
     IF(n_iter > config % conv_check_n_previous .AND.    &
        (J - J_last) > maxval(delta_J)) THEN  
       CALL message(msg_info,    &
           'Convergence assumed to be achieved as the cost function is ' //  &
           'increasing \n ')
       J = J_last
       x%state = state_last
       n_iter = n_iter - 1
       indic = 0
       RETURN
     ENDIF

     delta_J(i_pointer)     = J_last - J 
     J_last = J
     delta_state(i_pointer) = MAXVAL(ABS(state_last - x%state)/state_sigma) 
     state_last = x%state

     indic = 4

  ENDIF

! 2.9 Check for convergence
! -------------------------

  IF (config % conv_check_apply .AND. indic == 4) THEN

     IF (config%minropp%impres == 0) THEN
        WRITE(it_str, '(i4)') n_iter
        WRITE(ch_str, '(g15.5)') J                         
        ch_str = ADJUSTL(ch_str)
        IF (n_iter > 1) THEN
           WRITE(co_str, '(g15.5)') delta_state(i_pointer)  
           co_str = ADJUSTL(co_str)
        ELSE
           co_str = ' -             '
        ENDIF

        CALL message(msg_cont, &
             '   n_iter = ' // it_str // '   J = ' // ch_str //   &
             '   max(relative change in state) = ' // co_str)

     ENDIF

     IF (MAXVAL(delta_state) < config%conv_check_max_delta_state .AND.      &
        n_iter > config % conv_check_n_previous) THEN
        WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_state  
        ch_str = ADJUSTL(ch_str)
        WRITE(it_str, '(i2)')    config%conv_check_n_previous      
        it_str = ADJUSTL(it_str)
        CALL message(msg_cont, '')
        CALL message(msg_info,    &
           'Convergence assumed to be achieved as the state vector did ' // &
           'not change by more\n   ' // 'than ' // TRIM(ch_str) // ' ' //   &
           'relative to the assumed background errors for the last ' //     &
           TRIM(it_str) // ' iterations.\n')
        indic = 0
     ELSE IF (MAXVAL(ABS(delta_J)) < config%conv_check_max_delta_J .AND.    &
        n_iter > config % conv_check_n_previous) THEN
        WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_J 
        ch_str = ADJUSTL(ch_str)
        WRITE(it_str, '(i2)')    config%conv_check_n_previous  
        it_str = ADJUSTL(it_str)
        CALL message(msg_cont, '')
        CALL message(msg_info,    &
           'Convergence assumed to be achieved as the cost function did ' // &
           'not change by more\n   ' // 'than ' // TRIM(ch_str) //           &
           ' for the last ' // TRIM(it_str) // ' iterations.\n')
        indic = 0
     ELSE
        indic = 4
     ENDIF

  ENDIF

! 2.10 Clean up
! -------------

  CALL ropp_fm_free(y)
  CALL ropp_fm_free(x_ad)

  IF (config % conv_check_apply .AND. indic == 0) THEN
     DEALLOCATE(delta_state)
     DEALLOCATE(state_sigma)
     DEALLOCATE(state_last)
     DEALLOCATE(delta_J)
  ENDIF

  DEALLOCATE(delta_x)
  DEALLOCATE(delta_y)
  DEALLOCATE(y_ad)

  CALL message_set_routine(routine)
  
END SUBROUTINE ropp_1dvar_cost_refrac
