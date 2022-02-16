! $Id: ropp_1dvar_levmarq.f90 4452 2015-01-29 14:42:02Z idculv $
                                    
!****s* 1DVar/ropp_1dvar_levmarq
!
! NAME
!    ropp_1dvar_levmarq - Solve the 1DVar for background data using the
!                         Levenberg-Marquardt minimiser
!
! SYNOPSIS
!    CALL ropp_1dvar_levmarq(obs, bg, state, config, diag)
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
!    A solution for x is obtained by minimising J using the Levenberg-Marquardt
!    minimisation method.
!
! INPUTS
!    obs                Observation data structure.
!    bg                 Background data structure.
!    state              State vector structure.
!    config             Configuration structure.
!    diag               Diagnostics structure.
!
! OUTPUT
!    state 
!    diag
!
! REFERENCES
!   W.H. Press, S.A. Teukolsjy, W.T. Vetterling and B.P. Flannery,
!   Numerical Recipes in C - The Art of Scientific Computing.
!   2nd Ed., Cambridge University Press, 1992.
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
! 1. Bending angle
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_levmarq_bangle(obs, bg, state, config, diag)
  
! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_fm
! USE ropp_1dvar, not_this => ropp_1dvar_levmarq_bangle
  USE ropp_1dvar
  USE matrix

  IMPLICIT NONE

  TYPE(Obs1dBangle), INTENT(inout)      :: obs          ! Observation data
  TYPE(State1dFM),   INTENT(inout)      :: bg           ! Background data
  TYPE(State1dFM),   INTENT(inout)      :: state        ! State vector
  TYPE(VarConfig),   INTENT(in)         :: config       ! Configuration options
  TYPE(VarDiag),     INTENT(inout)      :: diag         ! Diagnostic output

  REAL(wp)                              :: J            ! Cost function value
  TYPE(State1dFM)                       :: x            ! Control vector    
  TYPE(State1dFM)                       :: xmin         ! Control at minimum J
  TYPE(Obs1dBangle)                     :: y            ! Forward model obs
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K            ! K-matrix
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_x      ! Change of state
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_y      ! Change of observation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: dJ_dx        ! Cost function gradient
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: diag_d2J     ! Diagonal d2J/dx2
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: d2J_dx2      ! 2nd derivative cost fn
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KO           ! K O^-1
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Bm1          ! B^-1
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Om1          ! O^-1
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_J      ! Change of cost fn
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: state_last   ! Previous state vector
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: state_sigma  ! State std deviation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_state  ! Change of state vector
  REAL(wp)                              :: J_last       ! Previous cost function
  REAL(wp)                              :: J_min        ! Minimum cost function
  INTEGER                               :: i_pointer    ! Value index
  INTEGER                               :: i            ! Counter
  INTEGER                               :: n_iter       ! Number of iterations
  INTEGER                               :: nobs         ! Number of observations
  INTEGER                               :: nstate       ! No. of state elements
  LOGICAL                               :: marq         ! Flag to L-M minimise
  REAL(wp)                              :: lambda       ! Iteration factor

  CHARACTER(len =   4)                  :: it_str
  CHARACTER(len =  15)                  :: ch_str, co_str
  CHARACTER(len = 256)                  :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_levmarq')

! 1.3 Initialise rolling buffers and pointer for convergence checks
! -----------------------------------------------------------------
      
  i_pointer = 0
  n_iter    = 0
  lambda    = 1.0E-4_wp
  IF (bg%direct_ion) lambda = 1.0E-2_wp
  marq      = .FALSE.
  J         = 0.0_wp
  J_last    = 1.0E30_wp
  J_min     = J_last

  nstate = SIZE(bg%state)
  nobs   = SIZE(obs%bangle)

  ALLOCATE(K(nobs, nstate))
  ALLOCATE(delta_x(nstate))
  ALLOCATE(delta_y(nobs))
  ALLOCATE(dJ_dx(nstate))
  ALLOCATE(diag_d2J(nstate))
  ALLOCATE(d2J_dx2(nstate,nstate))
  ALLOCATE(KO(nstate,nobs))
  ALLOCATE(Bm1(nstate,nstate))
  ALLOCATE(Om1(nobs,nobs))

  IF (ALLOCATED(delta_J)) DEALLOCATE(delta_J)
  ALLOCATE(delta_J(config%conv_check_n_previous))
  delta_J(:) = 0.0_wp
  
  IF (ALLOCATED(state_last)) DEALLOCATE(state_last)
  ALLOCATE(state_last(SIZE(bg%state)))
  state_last(:) = 0.0_wp
  
  IF (ALLOCATED(state_sigma)) DEALLOCATE(state_sigma)
  ALLOCATE(state_sigma(SIZE(bg%state)))
  DO i = 1, SIZE(state_sigma)
     state_sigma(i) = SQRT(bg%cov%d(i + i*(i-1)/2)) ! Direct read from matrix_pp
  ENDDO
  
  IF (ALLOCATED(delta_state)) DEALLOCATE(delta_state)
  ALLOCATE(delta_state(config%conv_check_n_previous))
  delta_state(:) = 0.0_wp

  Bm1 = matrix_invert(bg%cov)

  Om1 = matrix_invert(obs%cov)

 ! 1.4 First guess state
 ! ---------------------

   x = bg
   xmin = bg

! 1.5 Main minimisation iteration loop
! ------------------------------------

   DO WHILE(J <= J_min)

      ! 1.5.1 Update iteration counter

      n_iter = n_iter + 1

! 1.6 Compute cost function
! -------------------------
      
      ! 1.6.1 Calculate pseudo observations
      ! -----------------------------------

      ALLOCATE(y%bangle(SIZE(obs%bangle)))
      
      y % g_sfc        = obs % g_sfc
      y % r_earth      = obs % r_earth
      y % r_curve      = obs % r_curve
      y % undulation   = obs % undulation
      
      CALL copy_alloc(obs % impact,  y % impact)
      CALL copy_alloc(obs % weights, y % weights)
      
      ! 1.6.2 Forward model
      ! -------------------
      
      IF(ASSOCIATED(x%ak))THEN
         CALL ropp_fm_state2state_ecmwf(x)
      ELSE
         CALL ropp_fm_state2state_meto(x)
      ENDIF
      
      CALL ropp_fm_bangle_1d(x, y)

      ! 1.6.3 Compute cost function
      ! ---------------------------

      delta_x = (x%state - bg%state)

      delta_y = (y%bangle - obs%bangle) * y%weights

      J  = 0.5_wp * DOT_PRODUCT(delta_y, MATMUL(Om1, delta_y)) + &
           0.5_wp * DOT_PRODUCT(delta_x, MATMUL(Bm1, delta_x))

      IF (n_iter == 1) diag % J_init = J

      IF (J > J_last) marq = .TRUE.

! 1.7 Levenberg-Marquardt minimisation
! ------------------------------------

      IF(.NOT. marq)THEN

         ! 1.7.1 Normal Newtonian iteration
         ! --------------------------------

         lambda = 0.1_wp * lambda
         marq = .FALSE.

         ! 1.7.2 Evaluate K gradient matrix for current x
         ! ----------------------------------------------
         
         CALL ropp_fm_bangle_1d_grad(x, y, K)

         ! 1.7.3 Calculate -dJ_dx vector and d2J_dx2 matrix at x
         ! -----------------------------------------------------

         delta_x = x%state - bg%state

         delta_y = (obs%bangle - y%bangle) * y%weights

         WHERE(ABS(delta_y) > 50.0_wp)
            delta_y = 0.0_wp
         END WHERE

         KO = MATMUL(TRANSPOSE(K), Om1)
         dJ_dx = MATMUL(KO, delta_y) - MATMUL(Bm1, delta_x)
         d2J_dx2 = Bm1 + MATMUL(KO, K)

         ! 1.7.4 Store inverse of solution covariance matrix
         ! -------------------------------------------------

         DO i=1,SIZE(bg%state)
            diag_d2J(i) = d2J_dx2(i,i)
         ENDDO

      ELSE
        
        ! 1.7.5 Levenberg-Marquardt iteration
        ! ----------------------------------- 
        !   The previous increment increased the value of the penalty function. 
        !   Use previous values of -dJ_dx, d2J_dx2 and adjust value of lambda

         marq = .TRUE.
         lambda = 100.0_wp * lambda

      ENDIF

      ! 1.7.6 Levenberg-Marquardt adjustment to diagonal terms
      ! ------------------------------------------------------

      DO i=1, SIZE(bg%state)
         d2J_dx2(i,i) = diag_d2J(i) * (lambda+1.0_wp)
      ENDDO
      
      ! 1.7.7 Solve matrix equation d2J_dx2 . dx = -dJ_dx and update state
      ! ------------------------------------------------------------------

      delta_x = matrix_solve(d2J_dx2, dJ_dx)

      IF (x%use_logq) THEN
        x%state = x%state + SIGN(MIN(ABS(delta_x), ABS(x%state/2.0_wp)), &
                                 delta_x)
      ELSE
        x%state = x%state + delta_x
      END IF

      IF (bg%direct_ion) THEN

        i = nstate - 2
        IF (x%state(i) < ropp_ZERO) THEN
          CALL message(msg_warn, "Levenberg-Marquardt solver returns " // &
            "Ne_max < 0 ... suggest examining final value. \n")
        END IF

        i = nstate - 1
        IF (x%state(i) < 0.01_wp*bg%state(i)) THEN
          CALL message(msg_warn, "Levenberg-Marquardt solver returns " // &
            "H_peak < 1% of background ... resetting to background value. \n")
            x%state(i) = bg%state(i)
        END IF

        i = nstate
        IF (x%state(i) < 0.01_wp*bg%state(i)) THEN
          CALL message(msg_warn, "Levenberg-Marquardt solver returns " // &
            "H_width < 1% of background ... resetting to background value. \n")
            x%state(i) = bg%state(i)
        END IF

      END IF

! 1.8 Update state vector variables - NWP model-dependent
! -------------------------------------------------------

      IF(ASSOCIATED(x%ak))THEN
         CALL ropp_fm_state2state_ecmwf(x)
      ELSE
         CALL ropp_fm_state2state_meto(x)
      ENDIF

! 1.9 Keep track of current state
! -------------------------------
      
      IF( J_min > J ) THEN
         J_min = J
         xmin = x
      ENDIF

      IF (config % conv_check_apply) THEN
         
         i_pointer = MOD(i_pointer + 1, config % conv_check_n_previous)
         IF (i_pointer == 0) i_pointer = config % conv_check_n_previous

         delta_J(i_pointer)     = J_last - J
         J_last = J
         delta_state(i_pointer) = MAXVAL(ABS(state_last - x%state)/state_sigma)
         state_last = x%state

      ENDIF

! 1.10 Check for convergence
! --------------------------

      IF (config % conv_check_apply) THEN

         IF (config%minropp%impres == 0)THEN
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
                 '   n_iter = ' // it_str // '   J = ' // ch_str //  &
                 '   max(relative change in state) = ' // co_str)
         ENDIF

         IF (MAXVAL(delta_state) < config%conv_check_max_delta_state    &
              .AND. n_iter > config % conv_check_n_previous) THEN
            WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_state  
            ch_str = ADJUSTL(ch_str)
            WRITE(it_str, '(i2)')    config%conv_check_n_previous       
            it_str = ADJUSTL(it_str)
            CALL message(msg_cont, '')
            CALL message(msg_info,                                             &
              'Convergence assumed to be achieved as the state vector did ' // &
              'not change by more\n   ' // 'than ' // TRIM(ch_str) // ' ' //   &
              'relative to the assumed background errors for the last ' //     &
                 TRIM(it_str) // ' iterations.\n')
            EXIT
         ELSE IF (MAXVAL(ABS(delta_J)) < config%conv_check_max_delta_J   &
              .AND. n_iter > config % conv_check_n_previous) THEN
            WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_J  
            ch_str = ADJUSTL(ch_str)
            WRITE(it_str, '(i2)')    config%conv_check_n_previous   
            it_str = ADJUSTL(it_str)
            CALL message(msg_cont, '')
            CALL message(msg_info,                                             &
             'Convergence assumed to be achieved as the cost function did ' // &
             'not change by more\n   ' // 'than ' // TRIM(ch_str) //           &
             ' for the last ' // TRIM(it_str) // ' iterations.\n')
            EXIT

         ENDIF

      ENDIF

   ENDDO    ! end main iteration loop

! 1.11 Copy solution back to state
! --------------------------------
   
   state = xmin
   
! 1.12 Diagnostic data
! --------------------

   diag % J        = J_min

   IF (COUNT(obs % weights > 0.0_wp) > 0) THEN
     diag % J_scaled = 2.0_wp * J_min / REAL(COUNT(obs % weights > 0.0_wp), wp)
   ENDIF

   diag % n_iter   = n_iter

   ALLOCATE (diag%J_bgr(SIZE(state%state)))
   delta_x = state%state - bg%state
   diag % J_bgr = 0.5_wp * delta_x * matrix_solve(bg%cov, delta_x)

! 1.13 Clean up
! -------------

   DEALLOCATE(K)
   DEALLOCATE(delta_x)
   DEALLOCATE(delta_y)
   DEALLOCATE(dJ_dx)
   DEALLOCATE(diag_d2J)
   DEALLOCATE(d2J_dx2)
   DEALLOCATE(KO)
   DEALLOCATE(Bm1)
   DEALLOCATE(Om1)

   CALL message_set_routine(routine)


END SUBROUTINE ropp_1dvar_levmarq_bangle

!***

!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_levmarq_refrac(obs, bg, state, config, diag)
  
! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_fm
! USE ropp_1dvar, not_this => ropp_1dvar_levmarq_refrac
  USE ropp_1dvar
  USE matrix

  IMPLICIT NONE

  TYPE(Obs1dRefrac), INTENT(inout)      :: obs          ! Observation data
  TYPE(State1dFM),   INTENT(inout)      :: bg           ! Background data
  TYPE(State1dFM),   INTENT(inout)      :: state        ! State vector
  TYPE(VarConfig),   INTENT(in)         :: config       ! Configuration options
  TYPE(VarDiag),     INTENT(inout)      :: diag         ! Diagnostic output

  REAL(wp)                              :: J            ! Cost function value
  TYPE(State1dFM)                       :: x            ! Control vector
  TYPE(State1dFM)                       :: xmin         ! Control at minimum J
  TYPE(Obs1dRefrac)                     :: y            ! Forward model obs
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K            ! K-matrix
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_x      ! Change of state
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_y      ! Change of observation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: dJ_dx        ! Cost function gradient
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: diag_d2J     ! Diagonal d2J/dx2
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: d2J_dx2      ! 2nd derivative cost
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KO           ! K O^-1
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Bm1          ! B^-1
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Om1          ! O^-1
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_J      ! Change of cost fn
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: state_last   ! Previous state vector
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: state_sigma  ! State std deviation 
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: delta_state  ! Change of state vector
  REAL(wp)                              :: J_last       ! Previous cost function
  REAL(wp)                              :: J_min        ! Minimum cost function
  INTEGER                               :: i_pointer    ! Value index
  INTEGER                               :: i            ! Counter
  INTEGER                               :: n_iter       ! Number of iterations
  INTEGER                               :: nobs         ! Number of obserations
  INTEGER                               :: nstate       ! No. of state elements
  LOGICAL                               :: marq         ! Flag to L-M minimise
  REAL(wp)                              :: lambda       ! Iteration factor

  CHARACTER(len =   4)                  :: it_str
  CHARACTER(len =  15)                  :: ch_str, co_str
  CHARACTER(len = 256)                  :: routine

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_levmarq')

! 2.3 Initialise rolling buffers and pointer for convergence checks
! -----------------------------------------------------------------
      
  i_pointer = 0
  n_iter    = 0 
  lambda    = 1.0E-4_wp
  marq      = .FALSE.
  J         = 0.0_wp
  J_last    = 1.0E30_wp
  J_min     = J_last

  nstate = SIZE(bg%state)
  nobs   = SIZE(obs%refrac)
  
  ALLOCATE(K(nobs, nstate))
  ALLOCATE(delta_x(nstate))
  ALLOCATE(delta_y(nobs))
  ALLOCATE(dJ_dx(nstate))
  ALLOCATE(diag_d2J(nstate))
  ALLOCATE(d2J_dx2(nstate,nstate))
  ALLOCATE(KO(nstate,nobs))
  ALLOCATE(Bm1(nstate,nstate))
  ALLOCATE(Om1(nobs,nobs))

  Bm1 = matrix_invert(bg%cov)
  Om1 = matrix_invert(obs%cov)

  IF (ALLOCATED(delta_J)) DEALLOCATE(delta_J)
  ALLOCATE(delta_J(config%conv_check_n_previous))
  delta_J(:) = 0.0_wp
  
  IF (ALLOCATED(state_last)) DEALLOCATE(state_last)
  ALLOCATE(state_last(SIZE(bg%state)))
  state_last(:) = 0.0_wp
  
  IF (ALLOCATED(state_sigma)) DEALLOCATE(state_sigma)
  ALLOCATE(state_sigma(SIZE(bg%state)))
  DO i = 1, SIZE(state_sigma)
     state_sigma(i) = SQRT(bg%cov%d(i + i*(i-1)/2)) ! Direct read from matrix_pp
  ENDDO
  
  IF (ALLOCATED(delta_state)) DEALLOCATE(delta_state)
  ALLOCATE(delta_state(config%conv_check_n_previous))
  delta_state(:) = 0.0_wp

! 2.4 First guess state
! --------------------- 
   
   x = bg 
   xmin = bg
     
! 2.5 Main minimisation iteration loop
! ------------------------------------

   DO WHILE (J <= J_min)  

      ! 2.5.1 Update iteration counter

      n_iter = n_iter + 1
      
! 2.6 Compute cost function
! -------------------------
      
      ! 2.6.1 Calculate pseudo observations
      ! -----------------------------------

      ALLOCATE(y%refrac(SIZE(obs%refrac)))

      CALL copy_alloc(obs % geop,    y % geop)
      CALL copy_alloc(obs % weights, y % weights)
      
      ! 2.6.2 Forward model
      ! -------------------

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

      ! 2.6.3 Compute cost function
      ! ---------------------------

      delta_x = (x%state - bg%state)
      delta_y = (y%refrac - obs%refrac) * y%weights

      J  = 0.5_wp * DOT_PRODUCT(delta_y, MATMUL(Om1, delta_y)) &
           + 0.5_wp * DOT_PRODUCT(delta_x, MATMUL(Bm1, delta_x))
      
      IF (n_iter == 1) diag % J_init = J  
      
      IF (J > J_last) marq = .TRUE.

! 2.7 Levenberg-Marquardt minimisation
! ------------------------------------

      IF(.NOT. marq)THEN

         ! 2.7.1 Normal Newtonian iteration
         ! --------------------------------
         
         lambda = 0.1_wp * lambda
         marq = .FALSE.

         ! 2.7.2 Evaluate K gradient matrix for current x
         ! ----------------------------------------------
         
         CALL ropp_fm_refrac_1d_grad(x, y, K)

         ! 2.7.3 Calculate -dJ_dx vector and d2J_dx2 matrix at x
         ! -----------------------------------------------------

         delta_x = x%state - bg%state
         delta_y = (obs%refrac - y%refrac) * y%weights
         
         WHERE(ABS(delta_y) > 50.0_wp) 
            delta_y = 0.0_wp
         END WHERE

         KO = MATMUL(TRANSPOSE(K), Om1)
         dJ_dx = MATMUL(KO, delta_y) - MATMUL(Bm1, delta_x)
         d2J_dx2 = Bm1 + MATMUL(KO, K)
 
         ! 2.7.4 Store inverse of solution covariance matrix
         ! -------------------------------------------------
         
         DO i=1,SIZE(bg%state)
            diag_d2J(i) = d2J_dx2(i,i)
         ENDDO
                  
      ELSE
        
        ! 2.7.5 Levenberg-Marquardt iteration
        ! ----------------------------------- 
        !   The previous increment increased the value of the penalty function. 
        !   Use previous values of -dJ_dx, d2J_dx2 and adjust value of lambda

         marq = .TRUE.
         lambda = 100.0_wp * lambda
         
      ENDIF
      
      ! 2.7.6 Levenberg-Marquardt adjustment to diagonal terms
      ! ------------------------------------------------------

      DO i=1, SIZE(bg%state)
         d2J_dx2(i,i) = diag_d2J(i) * (lambda+1.0_wp)
      ENDDO
      
      ! 2.7.7 Solve matrix equation d2J_dx2 . dx = -dJ_dx and update state
      ! ------------------------------------------------------------------
      
      if (x%use_logq) then
        delta_x = matrix_solve(d2J_dx2, dJ_dx)
        x%state = x%state + SIGN( MIN(ABS(delta_x),ABS(x%state/2.0_wp)), &
                             delta_x)
      else
        x%state = x%state + matrix_solve(d2J_dx2, dJ_dx)
      endif

! 2.8 Update state vector variables - model dependent
! ---------------------------------

      IF(ASSOCIATED(x%ak))THEN
         CALL ropp_fm_state2state_ecmwf(x) 
      ELSE 
         CALL ropp_fm_state2state_meto(x) 
      ENDIF

! 2.9 Keep track of current state
! -------------------------------
      
      IF( J_min > J ) THEN
         J_min = J  
         xmin = x
      ENDIF 
      
      IF (config % conv_check_apply) THEN
         
         i_pointer = MOD(i_pointer + 1, config % conv_check_n_previous)
         IF (i_pointer == 0) i_pointer = config % conv_check_n_previous

         delta_J(i_pointer)     = J_last - J
         J_last = J
         delta_state(i_pointer) = MAXVAL(ABS(state_last - x%state)/state_sigma)
         state_last = x%state
         
      ENDIF

! 2.10 Check for convergence
! --------------------------

      IF (config % conv_check_apply) THEN
                  
         IF (config%minropp%impres == 0)THEN
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
                 '   n_iter = ' // it_str // '   J = ' // ch_str //  &
                 '   max(relative change in state) = ' // co_str)
         ENDIF
         
         
         IF (MAXVAL(delta_state) < config%conv_check_max_delta_state    &
              .AND. n_iter > config % conv_check_n_previous) THEN
            WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_state  
            ch_str = ADJUSTL(ch_str)
            WRITE(it_str, '(i2)')    config%conv_check_n_previous       
            it_str = ADJUSTL(it_str)
            CALL message(msg_cont, '')
            CALL message(msg_info,                                             &
              'Convergence assumed to be achieved as the state vector did ' // &
              'not change by more\n   ' // 'than ' // TRIM(ch_str) // ' ' //   &
              'relative to the assumed background errors for the last ' //     &
              TRIM(it_str) // ' iterations.\n')
            EXIT
         ELSE IF (MAXVAL(ABS(delta_J)) < config%conv_check_max_delta_J   &
              .AND. n_iter > config % conv_check_n_previous) THEN
            WRITE(ch_str, '(g15.5)') config%conv_check_max_delta_J  
            ch_str = ADJUSTL(ch_str)
            WRITE(it_str, '(i2)')    config%conv_check_n_previous   
            it_str = ADJUSTL(it_str)
            CALL message(msg_cont, '')
            CALL message(msg_info,                                             &
             'Convergence assumed to be achieved as the cost function did ' // &
             'not change by more\n   ' // 'than ' // TRIM(ch_str) //           &
             ' for the last ' // TRIM(it_str) // ' iterations.\n')
            EXIT
         ENDIF
         
      ENDIF

   ENDDO    ! end main iteration loop

! 2.11 Copy solution back to state
! --------------------------------
   
   state = xmin
   
! 2.12 Diagnostic data
! --------------------

   diag % J        = J_min

   IF (COUNT(obs % weights > 0.0_wp) > 0) THEN
     diag % J_scaled = 2.0_wp * J_min / REAL(COUNT(obs % weights > 0.0_wp), wp)
   ENDIF

   diag % n_iter   = n_iter

   ALLOCATE (diag%J_bgr(SIZE(state%state)))
   delta_x = state%state - bg%state
   diag % J_bgr = 0.5_wp * delta_x * matrix_solve(bg%cov, delta_x)

! 2.13 Clean up
! -------------

   DEALLOCATE(K)
   DEALLOCATE(delta_x)
   DEALLOCATE(delta_y)
   DEALLOCATE(dJ_dx)
   DEALLOCATE(diag_d2J)
   DEALLOCATE(d2J_dx2)
   DEALLOCATE(KO)
   DEALLOCATE(Bm1)
   DEALLOCATE(Om1)

   CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_levmarq_refrac
