! $Id: ropp_pp_fit_model_refraction.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_fit_model_refraction(impact_LC, bangle_LC, impact_model,   &
                                        bangle_model, config)

!****s* ModelRefraction/ropp_pp_fit_model_refraction *
!
! NAME
!    ropp_pp_fit_model_refraction - Fit model bending angle profile with 
!                                   observed bending angles
!
! SYNOPSIS
!    call ropp_pp_fit_model_refraction(impact_LC, bangle_LC, 
!                                      impact_model, bangle_model, config)
! 
! DESCRIPTION
!    This subroutine calculates a fitting factor of model bending angle profile
!    with observed bending angles by linear regression in height interval 
!    40-60 km.
!
! INPUTS
!    real(wp), dim(:) :: impact_LC     Observed impact parameters (m)   
!    real(wp), dim(:) :: bangle_LC     Observed bending angles (rad)
!    real(wp), dim(:) :: impact_model  Model impact parameter (m)
!    real(wp), dim(:) :: bangle_model  Model bending angles (rad)
!    type(ppConfig)   :: config        Configuration parameters
!
! OUTPUT
!    real(wp), dim(:) :: bangle_model  Fitted model bending angles (rad)
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE messages
  USE ropp_utils, ONLY: WHERE
  USE ropp_pp_utils, ONLY: ropp_pp_regression
! USE ropp_pp, not_this => ropp_pp_fit_model_refraction
  USE ropp_pp_types, ONLY: PPConfig
  
  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)    :: impact_LC    ! LC impact param (m)
  REAL(wp), DIMENSION(:), INTENT(in)    :: bangle_LC    ! LC bending angle (rad)
  REAL(wp), DIMENSION(:), INTENT(in)    :: impact_model ! Model level impact
  REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_model ! Model level bangle
  TYPE(PPConfig),         INTENT(in)    :: config       ! Configuration params

  REAL(wp), DIMENSION(:), ALLOCATABLE   :: bangle_LCM   ! LC bangle on model lvl
  REAL(wp), DIMENSION(2)                :: rf           ! Regression factor
  REAL(wp)                              :: delta        ! a posterior stddev(rf)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K            ! matrix for regression
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: Y            ! regression data

  INTEGER,  DIMENSION(:), POINTER :: idx => NULL()      ! Array indices
  INTEGER                         :: n_mod, n_reg       ! Counters

  
  INTEGER            :: IFmin    ! Start index of fitting area
  INTEGER            :: IFmax    ! End index of fitting area
  INTEGER            :: Imin     ! Lower index of fitting area
  INTEGER            :: Imax     ! Upper index of fitting area
  INTEGER            :: Jmax     ! Upper index of fitting area
  INTEGER            :: di       ! Scan direction
  INTEGER            :: i        ! Counter
  CHARACTER(len=247) :: outstr

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  n_mod = SIZE(bangle_model)
  ALLOCATE(bangle_LCM(n_mod))

!-------------------------------------------------------------------------------
! 3. Interpolate observations onto model levels
!-------------------------------------------------------------------------------

  CALL ropp_pp_interpol(impact_LC, impact_model, bangle_LC, bangle_LCM)
    
!-------------------------------------------------------------------------------
! 4. Linear regression
!-------------------------------------------------------------------------------

  idx => WHERE ( impact_model > config%r_curve + config%hmin_fit .AND.    &
                 impact_model < config%r_curve + config%hmax_fit, n_reg)

  ! 4.1 1-parameter fit

  IF (config%nparm_fit == 1) THEN     

    IF (n_reg > 100) THEN
      rf(1) = SUM(bangle_model(idx)*bangle_LCM(idx))/SUM(bangle_model(idx)**2)

      delta = SQRT(SUM((bangle_LCM(idx)/bangle_model(idx) - rf(1))**2)/n_reg)
      
      rf(1) = (rf(1) + (delta/config%omega_fit)**2) /           &
               (1.0_wp + (delta/config%omega_fit)**2)
    ELSE
      rf(1) = 1.0_wp
    ENDIF
    
    bangle_model(:) = rf(1) * bangle_model(:)

    WRITE(outstr,'(2X,2(A,F10.0),A,I6)')                        &
       '1-parameter model:ob fit: From ', impact_model(idx(1))-config%r_curve,&
       ' to ', impact_model(idx(n_reg))-config%r_curve, ' No. data = ',  n_reg
    CALL message(msg_diag, outstr)
    WRITE(outstr,'(2X,A,F7.4:",",1X)')   'RF        = ', RF(1)
    CALL message(msg_diag, outstr)


  ! 4.2 2-parameter fit

  ELSE IF (config%nparm_fit == 2) THEN

    imin = MINVAL(idx) !!SUM(MINLOC(impact_model(idx)))
    imax = MAXVAL(idx) !!SUM(MAXLOC(impact_model(idx)))
    di = SIGN(1, imax-imin)

    jmax = imin

    DO i=imin,imax,di
      IF (ABS(bangle_LCM(i) - bangle_model(i)) < 0.3*bangle_model(i)) THEN
        jmax = i
      ELSE
        EXIT
      END IF
    END DO

    imax = jmax

    ifmin = MIN(imin, imax)
    ifmax = MAX(imin, imax)
    n_reg = SIZE(impact_model(ifmin:ifmax))

    IF (n_reg > 100) THEN

      ALLOCATE(K(n_reg, 2))
      ALLOCATE(Y(n_reg))

      K(:,1) = 1.0_wp
      K(:,2) = LOG(bangle_model(ifmin:ifmax))
      Y(:)   = LOG(bangle_LCM(ifmin:ifmax))
      
      CALL ropp_pp_regression(K, Y, rf)
      
      DEALLOCATE(K)
      DEALLOCATE(Y)
      
    ELSE

      rf(:) = (/ 0.0_wp, 1.0_wp /)

    ENDIF

    bangle_model(:) = EXP(rf(1) + rf(2)*LOG(bangle_model(:)))
         
    WRITE(outstr,'(2X,2(A,F10.0),A,I6)')                        &
       '2-parameter model:ob fit: From ', impact_model(Imin)-config%r_curve, &
       ' to ', impact_model(Imax)-config%r_curve, ' No. data = ',  n_reg
    CALL message(msg_diag, outstr)
    WRITE(outstr,'(2X,A,2(F7.4:",",1X))')   'RF        = ', RF(:)
    CALL message(msg_diag, outstr)
    
    ELSE
      
    CALL message(msg_warn, "Number of parameters for model fit not recognised." //  &
                           " Check config%nparm_fit")

  ENDIF

!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(idx)
  DEALLOCATE(bangle_LCM)

END SUBROUTINE ropp_pp_fit_model_refraction
