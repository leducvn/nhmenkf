! $Id: ropp_pp_bangle2phase.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_bangle2phase(time, r_leo, r_gns, r_coc, impact, bangle,   &
                                phase, dphi, impact_LM, IL, IU) 

!****s* Preprocessing/ropp_pp_bangle2phase *
!
! NAME
!    ropp_pp_bangle2phase - Compute excess phase as a function of time from 
!                           bending angle as function of impact parameter
!                   
! SYNOPSIS
!    call ropp_pp_bangle2phase(time, r_leo, r_gns, r_coc, impact, bangle, 
!                              phase, dphi, impact_LM, IL, IU)
! 
! DESCRIPTION
!    This routine computes model excess phase from bending angle data
!      1. Compute t(p) for given satellite trajectories
!      2. Compute p(t) -> d(t) -> Phi(t)
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! cartesian LEO coordinates (ECF)
!    real(wp), dimension(:,:) :: r_gns     ! cartesian GPS coordinates (ECF)
!    real(wp), dimension(:)   :: r_coc     ! cartesian centre curvature (ECF)
!    real(wp), dimension(:)   :: impact    ! input impact parameters (m)
!    real(wp), dimension(:)   :: bangle    ! input bending angle profile (rad)
!
! OUTPUT
!    real(wp), dimension(:)   :: bangle    ! output bending angle profile (rad)
!    real(wp), dimension(:)   :: phase     ! excess phase as function time (m)
!    real(wp), optional, dimension(:) :: dphi      ! derivative phase with time
!    real(wp), optional, dimension(:) :: impact_LM ! model impact parameter (m)
!    integer,  optional               :: IL  ! lower index valid model interval
!    integer,  optional               :: IU  ! upper index valid model interval
!
! NOTES
!
! REFERENCES
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
! USE ropp_pp, not_this => ropp_pp_bangle2phase
  USE ropp_pp
  USE ropp_utils, ONLY: impact_parameter

  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)    :: time    ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_leo   ! LEO coordinates (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_gns   ! GPS coordinates (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: r_coc   ! Centre curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: impact  ! Impact parameter (m)
  REAL(wp), DIMENSION(:),   INTENT(inout) :: bangle  ! Bending angle (m)
  REAL(wp), DIMENSION(:),   INTENT(out)   :: phase   ! Excess phase (m)
  REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: dphi      ! dPhi/dt (m/s)
  REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: impact_LM ! Model IP (m)
  INTEGER,  OPTIONAL, INTENT(out)       :: IL ! Lower index valid model interval
  INTEGER,  OPTIONAL, INTENT(out)       :: IU ! Upper index valid model interval

  INTEGER, PARAMETER     :: nv = 5 ! Polynomial degree for calculating velocity
  INTEGER                :: nh     ! Number of data (bangle,impact)
  INTEGER                :: n      ! Number of data (time,phase,dPhi,impact_LM)
  INTEGER                :: i      ! Data index
  REAL(wp)               :: ps1    ! First straight-line impact parameter
  REAL(wp)               :: psN    ! Last straight-line impact parameter
  REAL(wp)               :: t_init ! Initial approximation for time
  REAL(wp), DIMENSION(3) :: xgnsT  ! GPS coordinates interpolated
  REAL(wp), DIMENSION(3) :: vgnsT  ! GPS velocity interpolated
  REAL(wp), DIMENSION(3) :: xleoT  ! LEO coordinates interpolated
  REAL(wp), DIMENSION(3) :: vleoT  ! LEO velocity interpolated
  INTEGER                :: iup    ! Index of upper point
  INTEGER                :: ilw    ! Index of lower point
  INTEGER                :: di     ! Index increment
  INTEGER                :: ih1    ! Lower index H of valid model interval
  INTEGER                :: ih2    ! Upper index H of valid model interval
  INTEGER                :: time1  ! Lower border of valid model interval
  INTEGER                :: time2  ! Upper border of valid model interval
  REAL(wp)               :: phase_min  ! Phase base

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: dphi_dt    ! Interpolated dphi_dt
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: t_norm     ! Normalised time
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KV         ! Regression matrix for v
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_vleo ! Regression coeffs vleo
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_vgns ! Regression coeffs vgns

  REAL(wp), DIMENSION(:), ALLOCATABLE   :: bangle_nh  ! Bangle on hi-res grid
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: time_nh    ! Time(impact)
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: doppler_nh ! Doppler(impact)
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: range_nh   ! Pseudo-range(impact)
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: phase_nh   ! ExcessPhase(impact)
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: dphi_nh    ! dPhi/dt
  REAL(wp) :: d0
  REAL(wp), DIMENSION(3) :: ugl

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

  n = SIZE(time)
  nh = SIZE(impact)

  ALLOCATE(t_norm(n))
  ALLOCATE(dphi_dt(n))

  ALLOCATE(bangle_nh(nh))
  ALLOCATE(time_nh(nh))
  ALLOCATE(doppler_nh(nh))
  ALLOCATE(range_nh(nh))
  ALLOCATE(phase_nh(nh))
  ALLOCATE(dphi_nh(nh))

  DO i=1,n
     t_norm(i) = (time(i) - time(1))/(time(n) - time(1))
  ENDDO

  ps1 = impact_parameter(r_leo(1,:) - r_coc(:), r_gns(1,:) - r_coc)
  psn = impact_parameter(r_leo(n,:) - r_coc(:), r_gns(n,:) - r_coc)

!-------------------------------------------------------------------------------
! 2. Calculate model phase from MSIS+humidity bending angle profile
!-------------------------------------------------------------------------------
  
  ! 2.1 Trajectory interpolation
  
  ALLOCATE(KV(n, 0:nv))
  
  CALL ropp_pp_init_polynomial(t_norm, KV)

  ! 2.2 Perform regression on positions
  DO i=1,3
     CALL ropp_pp_regression(KV, r_leo(:,i), coeff_vleo(:,i))
     CALL ropp_pp_regression(KV, r_gns(:,i), coeff_vgns(:,i))
  ENDDO

  ! 2.3 Perform regression on residual positions to gain higher accuracy
  DO i=1,3
    CALL ropp_pp_residual_regression(KV, t_norm, r_leo(:,i), coeff_vleo(:,i))
    CALL ropp_pp_residual_regression(KV, t_norm, r_gns(:,i), coeff_vgns(:,i))
  ENDDO

  DEALLOCATE(KV)

  ! 2.4 Compute doppler shift
  
  ih1 = nh
  ih2 = 1
  
  IF (impact(1) > impact(nh)) THEN
     iup = 1
     ilw = nh
     di  = 1
  ELSE
     iup = nh
     ilw = 1
     di  = -1
  ENDIF
  
  IF (ps1 > psN) THEN
     t_init  = time(1)
  ELSE
     t_init  = time(n)
  ENDIF

  DO i=iup,ilw,di

  ! 2.5 Find location of trajectory point for current impact parameter and 
  !     bending angle
     
    CALL locate_time(impact(i), bangle(i), time, coeff_vleo, coeff_vgns,    &
                     r_coc, t_init, time_nh(i))

     ih1 = MIN(i, ih1)
     ih2 = MAX(i, ih2)
     t_init = time_nh(i)

     CALL ropp_pp_interpolate_trajectory(time, coeff_vleo, coeff_vgns, r_coc,  &
                                         time_nh(i), xleoT, vleoT, xgnsT, vgnsT)
     
     CALL ropp_pp_impact2doppler(xleoT - r_coc, vleoT, xgnsT - r_coc, vgnsT,  &
                                 impact(i), doppler_nh(i), bangle_nh(i))

     range_nh(i) = SQRT(SUM((xgnsT - xleoT)**2))   
     
     ugl(:) = (xleoT(:) - xgnsT(:))/range_nh(i)
     d0 = DOT_PRODUCT(-(vleoT(:) - vgnsT(:)),ugl(:))/c_light   
     dphi_nh(i) = - c_light * (doppler_nh(i) - d0)

  ENDDO

  ! 2.6 Compute excess phase

  phase_nh(:) = 0.0_wp
  
  DO i=ih1+1,ih2
     phase_nh(i) = phase_nh(i-1) - c_light *   &
               (doppler_nh(i)+doppler_nh(i-1))*(time_nh(i)-time_nh(i-1))/2.0_wp
  ENDDO

  phase_nh(ih1:ih2) = phase_nh(ih1:ih2) - range_nh(ih1:ih2)

  phase_min    = MINVAL(phase_nh(:))

  phase_nh(ih1:ih2) = phase_nh(ih1:ih2) - phase_min

!-------------------------------------------------------------------------------
! 3. Interpolate 
!-------------------------------------------------------------------------------
  
  DO i=1,n
    CALL ropp_pp_interpol(time_nh(ih1+1:ih2-1),time(i),dphi_nh(ih1+1:ih2-1), dphi_dt(i))
  ENDDO
  
  IF ( PRESENT(impact_LM) ) THEN
    DO i=1,n
     CALL ropp_pp_interpol(time_nh(ih1+1:ih2-1), time(i),      &
                           impact(ih1+1:ih2-1),  impact_LM(i))
   ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! 4. Integrate to compute excess phase
!-------------------------------------------------------------------------------

  phase(1) = 0.0_wp
  DO i=2,n
     phase(i) = phase(i-1)+(dphi_dt(i)+dphi_dt(i-1))*(time(i)-time(i-1))/2.0_wp
  ENDDO
  
  phase_min = MINVAL(phase(:))
  phase(:) = phase(:) - phase_min

  IF ( PRESENT(dPhi) ) THEN
     dPhi = dphi_dt
  ENDIF

!-------------------------------------------------------------------------------
! 5. Determine valid model interval range
!-------------------------------------------------------------------------------
  
  IF (PRESENT(il) .AND. PRESENT(iu)) THEN
     
     time1 = MIN(time_nh(ih1+1),time_nh(ih2-1))
     time2 = MAX(time_nh(ih1+1),time_nh(ih2-1))
     
     il = n
     iu = 1
     
     DO i=1,n
        IF ((time1 <= time(i)) .AND. (time(i) <= time2)) THEN
           il = MIN(i,il)
           iu = MAX(i,iu)
        END IF
     END DO
     
  ENDIF

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(t_norm)
  DEALLOCATE(dphi_dt)  

  DEALLOCATE(bangle_nh)
  DEALLOCATE(time_nh)
  DEALLOCATE(doppler_nh)
  DEALLOCATE(range_nh)
  DEALLOCATE(phase_nh)
  DEALLOCATE(dphi_nh)

CONTAINS 

!-------------------------------------------------------------------------------
! 10. Location of trajectory point for given impact parameter and bending
!-------------------------------------------------------------------------------

  !! Newtonian-iterative solution of
  !!       XGPS^XLEO - acos(impact/|XGPS|) - acos(impact/|XLEO|) = bangle

  SUBROUTINE locate_time(impact,bangle,time,cleo,cgns,r_coc,t_init,t_final)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp, ONLY: ropp_pp_interpolate_trajectory

    IMPLICIT NONE

    ! 10.1 Declarations

    REAL(wp),                  INTENT(in)  :: impact  ! Impact parameter (m)
    REAL(wp),                  INTENT(in)  :: bangle  ! Bending angle (rad)
    REAL(wp), DIMENSION(:),    INTENT(in)  :: time    ! Time of samples (s)
    REAL(wp), DIMENSION(0:,:), INTENT(in)  :: cgns    ! Regression coeffs VGNS
    REAL(wp), DIMENSION(0:,:), INTENT(in)  :: cleo    ! Regression coeffs VLEO
    REAL(wp), DIMENSION(:),    INTENT(in)  :: r_coc   ! Centre of curvature
    REAL(wp),                  INTENT(in)  :: t_init  ! Initial approximation
    REAL(wp),                  INTENT(out) :: t_final ! Located time

    REAL(wp), PARAMETER    :: dl = 10.0_wp  ! Spatial scale for differentiation
    REAL(wp), PARAMETER    :: el = 1e-4_wp  ! Spatial accuracy (m)
    REAL(wp)               :: rg, rl        ! GPS and LEO radius (m)
    REAL(wp)               :: dt            ! Finite difference
    REAL(wp)               :: theta         ! Satellite-to-satellite angle
    REAL(wp)               :: dft           ! Function derivative
    REAL(wp)               :: st            ! Time increment
    INTEGER                :: i, k          ! Indices
    REAL(wp), DIMENSION(3) :: xgns, xleo    ! Cartesian GPS and LEO coordinates
    REAL(wp), DIMENSION(3) :: vgns, vleo    ! Cartesian GPS and LEO velocity  
    REAL(wp), DIMENSION(2) :: func          ! Cosequent values of function

    ! 10.2 Initialization

    t_final = t_init
    dt = 0.0_wp

    ! 10.3 Find initial interpolated coordinates
    
    CALL ropp_pp_interpolate_trajectory(time, cleo, cgns, r_coc, t_init,   &
                                        xleo, vleo, xgns, vgns, theta)  
    dt = dl / MAX( SQRT(SUM(vgns(:)**2)), SQRT(SUM(vleo(:)**2)))
    
    ! 10.4 Iterative solution

    DO k = 1,20

       DO i=1,2

          CALL ropp_pp_interpolate_trajectory(time, cleo, cgns, r_coc,       &
                                              t_final+(i-1)*dt,              &
                                              xleo, vleo, xgns, vgns, theta )
          
          rg = SQRT(SUM((xgns(:)-r_coc)**2))
          rl = SQRT(SUM((xleo(:)-r_coc)**2))
          
          func(i) = theta - ACOS(impact/rg) - ACOS(impact/rl) - bangle

       ENDDO

       dft = (func(2) - func(1))/dt
       st = -func(1)/dft

       IF (ABS(st) < 20*ABS(time(SIZE(time))-time(1))) THEN
          t_final = t_final + st
       ENDIF
       
       IF (ABS(st/dt) < el/dl) THEN
          EXIT 
       ENDIF
       
    ENDDO

  END SUBROUTINE locate_time

END SUBROUTINE ropp_pp_bangle2phase
