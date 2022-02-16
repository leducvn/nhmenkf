! $Id: ropp_pp_bending_angle_go.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_bending_angle_go(time, r_leo, r_gns, r_coc,            &
                                    phase_L1, phase_L2, w_smooth, filter, &
                                    impact_L1, bangle_L1, impact_L2, bangle_L2)

!****s* GeometricOptics/ropp_pp_bending_angle_go *
!
! NAME
!    ropp_pp_bending_angle_go - Calculate L1 and L2 bending angle profiles from
!                               occultation data by GEOMETRIC OPTICS
!                   
! SYNOPSIS
!    call ropp_pp_bending_angle_go(time, r_leo, r_gns, r_coc, phase_L1,
!                                  phase_L2, w_smooth, filter, impact, bangle)
! 
! DESCRIPTION
!    This routine calculates L1 and L2 bending angles.
!      1) Detrend excess phase using spline regression
!      2) Numerical differentiation of detrended excess phase
!      3) Calculation of impact parameter and bending angle from Doppler shift
!         using Snell's Law - GEOMETRIC OPTICS method
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! Relative time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! Cartesian LEO coordinates (m)
!    real(wp), dimension(:,:) :: r_gns     ! Cartesian GPS coordinates (m)
!    real(wp), dimension(:)   :: r_coc     ! Centre curvature coordinates (m)
!    integer                  :: w_smooth  ! Smoothing window (points)
!    character(len=*)         :: filter    ! Filter type (optest/slpoly)
!    real(wp), dimension(:)   :: phase_L1  ! L1 excess phase (m)
!    real(wp), dimension(:)   :: phase_L2  ! L2 excess phase (m)
!
! OUTPUT
!    real(wp), dimension(:)   :: impact_L1 ! L1 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L1 ! L1 bending angles (rad)
!    real(wp), dimension(:)   :: impact_L2 ! L2 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L2 ! L2 bending angles (rad)
! 
! NOTES
!    Method:
!        1. Detrending excess phase using spline regression.
!        2. Numerical differentiation of detrended phase using
!           optimal solution of integral equation in matrix form.
!        3. Calculation of impact parameter and bending
!           angle from Doppler shift using Snell's law.
!
! REFERENCES
!   Vorob'ev, V.V. and Krasil'nikova T.G. 1994, 
!   Estimation of the accuracy of the atmospheric refractive index recovery
!   from doppler Sshift measurements at frequencies used in the NAVSTAR system
!   Physics of the Atmosphere and Ocean (29) 602-609
!
!  Gorbunov, M.E. and Kornblueh, L. 2003,
!  Principles of variational assimilation of GNSS radio occultation data
!  Max Planck Institute Report 350 
!  http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/max_scirep_350.pdf
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
! USE ropp_pp, not_this => ropp_pp_bending_angle_go
  USE ropp_pp
  USE messages

  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! Centre curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L1  ! L1 excess phase (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L2  ! L2 excess phase (m)
  INTEGER,                  INTENT(in)  :: w_smooth  ! Smoothing window (points)
  CHARACTER(len=*),         INTENT(in)  :: filter    ! Filter type ('optest'/'slpoly')
  REAL(wp), DIMENSION(:),   INTENT(out) :: impact_L1 ! L1 impact parameters (m)
  REAL(wp), DIMENSION(:),   INTENT(out) :: bangle_L1 ! L1 bending angles (rad)
  REAL(wp), DIMENSION(:),   INTENT(out) :: impact_L2 ! L2 impact parameters (m)
  REAL(wp), DIMENSION(:),   INTENT(out) :: bangle_L2 ! L2 bending angles (rad)

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: tmp 
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: S0     ! Basic excess phase
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: trs    ! Delta spline nodes
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K      ! Basic matrix for regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KInv   ! Inverse matrix for regressn
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: b      ! Regression coefficients
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: d      ! Regression coefficients
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ST     ! Excess phase trend
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: SD     ! Detrended excess phase
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DST    ! Deriv of excess phase trend
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DSD    ! Deriv of detrended phase
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: FSD    ! Filtered detrended phase
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo   ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns   ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo   ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns   ! GPS velocity by regression

  INTEGER                :: npoints   ! Number of data points 
  INTEGER                :: i, ic     ! Counters
  INTEGER                :: nrg       ! Dimension of reduced grid
  INTEGER                :: nrf       ! Number of points of reduced grid
  REAL(wp)               :: dtf       ! Step size of reduced grid
  INTEGER                :: wf        ! Window on reduced grid
  REAL(wp), DIMENSION(3) :: U0        ! GPS-LEO direction
  REAL(wp)               :: d0        ! Vacuum doppler shift
  REAL(wp)               :: doppler   ! Doppler shift
  CHARACTER(len = 256)   :: routine

  INTEGER, PARAMETER     :: nc = 2    ! Number of channels
  INTEGER, PARAMETER     :: nr = 50   ! Default no. of nodes for spline regressn
  INTEGER, PARAMETER     :: nd = 11   ! No. of points for differentiation (odd)
  INTEGER, PARAMETER     :: np = 3    ! Polynomial degree for filtering



  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_bending_angle_go')

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  npoints = SIZE(time)
  nrf = CEILING(npoints/8000.0_wp)
  nrg = SIZE(time(1::nrf))

  ALLOCATE(S0(nc))
  ALLOCATE(trs(0:nr))
  ALLOCATE(K(nrg,0:nr))
  ALLOCATE(KInv(0:nr,nrg))
  ALLOCATE(b(nc,0:nr))
  ALLOCATE(d(nc,0:nr))
  ALLOCATE(ST(nc,npoints))
  ALLOCATE(SD(nc,npoints))
  ALLOCATE(DST(nc,npoints))  
  
!-------------------------------------------------------------------------------
! 3. Detrend excess phase
!-------------------------------------------------------------------------------

  ! 3.1 Define grid for spline regression

  DO i=0,nr
     trs(i) = time(1) + (i*1.0_wp)*(time(npoints)-time(1))/(nr*1.0_wp)
  ENDDO

  ! 3.2 Generate basic polynomials 

  CALL ropp_pp_basic_splines(time(1::nrf), trs, K)
  
  ! 3.3 Polynomial regression for L1 and L2 phase

  S0(1) = MINVAL(phase_L1) - 1.0_wp
  S0(2) = MINVAL(phase_L2) - 1.0_wp  

  ALLOCATE(tmp(nrg))
  tmp(:) = LOG(phase_L1(1::nrf)-S0(1))
  CALL ropp_pp_regression(K, tmp, b(1,:))
  tmp(:) = LOG(phase_L2(1::nrf)-S0(2))
  CALL ropp_pp_regression(K, tmp, b(2,:))
  DEALLOCATE(tmp)
  

  ! 3.4 Spline interpolation at each time point to obtain excess phase trend

  DO ic = 1,nc
     CALL ropp_pp_init_spline(trs, b(ic,:), d(ic,:))

     DO i=1,npoints
        CALL ropp_pp_interpol_spline(trs, b(ic,:), d(ic,:), time(i), &
                                       ST(ic,i), DST(ic,i))
        ST(ic,i) = EXP(ST(IC,i)) + S0(IC)          ! Trend of excess phase
     ENDDO
   ENDDO

  ! 3.5 Derivative excess phase trend

  DST(1,:) = DST(1,:) * (phase_L1(:) - S0(1))
  DST(2,:) = DST(2,:) * (phase_L2(:) - S0(2))
  
  ! 3.6 Detrended excess phase

  SD(1,:) = phase_L1(:) - ST(1,:)
  SD(2,:) = phase_L2(:) - ST(2,:)

  ! 3.7 Clean up

  DEALLOCATE(S0)
  DEALLOCATE(trs)
  DEALLOCATE(K)
  DEALLOCATE(KInv)
  DEALLOCATE(b)
  DEALLOCATE(d)

!-------------------------------------------------------------------------------
! 4. Differentiation of detrended excess phase
!-------------------------------------------------------------------------------

  ALLOCATE(DSD(nc,npoints)) 
  ALLOCATE(FSD(nc,npoints))
  
  ! 4.1 Determine data step

  dtf = nrf*(time(2)-time(1))
  wf  = CEILING(REAL(w_smooth, wp)/REAL(nrf, wp))

  IF (w_smooth*nd > npoints) THEN
     CALL message(msg_warn, &
          'Filter width larger than data array size - will not process')
     RETURN
   ENDIF
  
  ! 4.2 Filtering

  SELECT CASE(filter)
  CASE('optest')
    CALL ropp_pp_filter(dtf, SD(:,1::nrf), wf, nd, DS=DSD(:,1::nrf))
  CASE('slpoly')
    CALL ropp_pp_sliding_polynomial(time(1::nrf), SD(:,1::nrf), wf, np,   &
                                    FSD(:,1::nrf), DS=DSD(:,1::nrf))
  END SELECT

  ! 4.3 Interpolation

  IF (nrf > 1) THEN
     DO ic=1,nc
        CALL ropp_pp_interpol(time(1::nrf),   time(1:npoints),     &
                              DSD(ic,1::nrf), DSD(ic,1:npoints))
     ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! 5. Calculation of GPS and LEO velocities
!-------------------------------------------------------------------------------

  ALLOCATE(xleo(npoints,3))
  ALLOCATE(xgns(npoints,3))
  ALLOCATE(vleo(npoints,3))
  ALLOCATE(vgns(npoints,3))
  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo, xgns, vgns)
   
!-------------------------------------------------------------------------------
! 6. Calculation of bending angles - geometric optics
!-------------------------------------------------------------------------------

  DO i = 1, npoints

     ! 6.1 Straight line path

     U0(:) = (r_leo(i,:) - r_gns(i,:))/SQRT(SUM((r_leo(i,:) - r_gns(i,:))**2))

     d0 = (c_light - DOT_PRODUCT(vleo(i,:),U0)) /  &
               (c_light - DOT_PRODUCT(vgns(i,:),U0)) - 1.0_wp

     ! 6.2 L1 bending angles

     doppler = d0 - (DST(1,i) + DSD(1,i)) / c_light

     CALL ropp_pp_geometric_optics(xleo(i,:) - r_coc(:), vleo(i,:), &
                                   xgns(i,:) - r_coc(:), vgns(i,:), &
                                   doppler, impact_L1(i), bangle_L1(i))
     
     ! 6.3 L2 bending angles

     doppler = d0 - (DST(2,i) + DSD(2,i)) / c_light
     CALL ropp_pp_geometric_optics(xleo(i,:) - r_coc(:), vleo(i,:), &
                                   xgns(i,:) - r_coc(:), vgns(i,:), &
                                   doppler, impact_L2(i), bangle_L2(i))

   ENDDO

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

   DEALLOCATE(FSD)
   DEALLOCATE(SD)
   DEALLOCATE(ST)
   DEALLOCATE(DST)
   DEALLOCATE(DSD)
   DEALLOCATE(xleo, xgns, vleo, vgns)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_pp_bending_angle_go
