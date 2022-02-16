! $Id: ropp_pp_amplitude_go.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_amplitude_go(time, r_leo, r_gns, r_coc, roc, impact,   &
                                snr, w_smooth, snr_R)

!****s* Preprocess/ropp_pp_amplitude_go *
!
! NAME
!    ropp_pp_amplitude_go - Calculate RO signal amplitude in GO approximation.
!                           
! SYNOPSIS
!    call ropp_pp_amplitude_go(time, r_leo, r_gns, r_coc, roc, impact, snr, 
!                              w_smooth, snr_R)
! 
! DESCRIPTION
!    This routine calculates L2 amplitude in Geometric optics approximation
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! Relative time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! LEO coordinates (ECI or ECF) (m)
!    real(wp), dimension(:,:) :: r_gns     ! GPS coordinates (ECI or ECF) (m)
!    real(wp), dimension(:)   :: r_coc     ! centre curvature coordinates (m)
!    real(wp)                 :: roc       ! Radius curvature (m)
!    real(wp), dimension(:)   :: impact    ! Impact parameters (m)
!    real(wp), dimension(:)   :: snr       ! Observed amplitude
!    integer                  :: w_smooth  ! smoothing window (points)
!
! OUTPUT
!    real(wp), dimension(:)   :: snr_R     ! Refractive amplitude
! 
! NOTES
!
! REFERENCES
!   Gorbunov M.E. and Lauritsen K.B. 2004
!   Analysis of wave fields by Fourier integral operators and their application
!   for radio occultations
!   Radio Science (39) RS4010
!
!   Gorbunov M.E., Lauritsen K.B., Rodin A., Tomassini M., Kornblueh L. 2005 
!   Analysis of the CHAMP experimental data on radio-occultation sounding of
!   the Earth's atmosphere.
!   Izvestiya Atmospheric and Oceanic Physics (41) 726-740.
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
! USE ropp_pp, not_this => ropp_pp_amplitude_go
  USE ropp_utils, ONLY: vector_angle, ropp_MDFV

  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! Cartesian LEO coordinates
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! Cartesian GPS coordinates
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! Centre curvature coords
  REAL(wp),                 INTENT(in)  :: roc       ! Radius curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: impact    ! Impact parameters (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: snr       ! Observed amplitude
  INTEGER,                  INTENT(in)  :: w_smooth  ! Smoothing window (points)
  REAL(wp), DIMENSION(:),   INTENT(out) :: snr_R     ! Refractive amplitude

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo   ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns   ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo   ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns   ! GPS velocity by regression

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xrleo  ! LEO coordinates from r_coc
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xrgns  ! GPS coordinates from r_coc
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: rxleo  ! LEO radii from r_coc
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: rxgns  ! GPS radii from r_coc
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: rvleo  ! LEO radial velocity
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: rvgns  ! GPS radial velocity
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: theta  ! GPS-LEO angular separation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: vtheta ! Angular velocity
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ps     ! Smoothed impact parameter
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: pv     ! d(Impact parameter)/dt
  
  REAL(wp)                              :: dgns   ! GNS-limb distances 
  REAL(wp)                              :: dleo   ! LEO-limb distances 
  REAL(wp)                              :: arn    ! Normalizing constant
  INTEGER, PARAMETER                    :: np = 3 ! Polynomial degree for filter
  INTEGER                               :: npoints, i, nrf, wr, nv
  LOGICAL,  DIMENSION(:),   ALLOCATABLE :: m      ! Mask for array operations

  REAL(wp), PARAMETER  :: Hvac = 30000.0_wp       ! Height for computation of 
                                                  ! vacuum amplitude

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------
    
  npoints = SIZE(time)
  
  ALLOCATE(xleo(npoints,3))
  ALLOCATE(xgns(npoints,3))
  ALLOCATE(vleo(npoints,3))
  ALLOCATE(vgns(npoints,3))
  
  ALLOCATE(xrleo(npoints,3))
  ALLOCATE(xrgns(npoints,3))
  ALLOCATE(rxleo(npoints))
  ALLOCATE(rxgns(npoints))
  ALLOCATE(rvleo(npoints))
  ALLOCATE(rvgns(npoints))

  ALLOCATE(theta(npoints))
  ALLOCATE(vtheta(npoints))

  ALLOCATE(ps(npoints))
  ALLOCATE(pv(npoints))

!-------------------------------------------------------------------------------
! 3. Compute coordinates and velocities
!-------------------------------------------------------------------------------

  ! 3.1 Cartesian coordinates and velocities

  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo, xgns, vgns)

  DO i=1,npoints

  ! 3.2 Radii and radial velocities

     xrleo(i,:) = xleo(i,:) - r_coc(:)
     xrgns(i,:) = xgns(i,:) - r_coc(:)
     
     rxleo(i) = SQRT(SUM(xrleo(i,:)**2))
     rxgns(i) = SQRT(SUM(xrgns(i,:)**2))

     rvleo(i) = DOT_PRODUCT(vleo(i,:), (xrleo(i,:)/SQRT(SUM(xrleo(i,:)**2))))
     rvgns(i) = DOT_PRODUCT(vgns(i,:), (xrgns(i,:)/SQRT(SUM(xrgns(i,:)**2))))

  ! 3.3 Angle and angular velocity

     theta(i)  = vector_angle(xrgns(i,:), xrleo(i,:))
     vtheta(i) = (-1.0_wp/(SQRT((rxgns(i)*rxleo(i))**2 -             &
                    DOT_PRODUCT(xrgns(i,:),xrleo(i,:))**2)))*        &
                     (DOT_PRODUCT(vgns(i,:),xrleo(i,:)) +            &
                     DOT_PRODUCT(xrgns(i,:),vleo(i,:)) -             &
                     (DOT_PRODUCT(xrleo(i,:),xrgns(i,:))*            &
                     (rvgns(i)/rxgns(i) + rvleo(i)/rxleo(i))))

  ENDDO

  ! 3.4 Impact parameter rate of change

  nrf = CEILING(npoints/8000.0_wp)
  wr  = CEILING(REAL(w_smooth, wp)/REAL(nrf, wp)) 

  CALL ropp_pp_sliding_polynomial(time(1::nrf), impact(1::nrf), wr, np,     &
                                  ps(1::nrf), pv(1::nrf))

  IF ( nrf > 1 ) THEN
    CALL ropp_pp_interpol(time(1::nrf), time, ps(1::nrf), ps)
    CALL ropp_pp_interpol(time(1::nrf), time, pv(1::nrf), pv)
  ENDIF

!-------------------------------------------------------------------------------
! 4. Compute refractive amplitude ("a2")
!-------------------------------------------------------------------------------

  DO i=1,npoints
     dgns = SQRT(rxgns(i)**2 - ps(i)**2)
     dleo = SQRT(rxleo(i)**2 - ps(i)**2)
     
     snr_R(i) = SQRT(ABS( (1.0_wp/(dgns*dleo))*                        &
              (1.0_wp - (ps(i)/(rxgns(i)*dgns))*rvgns(i)/vtheta(i)     &
                  - (ps(i)/(rxleo(i)*dleo))*rvleo(i)/vtheta(i))*        &
                 pv(i)/vtheta(i)))

  ENDDO

  snr_R(:) = snr_R(:)*SQRT(ps(:)/(rxgns(:)*rxleo(:)*SIN(theta(:))))

  ALLOCATE(M(npoints))
  M(:) = (PS(:)-roc > Hvac) .AND. (PS(:)-roc < Hvac + 5000.0)
  nv = COUNT(Mask = M(:))

  IF (nv == 0) THEN
    snr_R(:) = ropp_MDFV
  ELSE
    arn = SUM(snr(:)/snr_R(:), Mask = M(:)) / nv 
    snr_R(:) = snr_R(:)*arn
  ENDIF
  
!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------
  DEALLOCATE(m)
  DEALLOCATE(pv)
  DEALLOCATE(ps)
  
  DEALLOCATE(vtheta)
  DEALLOCATE(theta)

  DEALLOCATE(rvgns)
  DEALLOCATE(rvleo)
  DEALLOCATE(rxgns)
  DEALLOCATE(rxleo)
  DEALLOCATE(xrgns)
  DEALLOCATE(xrleo)

  DEALLOCATE(vgns)
  DEALLOCATE(vleo)
  DEALLOCATE(xgns)
  DEALLOCATE(xleo)

END SUBROUTINE ropp_pp_amplitude_go


