! $Id: ropp_pp_modelphase.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_modelphase(month, lat, lon, time, r_leo, r_gns,    &
                              r_coc, roc, phase_LM, impact_LM, config) 

!****s* Preprocessing/ropp_pp_modelphase *
!
! NAME
!    ropp_pp_modelphase - Compute model excess phase from MSIS data
!                   
! SYNOPSIS
!    call ropp_pp_modelphase(month, lat, lon, time, r_leo, r_gns,
!                            r_coc, roc, phase_LM, impact_LM, config)
! 
! DESCRIPTION
!    This routine computes model excess phase
!      1. Compute MSIS model bending angle as a function of impact parameter p
!      2. Compute t(p) for given satellite trajectories
!      3. Compute p(t) -> d(t) -> phase(t)
!
! INPUTS
!    integer                  :: month     ! month of year
!    real(wp)                 :: lat       ! occultation point latitude
!    real(wp)                 :: lon       ! occultation point longitude
!    real(wp), dimension(:)   :: time      ! time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! cartesian LEO coordinates (ECF)
!    real(wp), dimension(:,:) :: r_gns     ! cartesian GPS coordinates (ECF)
!    real(wp), dimension(:)   :: r_coc     ! cartesian centre curvature (ECF)
!    real(wp)                 :: roc       ! radius curvature
!    type(PPConfig)           :: config    ! configuration options
!
! OUTPUT
!    real(wp), dimension(:)   :: phase_LM  ! model excess phase (m)
!    real(wp), dimension(:)   :: impact_LM ! model impact parameter (m)
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
! USE ropp_pp, not_this => ropp_pp_modelphase
  USE ropp_pp_types, ONLY: PPConfig
  USE ropp_pp_constants, ONLY: R_dry, R_vap, epsilon_water, kappa1, kappa2
  USE ropp_utils, ONLY: impact_parameter, ropp_MDFV

  IMPLICIT NONE

  INTEGER,                  INTENT(in)  :: month     ! Month of year
  REAL(wp),                 INTENT(in)  :: lat       ! Latitude
  REAL(wp),                 INTENT(in)  :: lon       ! Longitude
  REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates (ECF)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates (ECF)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! Centre curvature (ECF)
  REAL(wp),                 INTENT(in)  :: roc       ! Radius curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(out) :: phase_LM  ! Model excess phase (m)
  REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: impact_LM ! Model impact (m)
  TYPE(PPConfig),         INTENT(inout) :: config    ! Configuration options

  INTEGER, PARAMETER :: nh = 1200   ! Number of high-resolution grid points
  INTEGER, PARAMETER :: nr = 10     ! Polynomial degree for phase regression
  INTEGER            :: n, nl       ! Number of data
  INTEGER            :: i           ! Data index
  REAL(wp)           :: ps1         ! First straight-line impact parameter
  REAL(wp)           :: psN         ! Last straight-line impact parameter
  REAL(wp)           :: pmin        ! Minimum impact parameter
  REAL(wp)           :: bmax        ! Maximum model bending angle from trajectory
  REAL(wp)           :: msis_max    ! Maximum model bending angle from MSIS
  REAL(wp)           :: zmax        ! Maximum height
  REAL(wp)           :: scale       ! Vertical scale of refractivity
  INTEGER            :: i1          ! Lower index of valid model interval
  INTEGER            :: i2          ! Upper index of valid model interval
  REAL(wp)           :: phase_min   ! Excess phase base
  REAL(wp)           :: doppler     ! Relative Doppler frequency shift
  REAL(wp)           :: dpl         ! Extended grid resolution

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: dph_mod  ! Interpolated dph_dt

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo   ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns   ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo   ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns   ! GPS velocity by regression

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: alt    ! Altitude grid
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: refrac ! Refractive index  
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: grad   ! d(ln n)/dz
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: impact ! Grid of x = n*r
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: bangle ! Bending angle on hi-res grid
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: d2ba   ! 2nd derivative bangle
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: impact_nh ! Extended grid of x = n*r
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: bangle_nh ! Extended bending angle
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: shum   ! Specific humidity
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: temp   ! Temperature
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: pres   ! Pressure
  REAL(wp)                              :: pwvp   ! Partial water vapor pressure

  REAL(wp), PARAMETER :: DPmin = 100.0_wp         ! Step for locating Pmin
  REAL(wp), PARAMETER :: Qtrop = 90.0_wp          ! 90% tropospheric humidity  
  REAL(wp), PARAMETER :: aq =  R_dry/R_vap
  REAL(wp), PARAMETER :: bq =  1.0_wp - aq

!-------------------------------------------------------------------------------
! 2. Data allocation
!-------------------------------------------------------------------------------

  n = SIZE(time)

  ALLOCATE(dph_mod(n))

  ALLOCATE(alt(0:nh))
  ALLOCATE(refrac(0:nh))
  ALLOCATE(impact(0:nh))
  ALLOCATE(grad(0:nh))
  ALLOCATE(bangle(0:nh))
  ALLOCATE(d2ba(0:nh))
  ALLOCATE(shum(0:nh))
  ALLOCATE(temp(0:nh))
  ALLOCATE(pres(0:nh))

  ALLOCATE(xleo(n,3))
  ALLOCATE(xgns(n,3))
  ALLOCATE(vleo(n,3))
  ALLOCATE(vgns(n,3))

!-------------------------------------------------------------------------------
! 3. Calculate bending angles from MSIS data
!-------------------------------------------------------------------------------

  ! 3.1 Impact parameter limits - compute straight-line impact parameters

  ps1 = impact_parameter(r_leo(1,:) - r_coc(:), r_gns(1,:) - r_coc)
  psN = impact_parameter(r_leo(n,:) - r_coc(:), r_gns(n,:) - r_coc)

  ! 3.2 MSIS refractivity profile

  zmax = 20000.0_wp + MAX(120000.0_wp, MAX(ps1, psN) - roc)

  DO i=0,nh
     alt(i) = (i*zmax)/(1.0_wp*nh)
  ENDDO
  
  CALL ropp_pp_refrac_MSIS(config%mfile, month, lat, lon, alt, refrac, grad)

  ! 3.3 Compute dln(n)/dx

  grad(:) = grad(:)*1.e-6_wp
  grad(:) = grad(:)/((1.0_wp + 1.e-6_wp*refrac(:)) *  &
                        (1.0_wp + 1.e-6_wp*refrac(:) + (roc + alt(:))*grad(:)))

  ! 3.4 Compute impact parameter

  DO i=0,nh
     impact(i) = (roc + alt(i)) * (1.0_wp + 1e-6_wp*refrac(i))
  ENDDO

  ! 3.5 Adjust refractivity profile 
  !     (add tropospheric humidity to MSIS climatology)

  ! 3.5.1 Compute dry temperature and pressure
  shum(:) = 0.0_wp
  CALL ropp_pp_tdry(lat, alt, refrac, shum, temp, pres)

  ! 3.5.2 Add humidity and recompute refractivity ('MSIS+Q')
  DO i=0,nh
     IF(alt(i) < 15000.0_wp)THEN
        pwvp = Qtrop * 6.11_wp *      &
                EXP(17.27_wp*(temp(i)-273.16_wp)/(237.3_wp+(temp(i)-273.16_wp)))
        shum(i) = pwvp * aq / (pres(i) - pwvp*bq)
     ELSE
        shum(i) = 0.0_wp
     ENDIF
     
     pwvp = pres(i) * shum(i) / (epsilon_water + (1.0_wp-epsilon_water)*shum(i))
     refrac(i) = kappa1*pres(i) / temp(i) + kappa2* pwvp / temp(i)**2

  ENDDO
  
  ! 3.5.3 Re-compute impact parameter
  DO i=0,nh
     impact(i) = (roc + alt(i)) * (1.0_wp + 1e-6_wp*refrac(i))
  ENDDO

  ! 3.5.4 Re-compute profile gradient
  grad(0) = 1.e-6_wp*refrac(0) * LOG(refrac(1)/refrac(0))/(alt(1) - alt(0))
  grad(nh) = 1.e-6_wp*refrac(nh) * LOG(refrac(nh)/refrac(nh-1))/(alt(nh) - alt(nh-1)) 
  DO i=1,nh-1
     grad(i) = 1.e-6_wp*refrac(i) * LOG(refrac(i+1)/refrac(i-1))/(alt(i+1) - alt(i-1))
  ENDDO

  grad(:) = grad(:) / ((1.0_wp + 1.e-6*refrac(:)) *     &
                       (1.0_wp + 1.e-6_wp*refrac(:) + (roc + alt(:))*grad(:)))

  ! 3.6 Abel transform to obtain MSIS bending angle profile

  IF ( INDEX(config%abel, "EXP" ) == 1 ) THEN
     
     CALL ropp_pp_abel_EXP(impact(0:), refrac(0:), impact(0:), bangle(0:))
     
  ELSE
    
     scale = -1.e-6_wp*refrac(nh)/grad(nh)
     CALL ropp_pp_abel_LIN(impact(0:), refrac(0:), impact(0:), bangle(0:),  &
                            grad(0:), scale)

  ENDIF

!-------------------------------------------------------------------------------
! 4. Extend MSIS+Q bending angle profile
!-------------------------------------------------------------------------------

  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo, xgns, vgns)

  IF ( psN < ps1 ) THEN
    i = n
  ELSE
    i = 1
  ENDIF

  Pmin = impact(0)
  MSIS_max = bangle(0)

  ! 4.1 Find Pmin
  
  Locate_Pmin: DO 

    CALL ropp_pp_impact2doppler(xleo(i,:) - r_coc, vleo(i,:),    &
                                xgns(i,:) - r_coc, vgns(i,:),    &
                                Pmin, doppler, Bmax)
  
    IF ( Bmax < MSIS_max ) THEN
      EXIT Locate_Pmin
    
    ELSE
      Pmin = Pmin - DPmin
      CALL ropp_pp_interpol(impact, Pmin, LOG(bangle), MSIS_max)
      MSIS_max = EXP(MSIS_max)
    ENDIF
    
  ENDDO Locate_Pmin
    
  Pmin = Pmin - 5*(Zmax / nh)
  dpl = impact(1) - impact(0)

  nl = MAX(0, FLOOR((impact(0) - Pmin)/dpl))
  
  ! 4.2 Allocate arrays for extended grid

  ALLOCATE(impact_nh(-nl:nh))
  ALLOCATE(bangle_nh(-nl:nh))

  ! 4.3 Extend MSIS bending angle profile

  impact_nh(0:nh) = impact(0:nh)
  bangle_nh(0:nh) = bangle(0:nh)
  
  CALL ropp_pp_init_spline(impact(:), LOG(bangle(:)), d2ba(:))
  
  DO i=-nl,-1
    
    impact_nh(i) = impact(0) + REAL(i, wp) * dpl 
    CALL ropp_pp_interpol_spline(impact, LOG(bangle), d2ba,   &
                                 impact_nh(i), bangle_nh(i))
    bangle_nh(i) = EXP(bangle_nh(i))

  ENDDO
  
!-------------------------------------------------------------------------------
! 5. Calculate model phase from MSIS+Q bending angles
!-------------------------------------------------------------------------------

  CALL ropp_pp_bangle2phase(time(1:n), r_leo, r_gns, r_coc, impact_nh(-nl:nh), &
                            bangle_nh(-nl:nh), phase_LM(1:n), dph_mod(1:n),    &
                            impact_LM(1:n), i1, i2)

!-------------------------------------------------------------------------------
! 6. Integrate merged dPh/dt to compute excess phase
!-------------------------------------------------------------------------------

  phase_LM(:) = ropp_MDFV
  phase_LM(i1) = 0.0_wp
  DO i=i1+1,i2
     phase_LM(i) = phase_LM(i-1) +      &
                     (dph_mod(i)+dph_mod(i-1))*(time(i)-time(i-1))/2.0_wp
  END DO

  phase_min = MINVAL(phase_LM(i1:i2))
  phase_LM(i1:i2) = phase_LM(i1:i2) - phase_min

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(dph_mod)

  DEALLOCATE(xleo)
  DEALLOCATE(xgns)
  DEALLOCATE(vleo)
  DEALLOCATE(vgns)  

  DEALLOCATE(alt)
  DEALLOCATE(refrac)
  DEALLOCATE(impact)
  DEALLOCATE(grad)
  DEALLOCATE(bangle)
  DEALLOCATE(d2ba)
  DEALLOCATE(shum)
  DEALLOCATE(temp)
  DEALLOCATE(pres)

  DEALLOCATE(impact_nh)
  DEALLOCATE(bangle_nh)

END SUBROUTINE ropp_pp_modelphase
