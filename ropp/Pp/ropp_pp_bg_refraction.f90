! $Id: ropp_pp_model_refraction.f90 2369 2009-11-23 11:12:05Z frhl $

SUBROUTINE ropp_pp_bg_refraction(bfile, month, lat, lon, impact,   &
                                 bangle_BG, config)

!****s* ModelRefraction/ropp_pp_bg_refraction *
!
! NAME
!    ropp_pp_bg_refraction - Calculate bending angle profile for BG 
!                            refractivity
!
! SYNOPSIS
!    call ropp_pp_bg_refraction(bfile, month, lat, lon, impact,    &
!                               bangle_BG, config)
! 
! DESCRIPTION
!    This subroutine calculates the bending angle profile for a given location
!    from the background t,p,q fields.
!
! INPUT
!    character(len=*) :: bfile          Model coefficients file
!    integer,         :: month          Month of year
!    real(wp)         :: lat            Latitude  (deg)
!    real(wp)         :: lon            Longitude (deg)
!    real(wp), dim(:) :: impact         Impact parameter (m)
!    type(PPConfig)   :: config         Configuration parameters
!
! OUTPUT
!    real(wp), dim(:) :: bangle_BG      BG bending angles 
!                                       (on input impact parameter levels).
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
! USE ropp_pp, not_this => ropp_pp_bg_refraction
  USE ropp_pp_types, ONLY: PPConfig
  
  IMPLICIT NONE

  CHARACTER(len=*),       INTENT(in)  :: bfile       ! BG file path
  INTEGER,                INTENT(in)  :: month       ! Month of year
  REAL(wp),               INTENT(in)  :: lat         ! Latitude
  REAL(wp),               INTENT(in)  :: lon         ! Longitude
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact      ! Impact parameter (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_BG   ! BG bending angle (rad)
  TYPE(PPConfig),         INTENT(in)  :: config      ! Configuration options

  INTEGER                             :: i           ! Index
  INTEGER                             :: nh          ! No. of hi-res grid points
  REAL(wp)                            :: hmax
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: alt         ! Hi-res altitude scale (m)
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: impact_nh   ! Hi-res impact parameters
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: refrac_nh   ! Hi-res model refractivity
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: grad        ! Hi-res model ref gradient
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: ref_n       ! Scaled refractivity
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: bangle_nh   ! Hi-res model refractivity
  REAL(wp)                            :: scale       ! Asymptotic scale height
  REAL(wp), PARAMETER                 :: dz = 30.0   ! Integration step (m)

!-------------------------------------------------------------------------------
! 2. Define high-resolution vertical grid for BG refractivity
!-------------------------------------------------------------------------------

  hmax = MAX(config%ztop_invert, 40000.0_wp + (MAXVAL(impact - config%r_curve)))
  nh = CEILING(hmax/dz)
  
  ALLOCATE(alt(0:nh))
  ALLOCATE(impact_nh(0:nh))
  ALLOCATE(refrac_nh(0:nh))
  ALLOCATE(grad(0:nh))
  ALLOCATE(bangle_nh(0:nh))
  ALLOCATE(ref_n(0:nh))
  
  DO i=0,nh
     alt(i) = i*hmax/(nh*1.0_wp)
  ENDDO
  
!-------------------------------------------------------------------------------
! 3. Obtain BG refractivity values on high-resolution grid
!-------------------------------------------------------------------------------

  CALL ropp_pp_refrac_BG(bfile, month, lat, lon, alt(0:nh), refrac_nh(0:nh))

!-------------------------------------------------------------------------------
! 4. Compute BG bending angle on output impact parameter levels
!-------------------------------------------------------------------------------

  impact_nh(:)  = (1.0_wp + refrac_nh(:)*1.e-6_wp) * (alt(:) + config%r_curve)
  
  IF ( INDEX(config%abel, "EXP" ) == 1 ) THEN

     CALL ropp_pp_abel_EXP(impact_nh, refrac_nh, impact, bangle_BG)

  ELSE
    
    ref_n(:) = 1.e-6_wp*refrac_nh(:)
    grad(0) = ref_n(0) * LOG(ref_n(1)/ref_n(0))/(alt(1) - alt(0))
    grad(nh) = ref_n(nh) * LOG(ref_n(nh)/ref_n(nh-1))/(alt(nh) - alt(nh-1))
    DO i=1,nh-1
      grad(i) = ref_n(i) * LOG(ref_n(i+1)/ref_n(i-1))/(alt(i+1) - alt(i-1))
    ENDDO
    grad(:) = grad(:) / ((1.0_wp + 1.e-6*refrac_nh(:)) *     &
       (1.0_wp + 1.e-6_wp*refrac_nh(:) + (config%r_curve + alt(:))*grad(:)))
    
    scale = -1.e-6_wp*refrac_nh(nh)/grad(nh) 
    CALL ropp_pp_abel_LIN(impact_nh, refrac_nh, impact_nh, bangle_nh,  &
                          grad, scale)
    CALL ropp_pp_interpol(impact_nh, impact, bangle_nh, bangle_BG)
          
  ENDIF

  DEALLOCATE(alt)
  DEALLOCATE(impact_nh)
  DEALLOCATE(refrac_nh)
  DEALLOCATE(grad)
  DEALLOCATE(ref_n)
  DEALLOCATE(bangle_nh)

END SUBROUTINE ropp_pp_bg_refraction
