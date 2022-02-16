! $Id: ropp_pp_model_refraction.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_model_refraction(mfile, month, lat, lon, impact,   &
                                    bangle_MSIS, config)

!****s* ModelRefraction/ropp_pp_model_refraction *
!
! NAME
!    ropp_pp_model_refraction - Calculate bending angle profile for MSIS 
!                               refractivity
!
! SYNOPSIS
!    call ropp_pp_model_refraction(mfile, month, lat, lon, impact,    &
!                                  bangle_MSIS, config)
! 
! DESCRIPTION
!    This subroutine calculates the bending angle profile for a given location
!    from the MSIS refractivity field
!
! INPUT
!    character(len=*) :: mfile          Model coefficients file
!    integer,         :: month          Month of year
!    real(wp)         :: lat            Latitude  (deg)
!    real(wp)         :: lon            Longitude (deg)
!    real(wp), dim(:) :: impact         Impact parameter (m)
!    type(PPConfig)   :: config         Configuration parameters
!
! OUTPUT
!    real(wp), dim(:) :: bangle_MSIS    MSIS bending angles 
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
! USE ropp_pp, not_this => ropp_pp_model_refraction
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPConfig
  
  IMPLICIT NONE

  CHARACTER(len=*),    INTENT(inout)  :: mfile       ! MSIS file path
  INTEGER,                INTENT(in)  :: month       ! Month of year
  REAL(wp),               INTENT(in)  :: lat         ! Latitude
  REAL(wp),               INTENT(in)  :: lon         ! Longitude
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact      ! Impact parameter (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_MSIS ! MSIS bending angle (rad)
  TYPE(PPConfig),         INTENT(in)  :: config      ! Configuration options

  INTEGER                             :: i           ! Index
  INTEGER                             :: nh          ! No. of hi-res grid points
  REAL(wp)                            :: hmax
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: alt         ! Hi-res altitude scale (m)
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: impact_nh   ! Hi-res impact parameters
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: refrac_nh   ! Hi-res model refractivity
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: grad_refrac ! Hi-res model ref gradient
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: bangle_nh   ! Hi-res model refractivity
  REAL(wp)                            :: scale       ! Asymptotic scale height
  REAL(wp), PARAMETER                 :: dz = 30.0   ! Integration step (m)

!-------------------------------------------------------------------------------
! 2. Define high-resolution vertical grid for MSIS refractivity
!-------------------------------------------------------------------------------

  hmax = MAX(config%ztop_invert, 40000.0_wp + (MAXVAL(impact - config%r_curve)))
  nh = CEILING(hmax/dz)
  
  ALLOCATE(alt(0:nh))
  ALLOCATE(impact_nh(0:nh))
  ALLOCATE(refrac_nh(0:nh))
  ALLOCATE(grad_refrac(0:nh))
  ALLOCATE(bangle_nh(0:nh))
  
  DO i=0,nh
     alt(i) = i*hmax/(nh*1.0_wp)
  ENDDO
  
!-------------------------------------------------------------------------------
! 3. Obtain MSIS refractivity values on high-resolution grid
!-------------------------------------------------------------------------------

  CALL ropp_pp_refrac_MSIS(mfile, month, lat, lon, alt, refrac_nh, grad_refrac)

  ! 3.1 dln(n)/dx

  grad_refrac(:) = grad_refrac(:)*1.e-6_wp
  DO i=0,nh
     grad_refrac(i) = grad_refrac(i)/((1.0_wp + 1.e-6_wp*refrac_nh(i)) *   &
                        (1.0_wp + 1.e-6_wp*refrac_nh(i) +                  &
                        (config%r_curve + alt(i))*grad_refrac(i)))
  ENDDO
  
!-------------------------------------------------------------------------------
! 4. Compute MSIS bending angle on output impact parameter levels
!-------------------------------------------------------------------------------

  impact_nh(:)  = (1.0_wp + refrac_nh(:)*1.e-6_wp) * (alt(:) + config%r_curve)
  
  IF ( INDEX(config%abel, "EXP" ) == 1 ) THEN

     CALL ropp_pp_abel_EXP(impact_nh, refrac_nh, impact, bangle_MSIS)

  ELSE
    
    scale = -1.e-6_wp*refrac_nh(nh)/grad_refrac(nh) 
    CALL ropp_pp_abel_LIN(impact_nh, refrac_nh, impact_nh, bangle_nh,  &
                          grad_refrac, scale)
    CALL ropp_pp_interpol(impact_nh, impact, bangle_nh, bangle_MSIS)
          
  ENDIF

  DEALLOCATE(alt)
  DEALLOCATE(impact_nh)
  DEALLOCATE(refrac_nh)
  DEALLOCATE(grad_refrac)
  DEALLOCATE(bangle_nh)
  MSIS_read = .false.

END SUBROUTINE ropp_pp_model_refraction
