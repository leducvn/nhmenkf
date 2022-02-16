! $Id: ropp_pp_invert_refraction.f90 3551 2013-02-25 09:51:28Z idculv $


SUBROUTINE ropp_pp_invert_refraction(mfile, month, lat, lon, impact, bangle,  &
                                     geop, refrac, config)

!****s* IonosphericCorrection/ropp_pp_invert_refraction *
!
! NAME
!    ropp_pp_invert_refraction - Invert neutral atmosphere bending angle
!
! SYNOPSIS
!    call ropp_pp_invert_refraction(mfile, month, lat, lon, impact, bangle,
!                                   geop, refrac, config)
! 
! DESCRIPTION
!    This subroutine calculates the refractivity profile for a given location
!    from a bending angle profile combined with MSIS climatology
!
! INPUTS
!    character(len=*) :: mfile          Model coefficients file
!    integer          :: month          Month of year  
!    real(wp)         :: lat            Latitude  (deg)
!    real(wp)         :: lon            Longitude (deg)
!    real(wp), dim(:) :: impact         Impact parameter (m)
!    real(wp), dim(:) :: bangle         Bending angles (rad)
!    type(ppConfig)   :: config         Configuration parameters
!
! OUTPUT
!    real(wp), dim(:) :: geop           Geopotential height (m)
!    real(wp), dim(:) :: refrac         Refractivity (N-units)
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
  USE geodesy, ONLY: geometric2geopotential
! USE ropp_pp, not_this => ropp_pp_invert_refraction
  USE ropp_pp_types, ONLY: PPConfig

  IMPLICIT NONE

  CHARACTER(len=*),    INTENT(inout)  :: mfile      ! MSIS file path
  INTEGER,                INTENT(in)  :: month      ! Month of year
  REAL(wp),               INTENT(in)  :: lat        ! Latitude
  REAL(wp),               INTENT(in)  :: lon        ! Longitude
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact     ! Impact parameter (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle     ! Bending angle (rad)
  REAL(wp), DIMENSION(:), INTENT(out) :: geop       ! Geopotential height (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: refrac     ! Refractivity 
  TYPE(ppConfig),         INTENT(in)  :: config     ! Configuration options

  REAL(wp), DIMENSION(:), ALLOCATABLE :: impact_nh  ! Hi-res impact param grid
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_nh  ! Hi-res bending angle 
  REAL(wp), DIMENSION(:), ALLOCATABLE :: refrac_nh  ! Hi-res refractivity
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_nh    ! Hi-res geopotential height
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alt_nh     ! Hi-res altitude grid
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_nr  ! Model bending angle
  REAL(wp)                            :: scale      ! Estimated vertical scale

  INTEGER                             :: i          ! Index 
  REAL(wp)                            :: Pmax       ! Maximum impact parameter
  REAL(wp)                            :: rf         ! Regression factor
  INTEGER                             :: n, nh, nr, nu  
 
!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  Pmax = MAX(impact(1), impact(SIZE(impact)))
  n = SIZE(impact)
  nh = MAX(0, CEILING((config%r_curve + config%ztop_invert - Pmax) /   &
                        config%dzh_invert))

  ALLOCATE(impact_nh(n+nh))
  ALLOCATE(bangle_nh(n+nh))
  ALLOCATE(refrac_nh(n+nh))
  ALLOCATE(geop_nh(n+nh))
  ALLOCATE(alt_nh(n+nh))

  impact_nh(1:n) = impact
  bangle_nh(1:n) = bangle

!-------------------------------------------------------------------------------
! 3. Obtain MSIS refractivity values on high-resolution grid
!-------------------------------------------------------------------------------

  IF ( nh > 0 ) THEN

     ! 3.1 Calculate upper part of grid
     
     DO i=1,nh
        impact_nh(n+i) = Pmax + REAL(i,wp) *     &
                         (config%r_curve+config%ztop_invert-Pmax)/(1.0_wp*nh)
     ENDDO
     
     ! 3.2 Calculate model refraction
     CALL ropp_pp_model_refraction(mfile,month,lat,lon,impact_nh(n+1:n+nh),    &
                                   bangle_nh(n+1:n+nh), config) 

     ! 3.3 Regression

     nr = SUM(MINLOC(impact_nh, impact_nh(:) > Pmax-config%dzr_invert))
     nu = SUM(MAXLOC(impact_nh, impact_nh(:) < Pmax-3000.0_wp))

     ALLOCATE(bangle_nr(nr:nu))
     CALL ropp_pp_model_refraction(mfile, month, lat, lon, impact_nh(nr:nu),   &
                                   bangle_nr(nr:nu), config) 
     
     rf = SUM(bangle_nr(nr:nu)*bangle_nh(nr:nu))/SUM(bangle_nr(nr:nu)**2)

     bangle_nh(n+ 1:n+nh) = rf*bangle_nh(n+1:n+nh)
     
     DEALLOCATE(bangle_nr)
     
  ENDIF

!-------------------------------------------------------------------------------
! 4. Compute refractivity on output impact parameter levels
!-------------------------------------------------------------------------------
  
  IF ( INDEX(config%abel, "EXP" ) == 1 ) THEN

     CALL ropp_pp_invert_EXP(impact_nh,bangle_nh,impact_nh(1:n),refrac_nh(1:n))

  ELSE
     
     scale = -(impact_nh(n+nh) - impact_nh(3*(n+nh)/4)) /  &
          (LOG(bangle_nh(n+nh)) - LOG(bangle_nh(3*(n+nh)/4)))
     CALL ropp_pp_invert_LIN(impact_nh, bangle_nh, impact_nh(1:n),    &
                              refrac_nh(1:n), scale)

  ENDIF

  alt_nh = (impact_nh / (1.0_wp + refrac_nh * 1.e-6_wp)) -  config%r_curve
  geop_nh = geometric2geopotential(lat, alt_nh)

  refrac(1:n) = refrac_nh(1:n)
  geop(1:n)   = geop_nh(1:n)

!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(impact_nh)
  DEALLOCATE(bangle_nh)
  DEALLOCATE(refrac_nh)
  DEALLOCATE(geop_nh)
  DEALLOCATE(alt_nh)

END SUBROUTINE ropp_pp_invert_refraction
