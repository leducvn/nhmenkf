! $Id: ropp_pp_impact2doppler.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_impact2doppler(xleo,vleo,xgns,vgns,impact,doppler,bangle)

!****s* Preprocessing/ropp_pp_impact2doppler *
!
! NAME
!    ropp_pp_impact2doppler - Compute relative doppler shift from geometrical
!                             data
!                   
! SYNOPSIS
!    call ropp_pp_impact2doppler(xleo,vleo,xgns,vgns,impact,doppler,bangle)
! 
! DESCRIPTION
!
! INPUTS
!    real(wp), dimension(:) :: xleo     ! cartesian LEO coordinates
!    real(wp), dimension(:) :: vleo     ! cartesian LEO velocity
!    real(wp), dimension(:) :: xgns     ! cartesian GPS coordinates
!    real(wp), dimension(:) :: vgns     ! cartesian GPS velocity
!    real(wp)               :: impact   ! impact parameter
!
! OUTPUT
!    real(wp)               :: doppler  ! Relative doppler frequency shift
!    real(wp)               :: bangle   ! Bending angle
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
! USE ropp_pp, not_this => ropp_pp_impact2doppler
  USE ropp_pp_constants, ONLY: c_light
  USE ropp_utils, ONLY: rotate, vector_angle, vector_product
  
  IMPLICIT NONE

    REAL(wp), DIMENSION(3), INTENT(in)  :: xleo    ! Cartesian LEO coordinates
    REAL(wp), DIMENSION(3), INTENT(in)  :: vleo    ! Cartesian LEO velocity
    REAL(wp), DIMENSION(3), INTENT(in)  :: xgns    ! Cartesian GPS coordinates
    REAL(wp), DIMENSION(3), INTENT(in)  :: vgns    ! Cartesian GPS velocity
    REAL(wp),               INTENT(in)  :: impact  ! Impact parameter
    REAL(wp),               INTENT(out) :: doppler ! Doppler frequency shift
    REAL(wp),               INTENT(out) :: bangle  ! Bending angle

    REAL(wp)                            :: alpha   ! Ray elevation angle
    REAL(wp), DIMENSION(3)              :: Uleo    ! Ray direction at LEO
    REAL(wp), DIMENSION(3)              :: Ugns    ! Ray direction at GPS
    REAL(wp), DIMENSION(3)              :: work    ! Temporary storage

!-------------------------------------------------------------------------------
! 2. Calculate ray directions
!-------------------------------------------------------------------------------
    
    alpha = ASIN(impact/SQRT(SUM(xgns(:)**2)))
    work = rotate(-xgns, vector_product(xleo, xgns), alpha)
    Ugns(:) = work(:) / SQRT(SUM(work(:)**2)) 
    
    alpha = ASIN(impact/SQRT(SUM(xleo(:)**2)))
    work = rotate(-xleo, vector_product(xgns, xleo), alpha)
    Uleo(:) = (-1.0_wp) * work(:) / SQRT(SUM(work(:)**2)) 
    
!-------------------------------------------------------------------------------
! 3. Calculate bending angle and relative doppler shift
!-------------------------------------------------------------------------------

    doppler = (DOT_PRODUCT(vgns,Ugns) - DOT_PRODUCT(vleo,Uleo)) /  &
                    (c_light - DOT_PRODUCT(vgns, Ugns))
    
    bangle = vector_angle(Ugns, Uleo, vector_product(xgns,xleo))

  END SUBROUTINE ropp_pp_impact2doppler
