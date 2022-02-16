! $Id: ropp_pp_interpolate_trajectory.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_interpolate_trajectory(time, cleo, cgns, r_coc, t_init,  &
                                          xleo, vleo, xgns, vgns, theta)

!****s* Preprocessing/ropp_pp_interpolate_trajectory *
!
! NAME
!    ropp_pp_interpolate_trajectory - Calculate interpolated coordinates,
!                                     velocities, satellite-satellite angle
!                   
! SYNOPSIS
!     call ropp_pp_interpolate_trajectory(time, cleo, cgns, r_coc, t_init,  & 
!                                         xleo, vleo, xgns, vgns, theta)
! 
! DESCRIPTION
!    This routine calculates interpolated  coordinates, velocities
!    and satellite-satellite angle at time t_init by polynomial regression
!
! INPUTS
!    real(wp), dimension(:)   :: time    ! time of sample (s)
!    real(wp), dimension(:,:) :: cleo    ! LEO regression coefficients
!    real(wp), dimension(:,:) :: cgns    ! GPS regression coefficients
!    real(wp), dimension(:,:) :: r_coc   ! cartesian centre curvature
!    real(wp)                 :: t_init  ! interpolation time
!
! OUTPUTS
!    real(wp), dimension(:,:) :: xleo    ! LEO positions from regression
!    real(wp), dimension(:,:) :: vleo    ! LEO velocities from regression
!    real(wp), dimension(:,:) :: xgns    ! GPS positions from regression
!    real(wp), dimension(:,:) :: vgns    ! GPS velocities from regression
!    real(wp), optional       :: theta   ! Satellite-to-satellite angle
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
  USE ropp_pp, ONLY: ropp_pp_polynomial
  USE ropp_utils
  
  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:),   INTENT(in)  :: time   ! Relative time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: cleo   ! LEO regression coefficients
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: cgns   ! GPS regression coefficients
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc  ! Centre of curvature
  REAL(wp),                 INTENT(in)  :: t_init ! Interpolation time
  REAL(wp), DIMENSION(:),   INTENT(out) :: xleo   ! Cartesian LEO coordinates
  REAL(wp), DIMENSION(:),   INTENT(out) :: vleo   ! Cartesian LEO velocity
  REAL(wp), DIMENSION(:),   INTENT(out) :: xgns   ! Cartesian GPS coordinates
  REAL(wp), DIMENSION(:),   INTENT(out) :: vgns   ! Cartesian GPS velocity
  REAL(wp), OPTIONAL,       INTENT(out) :: theta  ! Satellite-to-satellite angle

  REAL(wp)                              :: t_norm ! Normalised time  
  INTEGER                               :: n, m   ! Indices      
  
!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------
  
  n = SIZE(time)
  
  t_norm = (t_init - time(1))/(time(n) - time(1))
  
!-------------------------------------------------------------------------------
! 3. Polynomial expression for position, and 1st derivative
!-------------------------------------------------------------------------------

  DO m=1,3
     CALL ropp_pp_polynomial(cleo(:,m), t_norm, xleo(m), vleo(m))
     CALL ropp_pp_polynomial(cgns(:,m), t_norm, xgns(m), vgns(m))
  ENDDO

!-------------------------------------------------------------------------------
! 4. Re-normalise velocity
!-------------------------------------------------------------------------------

  vleo = vleo / (time(n)-time(1))    
  vgns = vgns / (time(n)-time(1))
  
!-------------------------------------------------------------------------------
! 5. Compute satellite-to-satellite angle
!-------------------------------------------------------------------------------

  IF (PRESENT(theta)) THEN
     theta = vector_angle(xgns - r_coc, xleo - r_coc)
  ENDIF
  
END SUBROUTINE ropp_pp_interpolate_trajectory
