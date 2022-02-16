! $Id: ropp_pp_merge_profile.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_merge_profile(impact_L1, bangle_L1, impact_L2, bangle_L2,  &
                                 impact_I1, bangle_I1, impact_I2, bangle_I2,  &
                                 Pmin, Pmax)

!****s* IonosphericCorrection/ropp_pp_merge_profile *
!
! NAME
!    ropp_pp_merge_profile - 
!                   Merge and interpolate bending angles onto standard grid
!                   
! SYNOPSIS
!    call ropp_pp_merge_profile(impact_L1, bangle_L1, impact_L2, bangle_L2, 
!                               impact_I1, bangle_I1, impact_I2, bangle_I2,
!                               Pmin, Pmax)                               
! 
! DESCRIPTION
!    This routine calculates L1 and L2 bending angles on a standard vertical 
!    grid by linear interpolation
!
! INPUTS
!    real(wp), dimension(:) :: impact_L1   ! Impact parameters of channel L1
!    real(wp), dimension(:) :: bangle_L1   ! Bending angles for channel L1
!    real(wp), dimension(:) :: impact_L2   ! Impact parameters of channel L2
!    real(wp), dimension(:) :: bangle_L2   ! Bending angles for channel L2
!    real(wp), optional     :: Pmin        ! Minimum impact parameter
!    real(wp), optional     :: Pmax        ! Maximum impact parameter
!
! OUTPUT
!    real(wp), dimension(:) :: impact_I1   ! Standard impact parameter grid L1
!    real(wp), dimension(:) :: bangle_I1   ! L1 bending angles on impact_I1
!    real(wp), dimension(:) :: impact_I2   ! Standard impact parameter grid L1
!    real(wp), dimension(:) :: bangle_I2   ! L2 bending angles on impact_I2
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
! USE ropp_pp, not_this => ropp_pp_merge_profile
  USE ropp_pp

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L1  ! L1 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L1  ! L1 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L2  ! L2 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L2  ! L2 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(out) :: impact_I1  ! Standard IP grid (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_I1  ! L1 bangle on impact_I1
  REAL(wp), DIMENSION(:), INTENT(out) :: impact_I2  ! Standard IP grid (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_I2  ! L2 bangle on impact_I2
  REAL(wp), OPTIONAL, INTENT(in)      :: Pmin       ! Minimum impact parameter
  REAL(wp), OPTIONAL, INTENT(in)      :: Pmax       ! Maximum impact parameter

  REAL(wp), PARAMETER                 :: dhs = 5000.0_wp  ! Safety border (m)

  REAL(wp) :: Pmax1
  REAL(wp) :: Pmin1
  INTEGER  :: n_out
  INTEGER  :: i

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_out = SIZE(impact_I1)
  
  IF(PRESENT(Pmin))THEN 
     Pmin1 = Pmin
  ELSE
     Pmin1 = MINVAL(impact_L1(:))
  ENDIF

  IF(PRESENT(Pmax))THEN
     Pmax1 = Pmax
  ELSE
     Pmax1 = MAXVAL(impact_L1(:)) - dhs
  ENDIF
  
!-------------------------------------------------------------------------------
! 3. Determine impact parameter grids - monotonously increasing
!-------------------------------------------------------------------------------
  
  DO i = 1, n_out
     impact_I1(i) = Pmin1 + (i-1.0_wp)*(Pmax1 - Pmin1)/(n_out - 1.0_wp)
     impact_I2(i) = Pmin1 + (i-1.0_wp)*(Pmax1 - Pmin1)/(n_out - 1.0_wp)
  ENDDO
  
!-------------------------------------------------------------------------------
! 4. Interpolation
!-------------------------------------------------------------------------------

  CALL ropp_pp_interpol(impact_L1, impact_I1, bangle_L1, bangle_I1,Cext=.true.)
  CALL ropp_pp_interpol(impact_L2, impact_I2, bangle_L2, bangle_I2,Cext=.true.)

END SUBROUTINE ropp_pp_merge_profile

