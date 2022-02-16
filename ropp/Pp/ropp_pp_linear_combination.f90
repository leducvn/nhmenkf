! $Id: ropp_pp_linear_combination.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_linear_combination(impact_L1, bangle_L1, impact_L2,   &
                                      bangle_L2, impact_LC, bangle_LC)

!****s* IonosphericCorrection/ropp_pp_linear_combination *
!
! NAME
!    ropp_pp_linear_combination - 
!                   Calculate a one dimensional bending angle
!                   profile on L1 impact heights from L1/L2 bending angles 
!                   using linear combination
!                   
! SYNOPSIS
!    call ropp_pp_linear_combination(impact_L1, bangle_L1, 
!                                    impact_L2, bangle_L2,
!                                    impact_LC, bangle_LC)
! 
! DESCRIPTION
!    This routine calculates bending angles at a given set of impact parameters
!    from vertical profiles of bending angles at the two measurement 
!    frequencies (channels) L1 and L2.
!
! INPUTS
!    real(wp), dimension(:) :: impact_L1   ! Impact parameters of channel L1
!    real(wp), dimension(:) :: bangle_L1   ! Bending angles for channel L1
!    real(wp), dimension(:) :: impact_L2   ! Impact parameters of channel L2
!    real(wp), dimension(:) :: bangle_L2   ! Bending angles for channel L2
!
! OUTPUT
!    real(wp), dimension(:) :: impact_LC   ! Impact parameters of channel L1
!    real(wp), dimension(:) :: bangle_LC   ! Corrected bending angles 
!
! NOTES
!    Uses a linear combination of bending angles at the common impact parameter
!    (impact_L1)
!                        
!                          bangle1(impactL1)*f1*f1  - bangle2(impactL1)*f2*f2
!   bangle_LC(impact_L1) = --------------------------------------------------
!                                          f1*f1  - f2*f2
!
! REFERENCES
!     M.E. Gorbunov
!     Ionospheric correction and statistical optimization of radio 
!     occultation data
!     Radio Science, 37(5), 1084
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
  USE ropp_utils, ONLY: ropp_MDTV, ropp_MDFV
  USE ropp_pp_constants
! USE ropp_pp, not_this => ropp_pp_linear_combination

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L1  ! L1 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L1  ! L1 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L2  ! L2 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L2  ! L2 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(out) :: impact_LC  ! Output impact parameters
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_LC  ! Corrected bending angles 

  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_C2  ! Interpolated L2 bangle
  INTEGER                             :: n_obs      ! Number of observations
  INTEGER                             :: i          ! Index 

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_obs = SIZE(impact_L1)
  ALLOCATE(bangle_C2(n_obs))

!-------------------------------------------------------------------------------
! 3. Linear combination
!-------------------------------------------------------------------------------
  
  CALL ropp_pp_interpol(impact_L2, impact_L1, bangle_L2, bangle_C2)

  DO i=1,n_obs

     IF (bangle_C2(i) > ropp_MDTV) THEN
        
        IF (bangle_L1(i) > ropp_MDTV) THEN 
           bangle_LC(i) = (bangle_L1(i) * f_L1**2 - bangle_C2(i) * f_L2**2) / &
                              (f_L1**2 - f_L2**2)
        ELSE
           bangle_LC(i) = bangle_L1(i)
        ENDIF
        
     ELSE
        
        bangle_LC(i) = ropp_MDFV

     ENDIF

     impact_LC(i) = impact_L1(i)

  ENDDO

  DEALLOCATE(bangle_C2)

END SUBROUTINE ropp_pp_linear_combination
