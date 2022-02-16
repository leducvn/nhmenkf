! $Id: ropp_pp_smooth_profile.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_smooth_profile(impact, bangle, smooth, config)

!****s* IonosphericCorrection/ropp_pp_smooth_profile *
!
! NAME
!    ropp_pp_smooth_profile - Filter bending angle profile 
!                   
!                   
! SYNOPSIS
!    call ropp_pp_smooth_profile(impact, bangle, smooth, config)
! 
! DESCRIPTION
!    This routine filters a signal by least-square fitting a polynomial
!    in sliding windows
!
! INPUTS
!    real(wp), dimension(:) :: impact      ! Impact parameters 
!    real(wp), dimension(:) :: bangle      ! Bending angles
!    type(ppConfig)         :: config      ! Configuration parameters 
!
! OUTPUT
!    real(wp), dimension(:) :: smooth      ! Smoothed bending angles 
!    type(ppConfig)         :: config      ! Configuration parameters 
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
! USE ropp_pp, not_this => ropp_pp_smooth_profile
  USE ropp_pp_utils
  USE ropp_pp_types, ONLY: PPConfig

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), TARGET, INTENT(in)  :: impact  ! Impact parameters 
  REAL(wp), DIMENSION(:), TARGET, INTENT(in)  :: bangle  ! Bending angles
  REAL(wp), DIMENSION(:),         INTENT(out) :: smooth  ! Smoothed bending
  TYPE(ppConfig),               INTENT(inout) :: config  ! Config parameters

  INTEGER :: n, ws, w_model
  INTEGER :: i, imin, imax
  REAL(wp) :: pmax, pmin

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: A

  REAL(wp), DIMENSION(:), POINTER :: window_bangle => NULL()
  REAL(wp), DIMENSION(:), POINTER :: window_impact => NULL()

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------
  
  n = SIZE(impact)
  pmax = MAXVAL(impact)
  pmin = MINVAL(impact)
  
  w_model = CEILING(config%fw_smooth*(n-1.0_wp)/ABS( pmax - pmin ))

  ! 2.1 Smoothing window size
  
  ws = MIN(2*(w_model/2) + 1, n)

  ! 2.2 Memory allocation
  
  ALLOCATE(K(ws,0:config%np_smooth))
  ALLOCATE(A(0:config%np_smooth))

!-------------------------------------------------------------------------------
! 3. Sliding polynomial regression
!-------------------------------------------------------------------------------

  DO i=1,n

     ! 3.1 Positioning sliding window
     
     Imin = MAX(1, i-ws/2)
     Imax = MIN(n, i+ws/2)
     IF(Imin == 1) Imax = ws
     IF(Imax == n) Imin = n - ws + 1
     
     window_impact => impact(imin:imax)

     ! 3.2 Computation of basic polynomials
     
     CALL ropp_pp_init_polynomial(window_impact-impact(i), K)
     
     ! 3.3 Sliding polynomial regression

     window_bangle => bangle(imin:imax)
     
     CALL ropp_pp_regression(K(:,:), window_bangle(:), A(:))

     smooth(i) = A(0)

  ENDDO

  DEALLOCATE(K)
  DEALLOCATE(A)

END SUBROUTINE ropp_pp_smooth_profile

