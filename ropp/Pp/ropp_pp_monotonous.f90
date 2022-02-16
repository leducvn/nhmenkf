! $Id: ropp_pp_monotonous.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_monotonous(x, d)

!****s* Monotonous/ropp_pp_monotonous *
!
! NAME
!    ropp_pp_monotonous - Find monotonous sequence.
!
! SYNOPSIS
!    call ropp_pp_monotonous(x, d)
! 
! DESCRIPTION
!    This subroutine finds a monotonous sequence by an iterative orthogonal
!    projections and crawling monotonoization method.
!
! INPUTS
!    real(wp), dim(:)  :: x      Coordinate values to be transformed.
!    integer, optional :: d      Sort direction flag 
!                                             > 0 - rising
!                                             < 0 - setting
!
! OUTPUT
!    real(wp), dim(:)  :: x      Transformed sequence
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

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), TARGET, INTENT(inout) :: x       ! Data array 
  INTEGER, OPTIONAL                             :: d       ! Sort direction flag

  REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET   :: y       ! Work array
  REAL(wp), DIMENSION(:), ALLOCATABLE           :: dx_max  ! Upper limit 
  REAL(wp), DIMENSION(:), ALLOCATABLE           :: dx_min  ! Lower limit
  
  REAL(wp), DIMENSION(:), POINTER  :: x0 => null()  ! Pointer for x(1:N-1)
  REAL(wp), DIMENSION(:), POINTER  :: x1 => null()  ! Pointer for x(2:N)
  REAL(wp), DIMENSION(:), POINTER  :: y0 => null()  ! Pointer for y(1:N-1)
  REAL(wp), DIMENSION(:), POINTER  :: y1 => null()  ! Pointer for y(2:N)
  REAL(wp), DIMENSION(:), POINTER  :: xp => null()  ! Pointer for X

  REAL(wp)                         :: xmax          ! Maximum x value         
  REAL(wp)                         :: dx            ! Lower limit array incremnt
  INTEGER                          :: dir           ! Sort direction
  INTEGER                          :: n             ! Number of data points
  INTEGER                          :: i     ! Index
  
  REAL(wp), PARAMETER  :: db = 500.0_wp     ! Draw-back limit
  REAL(wp), PARAMETER  :: eps = 1000.0_wp   ! Accuracy for monotonization
  
!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------
  
  ! 2.1 Sorting direction (default increasing with index)

  IF(PRESENT(d))THEN 
     dir = d
  ELSE
     dir = 1
  ENDIF

  ! 2.2 Allocate arrays
  
  n = SIZE(x)

  ALLOCATE(y(n))
  ALLOCATE(dx_max(n))
  ALLOCATE(dx_min(n))

  IF ( x(1) > x(n) ) THEN
     xp => x(:)
  ELSE
     xp => x(n:1:-1)
  ENDIF
  
  xmax = dir*xp(1)

!-------------------------------------------------------------------------------
! 3. Raw pre-filtering
!-------------------------------------------------------------------------------

  DO i=2, n
     xp(i) = MAX(xmax - db, dir*xp(i))/dir
     xmax  = MAX(xmax, dir*xp(i))
  ENDDO

  y(:) = xp(:)
  
  x0 => xp(1:n-1)
  x1 => xp(2:n)
  y0 => y(1:n-1)
  y1 => y(2:n)
  
!-------------------------------------------------------------------------------
! 4. Orthogonal projections
!-------------------------------------------------------------------------------

  DO 

     IF ( MAXVAL( dir*x0(:) - dir*x1(:)) < eps ) THEN
        EXIT
     ENDIF

     WHERE ( dir*x1(:) < dir*x0(:) )
        y0(:) = (x0(:) + x1(:))/2.0_wp
        y1(:) = (x0(:) + x1(:))/2.0_wp
     endwhere

     xp(:) = y(:)

  ENDDO

!-------------------------------------------------------------------------------
! 5. Crawling monotonization
!-------------------------------------------------------------------------------

  ! 5.1 Estimation of increment

  dx = 0.01_wp * (dir*xp(n) - dir*xp(1))/n

  ! 5.2 Forward and back passes

  dx_max(1) = dir*xp(1)
  DO i=2,n
     dx_max(i) = MAX(dx_max(i-1) + dx, dir*xp(i))
  ENDDO
  
  dx_min(n) = dir*xp(n)
  DO i=n-1,1,-1
     dx_min(i) = MIN(dx_min(i+1) - dx, dir*xp(i))
  ENDDO

  ! 5.3 Monotonization

  xp(:) = (dx_min(:) + dx_max(:))/(2.0_wp*dir)

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(y)
  DEALLOCATE(dx_max)
  DEALLOCATE(dx_min)

END SUBROUTINE ropp_pp_monotonous
