! $Id: ropp_pp_spline.f90 1960 2008-11-13 16:28:24Z frhl $

MODULE ropp_pp_spline

!****m* Modules/ropp_pp_spline *
!
! NAME
!    ropp_pp_spline - Spline interpolation routines
!
! SYNOPSIS
!    use ropp_pp_spline
!
! DESCRIPTION
!    This module provides spline interpolation routines used by
!    the ROPP pre-processor package.
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

  USE typesizes, ONLY: wp => EightByteReal
 
CONTAINS

!****s* PPSpline/ropp_pp_basic_splines *
!
! NAME
!    ropp_pp_basic_splines - Generate a matrix of basic polynomials for
!                            polynomial regression
!
! SYNOPSIS
!    call ropp_pp_basic_splines(X, Xs, K)
!
! DESCRIPTION
!    Generation of matrix of basic polynomials for polynomial regression
!
!****
!
  SUBROUTINE ropp_pp_basic_splines(X, Xs, K)

    IMPLICIT NONE

    ! 1.1 Declarations

    REAL(wp), DIMENSION(:), INTENT(in)    :: X  ! Grid of argument Z
    REAL(wp), DIMENSION(:), INTENT(in)    :: Xs ! X-grid of delta-splines
    REAL(wp), DIMENSION(:,:), INTENT(out) :: K  ! Matrix of basic polynomials

    INTEGER                       :: i    ! x index
    INTEGER                       :: j    ! Delta-spline number
    REAL(wp), DIMENSION(:), ALLOCATABLE :: S    ! Delta-spline
    REAL(wp), DIMENSION(:), ALLOCATABLE :: D2S  ! Delta-spline 2nd derivative

    ALLOCATE(S(SIZE(xs)))
    ALLOCATE(D2S(SIZE(xs)))

    ! 1.2 Generate basic function matrix K

    S(:) = 0.0_wp

    DO j=1, SIZE(K,2)
       S(j) = 1.0_wp
       CALL ropp_pp_init_spline(Xs, S, D2S)
       DO i=1, SIZE(X)
          CALL ropp_pp_interpol_spline(Xs, S, D2S, X(i), K(i,j))
       ENDDO
       S(j) = 0.0_wp
    ENDDO

    DEALLOCATE(S)
    DEALLOCATE(D2S)

  END SUBROUTINE ropp_pp_basic_splines

!****s* PPSpline/ropp_pp_init_spline *
!
! NAME
!    ropp_pp_init_spline - Calculation of 2nd derivative of a spline
!
! SYNOPSIS
!    call ropp_pp_init_spline(x, f, d2)
!
! DESCRIPTION
!    Calculate 2nd derivative of a spline, drive-through solution
!
!****
!
  SUBROUTINE ropp_pp_init_spline(x, f, d2)

    IMPLICIT NONE

    ! 2.1 Declarations
    
    REAL(wp), DIMENSION(:), INTENT(in)  :: x  ! Argument grid (monotonous)
    REAL(wp), DIMENSION(:), INTENT(in)  :: f  ! Gridded function
    REAL(wp), DIMENSION(:), INTENT(out) :: d2 ! 2nd derivative of spline

    REAL(wp), DIMENSION(:), ALLOCATABLE :: d1
    REAL(wp)                     :: dfl, dfr, df, a, b, c
    INTEGER                      :: i, N

    ! 2.2 Initialisation

    N = SIZE(x) 
    ALLOCATE(d1(N))
    d2(:) = 0.0_wp
    
    ! 2.3 Drive-through calculation of spline coefficients
    
    d1(1) = 0.0_wp
    d2(N) = 0.0_wp
    d1(1) = 0.0_wp

    DO i = 2, N-1

       dfl   = (f(i) - f(i-1))/(x(i) - x(i-1))
       dfr   = (f(i+1) - f(i))/(x(i+1) - x(i))
       df    = (dfr - dfl)/(x(i+1) - x(i-1))
       a     = (x(i) - x(i-1))/(2*(x(i+1) - x(i-1)))
       b     = (x(i+1) - x(i))/(2*(x(i+1) - x(i-1)))
       c     = 1 + a*d1(i-1)
       d1(i) = -b/c
       d2(i) = (3*df - a*d2(i-1))/c
    ENDDO
    
    DO i = N-1, 2, -1
       d2(i) = d1(i)*d2(i+1) + d2(i)
    ENDDO

    DEALLOCATE(d1)
    
  END SUBROUTINE ropp_pp_init_spline

!****s* PPSpline/ropp_pp_interpol_spline *
!
! NAME
!    ropp_pp_interpol_spline - Spline interpolation of a gridded function with
!                              linear extrapolation outside grid extent
!
! SYNOPSIS  
!    call ropp_pp_interpol_spline(x, f, d2, x_int, f_int, fd_int, fd2_int)
!
! DESCRIPTION
!    Search for grid interval containing the interpolation point and summation
!    of polynomail with given spline coefficients
!
!****
!
  SUBROUTINE ropp_pp_interpol_spline(x, f, d2, x_int, f_int, fd_int, fd2_int)
    
    IMPLICIT NONE

    ! 3.1 Declarations

    REAL(wp), DIMENSION(:), INTENT(in)  :: x  ! Argument grid (monotonous)
    REAL(wp), DIMENSION(:), INTENT(in)  :: f  ! Gridded function
    REAL(wp), DIMENSION(:), INTENT(in)  :: d2 ! 2nd derivative of spline
    REAL(wp),               INTENT(in)  :: x_int  ! Interpolation point
    REAL(wp),               INTENT(out) :: f_int  ! Interpolated function value
    REAL(wp), OPTIONAL,     INTENT(out) :: fd_int ! Interpolated 1st derivative
    REAL(wp), OPTIONAL,     INTENT(out) :: fd2_int ! Interpolated 2nd deriv
    
    REAL(wp) :: a1, a2, a3   ! Polynomial coefficients
    REAL(wp) :: dx           ! Grid interval
    REAL(wp) :: dx_t         ! Grid-point to interpolation-point distance
    REAL(wp) :: x_t          ! Interpolation point projected to grid extent
    REAL(wp) :: fd           ! Interpolated derivative
    INTEGER  :: i            ! Array index
    INTEGER  :: N            ! Number of data
    INTEGER  :: i_int        ! Interpolation interval index

    ! 3.2 Location of interpolation point inside grid

    N     = SIZE(x)
    
    x_t   = MIN(MAX(x_int, MIN(x(1),x(N))), MAX(x(1),x(N)))
    
    i_int = ropp_pp_seek_index(x, x_t)
    i     = MAX(i_int, 1)
    
    ! 3.3 Calculation of interpolation coefficients

    dx    = x(i+1) - x(i)
    a2    = d2(i)/2.0_wp
    a3    = (d2(i+1) - d2(i))/(6.0_wp*dx)
    a1    = (f(i+1) - f(i))/dx - dx*(a2 + dx*a3)
    
    ! 3.4 Calculated interpolated value

    dx_t  = x_t - x(i)
    fd    = a1 + dx_t*(2*a2 + dx_t*3*a3)
    f_int = f(i) + dx_t*(a1 + dx_t*(a2 + dx_t*a3)) + &
         fd*(x_int - x_t)
    
    IF (PRESENT(fd_int)) THEN
       fd_int = fd
    ENDIF
    
    IF (PRESENT(fd2_int)) THEN
       fd2_int  = 2.0_wp*a2 + dx_t*6.0_wp*a3
    ENDIF

  END SUBROUTINE ropp_pp_interpol_spline
    
!****s* PPUtils/ropp_pp_seek_index *
!
! NAME
!    ropp_pp_seek_index - Find grid interval containing given point
!
! SYNOPSIS  
!    ip = ropp_pp_seek_index(x, xp)
!
! DESCRIPTION
!    Search for grid interval containing the interpolation point.
!    Combination of Newton-like and dichotomic iterative search.
!     Result:
!         >0   - index of grid interval containing xp (x(ip) <= xp <= x(ip+1))
!          0   - xp is not inside grid
!         -1   - iterations do not converge
!
!****

  FUNCTION ropp_pp_seek_index(x, xp) RESULT(ip)
    
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(in) :: x   ! x-grid (homogeneous)
    REAL(wp),               INTENT(in) :: xp  ! point to locate inside grid
    INTEGER                            :: ip
    
    INTEGER  :: N      ! Number of grid points
    INTEGER  :: Imin   ! Upper estimate of index
    INTEGER  :: Imax   ! Lower estimate of index
    INTEGER  :: Dir    ! Direction of argument change
    INTEGER  :: It     ! Iteration count
    INTEGER  :: di     ! Index increment in iterations
    INTEGER  :: is     ! Step direction count

    ! 4.1 Grid size and direction calculation
   
    N    = SIZE(x)
    Dir  = NINT(SIGN(1.0_wp, x(N)-x(1)))
        
    ! 4.2 Checking if point is inside grid

    IF ((Dir*xp < Dir*x(1)) .OR. (Dir*xp > Dir*x(N))) THEN
       ip = 0
       RETURN
    END IF

    ! 4.3 Initial approximation
    
    Imin = 1
    Imax = N
    ip   = Imin+FLOOR(REAL(Imax - Imin, wp)*(xp - x(Imin))/(x(Imax) - x(Imin)))
    ip   = MAX(1, MIN(ip, N-1))
    
    ! 4.4 Iterative index search

    It = 0
    is = 0
    
    Search: DO
       IF ((Dir*x(ip) <= Dir*xp) .AND. (Dir*xp <= Dir*x(ip+1))) THEN
          EXIT Search
       END IF
       IF (ABS(is) > 1) THEN
          ip = (Imax + Imin)/2
          is = 0
       END IF
       IF (Dir*x(ip+1) < Dir*xp) THEN
          Imin = ip + 1
          di   = FLOOR(REAL(Imax - Imin, wp)*(xp - x(Imin))/(x(Imax) - x(Imin)))
          ip   = Imin + di
          is   = is + 1
       ELSE IF (Dir*xp < Dir*x(ip)) THEN
          Imax = ip
          di   = FLOOR(REAL(Imin - Imax, wp)*(xp - x(Imax))/(x(Imin) - x(Imax)))
          ip   = Imax + di
          is   = is - 1
       END IF
       ip = MAX(1, MIN(ip, N-1))
       It = It + 1
       IF (It > N) THEN
          ip = -1
          EXIT Search
       END IF
    END DO Search
    
  END FUNCTION ropp_pp_seek_index

END MODULE ropp_pp_spline
