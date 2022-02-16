! $Id: ropp_pp_interpol.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Interpolation/ropp_pp_interpol *
!
! NAME
!    ropp_pp_interpol - Interpolate linearly.
!
! SYNOPSIS
!    call ropp_pp_interpol(x, newx, array, interp)
!
! DESCRIPTION
!    This subroutine interpolates an array assuming it varies linearly in x.
!
! INPUTS
!    real(wp), dim(:) :: x       Coordinate values.
!    real(wp), dim(:) :: newx    New coordinate values.
!    real(wp), dim(:) :: array   Data to be interpolated (lives on x).
!    logical,optional :: Cext    Constant or linear (default) extrapolation

!
! OUTPUT
!    real(wp), dim(:) :: interp  Interpolated data (lives on newx).
!
! NOTES
!    The coordinate array x must be strictly monotonically increasing. If
!    elements of newx are outside the range of x, data will be extrapolated.
!
!    None of the above conditions are checked for, but wrong or unexpected
!    results will be obtained if thone of them is not met.
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
! 1. Linear interpolation - Scalar output
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_interpol_scl(x, newx, array, interp, Cext)

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_pp_spline

  IMPLICIT NONE

  ! 1.1 Declarations

  REAL(wp), DIMENSION(:), INTENT(in)  :: x
  REAL(wp),               INTENT(in)  :: newx
  REAL(wp), DIMENSION(:), INTENT(in)  :: array
  REAL(wp),               INTENT(out) :: interp
  LOGICAL,  OPTIONAL,     INTENT(in)  :: Cext    ! constant/linear (default)

  REAL(wp) :: x_t                 ! Interpolation point projected to grid
  REAL(wp) :: a1                  ! Slope
  INTEGER  :: i, i_int            ! Indices
  LOGICAL  :: cx                  ! Constant/linear extrapolation

  cx = .FALSE.

  ! 1.2 Linear interpolation

  x_t = MIN(MAX(newx, MIN(x(1),x(SIZE(x)))), MAX(x(1),x(SIZE(x))))
  i_int = ropp_pp_seek_index(x, x_t)
  i = MAX(i_int, 1)

  a1 = (array(i+1) - array(i)) / (x(i+1) - x(i))
  IF (PRESENT(Cext)) cx = Cext

  IF (cx) THEN
    interp = array(i) + ((x_t - x(i))*a1)
  ELSE
    interp = array(i) + ((newx - x(i))*a1)
  ENDIF

END SUBROUTINE ropp_pp_interpol_scl

!-------------------------------------------------------------------------------
! 2. Linear interpolation - Array output
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_interpol_arr(x, newx, array, interp, Cext)

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_pp_spline

  IMPLICIT NONE

  ! 2.1 Declarations

  REAL(wp), DIMENSION(:), INTENT(in)  :: x
  REAL(wp), DIMENSION(:), INTENT(in)  :: newx
  REAL(wp), DIMENSION(:), INTENT(in)  :: array
  REAL(wp), DIMENSION(:), INTENT(out) :: interp
  LOGICAL,  OPTIONAL,     INTENT(in)  :: Cext    ! constant/linear (default)

  REAL(wp) :: x_t                 ! Interpolation point projected to grid
  INTEGER  :: i, k, i_int         ! Indices
  LOGICAL  :: cx                  ! Constant/linear extrapolation

  cx = .FALSE.
  IF (PRESENT(Cext)) cx = Cext

  ! 2.2 Interpolation

  DO k = 1, SIZE(newx)

     x_t = MIN(MAX(newx(k), MIN(x(1),x(SIZE(x)))), MAX(x(1),x(SIZE(x))))
     i_int = ropp_pp_seek_index(x, x_t)
     i = MAX(i_int, 1)

     IF (cx) THEN
       interp(k) = array(i) + &
          ( (x_t - x(i)) / (x(i+1) - x(i)) * (array(i+1) - array(i)) )
     ELSE
        interp(k) = array(i) + &
          ( (newx(k) - x(i)) / (x(i+1) - x(i)) * (array(i+1) - array(i)) )
      ENDIF

  ENDDO

END SUBROUTINE ropp_pp_interpol_arr

!-------------------------------------------------------------------------------
! 3. Linear interpolation - Integer array output
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_interpol_int(x, newx, array, interp)

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_pp_spline

  IMPLICIT NONE

  ! 2.1 Declarations

  REAL(wp), DIMENSION(:), INTENT(in)  :: x
  REAL(wp),               INTENT(in)  :: newx
  INTEGER,  DIMENSION(:), INTENT(in)  :: array
  INTEGER,                INTENT(out) :: interp

  INTEGER  :: i, i_int            ! Indices

  ! 2.2 Interpolation

  i_int = ropp_pp_seek_index(x, newx)
  i = MAX(i_int, 1)

  IF (i == 0) THEN
    interp = array(1)
  ELSE IF (i == SIZE(x)) THEN
    interp = array(SIZE(x))
  ELSE
    IF (ABS(newx - x(i)) < ABS(newx - x(i+1))) THEN
      interp = array(i)
    ELSE
      interp = array(i+1)
    ENDIF
  ENDIF

END SUBROUTINE ropp_pp_interpol_int
