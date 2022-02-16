! $Id: ropp_pp_interpol_log.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_interpol_log(x, newx, array, interp)

!****s* Interpolation/ropp_pp_interpol_log *
!
! NAME
!    ropp_pp_interpol - Interpolate logarithmically.
!
! SYNOPSIS
!    call ropp_pp_interpol_log(x, newx, array, interp)
! 
! DESCRIPTION
!    This subroutine interpolates an array assuming it varies exponentially
!    in x i.e. assuming its log varies linearly as a function of x.  
!
! INPUTS
!    real(wp), dim(:) :: x       Coordinate values.
!    real(wp), dim(:) :: newx    New coordinate values.
!    real(wp), dim(:) :: array   Data to be interpolated (lives on x).
!
! OUTPUT
!    real(wp), dim(:) :: interp  Interpolated data (lives on newx).
!
! NOTES
!    Array must be strictly positive.
!
!    The coordinate array x must be strictly monotonically increasing. If
!    elements of newx are outside the range of x, data will be extrapolated.
!
!    None of the above conditions are checked for, but wrong or unexpected
!    results will be obtained if thone of them is not met. 
!
! SEE ALSO
!    ropp_pp_interpol
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

  REAL(wp), DIMENSION(:), INTENT(in)  :: x
  REAL(wp), DIMENSION(:), INTENT(in)  :: newx
  REAL(wp), DIMENSION(:), INTENT(in)  :: array
  REAL(wp), DIMENSION(:), INTENT(out) :: interp

  INTEGER                             :: i, j, k

!-------------------------------------------------------------------------------
! 2. Do the interpolation
!-------------------------------------------------------------------------------

  DO k = 1, SIZE(newx)
     j = 2
     DO WHILE (j < SIZE(x) .AND. x(j) < newx(k))
        j = j + 1
     ENDDO
     i = MAX(j - 1,1)

     ! Logarithmic interpolation for +ve data
     IF(i > 0 .AND. array(i) > 0.0_wp .AND. array(j) > 0.0_wp)THEN

        interp(k) = LOG(array(j)) + &
          ( (newx(k) - x(j)) / (x(i) - x(j)) * (LOG(array(i)) - LOG(array(j))))
        interp(k) = EXP(interp(k))
        
     ! Linear interpolation for -ve data   
     ELSE
       
        interp(k) = array(j) + &
          ( (newx(k) - x(j)) / (x(i) - x(j)) * (array(i) - array(j)))

     ENDIF

  ENDDO

END SUBROUTINE ropp_pp_interpol_log

