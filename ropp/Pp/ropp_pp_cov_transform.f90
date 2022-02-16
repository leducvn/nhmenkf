! $Id: ropp_pp_cov_transform.f90 3491 2013-02-06 12:43:43Z idculv $

SUBROUTINE ropp_pp_cov_transform(y, x, a, cov)

!****s* TropoPauseHeight/ropp_pp_cov_transform *
!
! NAME
!   ropp_pp_cov_transform
!
! SYNOPSIS
!   Calculate covariance transform
!
!   CALL ropp_pp_cov_transform(y, x, a, cov)
!
! DESCRIPTION
!   Evaluate cov = integral from max(y-a/2, min(y)) to min(y+a/2, max(y))
!   of x(y') [x(y')-x(y)] dy', divided by a.
!
! INPUTS
!    REAL(wp), DIMENSION(:), INTENT(IN)     :: y, x
!    REAL(wp),               INTENT(IN)     :: a
!
! OUTPUTS
!    REAL(wp), DIMENSION(:), INTENT(OUT)    :: cov
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes,  ONLY: wp => EightByteReal

  IMPLICIT NONE

! 1.1 Declarations
! ----------------

  REAL(wp), DIMENSION(:), INTENT(IN)        :: y, x
  REAL(wp),               INTENT(IN)        :: a
  REAL(wp), DIMENSION(:), INTENT(INOUT)     :: cov

! Local variables
  REAL(wp), DIMENSION(:), ALLOCATABLE       :: h, xh
  REAL(wp)                                  :: xi, yi, aby2, y_min, y_max, y_l, y_u
  REAL(wp)                                  :: m, delta_y, mdelta_y, corr_l, corr_u
  INTEGER                                   :: i, j, n, i_l, i_u

! 1.2 Compute the transform
! -------------------------

  aby2 = a / 2.0_wp

  y_min = MINVAL(y)  ;  y_max = MAXVAL(y)

  n = SIZE(y)  ! Can safely assume to be the common size

  ALLOCATE(h(n), xh(n))

  DO i=1,n

    xi = x(i)  ;  yi = y(i)

    h = x - xi

    xh = x*h


! Locate lower index
    y_l = MAX(yi - aby2, y_min)
    i_l = SUM(MINLOC(y, mask=y>y_l))
    IF (i_l > 1) h(1:i_l-1) = 0.0_wp


! Locate upper index
    y_u = MIN(yi + aby2, y_max)
    i_u = SUM(MAXLOC(y, mask=y<y_u))
    IF (i_u < n) h(i_u+1:n) = 0.0_wp


! Body of integral (trapezoidal rule)

    cov(i) = 0.0_wp

    IF (i_l <= (i_u-1)) THEN

      DO j=i_l,i_u-1
        cov(i) = cov(i) + &
                 0.5_wp*(xh(j) + xh(j+1))*(y(j+1) - y(j))
      END DO

    END IF


! Correction at lower limit (linearly extrapolate)

    IF (i_l < (n-1)) THEN

      m = (x(i_l+1) - x(i_l)) / ((y(i_l+1) - y(i_l)) + TINY(1.0_wp))

      delta_y = y(i_l) - y_l

      mdelta_y = m * delta_y

      corr_l = (x(i_l)*(x(i_l) - xi))          - &
               ((x(i_l) - 0.5_wp*xi)*mdelta_y) + &
               ((mdelta_y*mdelta_y)/3.0_wp)

      corr_l = corr_l * delta_y

      cov(i) = cov(i) + corr_l

    END IF


! Correction at upper limit (linearly extrapolate)

    IF (i_u > 1) THEN

      m = (x(i_u) - x(i_u-1)) / ((y(i_u) - y(i_u-1)) + TINY(1.0_wp))

      delta_y = y_u - y(i_u)

      mdelta_y = m * delta_y

      corr_u = (x(i_u)*(x(i_u) - xi))          + &
               ((x(i_u) - 0.5_wp*xi)*mdelta_y) + &
               ((mdelta_y*mdelta_y)/3.0_wp)

      corr_u = corr_u * delta_y

      cov(i) = cov(i) + corr_u

    END IF


  END DO

! Normalise

  cov = cov / a

! Clean up

  DEALLOCATE(h, xh)


END SUBROUTINE ropp_pp_cov_transform
