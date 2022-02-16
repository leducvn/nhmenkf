! $Id: ropp_io_thin_sg.f90 3551 2013-02-25 09:51:28Z idculv $
!
!****s* Thin/ropp_io_thin_sg *
!
! NAME
!   ropp_io_thin_sg
!
! SYNOPSIS
!   call ropp_io_thin_sg (nLev, Val[,npoints=npoints][,order=order])
!
! DESCRIPTION
!   This routine calculates coefficients of a linear filter following
!   the Savitzky-Golay procedure as described in Numerical Recipes
!   (see Ref.) and applies the smoothing kernel to the given data array.
!   Optional arguments can be used to control the smoothing.
!
! INPUTS
!   nLev      int   Number of levels in Val
!   Val       dflt  Original data values
!   order     int   Order of fitted polynomial (optional, default: 1)
!   npoints   int   Number of points left & right of centre
!                    (optional, default: 1)
!   order     int   Order of fitted polynomial. Must be <= npoints
!                    (optional, default: 1)
!
! OUTPUTS
!   Val       dflt  Smoothed data values
!
! CALLS
!
! USES
!   typesizes
!
! NOTES
!   This is a re-implementation of the Numerical Recipes routine
!   savgol to solve the least squares problem. In contrast to the 
!   Numerical Recipes routine, the full matrix inversion is done, allowing 
!   for the simultaneous calculation of the smoothing polynomial and all its 
!   derivatives at once. The smoothing filter is then applied to the given 
!   array of data points.
!
! REFERENCES
!   W.H. Press, S.A. Teukolsjy, W.T. Vetterling and B.P. Flannery,
!   Numerical Recipes in C - The Art of Scientific Computing.
!   2nd Ed., Cambridge University Press, 1992.
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

SUBROUTINE ropp_io_thin_sg ( nLev,    & ! (in)
                             Val,     & ! (inout)
                             npoints, & ! (optional in)
                             order )    ! (optional in)

! Declarations
! ------------

  USE typeSizes,  ONLY: wp => EightByteReal
  USE messages
  USE ropp_utils, ONLY: ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER,  PARAMETER  :: lder = 0  ! SG: Order of derivative

! Argument list parameters

  INTEGER,                     INTENT(IN)    :: nLev
  REAL(wp), DIMENSION(1:nLev), INTENT(INOUT) :: Val
  INTEGER,  OPTIONAL,          INTENT(IN)    :: npoints
  INTEGER,  OPTIONAL,          INTENT(IN)    :: order

! Local variables

  REAL(wp), DIMENSION(1:nLev)                :: arr
  REAL(wp)                                   :: sum, factor
  INTEGER                                    :: i, j, k, ii, jj, kk, jk, iy, jx, jy
  INTEGER                                    :: mm, nn, np, mmp1
  INTEGER                                    :: shift
  INTEGER                                    :: m, nl, nr
  CHARACTER(len = 256)                       :: routine


! Local allocatable variables

  REAL(wp), DIMENSION(:),   ALLOCATABLE      :: coeff, b, kernel2
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: matrx, design 

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_thin_sg')

! 1. Initialisation
! -----------------

!   1.1. Set optional values
!   ------------------------

  IF (PRESENT(npoints)) THEN
    nl = npoints
    nr = npoints
  ELSE
    nl = 1
    nr = 1
  ENDIF

  IF (PRESENT(order)) THEN
    m = MIN(order,nl,nr)
  ELSE
    m = 1
  ENDIF

!   1.2 Useful values
!   ----------------

  np = nl + nr + 1

!   1.3 Allocate array variables
!   ----------------------------

  ALLOCATE(coeff(np))
  ALLOCATE(kernel2(1:np))
  ALLOCATE(matrx(1:np,1:nLev))
  ALLOCATE(design(m+1,m+1))
  ALLOCATE(b(m+1))

! 2. Set up the design matrix and right hand side
! -----------------------------------------------

  b = 0.0_wp
  b(lder+1) = 1.0_wp

  DO i=0,2*m
     sum = 0.0_wp
     IF(i == 0) sum = 1.0_wp
     
     DO k=1,nr
        sum = sum + REAL(k, wp)**i
     ENDDO
     DO k=1,nl
        sum = sum + REAL(-k, wp)**i
     ENDDO

     mm = MIN(i, 2*m - i)

     DO j = -mm, mm, 2
        design(1+(i+j)/2, 1+(i-j)/2) = sum
     ENDDO

  ENDDO

! 3. Solve the minimisation problem
! ---------------------------------

  CALL solve(design,b)

! 4. Calculate Savitzky-Golay coefficients
! ----------------------------------------

  coeff = 0.0_wp
  DO k=-nl,nr
     sum=b(1)
     factor=1.0_wp
     DO mm=1,m
        factor=factor*k
        sum=sum+b(mm+1)*factor
     ENDDO
     kk=k+nl+1
     coeff(kk)=sum
  ENDDO

! 5. Initialisation of the matrix
! -------------------------------

  mm    = SIZE(coeff)/2        ! Note integer division
  nn    = SIZE(Val)
  matrx = 0.0_wp

!   5.1 Lower index edge
!   --------------------

  DO i = 1, mm
    shift = mm - i + 1
    kernel2 = EOSHIFT(coeff, shift)
    DO jk = 1, SIZE(coeff) - shift
      jj = jk
      ii = mm + 1 + i - jj
      matrx(ii,jj) = kernel2(jk)
    END DO
  END DO

!   5.2 Upper index edge
!   --------------------

  DO i = nn-mm+1, nn
    shift = i - (nn-mm+1) + 1
    kernel2 = EOSHIFT(coeff, -shift)
    DO jk = shift+1, SIZE(kernel2)
      jj = i - mm + jk - shift - 1
      ii = mm + 1 + i - jj
      matrx(ii,jj) = kernel2(jk)
    END DO
  END DO

!   5.3 In between edges
!   --------------------

  DO i = mm+1, nn-mm
    DO jk = 1, SIZE(coeff)
      jj = i - mm + jk - 1
      ii = mm + 1 + i - jj
      matrx(ii,jj) = coeff(jk)
    END DO
  END DO

! 6. Do the multiplication  y := A*x
! -----------------------------------

  arr  = 0.0_wp
  jx   = 1
  jy   = 1
  mmp1 = mm + 1
  DO j = 1, nn
    IF ( val(jx) > ropp_mdtv ) THEN
      iy = jy
      k  = mmp1 - j
      DO i = MAX(1,j-mm), MIN(nn,j+mm)
        arr(iy) = arr(iy) + matrx(k+i,j) * Val(jx)
        iy      = iy  + 1
      END DO
    END IF
    jx = jx + 1
    IF ( j > mm ) jy = jy + 1
   END DO

! 7. Deallocate working arrays
! ----------------------------

  DEALLOCATE(coeff, kernel2, matrx, b, design)

! 8. Mask off end-effects and points contaminated by missing data
! ---------------------------------------------------------------

  arr(1:1+nl)   = ropp_MDFV
  arr(nn-nr:nn) = ropp_MDFV

  DO i = 1, nn
    IF ( Val(i) < ropp_MDTV ) THEN
      ii = MAX ( 1,  i-nl )
      jj = MIN ( nn, i+nr )
      arr(ii:jj) = ropp_MDFV
    END IF
  END DO

! 9. Return smoothed array
! ------------------------

  Val = arr

  CALL message_set_routine(routine)

CONTAINS


! 10. Solve linear matrix equation A x = b using Cholesky decomposition
! ---------------------------------------------------------------------
 
  SUBROUTINE solve(a,b)

    USE messages

    IMPLICIT NONE
    
    REAL(wp), DIMENSION(:,:) :: A
    REAL(wp), DIMENSION(:)   :: b
    
    REAL(wp), DIMENSION(SIZE(b))          :: x
    REAL(wp), DIMENSION(SIZE(b), SIZE(b)) :: G
    REAL(wp), DIMENSION(SIZE(b))          :: temp
    REAL(wp), DIMENSION(SIZE(b))          :: y
    INTEGER                               :: i, j, n
    
    ! 10.1 Check that matrix is a square matrix
    
    n = SIZE(A,1)
    
    IF(n /= SIZE(A,2))THEN
      CALL message(msg_fatal, "LSH matrix is not square - aborting")
    ENDIF
    
    ! 10.2 Determine Cholesky matrix G
    
    G(:,:) = 0.0_wp

    DO i = 1, n
       temp(i:n) = A(i:n, i)      ! store lower triangle matrix
       
       IF(i /= 1)THEN
          DO j = 1, i-1
             temp(i:n) = temp(i:n) - G(i,j)*G(i:n, j)
          ENDDO
       ENDIF
       
       !  10.3 Compute elements of Cholesky matrix G
       G(i:n, i) = temp(i:n) / SQRT(temp(i))
    ENDDO
    
    !  10.4 Solve G.y = b for y by forward substitution
    
    y = b
    y(1) = y(1) / G(1,1)
    DO j = 2, n
       y(j) = (y(j) - DOT_PRODUCT( G(j,1:j-1),y(1:j-1))) / G(j,j)
    ENDDO
    
    !  10.5 Solve G^T.x = y for x by backward substitution
    
    x = y
    x(n) = x(n) / G(n,n)
    DO j = n-1, 1, -1
       x(j) = (x(j) - DOT_PRODUCT(G(j+1:n,j),x(j+1:n))) / G(j,j)
    ENDDO
    
    b = x
    
  END SUBROUTINE solve
  
END SUBROUTINE ropp_io_thin_sg
