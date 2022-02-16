! $Id: ropp_pp_utils.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE ropp_pp_utils

!****m* Modules/ropp_pp_utils *
!
! NAME
!    ropp_pp_utils - Utility routines for ionospheric correction 
!
! SYNOPSIS
!    use ropp_pp_utils
!
! DESCRIPTION
!    This module provides signal and matrix processing routines used by
!    the ROPP pre-processor ionospheric correction package.
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
  
  INTERFACE ropp_pp_matmul
     MODULE PROCEDURE ropp_pp_matmul_arr
     MODULE PROCEDURE ropp_pp_matmul_vecl
     MODULE PROCEDURE ropp_pp_matmul_vecr
  END INTERFACE

CONTAINS

!****s* PPUtils/ropp_pp_matmul *
!
! NAME
!    ropp_pp_matmul - Matrix multiplication 
!
! SYNOPSIS
!    C = ropp_pp_matmul(A, B)
!
! DESCRIPTION
!    Replicates intrinsic 'MATMUL' function
!
!****
!
SUBROUTINE ropp_pp_matmul_sub(A, B, C)

  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: A
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: B
  REAL(wp), DIMENSION(:,:), INTENT(out) :: C
  INTEGER :: i, j, k

  DO i=1,SIZE(A,1)
    DO j=1, SIZE(B,2)
      C(i,j) = 0.0_wp
      DO k=1, SIZE(B,1)
        C(i,j) = C(i,j) + A(i,k)*B(k,j)
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE ropp_pp_matmul_sub

FUNCTION ropp_pp_matmul_arr(A, B) RESULT (C)

  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: A
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: B
  REAL(wp), DIMENSION(SIZE(A,1),SIZE(B,2)) :: C
  INTEGER :: i, j, k

  DO i=1,SIZE(A,1)
    DO j=1, SIZE(B,2)
      C(i,j) = 0.0_wp
      DO k=1, SIZE(B,1)
        C(i,j) = C(i,j) + A(i,k)*B(k,j)
      ENDDO
    ENDDO
  ENDDO
  
END FUNCTION ropp_pp_matmul_arr

FUNCTION ropp_pp_matmul_vecr(A, B) RESULT (C)

  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: A
  REAL(wp), DIMENSION(:), INTENT(in)  :: B
  REAL(wp), DIMENSION(SIZE(A,1)) :: C
  INTEGER :: j, k

  DO j=1, SIZE(A,1)
    C(j) = 0.0_wp
    DO k=1, SIZE(B)
      C(j) = C(j) + A(j,k)*B(k)
    ENDDO
  ENDDO
  
END FUNCTION ropp_pp_matmul_vecr

FUNCTION ropp_pp_matmul_vecl(A, B) RESULT (C)

  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:), INTENT(in)  :: A
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: B
  REAL(wp), DIMENSION(SIZE(B,2)) :: C
  INTEGER :: j, k

  DO j=1, SIZE(B,2)
    C(j) = 0.0_wp
    DO k=1, SIZE(B,1)
      C(j) = C(j) + A(k)*B(k,j)
    ENDDO
  ENDDO
  
END FUNCTION ropp_pp_matmul_vecl


!****s* PPUtils/ropp_pp_invert_matrix *
!
! NAME
!    ropp_pp_invert_matrix - Invert a matrix A(dimY, dimX) 
!
! SYNOPSIS
!    call ropp_pp_invert_matrix(A, B)
!
! DESCRIPTION
!    Gauss elimination
!
!****
!
SUBROUTINE ropp_pp_invert_matrix(A, B)

  IMPLICIT NONE
  
  ! 2.1 Declarations
  
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: A     ! Matrix to invert
  REAL(wp), DIMENSION(:,:), INTENT(out)   :: B     ! Inverted matrix
  
  INTEGER  :: N             ! Matrix dimension
  INTEGER  :: i, k          ! Matrix indeces
  INTEGER  :: m(1)          ! Search index
  REAL(wp) :: Alpha         ! Row combination coefficien
  REAL(wp) :: Det           ! Matrix determinant
  REAL(wp) :: Amin, Amax    ! Max and min diagonal elements
  REAL(wp), DIMENSION(:), ALLOCATABLE :: F  ! Exchange work space
  
  ! 2.2 Check array shape

  IF (SIZE(A,1) /= SIZE(A,2)) THEN
     RETURN
  ENDIF
  IF (ANY(SHAPE(B) /= SHAPE(B))) THEN
     RETURN
  ENDIF

  ALLOCATE(F(SIZE(A,1)))
    
  ! 2.3 Initialization

  N      = SIZE(A,1)
  B(:,:) = 0.0_wp
  
  DO i=1,N
     B(i,i) = 1.0_wp
  ENDDO
  
  ! 2.4 Elimination of sub-diagonal elements

  Lower: DO k=1,N-1
     
     IF (A(k,k) == 0.0_wp) THEN
        m(:) = MAXLOC(A(k+1:N,k), A(k+1:N,k) /= 0.0_wp)
        IF (m(1) == 0) THEN
           EXIT Lower
        END IF
        i      = k + m(1)
        F(:)   = A(k,:)
        A(k,:) = A(i,:)
        A(i,:) = F(:)
        F(:)   = B(k,:)
        B(k,:) = B(i,:)
        B(i,:) = F(:)
     END IF
     
     DO i=k+1,N
        Alpha  = A(i,k)/A(k,k)
        A(i,:) = A(i,:) - Alpha*A(k,:)
        B(i,:) = B(i,:) - Alpha*B(k,:)
     END DO
     
  END DO Lower

  ! 2.5 Checking for degenerated matrix
  
  Det  = 1.0_wp
  Amin = ABS (A(1,1))
  Amax = ABS (A(1,1))
  Diagonal: DO i=1,N
     Amin = MIN(ABS (A(i,i)), Amin)
     Amax = MAX(ABS (A(i,i)), Amax)
     Det  = Det*A(i,i)
  END DO Diagonal
  
  IF (Det == 0.0_wp) THEN
     RETURN
  END IF

  ! 2.6 Elimination of super-diagonal elements

  Upper: DO k=N,2,-1
     DO i=k-1,1,-1
        Alpha  = A(i,k)/A(k,k)
        A(i,:) = A(i,:) - Alpha*A(k,:)
        B(i,:) = B(i,:) - Alpha*B(k,:)
     END DO
  END DO Upper
  
  DO i=1,N
     B(i,:) = B(i,:)/A(i,i)
  END DO
  
  DEALLOCATE(F)

END SUBROUTINE ropp_pp_invert_matrix

!****s* PPUtils/ropp_pp_quasi_invert *
!
! NAME
!    ropp_pp_quasi_invert - Quasi-inverse of a matrix K(dimY, dimX) 
!
! SYNOPSIS
!    call ropp_pp_quasi_invert(K, Q)
!
! DESCRIPTION
!    1. dimY >= dimX:
!       x = Qy - vector minimizing ||Kx - y||
!       QK = E; Q is left inverse operator.
!       Q = (K^T K)^-1 K^T
!    2. dimY <= dimX:
!       x = Qy - solution minimizing ||x||
!       KQ = E; Q is right inverse operator.
!       Q = K^T (KK^T)^-1
!
!****
!
  SUBROUTINE ropp_pp_quasi_invert(K, Q)

    IMPLICIT NONE

    ! 3.1 Declarations

    REAL(wp), DIMENSION(:,:), INTENT(in)  :: K  ! Matrix to quasi-invert
    REAL(wp), DIMENSION(:,:), INTENT(out) :: Q  ! Quasi-inverse
    
    INTEGER :: DimX, DimY  ! Dimension of x and y spaces
    INTEGER                               :: N  ! Work matrix dimension
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: W  ! Work arrays for inversion
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: WI ! Work arrays for inversion
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KT

    ! 3.2 Initialization

    DimX = SIZE(K,2)
    DimY = SIZE(K,1)
    
    IF (ANY(SHAPE(Q) /= (/ DimX, DimY /))) THEN
       RETURN
    END IF
    
    N = MIN(DimX, DimY)
    ALLOCATE (W(N,N))
    ALLOCATE (WI(N,N))
    ALLOCATE (KT(DimX,DimY))
     
    KT = TRANSPOSE(K)

    ! 3.3 Quasi-Inversion
    
    IF (DimY >= DimX) THEN
       
       ! 3.3.1 Left inversion   W=K^TK, WI=(K^TK)^-1, Q=(K^TK)^-1 K^T

       W = ropp_pp_matmul(KT, K)
       CALL ropp_pp_invert_matrix(W, WI)
!       Q = ropp_pp_matmul(WI, KT)
       CALL ropp_pp_matmul_sub(WI, KT, Q)

    ELSE
       
       ! 3.3.2. Right inversion W=KK^T, WI=(KK^T)^-1, Q=K^T(KK^T)^-1
            
       W = ropp_pp_matmul(K, KT)
       CALL ropp_pp_invert_matrix(W, WI)
!       Q = ropp_pp_matmul(KT, WI)
       CALL ropp_pp_matmul_sub(KT, WI, Q)

    END IF
    
    ! 3.4 Clean up
    
    DEALLOCATE (W)
    DEALLOCATE (WI)
    DEALLOCATE (KT)

  END SUBROUTINE ropp_pp_quasi_invert
  
!****s* PPUtils/ropp_pp_polynomial *
!
! NAME
!    ropp_pp_polynomial - Calculation of a polynomial and its derivative
!
! SYNOPSIS
!    call ropp_pp_polynomial(c, x, P, DP)
!
! DESCRIPTION
!    Horner scheme
!
!****
!
  SUBROUTINE ropp_pp_polynomial(c, x, P, DP)

    USE typesizes, ONLY: wp => EightByteReal
    
    ! 4.1 Declarations
    
    IMPLICIT NONE
    
    REAL(wp), DIMENSION(0:), INTENT(in)  :: c     ! polynomial coefficients
    REAL(wp),                INTENT(in)  :: x     ! polynomial argument
    REAL(wp),                INTENT(out) :: P     ! polynomial value
    REAL(wp), OPTIONAL,      INTENT(out) :: DP    ! polynomial derivative
    
    INTEGER :: n, i
    
    ! 4.2 Calculate polynomial value

    n = UBOUND(c,1)
    
    P = c(n)
    DO i = n-1, 0, -1
       P = c(i) + x*P
    ENDDO
    
    ! 4.3 Calculate polynomial derivative
    
    IF ( PRESENT(DP) ) THEN
       DP = c(n)*REAL(n,wp)
       DO i = n-1, 1, -1
          DP = c(i)*REAL(i,wp) + x*DP
       ENDDO
    ENDIF
    
  END SUBROUTINE ropp_pp_polynomial

!****s* PPUtils/ropp_pp_init_polynomial *
!
! NAME
!    ropp_pp_init_polynomial - Generate matrix of basic polynomials 
!                              for regression
!
! SYNOPSIS
!    call ropp_pp_init_polynomial(x, K)
!
! DESCRIPTION
!    Compute basic polynomials as
!        K_ij = f_j(x_i) = (x_i)^j     for j=0..UBound(K,2)
!    
!****
!
  SUBROUTINE ropp_pp_init_polynomial(x, K)

    USE typesizes, ONLY: wp => EightByteReal
    
    ! 5.1 Declarations
    
    IMPLICIT NONE
    
    REAL(wp), DIMENSION(1:),    INTENT(in)  :: x  ! grid of argument x
    REAL(wp), DIMENSION(1:,0:), INTENT(out) :: K  ! matrix of polynomials
    INTEGER                                 :: i,j
    
    ! 5.2 Generate polynomials
        
    K(:,0) = 1.0_wp
    DO j=1,UBOUND(K,2)
      DO i=1,SIZE(x)
        IF (x(i) > 0) THEN
          K(i, j) = x(i)**REAL(j,wp)
        ELSE
          K(i, j) = x(i)**INT(REAL(j,wp))
        ENDIF
      ENDDO
    ENDDO
    
  END SUBROUTINE ropp_pp_init_polynomial

!****s* PPUtils/ropp_pp_regression *
!
! NAME
!    ropp_pp_regression - Linear regression
!
! SYNOPSIS
!    call ropp_pp_regression(K, y, a)
!
! DESCRIPTION
!    Quasi-inversion of matrix of basic functions:
!       || Sum_j a_j f_j(x_i) - y_i || = min
!    Solution is a = Q y where Q is left inverse of K_ij = f_j(x_i) 
!
!****
!
  SUBROUTINE ropp_pp_regression(K, y, a)

    USE typesizes, ONLY: wp => EightByteReal
    
    ! 6.1 Declarations
    
    IMPLICIT NONE
    
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: K  ! matrix of functions
    REAL(wp), DIMENSION(:),   INTENT(in)  :: y  ! regression data
    REAL(wp), DIMENSION(:),   INTENT(out) :: a  ! regression coefficients
    
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Q ! quasi-inverse matrix
    
    ALLOCATE(Q(SIZE(K,2), SIZE(K,1)))

    ! 6.2 Invert matrix of functions
        
    CALL ropp_pp_quasi_invert(K, Q)

    ! 6.3 Solve to find coefficients
    
    
    a(:) = ropp_pp_matmul(Q, y(:))

    DEALLOCATE(Q)
    
  END SUBROUTINE ropp_pp_regression

!****s* PPUtils/ropp_pp_residual_regression *
!
! NAME
!    ropp_pp_residual_regression - linear regression on residual
!
! SYNOPSIS
!    call ropp_pp_residual_regression(K, x, y, a)
!
! DESCRIPTION
!    Updates the regression coefficients using
!    the residual from a previous regression estimate
!
!****
!

  SUBROUTINE ropp_pp_residual_regression(K, x, y, a)

    USE typesizes, ONLY: wp => EightByteReal
    
    ! 7.1 Declarations
    
    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:),INTENT(in)   :: K  ! matrix of functions
    REAL(wp), DIMENSION(:),  INTENT(in)   :: x  ! independent variable
    REAL(wp), DIMENSION(:),  INTENT(in)   :: y  ! regression data
    REAL(wp), DIMENSION(:),  INTENT(inout):: a  ! regression coefficients

    REAL(wp), DIMENSION(:),   ALLOCATABLE :: y_a   ! polynomial data
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: y_res ! residual data
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: a_res ! residual coefficients
    INTEGER :: i

    ALLOCATE(y_a(SIZE(y)))
    ALLOCATE(y_res(SIZE(y)))
    ALLOCATE(a_res(SIZE(a)))

    ! 7.2 Calculate polynomial data from input regression coefficients

    DO i=1,SIZE(y)
      CALL ropp_pp_polynomial(a, x(i), y_a(i))
    END DO

    ! 7.3 Perform regression on the residuals

    y_res(:) = y(:) - y_a(:)

    CALL ropp_pp_regression(K, y_res, a_res)

    ! 7.4 Update the regression coefficients

    a(:) = a(:) + a_res(:)

    DEALLOCATE(y_a)
    DEALLOCATE(y_res)
    DEALLOCATE(a_res)

  END SUBROUTINE ropp_pp_residual_regression

!****s* PPUtils/ropp_pp_nearest_power2 *
!
! NAME
!    ropp_pp_nearest_power2 - Find power of 2 nearest input integer
!
! SYNOPSIS
!    n2 = ropp_pp_nearest_power2(n)
!
!**** 

  FUNCTION ropp_pp_nearest_power2(N) RESULT(nearest_power2)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: N                  ! given integer number
    INTEGER             :: nearest_power2     ! power of 2 right-nearest to N
    INTEGER             :: i

    i = 1

    DO 
       IF (i >= N) THEN
          Nearest_Power2 = i
          RETURN
       ENDIF
       i = 2*i
    ENDDO
    
  END FUNCTION ropp_pp_nearest_power2


!****s* PPUtils/ropp_pp_isnan *
!
! NAME
!    ropp_pp_isnan - Test if variable value is NaN
!
! SYNOPSIS
!    true_or_false = ropp_pp_isnan(x)
!
!**** 

  FUNCTION ropp_pp_isnan(x) RESULT(it_is)
    
    USE typesizes, ONLY: wp => EightByteReal
    IMPLICIT NONE
    
    REAL(wp), INTENT(in) :: x
    LOGICAL              :: it_is
    
    it_is = .FALSE.

    IF ( x /= x ) it_is = .TRUE.
    IF ( x + 1.0_wp == x ) it_is = .TRUE.
    IF ((x > 0) .EQV. (x <= 0)) it_is=.true.

  END FUNCTION ropp_pp_isnan

END MODULE ropp_pp_utils
