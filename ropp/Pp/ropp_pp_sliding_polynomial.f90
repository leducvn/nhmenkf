! $Id: ropp_pp_sliding_polynomial.f90 2228 2009-09-01 15:36:10Z frhl $

!****s* FFT/ropp_pp_sliding_polynomial *
!
! NAME
!    ropp_pp_sliding_poly - Least-square fitting polynomial in sliding
!                           windows
!
! SYNOPSIS
!    call ropp_pp_sliding_polynomial(t, s, w, np, fs, ds)
!
! DESCRIPTION
!    Least-square fitting polynomial in sliding windows 
!
! INPUTS
!    real(wp)             :: t       Time
!    real(wp), dim([:],:) :: s       Signal samples ([channel],time)
!    integer, [dim(:)],   :: w       Window width [npoints]
!    integer              :: np      Polynomial degree
!
! OUTPUT
!    real(wp), dim([:],:), optional  :: fs      Filtered signal
!    real(wp), dim([:],:), optional  :: ds      Signal derivative
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
!
!-------------------------------------------------------------------------------
! 1. 1-dimensional signal
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_sliding_poly_1d(t, s, w, np, fs, ds)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_utils

    ! 7.1 Declarations

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: T  ! Time
    REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: S  ! Signal samples (time)
    INTEGER,                  INTENT(in)  :: W  ! Window width [samples]
    INTEGER,                  INTENT(in)  :: NP ! Polynomial degree
    REAL(wp), DIMENSION(:),   INTENT(out) :: FS ! Filtered signal 
    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal derivative 

    INTEGER     :: NS       ! Number of signal samples
    INTEGER     :: WS       ! Odd sliding window
    INTEGER     :: i        ! Sample index
    INTEGER     :: Imin     ! Lower limit of sliding window
    INTEGER     :: Imax     ! Upper limit of sliding window

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K   ! Basic polynomials
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: A   ! Polynomial coefficients
    
    ! 7.2 Initialization

    ! 7.2.1 Determination of dimensions

    NS = SIZE(S)

    ! 7.2.2 Adjustment of odd filter width
    
    WS = MIN(2*(W/2) + 1, NS)

    ! 7.2.3. Memory allocation
    
    ALLOCATE(K(WS,0:NP))
    ALLOCATE(A(0:NP))

    ! 7.3 Sliding polynomial regression
    
    DO i = 1, NS

       ! 7.3.1 Positioning sliding window
       
       Imin = MAX(1,  i - WS/2)
       Imax = MIN(NS, i + WS/2)
       IF (Imin == 1) THEN
          Imax = WS
       END IF
       IF (Imax == NS) THEN
          Imin = NS - WS + 1
       END IF
       
       ! 7.3.2 Computation of basic polynomials

       CALL ropp_pp_init_polynomial(T(imin:imax)-T(i), K)

       ! 7.3.3 Sliding polynomial regression
       
       CALL ropp_pp_regression(K(:,:), S(imin:imax), A(:))
       
       FS(i) = A(0)
       
       IF (PRESENT(DS)) THEN
          DS(i) = A(1)
       END IF
              
    END DO
    
    ! 7.4 Clean up

    DEALLOCATE(K)
    DEALLOCATE(A)
  
  END SUBROUTINE ropp_pp_sliding_poly_1d

!-------------------------------------------------------------------------------
! 2. 2-dimensional signal
!-------------------------------------------------------------------------------

  SUBROUTINE ropp_pp_sliding_poly_2d(t, s, w, np, fs, ds)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_utils

    ! 8.1 Declarations

    IMPLICIT NONE

    REAL(wp), DIMENSION(:),   TARGET, INTENT(in)  :: T  ! Time
    REAL(wp), DIMENSION(:,:), TARGET, INTENT(in)  :: S  ! Signal samples [ch,t]
    INTEGER,                  INTENT(in)  :: W  ! Window width [samples]
    INTEGER,                  INTENT(in)  :: NP ! Polynomial degree
    REAL(wp), DIMENSION(:,:), INTENT(out) :: FS ! Filtered signal 
    REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: DS ! Signal derivative 

    INTEGER     :: NS       ! Number of signal samples
    INTEGER     :: NC       ! Number of channels
    INTEGER     :: WS       ! Odd sliding window
    INTEGER     :: IC       ! Channel index
    INTEGER     :: i        ! Sample index
    INTEGER     :: Imin     ! Lower limit of sliding window
    INTEGER     :: Imax     ! Upper limit of sliding window

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K   ! Basic polynomials
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: A   ! Polynomial coefficients

    ! 8.2 Initialization
    
    ! 8.2.1 Determination of dimensions

    NC = SIZE(S, 1)
    NS = SIZE(S, 2)
    
    ! 8.2.2 Adjustment of odd filter width
    
    WS = MIN(2*(W/2) + 1, NS)

    ! 8.2.3. Memory allocation
    
    ALLOCATE(K(WS,0:NP))
    ALLOCATE(A(NC,0:NP))

    ! 8.3 Sliding polynomial regression

    DO i = 1, NS

       ! 8.3.1 Positioning sliding window

       Imin = MAX(1,  i - WS/2)
       Imax = MIN(NS, i + WS/2)

       IF (Imin == 1) THEN
          Imax = WS
       END IF
       IF (Imax == NS) THEN
          Imin = NS - WS + 1
       END IF

       ! 8.3.2 Computation of basic polynomials

       CALL ropp_pp_init_polynomial(T(imin:imax)-T(i), K)

       ! 8.3.3 Sliding polynomial regression

       DO IC = 1, NC
          
          CALL ropp_pp_regression(K(:,:), S(IC,Imin:Imax), A(IC,:))
          
          FS(IC,i) = A(IC,0)

          IF (PRESENT(DS)) THEN
             DS(IC,i) = A(IC,1)
          END IF
          
       END DO
    END DO
     
    ! 8.4 Clean up
    
    DEALLOCATE(K)
    DEALLOCATE(A)

  END SUBROUTINE ropp_pp_sliding_poly_2d

!-------------------------------------------------------------------------------
! 2. 1-dimensional signal (vector window widths)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_sliding_poly_vec1d(t, s, w, np, fs, ds)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_utils
  
    ! 9.1 Declarations

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: T  ! Time
    REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: S  ! Signal samples (time)
    INTEGER,  DIMENSION(:),         INTENT(in) :: W  ! Window width (time)
    INTEGER,                  INTENT(in)  :: NP ! Polynomial degree
    REAL(wp), DIMENSION(:),           INTENT(out) :: FS ! Filtered signal 
    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal derivative 

    INTEGER     :: NS       ! Number of signal samples
    INTEGER     :: i        ! Sample index
    INTEGER     :: Imin     ! Lower limit of sliding window
    INTEGER     :: Imax     ! Upper limit of sliding window

    INTEGER,  DIMENSION(:), ALLOCATABLE   :: WS  ! Sliding window width
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K   ! Basic polynomials
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: A   ! Polynomial coefficients
    
    ! 9.2 Initialization
    
    ! 9.2.1 Determination of dimensions
    
    NS = SIZE(S)
    
    ! 9.2.3. Memory allocation
   
    ALLOCATE(A(0:NP))
    ALLOCATE(WS(NS))

    ! 9.2.4 Adjustment of odd filter width
    
    DO i = 1, NS
      WS(i) = MIN(2*(W(i)/2) + 1, NS)
    ENDDO
    
    ! 9.3 Sliding polynomial regression
    
    DO i = 1, NS

       ! 9.3.1 Positioning sliding window
       
       Imin = MAX(1,  i - WS(i)/2)
       Imax = MIN(NS, i + WS(i)/2)
       IF (Imin == 1) THEN
          Imax = WS(i)
       END IF
       IF (Imax == NS) THEN
          Imin = NS - WS(i) + 1
       END IF
       
       ! 7.3.2 Computation of basic polynomials
       ALLOCATE(K(WS(i),0:NP))

       CALL ropp_pp_init_polynomial(T(imin:imax)-T(i), K)

       ! 7.3.3 Sliding polynomial regression
       
       CALL ropp_pp_regression(K(:,:), S(imin:imax), A(:))
       
       FS(i) = A(0)
       
       IF (PRESENT(DS)) THEN
          DS(i) = A(1)
       END IF

       DEALLOCATE(K)
              
    END DO

    ! 7.4 Clean up

    DEALLOCATE(A)
    DEALLOCATE(WS)
  
  END SUBROUTINE ropp_pp_sliding_poly_vec1d
!
