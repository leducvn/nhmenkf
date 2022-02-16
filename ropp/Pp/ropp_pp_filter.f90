! $Id: ropp_pp_filter.f90 2228 2009-09-01 15:36:10Z frhl $

!****s* FFT/ropp_pp_filter *
!
! NAME
!    ropp_pp_filter - Filtering and differentitation of a signal
!
! SYNOPSIS
!    call ropp_pp_filter(dt, s, w, nd, fs, ds)
!
! DESCRIPTION
!    Optimal solution of integral equation
!
! INPUTS
!    real(wp)             :: dt      Time step
!    real(wp), dim([:],:) :: s       Signal samples ([channel],time)
!    integer              :: w       Window width [npoints]
!    integer              :: nd      Number of points for differentiation 
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

!-------------------------------------------------------------------------------
! 1. 1d Signal Filter
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_filter_1d(DT, S, W, ND, FS, DS) 
  
  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_pp_utils, ONLY: ropp_pp_quasi_invert, ropp_pp_matmul

  IMPLICIT NONE
  
    ! 1.1 Declarations
 
    REAL(wp),                 INTENT(in)  :: DT ! Time step
    REAL(wp), DIMENSION(:),   INTENT(in)  :: S  ! Signal samples (time)
    INTEGER,                  INTENT(in)  :: W  ! Window width [points]
    INTEGER,                  INTENT(in)  :: ND ! Number of points for 
                                                ! differentiation (must be odd)
    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: FS ! Filtered signal 
    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal derivative 

    INTEGER         :: N           ! Number of samples
    INTEGER         :: i           ! Sample index
    INTEGER         :: m, j        ! Dimension indices
    INTEGER         :: Imin, Imax  ! Differentiation index limits
    INTEGER         :: I0          ! Index of middle point

    INTEGER,  DIMENSION(nd)     :: ID  ! Point indices for differentiation
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: GL  ! Lower-border mean and deriv
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: GU  ! Upper-border mean and deriv
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KB  ! Border matrix of integral eq
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QB  ! Quasi-inverse of KB
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: SX  ! Extended data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KD  ! Matrix of integral equation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QD  ! Quasi-inverse of KD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QI  ! Internal part of QD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QO  ! Outside part of QD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KO  ! Inverse of QO

    ! 1.1 Initialization

    N  = SIZE(S)

    ALLOCATE(SX(1-(ND/2)*W:N+(ND/2)*W))
    ALLOCATE(KB(0:(ND/2)*W, 0:1+ND/2))
    ALLOCATE(QB(0:1+ND/2, 0:(ND/2)*W))
    ALLOCATE(GL(0:1+ND/2))
    ALLOCATE(GU(0:1+ND/2))
    ALLOCATE(KD(0:2*(ND/2)*W, 0:ND))
    ALLOCATE(QD(0:ND, 0:2*(ND/2)*W))
    ALLOCATE(QI(0:1, 0:(ND/2)*W))
    ALLOCATE(QO(0:1, (ND/2)*W))
    ALLOCATE(KO((ND/2)*W, 0:1))

    ! 1.2 Lower border filtering

    ! 1.2.1 Calculation of filtering grid
    
    i = 1
    DO j=1,1+ND/2
       ID(j) = i + (j-1)*W
    END DO

    Imin = ID(1)
    Imax = ID(1+ND/2)
    I0   = ID(1)
    
    ! 1.2.2. Forming integral equation matrix
    
    KB(:,0) = 1.0_wp
    
    DO j=1,1+ND/2 
       DO m=Imin,Imax
          KB(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO
    
    ! 1.2.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KB, QB)
    
    GL(:) = ropp_pp_MATMUL(QB, S(Imin:Imax))
        
    ! 1.3 Upper border filtering

    ! 1.3.1. Calculation of filtering grid
    
    i = N
    DO j=1,1+ND/2
       ID(j) = i + (j-(1+ND/2))*W
    END DO
    
    Imin = ID(1)
    Imax = ID(1+ND/2)
    I0   = ID(1+ND/2)
        
    ! 1.3.2. Forming integral equation matrix
    
    KB(:,0) = 1

    DO j=1,1+ND/2
       DO m=Imin,Imax
          KB(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO
    
    !1.3.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KB, QB)
    
    GU(:) = MATMUL(QB, S(Imin:Imax))
    
    ! 1.4 Calculation of filtering matrix

    ! 1.4.1. Calculation of filtering grid

    i = N/2
    DO j=1,ND
       ID(j) = i + (j-(1+ND/2))*W
    END DO
    
    Imin = ID(1)
    Imax = ID(ND)
    I0   = ID(1+ND/2)
    
    ! 1.4.2. Forming integral equation matrix
    
    KD(:,0) = 1
    
    DO j=1,ND
       DO m=Imin,Imax
          KD(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO

    ! 1.4.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KD, QD)

    ! 1.5 Extrapolation

    ! 1.5.1. Copying the internal part
    
    SX(1:N) = S(1:N)
    
    ! 1.5.2. Lower border
    
    QI(0,:) = QD(0,      (ND/2)*W:2*(ND/2)*W)
    QI(1,:) = QD(1+ND/2, (ND/2)*W:2*(ND/2)*W)
    
    QO(0,:) = QD(0,      0:(ND/2)*W-1)
    QO(1,:) = QD(1+ND/2, 0:(ND/2)*W-1)
    
    CALL ropp_pp_quasi_invert(QO, KO)

    SX(1-(ND/2)*W:0) = &
            MATMUL(KO, GL(0:1) - MATMUL(QI, S(1:1+(ND/2)*W)))
    
    ! 1.5.3. Upper border
    
    QI(0,:) = QD(0,      0:(ND/2)*W)
    QI(1,:) = QD(1+ND/2, 0:(ND/2)*W)
    
    QO(0,:) = QD(0,      (ND/2)*W+1:2*(ND/2)*W)
    QO(1,:) = QD(1+ND/2, (ND/2)*W+1:2*(ND/2)*W)
    
    CALL ropp_pp_quasi_invert(QO, KO)
    
    GU(1) = GU(1+ND/2)
    
    SX(N+1:N+(ND/2)*W) = &
            MATMUL(KO, GU(0:1) - MATMUL(QI, S(N-(ND/2)*W:N)))
    
    ! 1.6 Filtering

    DO i=1,N
       
       ! 1.6.1. Determination of point interval for filtering
       
       Imin = i - (ND/2)*W
       Imax = i + (ND/2)*W
       
       ! 1.6.2. Filtering signal
       
       IF (PRESENT(FS)) THEN
          FS(i) = SUM(QD(0,:)*SX(Imin:Imax))
       END IF
       
       ! 1.6.3. Differentiating signal

       IF (PRESENT(DS)) THEN
          DS(i) = SUM(QD(1+ND/2,:)*SX(Imin:Imax))
       END IF
    END DO
    
    DEALLOCATE(SX, KB, QB, GL, GU, KD, QD, QI, QO, KO)
    
  CONTAINS

    ! Integral of cap function of width 2*TW
    ! located at T0, from T0 to T

    REAL(wp) FUNCTION Cap_Integral(T, T0, TW)
      
      REAL(wp) :: T   ! Time
      REAL(wp) :: T0  ! Cap function location
      REAL(wp) :: TW  ! Cap function half width
      REAL(wp) :: TM  ! Absolute time interval limited
      
      TM = MIN(TW, ABS(T-T0))
      
      Cap_Integral = SIGN(1.0_wp, T-T0)*  &
           TM*(1-TM/(2*TW))
      
    END FUNCTION Cap_Integral
    
  END SUBROUTINE ropp_pp_filter_1d


!-------------------------------------------------------------------------------
! 2. 2d Signal Filter
!-------------------------------------------------------------------------------

  SUBROUTINE ropp_pp_filter_2d(DT, S, W, ND, FS, DS)        

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_utils, ONLY: ropp_pp_quasi_invert, ropp_pp_matmul

    IMPLICIT NONE

    ! 1.1 Declarations
 
    REAL(wp),                 INTENT(in)  :: DT ! Time step
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: S  ! Signal samples (channel,time)
    INTEGER,                  INTENT(in)  :: W  ! Window width [points]
    INTEGER,                  INTENT(in)  :: ND ! Number of points for 
                                                ! differentiation (must be odd)
    
    REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: FS ! Filtered signal 
    REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: DS ! Signal derivative 

    INTEGER         :: N           ! Number of samples
    INTEGER         :: NC          ! Nunber of channels
    INTEGER         :: IC          ! Channel number
    INTEGER         :: i           ! Sample index
    INTEGER         :: m, j        ! Dimension indices
    INTEGER         :: Imin, Imax  ! Differentiation index limits
    INTEGER         :: I0          ! Index of middle point

    INTEGER,  DIMENSION(nd)     :: ID  ! Point indices for differentiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: GL  ! Lower-border mean and deriv
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: GU  ! Upper-border mean and deriv
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KB  ! Border matrix of integral eq
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QB  ! Quasi-inverse of KB
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: SX  ! Extended data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KD  ! Matrix of integral equation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QD  ! Quasi-inverse of KD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QI  ! Internal part of QD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: QO  ! Outside part of QD
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KO  ! Inverse of QO

    ! 1.1 Initialization
    
    NC = SIZE(S, 1)
    N  = SIZE(S, 2)

    ALLOCATE(SX(NC, 1-(ND/2)*W:N+(ND/2)*W))
    ALLOCATE(KB(0:(ND/2)*W, 0:1+ND/2))
    ALLOCATE(QB(0:1+ND/2, 0:(ND/2)*W))
    ALLOCATE(GL(NC, 0:1+ND/2))
    ALLOCATE(GU(NC, 0:1+ND/2))
    ALLOCATE(KD(0:2*(ND/2)*W, 0:ND))
    ALLOCATE(QD(0:ND, 0:2*(ND/2)*W))
    ALLOCATE(QI(0:1, 0:(ND/2)*W))
    ALLOCATE(QO(0:1, (ND/2)*W))
    ALLOCATE(KO((ND/2)*W, 0:1))

    ! 1.2 Lower border filtering

    ! 1.2.1 Calculation of filtering grid
    
    i = 1
    DO j=1,1+ND/2
       ID(j) = i + (j-1)*W
    END DO

    Imin = ID(1)
    Imax = ID(1+ND/2)
    I0   = ID(1)
    
    ! 1.2.2. Forming integral equation matrix
    
    KB(:,0) = 1.0_wp
    
    DO j=1,1+ND/2 
       DO m=Imin,Imax
          KB(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO
    
    ! 1.2.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KB, QB)
    
    DO IC=1,NC
       GL(IC,:) = ropp_pp_MATMUL(QB, S(IC,Imin:Imax))
    END DO
    
    ! 1.3 Upper border filtering

    ! 1.3.1. Calculation of filtering grid
    
    i = N
    DO j=1,1+ND/2
       ID(j) = i + (j-(1+ND/2))*W
    END DO
    
    Imin = ID(1)
    Imax = ID(1+ND/2)
    I0   = ID(1+ND/2)
        
    ! 1.3.2. Forming integral equation matrix
    
    KB(:,0) = 1

    DO j=1,1+ND/2
       DO m=Imin,Imax
          KB(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO
    
    !1.3.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KB, QB)
    
    DO IC=1,NC
       GU(IC,:) = ropp_pp_MATMUL(QB, S(IC,Imin:Imax))
    END DO

    ! 1.4 Calculation of filtering matrix

    ! 1.4.1. Calculation of filtering grid

    i = N/2
    DO j=1,ND
       ID(j) = i + (j-(1+ND/2))*W
    END DO
    
    Imin = ID(1)
    Imax = ID(ND)
    I0   = ID(1+ND/2)
    
    ! 1.4.2. Forming integral equation matrix
    
    KD(:,0) = 1
    
    DO j=1,ND
       DO m=Imin,Imax
          KD(m-Imin, j) = &
               Cap_Integral(m*DT,  ID(j)*DT, W*DT) -  &
               Cap_Integral(I0*DT, ID(j)*DT, W*DT)
       END DO
    END DO

    ! 1.4.3. Quasi-solution of integral equation
    
    CALL ropp_pp_quasi_invert(KD, QD)

    ! 1.5 Extrapolation

    ! 1.5.1. Copying the internal part
    
    SX(:,1:N) = S(:,1:N)

    ! 1.5.2. Lower border
    
    QI(0,:) = QD(0,      (ND/2)*W:2*(ND/2)*W)
    QI(1,:) = QD(1+ND/2, (ND/2)*W:2*(ND/2)*W)
    
    QO(0,:) = QD(0,      0:(ND/2)*W-1)
    QO(1,:) = QD(1+ND/2, 0:(ND/2)*W-1)
    
    CALL ropp_pp_quasi_invert(QO, KO)

    DO IC=1,NC
       SX(IC, 1-(ND/2)*W:0) = &
            MATMUL(KO, GL(IC,0:1) - MATMUL(QI, S(IC,1:1+(ND/2)*W)))
    END DO

    ! 1.5.3. Upper border
    
    QI(0,:) = QD(0,      0:(ND/2)*W)
    QI(1,:) = QD(1+ND/2, 0:(ND/2)*W)
    
    QO(0,:) = QD(0,      (ND/2)*W+1:2*(ND/2)*W)
    QO(1,:) = QD(1+ND/2, (ND/2)*W+1:2*(ND/2)*W)

    CALL ropp_pp_quasi_invert(QO, KO)

    GU(:,1) = GU(:,1+ND/2)
    
    DO IC=1,NC
       SX(IC, N+1:N+(ND/2)*W) = &
            MATMUL(KO, GU(IC,0:1) - MATMUL(QI, S(IC,N-(ND/2)*W:N)))
    END DO

    ! 1.6 Filtering

    DO i=1,N
       
       ! 1.6.1. Determination of point interval for filtering
       
       Imin = i - (ND/2)*W
       Imax = i + (ND/2)*W
       
       ! 1.6.2. Filtering signal
       
       IF (PRESENT(FS)) THEN
          DO IC=1,NC
             FS(IC,i) = SUM(QD(0,:)*SX(IC,Imin:Imax))
          END DO
       END IF
       
       ! 1.6.3. Differentiating signal

       IF (PRESENT(DS)) THEN
          DO IC=1,NC
            if (maxval(SX(IC,Imin:imax)) > 100000000000.0_wp) STOP
             DS(IC,i) = SUM(QD(1+ND/2,:)*SX(IC,Imin:Imax))
          END DO
       END IF
    END DO
    
    DEALLOCATE(SX, KB, QB, GL, GU, KD, QD, QI, QO, KO)
    
  CONTAINS

    ! Integral of cap function of width 2*TW
    ! located at T0, from T0 to T

    REAL(wp) FUNCTION Cap_Integral(T, T0, TW)
      
      REAL(wp) :: T   ! Time
      REAL(wp) :: T0  ! Cap function location
      REAL(wp) :: TW  ! Cap function half width
      REAL(wp) :: TM  ! Absolute time interval limited
      
      TM = MIN(TW, ABS(T-T0))
      
      Cap_Integral = SIGN(1.0_wp, T-T0)*  &
           TM*(1-TM/(2*TW))
      
    END FUNCTION Cap_Integral
    
  END SUBROUTINE ropp_pp_filter_2d
