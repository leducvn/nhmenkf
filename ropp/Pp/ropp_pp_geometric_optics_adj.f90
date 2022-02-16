! $Id: ropp_pp_geometric_optics.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_geometric_optics_adj(r_leo, v_leo, r_gns, v_gns, doppler,  &
                                        impact, bangle, impact_dd, impact_dr, &
                                        bangle_dd, bangle_dr)

!****s* GeometricOptics/ropp_pp_geometric_optics_adj *
!
! NAME
!    ropp_pp_geometric_optics_adj - Calculate bending angle and impact 
!                                   parameter from relative Doppler frequency 
!                                   shift. ADJOINT VERSION.
!                   
! SYNOPSIS
!    call ropp_pp_geometric_optics_adj(r_leo, v_leo, r_gns, v_gns, doppler,
!                                      impact, bangle, impact_dd, impact_dr, 
!                                      bangle_dd, bangle_dr)
! 
! DESCRIPTION
!    This routine is the adjoint of ropp_pp_geometric_optics, which calculates 
!    bending angle and impact parameter from relative Doppler frequency shift.
!    Iterative solution of the system of equations:
!               (c - (v_leo, U_leo))/(c - (v_gns, U_gns)) - 1 = doppler
!               [r_leo, U_leo] - [r_gns, U_gns] = 0
!               (U_leo, U_leo) = 1
!               (U_gns, U_gns) = 1
!   where U_leo and U_gns are the ray directions at the receiver and 
!   transmitter respectively.
!
! INPUTS
!    real(wp), dimension(:) :: r_leo     ! relative LEO position (ECI)
!    real(wp), dimension(:) :: v_leo     ! LEO velocity (ECI)  
!    real(wp), dimension(:) :: r_gns     ! relative GPS position (ECI)
!    real(wp), dimension(:) :: v_gns     ! GPS velocity (ECI)  
!    real(wp)               :: doppler   ! relative Doppler frequency shift
!
! OUTPUT
!    real(wp)               :: impact    ! impact parameter (m)
!    real(wp)               :: bangle    ! bending angle (rad)
!    real(wp)               :: impact_dd ! d(IP)/d(d) 
!    real(wp), dimension(:) :: impact_dr ! d(IP)/d(r_leo,r_gns)
!    real(wp)               :: bangle_dd ! d(BA)/d(d) 
!    real(wp), dimension(:) :: bangle_dr ! d(BA)/d(r_leo,r_gns) 
!
! NOTES
!
! REFERENCES
!   Vorob'ev and Krasil'nikova 1994, 
!   Estimation of the accuracy of the atmospheric refractive index recovery
!   from doppler Sshift measurements at frequencies used in the NAVSTAR system
!   Physics of the Atmosphere and Ocean (29) 602-609
!
!  Gorbunov, M.E. and Kornblueh, L. 2003,
!  Principles of variational assimilation of GNSS radio occultation data
!  Max PlancK Institute Report 350 
!  http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/max_scirep_350.pdf
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
! USE ropp_pp, not_this => ropp_pp_geometric_optics_adj
  USE ropp_pp
  USE ropp_utils

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo      ! LEO position (m) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo      ! LEO velocity (m/s) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns      ! GPS position (m) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns      ! LEO velocity (m/s) [ECI]
  REAL(wp),               INTENT(in)  :: doppler    ! Relative Doppler shift
  REAL(wp),               INTENT(out) :: impact     ! Impact parameter (m)
  REAL(wp),               INTENT(out) :: bangle     ! Bending angle (rad)
  REAL(wp),               INTENT(out) :: impact_dd  ! d(IP)/d(d) 
  REAL(wp), DIMENSION(6), INTENT(out) :: impact_dr  ! d(IP)/d(r_leo,r_gns)
  REAL(wp),               INTENT(out) :: bangle_dd  ! d(BA)/d(d)
  REAL(wp), DIMENSION(6), INTENT(out) :: bangle_dr  ! d(BA)/d(r_leo,r_gns)

  REAL(wp), DIMENSION(:), ALLOCATABLE :: U_leo   ! Ray direction at LEO
  REAL(wp), DIMENSION(:), ALLOCATABLE :: U_gns   ! Ray direction at GPS
  REAL(wp), DIMENSION(:), ALLOCATABLE :: TC      ! Temporary storage
  REAL(wp), DIMENSION(6)              :: DZ      ! Perturbation of (U_leo,U_gns)
  REAL(wp), DIMENSION(6)              :: DF      ! Discrepancy vector
  REAL(wp), DIMENSION(6,6)            :: A       ! Matrix system A*DZ = DF
  REAL(wp), DIMENSION(6,6)            :: AI      ! Inverted system matrix
  REAL(wp), DIMENSION(6)              :: E_URT   ! d(BA)/d(U_leo,U_gns)
  REAL(wp), DIMENSION(3)              :: P_UR    ! d(IP)/d(U_leo)
  REAL(wp), DIMENSION(3)              :: P_XR    ! d(IP)/d(r_leo)
  REAL(wp), DIMENSION(6)              :: URT_D   ! d(U_leo)/d(D)
  REAL(wp), DIMENSION(3,3)            :: BXR     ! Influence of X_leo
  REAL(wp), DIMENSION(3,3)            :: BXT     ! Influence of X_gns
  REAL(wp), DIMENSION(6,6)            :: B       ! Influence of X_leo and X_gns
  REAL(wp), DIMENSION(6,6)            :: URT_XRT ! d(U_leo,U_gns)/d(X_leo,X_gns)
  INTEGER                             :: i, j, k ! Counters
  
  INTEGER,  PARAMETER :: nit = 10                ! Number of iterations
  REAL(wp), PARAMETER :: tensor(3,3,3) =   &     ! Antisymmetrical tensor
                             RESHAPE((/0,  0,  0,  0,  0,  1,  0, -1,  0,   &
                                       0,  0, -1,  0,  0,  0,  1,  0,  0,   &
                                       0,  1,  0, -1,  0,  0,  0,  0,  0/), &
                                       Shape = (/3,3,3/), Order = (/3,2,1/))

!-------------------------------------------------------------------------------
! 2. Array allocation
!-------------------------------------------------------------------------------

  ALLOCATE(U_leo(SIZE(r_leo)))
  ALLOCATE(U_gns(SIZE(r_gns)))
  ALLOCATE(TC(SIZE(r_gns)))

!-------------------------------------------------------------------------------
! 3. Initial approximation - straight line r_gns --> r_leo
!-------------------------------------------------------------------------------

  U_leo(:) = (r_leo(:) - r_gns(:))/SQRT(SUM((r_leo(:) - r_gns(:))**2))
  U_gns(:) = U_leo(:)

!-------------------------------------------------------------------------------
! 4. Iterative solution
!-------------------------------------------------------------------------------

  DO k = 1, Nit
     
     ! 4.1 Calculation of discrepancy

     DF(1)   = doppler - (DOT_PRODUCT(v_gns,U_gns)-DOT_PRODUCT(v_leo,U_leo)) / &
                             (c_light - DOT_PRODUCT(v_gns,U_gns))
     DF(2:4) = vector_product(r_gns, U_gns) - vector_product(r_leo, U_leo)
     DF(5)   = 1.0_wp - DOT_PRODUCT(U_leo,U_leo)
     DF(6)   = 1.0_wp - DOT_PRODUCT(U_gns,U_gns)
     
     ! 4.2 Calculation of matrix of linearized system
     
     A(:,:)   = 0.0_wp
     A(1,1:3) = -v_leo/(c_light - DOT_PRODUCT(v_gns,U_gns))
     A(1,4:6) =  v_gns * (c_light - DOT_PRODUCT(v_leo,U_leo)) /      &
                            (C_Light - DOT_PRODUCT(v_gns,U_gns))**2
     DO i=1,3
        DO j=1,3
           A(i+1,j)   =  SUM(tensor(i,:,j)*r_leo(:))
           A(i+1,j+3) = -SUM(tensor(i,:,j)*r_gns(:))
        END DO
     END DO
     A(5,1:3) = 2.0_wp*U_leo
     A(6,4:6) = 2.0_wp*U_gns
     
     ! 4.3 Solve liearized system

     CALL ropp_pp_invert_matrix(A, AI)
     DZ = MATMUL(AI, DF)

     ! 4.4 Calculation of next approximation
     
     TC = DZ(1:3)
     U_leo = U_leo + TC
     TC = DZ(4:6)
     U_gns = U_gns + TC

  ENDDO

!-------------------------------------------------------------------------------
! 5. Calculate bending angle and impact parameter
!-------------------------------------------------------------------------------

  ! 5.1 Bending angle 

!  bangle = vector_angle(U_gns, U_leo, vector_product(r_gns, r_leo))
  
  CALL vector_angle_adj(U_gns, U_leo, vector_product(r_gns, r_leo), bangle,   &
                        E_URT(4:6), E_URT(1:3))

  ! 5.2 impact = r_L sin(phi_l) = r_G sin(phi_G)
  
  impact = SQRT(DOT_PRODUCT(r_leo, r_leo)) * SIN(vector_angle(r_leo, U_leo))

!-------------------------------------------------------------------------------
! 6. Adjoint calculations
!-------------------------------------------------------------------------------

  ! 6.1 Derivatives of U_leo and U_gns

  URT_D(:) = AI(:,1)

  DO i=1,3
     DO j=1,3
        BXR(i,j) = -SUM(tensor(i,j,:)*U_leo(:))
        BXT(i,j) =  SUM(tensor(i,j,:)*U_gns(:))
     ENDDO
  ENDDO
  
  B(:,:)       = 0
  B(2:4,1:3)   = BXR(:,:)
  B(2:4,4:6)   = BXT(:,:)
  URT_XRT(:,:) = MATMUL(AI(:,:),B(:,:))
  
  ! 6.2 Derivatives of bending angle

  bangle_dd = SUM(E_URT(:)*URT_D(:))
  bangle_dr = MATMUL(E_URT(:),URT_XRT(:,:))

  ! 6.3 Derivatives of impact parameter
  
  P_UR(:)    = -(vector_product(r_leo,(vector_product(r_leo,U_leo))))/impact
  P_XR(:)    = -(vector_product(U_leo,(vector_product(U_leo,r_leo))))/impact
  
  impact_dd  = SUM(P_UR(:)*URT_D(1:3))
  
  impact_dr(:)   = MATMUL(P_UR(:),URT_XRT(1:3,:))
  impact_dr(1:3) = impact_dr(1:3) + P_XR(:)

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(U_leo)
  DEALLOCATE(U_gns)
  DEALLOCATE(TC)

CONTAINS
  
!-------------------------------------------------------------------------------
! 8. Adjoint code for angle between cartesian vectors
!-------------------------------------------------------------------------------

  SUBROUTINE vector_angle_adj(X, Y, A, angle,AXY_X, AXY_Y)

    USE typesizes,  ONLY: wp => EightByteReal
    USE ropp_utils, ONLY: vector_product
    
    REAL(wp), DIMENSION(3), INTENT(in)  :: X          ! Cartesian vector
    REAL(wp), DIMENSION(3), INTENT(in)  :: Y          ! Cartesian vector
    REAL(wp), DIMENSION(3), INTENT(in)  :: A          ! Orientation axis
    REAL(wp),               INTENT(out) :: angle      ! Angle between X and Y
    REAL(wp), DIMENSION(3), INTENT(out) :: AXY_X      ! d(AXY)/d(X)
    REAL(wp), DIMENSION(3), INTENT(out) :: AXY_Y      ! d(AXY)/d(Y)
    
    REAL(wp), DIMENSION(3) :: n, alpha, beta, gamma
    REAL(wp)               :: nn
    REAL(wp)               :: ag, bg

    nn = DOT_PRODUCT(A, A)
    
    IF (nn == 0) THEN
       angle = 0.0_wp
       axy_x = 0.0_wp
       axy_y = 0.0_wp
    ELSE
       n = A/SQRT(nn)
       alpha = vector_product(n, X)
       
       beta = X - DOT_PRODUCT(n, X) * n
       gamma = Y - DOT_PRODUCT(n, Y) * n
       
       ag = DOT_PRODUCT(alpha,gamma)
       bg = DOT_PRODUCT(beta,gamma)

       angle = ATAN2(ag, bg)

       AXY_X = -(ag*gamma - bg*(vector_product(Y,n)))/(ag**2 + bg**2)
       AXY_Y = -(ag*beta  - bg*alpha) / (ag**2 + bg**2)

    ENDIF
    
  END SUBROUTINE vector_angle_adj

END SUBROUTINE ropp_pp_geometric_optics_adj
