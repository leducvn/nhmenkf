! $Id: ropp_pp_geometric_optics.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_geometric_optics(r_leo, v_leo, r_gns, v_gns, doppler,    &
                                    impact, bangle)

!****s* GeometricOptics/ropp_pp_geometric_optics *
!
! NAME
!    ropp_pp_geometric_optics - Calculate bending angle and impact parameter 
!                               from relative Doppler frequency shift.
!                   
! SYNOPSIS
!    call ropp_pp_geometric_optics(r_leo, v_leo, r_gns, v_gns, doppler,
!                                  impact, bangle)
! 
! DESCRIPTION
!    This routine calculates bending angle and impact parameter from relative
!    Doppler frequency shift.
!    Iterative solution of the system of equations:
!               (c - (v_leo, U_leo))/(c - (v_gns, U_gns)) - 1 = doppler
!               [r_leo, U_leo] - [r_gns, U_gns] = 0
!               (U_leo, U_leo) = 1
!               (U_gns, U_gns) = 1
!   where U_leo and U_gns are the ray directions at the receiver and 
!   transmitter respectively.
!
! INPUTS
!    real(wp), dimension(:) :: r_leo     ! relative LEO position (m) [ECI]
!    real(wp), dimension(:) :: v_leo     ! LEO velocity (m/s) [ECI]  
!    real(wp), dimension(:) :: r_gns     ! relative GPS position (m) [ECI]
!    real(wp), dimension(:) :: v_gns     ! GPS velocity (m/s) [ECI]  
!    real(wp)               :: doppler   ! relative Doppler frequency shift
!
! OUTPUT
!    real(wp)               :: impact    ! impact parameter (m)
!    real(wp)               :: bangle    ! bending angle (rad)
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
! USE ropp_pp, not_this => ropp_pp_geometric_optics
  USE ropp_pp
  USE ropp_utils

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo   ! LEO position (m) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo   ! LEO velocity (m/s) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns   ! GPS position (m) [ECI]
  REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns   ! GPS velocity (m/s) [ECI]
  REAL(wp),               INTENT(in)  :: doppler ! Relative Doppler shift
  REAL(wp),               INTENT(out) :: impact  ! Impact parameter (m)
  REAL(wp),               INTENT(out) :: bangle  ! Bending angle (rad)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: U_leo   ! Ray direction at receiver
  REAL(wp), DIMENSION(:), ALLOCATABLE :: U_gns   ! Ray direction at transmitter
  REAL(wp), DIMENSION(:), ALLOCATABLE :: TC      ! Temporary storage
  REAL(wp), DIMENSION(6)              :: DZ      ! Perturbation of (U_leo,U_gns)
  REAL(wp), DIMENSION(6)              :: DF      ! Discrepancy vector
  REAL(wp), DIMENSION(6,6)            :: A       ! Matrix system A*DZ = DF
  REAL(wp), DIMENSION(6,6)            :: AI      ! Inverted system matrix
  INTEGER                             :: i, j, k ! Counters
  
  INTEGER,  PARAMETER :: Nit = 5                 ! Number of iterations
  REAL(wp), PARAMETER ::       &                 ! Antisymmetrical tensor
              tensor(3,3,3) = RESHAPE((/0,  0,  0,  0,  0,  1,  0, -1,  0,   &
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
     A(1,4:6) =  v_gns * (c_light - DOT_PRODUCT(v_leo,U_leo)) /        &
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

  bangle = vector_angle(U_gns, U_leo, vector_product(r_gns, r_leo))
  
  ! 5.2 impact = r_L sin(phi_l) = r_G sin(phi_G)
  
  impact = SQRT(DOT_PRODUCT(r_leo, r_leo)) * SIN(vector_angle(r_leo, U_leo))

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(U_leo)
  DEALLOCATE(U_gns)
  DEALLOCATE(TC)

END SUBROUTINE ropp_pp_geometric_optics
