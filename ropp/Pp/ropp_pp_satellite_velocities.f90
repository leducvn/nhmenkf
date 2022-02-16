! $Id: ropp_pp_satellite_velocities.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,   &
                                        xgns, vgns, abl, abg)

!****s* Preprocessing/ropp_pp_satellite_velocities *
!
! NAME
!    ropp_pp_satellite_velocities - Calculate satellite velocities by
!                                   polynomial regression
!                   
! SYNOPSIS
!    call ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,    &
!                                      xgns, vgns, abl, abg)
! 
! DESCRIPTION
!    This routine calculates satellite velocities by polynomial regression
!
! INPUTS
!    real(wp), dimension(:)   :: time    ! time of sample (s)
!    real(wp), dimension(:,:) :: r_leo   ! cartesian LEO position (ECI or ECF)
!    real(wp), dimension(:,:) :: r_gns   ! cartesian GPS position (ECI or ECF)
!
! OUTPUTS
!    real(wp), dimension(:,:) :: xleo    ! LEO positions from regression
!    real(wp), dimension(:,:) :: vleo    ! LEO velocities from regression
!    real(wp), dimension(:,:) :: xgns    ! GPS positions from regression
!    real(wp), dimension(:,:) :: vgns    ! GPS velocities from regression
!    real(wp), dimension(:,:) :: abl     ! LEO regression coefficients (opt)
!    real(wp), dimension(:,:) :: abg     ! GPS regression coefficients (opt)
!
! NOTES
!
! REFERENCES
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
! USE ropp_pp, not_this => ropp_pp_satellite_velocities

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)    :: time    ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo   ! Cartesian LEO position (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns   ! Cartesian GPS position (m)
  REAL(wp), DIMENSION(:,:), INTENT(out) :: xleo    ! LEO position vector (m)
  REAL(wp), DIMENSION(:,:), INTENT(out) :: vleo    ! LEO velocity vector (m)
  REAL(wp), DIMENSION(:,:), INTENT(out) :: xgns    ! GPS position vector (m)
  REAL(wp), DIMENSION(:,:), INTENT(out) :: vgns    ! GPS velocity vector (m)
  REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abl ! LEO regression coeffs
  REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abg ! GPS regression coeffs

  INTEGER, PARAMETER :: nv = 5  ! polynomial degree for calculation of velocity
  INTEGER :: n, i, m
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: t_norm      ! Normalised time
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KV          ! Regress matrix velocity
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_rleo  ! Coefficients for RLEO
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_rgns  ! Coefficients for RGNS

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  n = SIZE(time)

  ALLOCATE(t_norm(n))
  ALLOCATE(KV(n, 0:nv))

!-------------------------------------------------------------------------------
! 3. Re-normalise time
!-------------------------------------------------------------------------------

  t_norm(:) = (time(:) - time(1))/(time(n) - time(1))

!-------------------------------------------------------------------------------
! 4. Calculate regression coefficients
!-------------------------------------------------------------------------------

  CALL ropp_pp_init_polynomial(t_norm, KV)

  ! 4.1 Perform regression on positions
  DO m=1,3
     CALL ropp_pp_regression(KV, r_leo(:,m), coeff_rleo(:,m))
     CALL ropp_pp_regression(KV, r_gns(:,m), coeff_rgns(:,m))
  ENDDO

  ! 4.2 Perform regression on residual positions to gain higher accuracy
  DO m=1,3
    CALL ropp_pp_residual_regression(KV, t_norm, r_leo(:,m), coeff_rleo(:,m))
    CALL ropp_pp_residual_regression(KV, t_norm, r_gns(:,m), coeff_rgns(:,m))
  ENDDO

!-------------------------------------------------------------------------------
! 5. Regression
!-------------------------------------------------------------------------------

  DO i=1,n
     DO m=1,3
        CALL ropp_pp_polynomial(coeff_rleo(:,m),t_norm(i),xleo(i,m),vleo(i,m))
        CALL ropp_pp_polynomial(coeff_rgns(:,m),t_norm(i),xgns(i,m),vgns(i,m))
     ENDDO
     
     vleo(i,:) = vleo(i,:)/ (time(n) - time(1))
     vgns(i,:) = vgns(i,:)/ (time(n) - time(1))
     
  ENDDO

  IF (PRESENT(abl)) THEN
     abl(:,:) = coeff_rleo(:,:)
  ENDIF
  
  IF (PRESENT(abg)) THEN 
     abg(:,:) = coeff_rgns(:,:)
  ENDIF
  
!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(t_norm)
  DEALLOCATE(KV)

END SUBROUTINE ropp_pp_satellite_velocities
