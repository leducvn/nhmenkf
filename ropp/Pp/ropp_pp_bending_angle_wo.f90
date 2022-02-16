! $Id: ropp_pp_bending_angle_wo.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_bending_angle_wo(time, r_leo, r_gns, r_coc, roc,         &
                                    phase_L1, phase_L2, snr_L1, snr_L2, w_ls, &
                                    w_smooth, w_low, hmax, filter, opt_DL2,   &
                                    cff, dsh,                                 &
                                    impact_L1, bangle_L1, ba_sigma_L1,        &
                                    impact_L2, bangle_L2, ba_sigma_L2, diag) 

!****s* WaveOptics/ropp_pp_bending_angle_wo *
!
! NAME
!    ropp_pp_bending_angle_wo - Calculate L1 and L2 bending angle profiles from
!                               occultation data by WAVE OPTICS
!                   
! SYNOPSIS
!    call ropp_pp_bending_angle_wo(time, r_leo, r_gns, r_coc, phase_L1,
!                                  phase_L2, snr_L1, snr_L2, w_ls,
!                                  w_smooth, w_low, hmax, filter, opt_DL2,
!                                  cff, dsh, 
!                                  impact_L1, bangle_L1, ba_sigma_L1,
!                                  impact_L2, bangle_L2, ba_sigma_L2, diag)
! 
! DESCRIPTION
!    This routine calculates L1 and L2 bending angles using the CT2 algorithm.
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! Relative time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! LEO coordinates (m) (ECI or ECF)
!    real(wp), dimension(:,:) :: r_gns     ! GPS coordinates (m) (ECI or ECF)
!    real(wp), dimension(:)   :: r_coc     ! Centre curvature coords (m)
!    real(wp)                 :: roc       ! Radius curvature (m)
!    integer                  :: w_ls      ! Large-scale smoothing (points)
!    integer                  :: w_smooth  ! Smoothing window above 7km (point)
!    integer                  :: w_low     ! Smoothing window below 7km (point)
!    real(wp)                 :: hmax      ! Maximum height for WO (m)
!    character(len=*)         :: filter    ! Filter method ('optest'/'slpoly')
!    logical                  :: opt_DL2   ! Degraded L2 flag
!    integer                  :: cff       ! Complex filtering flag
!    real(wp)                 :: dsh       ! Shadow border width (m)
!    real(wp), dimension(:)   :: phase_L1  ! L1 excess phase (m)
!    real(wp), dimension(:)   :: phase_L2  ! L2 excess phase (m)
!    real(wp), dimension(:)   :: snr_L1    ! L1 amplitude
!    real(wp), dimension(:)   :: snr_L2    ! L2 amplitude
!    real(wp), dimension(:)   :: impact_L1 ! L1 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L1 ! L1 bending angles (rad)
!    real(wp), dimension(:)   :: impact_L2 ! L2 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L2 ! L2 bending angles (rad)

!
! OUTPUT
!    real(wp), dimension(:)   :: impact_L1 ! L1 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L1 ! L1 bending angles (rad)
!    real(wp), dimension(:)   :: ba_sigma_L1 ! L1 bending angle error (rad)
!    real(wp), dimension(:)   :: impact_L2 ! L2 impact parameters (m)
!    real(wp), dimension(:)   :: bangle_L2 ! L2 bending angles (rad)
!    real(wp), dimension(:)   :: ba_sigma_L2 ! L2 bending angle error (rad)
!    type(PPdiag), optional   :: diag      ! Additional output diagnostic strt
! 
! NOTES
!
! REFERENCES
!   Gorbunov M.E. and Lauritsen K.B. 2004
!   Analysis of wave fields by Fourier integral operators and their application
!   for radio occultations
!   Radio Science (39) RS4010
!
!   Gorbunov M.E., Lauritsen K.B., Rodin A., Tomassini M., Kornblueh L. 2005 
!   Analysis of the CHAMP experimental data on radio-occultation sounding of
!   the Earth's atmosphere.
!   Izvestiya Atmospheric and Oceanic Physics (41) 726-740.
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
! USE ropp_pp, not_this => ropp_pp_bending_angle_wo
  USE ropp_pp_types, ONLY: PPDiag
  USE messages

  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)    :: time      ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_leo     ! LEO coordinates (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_gns     ! GPS coordinates (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: r_coc     ! Centre curvature (m)
  REAL(wp),                 INTENT(in)    :: roc       ! Radius of curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: phase_L1  ! L1 excess phase (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: phase_L2  ! L2 excess phase (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: snr_L1    ! L1 amplitude
  REAL(wp), DIMENSION(:),   INTENT(in)    :: snr_L2    ! L2 amplitude
  INTEGER,                  INTENT(in)    :: w_ls      ! Large-scale smoothing
  INTEGER,                  INTENT(in)    :: w_smooth  ! Smooth window > 7km
  INTEGER,                  INTENT(in)    :: w_low     ! Smooth window < 7km
  REAL(wp),                 INTENT(in)    :: hmax      ! Max height for WO (m)
  CHARACTER(len=*),         INTENT(in)    :: filter    ! Filter method
  LOGICAL,                  INTENT(in)    :: opt_DL2   ! Degraded L2 flag
  INTEGER,                  INTENT(in)    :: cff       ! Complex filter flag
  REAL(wp),                 INTENT(in)    :: dsh       ! Shadow border width (m)
  REAL(wp), DIMENSION(:),   INTENT(inout) :: impact_L1 ! L1 impact parameter (m)
  REAL(wp), DIMENSION(:),   INTENT(inout) :: bangle_L1 ! L1 bending angles (rad)
  REAL(wp), DIMENSION(:),   INTENT(out)   :: ba_sigma_L1 ! L1 bangle std dev
  REAL(wp), DIMENSION(:),   INTENT(inout) :: impact_L2 ! L2 impact parameter (m)
  REAL(wp), DIMENSION(:),   INTENT(inout) :: bangle_L2 ! L2 bending angles (rad)
  REAL(wp), DIMENSION(:),   INTENT(out)   :: ba_sigma_L2 ! L2 bangle std dev
  TYPE(PPdiag), OPTIONAL,   INTENT(inout) :: diag      ! Additional diagnostics

  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: AP      
  REAL(wp), DIMENSION(:),   ALLOCATABLE   :: A0
  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: snr       ! Amplitude array [ch,t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: phase     ! Phase array [ch,t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: impact    ! Impact array [ch, t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: bangle    ! Bangle array [ch, t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: ba_cov    ! Bangle covariance
  INTEGER                                 :: ocd       ! Occultation direction
  INTEGER                                 :: n         ! Number of data points
  CHARACTER(len = 256)                    :: routine
  
  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_bending_angle_wo')

!-------------------------------------------------------------------------------
! 2. Ensure monotonous impact parameter grid
!-------------------------------------------------------------------------------

  n = SIZE(impact_L1)
  ocd = NINT(SIGN(1.0_wp, impact_L1(n)-impact_L1(1)))
  CALL ropp_pp_monotonous(impact_L1, -1)
  CALL ropp_pp_monotonous(impact_L2, -1)

!-------------------------------------------------------------------------------
! 3. Create combined data arrays [channel, time]
!-------------------------------------------------------------------------------

  ALLOCATE(snr(2,n))
  ALLOCATE(phase(2,n))
  ALLOCATE(impact(2,n))
  ALLOCATE(bangle(2,n))
  ALLOCATE(ba_cov(2,n))
  ALLOCATE(AP(2,n))
  ALLOCATE(A0(2))

  snr(1,:) = snr_L1(:)
  snr(2,:) = snr_L2(:)

  phase(1,:) = phase_L1(:)
  phase(2,:) = phase_L2(:)

  impact(1,:) = impact_L1(:)
  impact(2,:) = impact_L2(:)
  
  bangle(1,:) = bangle_L1(:)
  bangle(2,:) = bangle_L2(:)

!-------------------------------------------------------------------------------
! 4. Scale lowest amplitude data point
!-------------------------------------------------------------------------------
 
  IF (ocd == 1) THEN
     snr(:,1) = 1e-5_wp*MAXVAL(snr(:,:))
  ELSE
     snr(:,n) = 1e-5_wp*MAXVAL(snr(:,:))
  END IF

!-------------------------------------------------------------------------------
! 5. Perform canonical transform (CT2) 
!-------------------------------------------------------------------------------
  
  IF (PRESENT(diag)) THEN
    CALL ropp_pp_DCT(time, snr, phase, r_leo, r_gns, r_coc, roc, w_ls,   &
                     w_smooth, w_low, hmax, filter, opt_DL2, cff, dsh,   &
                     impact, bangle, ba_cov, diag)
  ELSE
    CALL ropp_pp_DCT(time, snr, phase, r_leo, r_gns, r_coc, roc, w_ls,   &
                     w_smooth, w_low, hmax, filter, opt_DL2, cff, dsh,   &
                     impact, bangle, ba_cov) 
  ENDIF
  
!-------------------------------------------------------------------------------
! 6. Update output variables with computed bending angle and impact parameter
!-------------------------------------------------------------------------------

  impact_L1(:) = impact(1,:)
  impact_L2(:) = impact(2,:)
  bangle_L1(:) = bangle(1,:)
  bangle_L2(:) = bangle(2,:)
  ba_sigma_L1(:) = SQRT(ba_cov(1,:))
  ba_sigma_L2(:) = SQRT(ba_cov(2,:))

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(ba_cov)
  DEALLOCATE(bangle)
  DEALLOCATE(impact)
  DEALLOCATE(phase)
  DEALLOCATE(snr)
  DEALLOCATE(AP)
  DEALLOCATE(A0)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_pp_bending_angle_wo
