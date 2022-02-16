! $Id: ropp_pp_correct_L2.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_correct_L2(time, r_leo, r_gns, r_coc, roc, impact_LM,    &
                              phase_LM, phase_L1, phase_L2, snr_L1, snr_L2, &
                              lcf, Hmid, L2Q)

!****s* Preprocess/ropp_pp_correct_L2 *
!
! NAME
!    ropp_pp_correct_L2 - Correction of unusable L2 data.
!                           
! SYNOPSIS
!    call ropp_pp_correct_L2(time, r_leo, r_gns, r_coc, roc, impact_LM, 
!                             phase_LM, phase_L1, phase_L2, snr_L1, snr_L2,
!                             lcf, Hmid, L2Q)
! 
! DESCRIPTION
!    This routine corrects unusable L2 data.
!       1) Computation of penalty function for L2 channel from radio optical 
!          analysis.
!       2) Lower region correction:
!             A2 is multiplied with penalty function;
!             Corrected S2 is linear combination of S1 and S2.
!       3) Upper region correction:
!             Corrected S2 is linear combination of S2 and and smoothed S2.
!       4) Click removal in L2 phase.
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! Relative time of samples (s)
!    real(wp), dimension(:,:) :: r_leo     ! LEO coordinates (ECI or ECF) (m)
!    real(wp), dimension(:,:) :: r_gns     ! GPS coordinates (ECI or ECF) (m)
!    real(wp), dimension(:)   :: r_coc     ! centre curvature coordinates (m)
!    real(wp)                 :: roc       ! Radius curvature (m)
!    real(wp), dimension(:)   :: impact_LM ! Model impact parameters (m)
!    real(wp), dimension(:)   :: phase_LM  ! Model excess phase (m)
!    real(wp), dimension(:)   :: phase_L1  ! L1 excess phase (m)
!    real(wp), dimension(:)   :: phase_L2  ! L2 excess phase (m)
!    real(wp), dimension(:)   :: snr_L1    ! L1 amplitude
!    real(wp), dimension(:)   :: snr_L2    ! L2 amplitude
!    integer,  dimension(:)   :: lcf       ! Lost carrier flag  
!    real(wp)                 :: hmid      ! Border between upper/lower region 
!
! OUTPUT
!    real(wp), dimension(:)   :: phase_L1  ! L1 excess phase (m)
!    real(wp), dimension(:)   :: phase_L2  ! L2 excess phase (m)
!    real(wp), dimension(:)   :: snr_L1    ! L1 amplitude
!    real(wp), dimension(:)   :: snr_L2    ! L2 amplitude
!    real(wp)                 :: L2Q       ! L2 badness score
! 
! NOTES
!
! REFERENCES
!   Gorbunov M.E., Lauritsen K.B., Rhodin A., Tomassini M. and Kornblueh L. 
!   2006
!   Radio holographic filtering, error estimation, and quality control of 
!   radio occultation data
!   Journal of Geophysical Research (111) D10105
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
  USE messages
! USE ropp_pp, not_this => ropp_pp_correct_L2
  USE ropp_pp

  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! Cartesian LEO coordinates
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! Cartesian GPS coordinates
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! Centre curvature (m)
  REAL(wp),                 INTENT(in)  :: roc       ! Radius curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: impact_LM ! Model impact parameter
  REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM  ! Model excess phase (m)
  REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L1  ! L1 excess phase (m)
  REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L2  ! L2 excess phase (m)
  REAL(wp), DIMENSION(:), INTENT(inout) :: snr_L1    ! L1 amplitude
  REAL(wp), DIMENSION(:), INTENT(inout) :: snr_L2    ! L2 amplitude
  INTEGER,  DIMENSION(:), INTENT(inout) :: lcf       ! Lost carrier flag
  REAL(wp),                 INTENT(in)  :: hmid   ! Border upper/lower region
  REAL(wp),                 INTENT(out) :: L2Q       ! L2 badness score

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo   ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns   ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo   ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns   ! GPS velocity by regression
  INTEGER, PARAMETER                    :: np = 3 ! Polynomial degree for filter
  INTEGER                               :: npoints ! Number of data points
  INTEGER                               :: i, ic  ! Indices
! INTEGER                               :: ocd    ! Occultation direction ! Commented at 20 July, 2016

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: A   ! L1 and L2 amplitude
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: S   ! L1 and L2 excess phase
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: PA  ! Average impact parameter
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: PD  ! RMS deviation of impact
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DA  ! Smoothed doppler
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: EA  ! Smoothed refraction angles
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: SA  ! Smoothed phase excess

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Q   ! Badness [channel,height]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: WL  ! Lower L2 weighting function
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: WU  ! Upper L2 weighting function
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: D1S ! Finite differences of phase
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: DSI ! S2 - S1 from regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DS  ! S - SM
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: D2S ! 2nd derivative of phase
  LOGICAL,  DIMENSION(:),   ALLOCATABLE :: M   ! Mask for array operations

  INTEGER  :: WS        ! Window width for filtering S
  INTEGER  :: NR        ! Number of points of reduced grid
  INTEGER  :: mc        ! Mask counter
  REAL(wp) :: RPA       ! Badness term for PA1-PA2
  REAL(wp) :: RPD       ! Badness term for PD2
  REAL(wp) :: DEI       ! Estimate of ionospheric E2-E1
  REAL(wp) :: DPI       ! Estimate of ionospheric P2-P1
  REAL(wp) :: PC        ! Impact parameter for A2 cut-off
  
  REAL(wp), PARAMETER :: CPA = 200.0_wp   ! Weighting parameter for PA1-PA2 (m)
  REAL(wp), PARAMETER :: CPD = 150.0_wp   ! Weighting parameter for PD2 (m)
  REAL(wp), PARAMETER :: QC = 15.0_wp     ! Limiting badness for L2 cut-off
  REAL(wp), PARAMETER :: WSS = 2000.0_wp  ! Smoothing step width (m)

  REAL(wp), PARAMETER :: HQmax = 50000.0_wp ! Upper height for L2Q estimate
  REAL(wp), PARAMETER :: HQmin = 15000.0_wp ! Lower height for L2Q estimate

  INTEGER :: Icl     ! Click index
  INTEGER :: nl2q    ! Number of data for L2Q estimate
  CHARACTER(len = 247) :: outstr
  CHARACTER(len = 256)                :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_correct_L2')


!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

  npoints = SIZE(time)

! ocd = NINT(SIGN(1.0_wp, impact_LM(npoints)-impact_LM(1))) ! Commented at 20 July, 2016

  ALLOCATE(A(2,npoints))
  ALLOCATE(S(2,npoints))
  ALLOCATE(PA(2,npoints))
  ALLOCATE(PD(2,npoints))
  ALLOCATE(SA(2,npoints))
  
  A(1,:) = snr_L1(:)
  A(2,:) = snr_L2(:)
  S(1,:) = phase_L1(:)
  S(2,:) = phase_L2(:)
    
!-------------------------------------------------------------------------------
! 3. Analyse L2
!-------------------------------------------------------------------------------

  ! 3.1 Re-accumulation of L2

  ALLOCATE(DS(2,npoints))
  
  DS(2,:) = (2.0_wp * Pi * f_L2 / C_Light) * (S(2,:) - phase_LM(:))
  
  CALL Accumulate_Phase(DS(2,:))
  
  S(2,:) = phase_LM(:) + DS(2,:)/(2.0_wp * Pi * f_L2 / C_Light)

  DEALLOCATE(DS)
  

  ! 3.2 Radioholographic filter

  CALL ropp_pp_radioholographic_filter(time, r_leo, r_gns, r_coc, roc,        &
                                       phase_LM, S(2:2,:), A(2:2,:))

  ! 3.3 Radio holographic analysis

  CALL ropp_pp_radiooptic_analysis(time, r_leo, r_gns, r_coc, roc, phase_LM,  &
                                   S, A, PA, PD)

  snr_L1(:)   = A(1,:)
  snr_L2(:)   = A(2,:)
  phase_L1(:) = S(1,:)
  phase_L2(:) = S(2,:)

  ! 3.4 Excess phase smoothing
  
  WS = CEILING(250.0_wp*(npoints-1)/ABS(impact_LM(npoints)-impact_LM(1)))

  ALLOCATE(DS(2,npoints))

  DS(1,:) = phase_L1(:) - phase_LM(:)
  DS(2,:) = phase_L2(:) - phase_LM(:)

  nr = CEILING(REAL(npoints)/8000.0_wp)

  CALL ropp_pp_sliding_polynomial(time(1::nr),DS(:,1::nr),WS,np,SA(:,1::nr))

  IF (nr > 1) THEN
    DO ic=1,2 
      CALL ropp_pp_interpol(time(1::nr),time(:),SA(ic,1::nr),SA(ic,:))
    ENDDO
  ENDIF
  
  DO IC=1,2
     SA(IC,:) = SA(IC,:) + phase_LM(:)
  ENDDO

  phase_L2(:) = SA(2,:)
  
  DEALLOCATE(DS)

!-------------------------------------------------------------------------------
! 4. Weighting functions and badness score evaluation
!-------------------------------------------------------------------------------

  ! 4.1 Calculate weighting functions

  ALLOCATE(Q(2,npoints))
  ALLOCATE(WL(npoints))
  ALLOCATE(WU(npoints))

  DO i=1,npoints
     RPA    = ABS(PA(1,i)-PA(2,i))/CPA
     RPD    = PD(2,i)/CPD
     Q(2,i) = (RPA + RPD)**2
     Q(1,i) = (PD(1,i)/CPD)**2
     WL(i)=(1.0_wp-Smooth_Step((impact_LM(i)-roc-Hmid)/WSS))/2.0_wp
     WU(i)=1.0_wp - WL(i)
  ENDDO

  PC = MAXVAL(impact_LM(:), Mask = (Q(2,:)*WL(:) > QC .OR. BTEST(LCF(:),0)))
  PC = MIN(roc + Hmid - WSS, MAX(MINVAL(impact_LM(:)) - 2.0_wp*WSS, PC)) + WSS

  DO i=1,npoints
    WL(i) = (1.0_wp - Smooth_Step((impact_LM(i) - PC)/WSS))/2.0_wp
    WU(i) = 1.0_wp - WL(i)
  ENDDO

  ! 4.2 Estimation of badness

  ALLOCATE(M(npoints))
  
  M(:) = (impact_LM(:)-roc > HQmin) .AND. (impact_LM(:)-roc < HQmax)
  nl2q = COUNT(Mask = M(:))

  IF (nl2q > 0) THEN
    l2q = MAXVAL(Q(2,:), Mask = M(:))
  ELSE
    l2q = HUGE(l2q)
  ENDIF
  
  WRITE(outstr, '(2X,A,F5.1,A,F5.1,A,F9.3)') &
        'L2 badness between ', HQmin/1000.0,' and ', HQmax/1000.0,' km: ', L2Q
  CALL message(msg_diag, outstr)

  DEALLOCATE(PD)
  DEALLOCATE(Q)

!-------------------------------------------------------------------------------
! 5. Compute coordinates and velocities
!-------------------------------------------------------------------------------

  ALLOCATE(DA(2,npoints))
  ALLOCATE(EA(2,npoints))
  ALLOCATE(xleo(npoints,3))
  ALLOCATE(xgns(npoints,3))
  ALLOCATE(vleo(npoints,3))
  ALLOCATE(vgns(npoints,3))

  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo, xgns, vgns)

!-------------------------------------------------------------------------------
! 6. Estimate L2 - L1 ionospheric difference
!-------------------------------------------------------------------------------

 ! 6.1 Calculate smoothed doppler and bending angles

  DO ic=1,2
     DO i=1,npoints
        CALL ropp_pp_impact2doppler(xleo(i,:)-r_coc(:), vleo(i,:),    &
                                    xgns(i,:)-r_coc(:), vgns(i,:),    &
                                    PA(ic,i), DA(ic,i), EA(ic,i))
     ENDDO
  ENDDO

  DEALLOCATE(DA)
  DEALLOCATE(xleo)
  DEALLOCATE(xgns)
  DEALLOCATE(vleo)
  DEALLOCATE(vgns)
  
  ! 6.2 Estimate ionospheric difference of bending angles

  M(:) = (PC <= impact_LM(:)) .AND. (impact_LM(:) <= PC + 2000.0_wp)
  mc = COUNT(Mask = M(:))

  IF (MC > 0) THEN
    DEI = SUM(EA(2,:) - EA(1,:), Mask = M(:)) / MC
    DPI = SUM(PA(2,:) - PA(1,:), Mask = M(:)) / MC
  ELSE
    DEI = 0.0_wp
    DPI = 0.0_wp
  ENDIF
  
  ! 6.3 Calculate ionospheric excess phase difference

  EA(2,:) = EA(1,:) + DEI
  PA(2,:) = PA(1,:) + DPI

  DO ic=1,2
     CALL ropp_pp_bangle2phase(time, r_leo, r_gns, r_coc,   &
                               PA(ic,:), EA(ic,:), SA(ic,:))
  ENDDO

  DEALLOCATE(M)
  DEALLOCATE(PA)
  DEALLOCATE(EA)
  ALLOCATE(DSI(npoints))

  DSI(:) = SA(2,:) - SA(1,:)
  
!-------------------------------------------------------------------------------
! 7. Correction of L2 phase and amplitude data
!-------------------------------------------------------------------------------

  ! 7.1 Amplitude correction

  DO i=1,npoints
    snr_L2(i)   = snr_L2(i)*(0.1_wp + 0.9_wp*WU(i))
  ENDDO
  
  ! 7.2 Excess phase correction

  ALLOCATE(D1S(2,npoints))

  D1S(1,1) = 0.0_wp
  D1S(2,1) = 0.0_wp
  DO i=2,npoints
     D1S(1,i) = phase_L1(i) - phase_L1(i-1) + DSI(i) - DSI(i-1)
     D1S(2,i) = phase_L2(i) - phase_L2(i-1)
  ENDDO
  
  DO i=1,npoints
     D1S(2,i) = D1S(2,i)*WU(i) + D1S(1,i)*WL(i)
  ENDDO

  DO i=2,npoints
     phase_L2(i) = phase_L2(i-1) + D1S(2,i)
  ENDDO
  
  ! 7.3 Click removal
  
  ALLOCATE(D2S(2:npoints-1))

  DO i=1,3
     D2S(:) = -(phase_L2(2:npoints-1) - phase_L2(1:npoints-2)) *     &
                   (phase_L2(3:npoints) - phase_L2(2:npoints-1))
     
     Icl = INT(SUM(MAXLOC(D2S(:)) + LBOUND(D2S)-1.0_wp))

     phase_L2(Icl) = (phase_L2(Icl+1) + phase_L2(Icl-1))/2.0_wp
     
     WRITE(outstr,'(2X,A,I5,A,F8.3)') &
       'Detected click in L2 phase at index ', Icl, ', time = ', time(Icl)
     CALL message(msg_diag, outstr)

  ENDDO

  DEALLOCATE(D2S)

!-------------------------------------------------------------------------------
! 8. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(S)
  DEALLOCATE(A)
  DEALLOCATE(SA)

  DEALLOCATE(WL)
  DEALLOCATE(WU)
  DEALLOCATE(DSI)
  DEALLOCATE(D1S)

  CALL message_set_routine(routine)

CONTAINS

!-------------------------------------------------------------------------------
! 7. Smooth step from -1 to 1 for X changing from -1 to 1
!-------------------------------------------------------------------------------

  FUNCTION smooth_step(x) RESULT(f)

    USE typesizes, ONLY: wp => EightByteReal
    IMPLICIT NONE

    REAL(wp), INTENT(in) :: x ! function argument
    REAL(wp)             :: f ! smooth step from -1 to 1
    REAL(wp)             :: a

! 7.1 Result constant +1 for x above 1 
    IF ( x >= 0.9999999_wp) THEN 
       f = 1.0_wp
! 7.2 Result constant -1 for x below -1 
    ELSE  IF ( x <= -0.9999999_wp) THEN  
       f = -1.0_wp
! 7.3 Result smoothly varying from -1 to 1 as x changes from -1 to 1
    ELSE
       a = 1.0_wp - 1.0_wp/(1.0_wp-MIN(1.0_wp,ABS(x)))**2
       IF ( ABS(a) < 500.0_wp ) THEN 
          f = SIGN(1.0_wp,x)*(1.0_wp - EXP(a))
       ELSE
          f = SIGN(1.0_wp,x)
       ENDIF
    ENDIF

  END FUNCTION smooth_step

  SUBROUTINE Accumulate_Phase(Ph, Sign)   ! (Array of (accumulated) phase, dir)
  
! Method:
!   Sign = 0 or no Sign:
!      Adding +-2*Pi where phase jumps from
!      +-Pi to -+Pi,
!   Sign > 0:
!      Adding +2*Pi where phase jumps from
!      - to +
!   Sign < 0
!      Adding -2*Pi where phase jumps from
!      + to -
  
    ! 15.1 Declarations
    
    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_constants, ONLY: pi
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(inout) :: Ph   ! Phase --> accumulated phase
    INTEGER, OPTIONAL,      INTENT(in)    :: Sign ! Phase change sign

    INTEGER  :: i     ! Array index
    INTEGER  :: PSign ! Phase change sign

    ! 15.2 Determine phase change sign
    
    IF (.NOT. PRESENT(Sign)) THEN
      PSign = 0
    ELSE 
      PSign = Sign
    ENDIF
    
    ! 15.3 Accumulate phase

    IF (PSign == 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1)+pi, 2*pi) - pi
      ENDDO
    ELSEIF (PSign > 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1), 2*pi)
      ENDDO
    ELSEIF (PSign < 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1)+2*pi, 2*pi) - 2*pi
      ENDDO
    ENDIF
    
  END SUBROUTINE Accumulate_Phase


END SUBROUTINE ropp_pp_correct_L2


