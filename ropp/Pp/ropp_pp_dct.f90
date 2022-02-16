! $Id: ropp_pp_DCT.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_DCT(time, snr, phase, r_leo, r_gns, r_coc, roc, w_ls,    &
                       w_smooth, w_low, hmax, filter, opt_DL2, cff, dsh,    &
                       impact, bangle, ba_cov, diag)

!****s* WaveOptics/ropp_pp_DCT *
!
! NAME
!    ropp_pp_DCT - Calculate L1 and L2 bending angle profiles using
!                  Canonical Transform (CT2).
!
! SYNOPSIS
!    call ropp_pp_DCT(time, snr, phase, r_leo, r_gns, r_coc, roc, w_ls,
!                     w_smooth, hmax, filter, opt_DL2, cff, dsh,
!                     impact, bangle, ba_cov, diag)
!
! DESCRIPTION
!    This routine calculates L1 and L2 bending angles using a CT2 algorithm.
!
! INPUTS
!    real(wp), dimension(:)   :: time     ! Relative time of samples (s)
!    real(wp), dimension(:,:) :: snr      ! L1,L2 amplitudes [channel, time]
!    real(wp), dimension(:,:) :: phase    ! L1,L2 excess phase (m) [ch, time]
!    real(wp), dimension(:,:) :: r_leo    ! LEO coordinates (m) (ECI or ECF)
!    real(wp), dimension(:,:) :: r_gns    ! GPS coordinates (m) (ECI or ECF)
!    real(wp), dimension(:)   :: r_coc    ! Centre curvature coordinates (m)
!    real(wp)                 :: roc      ! Radius of curvature (m)
!    integer                  :: w_ls     ! Large-scale smoothing window
!    integer                  :: w_smooth ! Smoothing window above 7km (points)
!    integer                  :: w_low    ! Smoothing window below 7km (points)
!    real(wp)                 :: hmax     ! Maximum height for WO processing
!    character(len=*)         :: filter   ! Filter method ('optest','slpoly')
!    logical                  :: opt_DL2  ! Degraded L2 flag
!    integer                  :: cff      ! Complex filtering flag
!    real(wp)                 :: dsh      ! Shadow border width (m)
!    real(wp), dimension(:,:) :: impact   ! L1,L2 impact parameters (m) [ch,t]
!    real(wp), dimension(:,:) :: bangle   ! L1,L2 bending angles (rad) [ch,t]
!
! OUTPUT
!    real(wp), dimension(:,:) :: impact   ! L1,L2 impact parameters (m) [ch, t]
!    real(wp), dimension(:,:) :: bangle   ! L1,L2 bending angles (rad) [ch, ip]
!    real(wp), dimension(:,:) :: ba_cov   ! Estimate of bangle covariance
!    type(PPdiag),   optional :: diag     ! Additional output diagnostics structure
!
! NOTES
!   Impact parameters are calculated with respect to local centre of curvature
!   Bending angles are calculated from Doppler shift in the approximation of
!   local spherical symmetry.
!   Variable names follow Gorbunov and Lauritsen (2004)
!                   P - impact parameter
!                   E - bending angle (epsilon)
!                   Y - new coordinate
!
! REFERENCES
!   Gorbunov M.E. and Lauritsen K.B. 2004
!   Analysis of wave fields by Fourier integral operators and their application
!   for radio occultations
!   Radio Science (39) RS4010
!
!   Gorbunov M.E., Lauritsen K.B., Rhodin A., Tomassini M. and Kornblueh L.
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
!   M. Gorbunov, Russian Academy of Sciences, Russia.
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
  USE ropp_pp, ONLY: ropp_pp_geometric_optics,      &
                     ropp_pp_geometric_optics_adj,  &
                     ropp_pp_satellite_velocities,  &
                     ropp_pp_interpol,              &
                     ropp_pp_FFT,                   &
                     ropp_pp_fourier_filter,        &
                     ropp_pp_filter,                &
                     ropp_pp_sliding_polynomial
  USE ropp_pp_constants, ONLY: c_light, pi, f_L1, f_L2
  USE ropp_pp_spline
  USE ropp_pp_utils
  USE ropp_pp_types
  USE ropp_utils, ONLY: vector_angle
  USE messages

  IMPLICIT NONE

  ! 1.1 Subroutine arguments

  REAL(wp), DIMENSION(:),   INTENT(in)    :: time     ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: snr      ! Amplitudes [ch, time]
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: phase    ! Excess phase (m) [ch, t]
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_leo    ! LEO coordinates (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_gns    ! GPS coordinates (m)
  REAL(wp), DIMENSION(:),   INTENT(in)    :: r_coc    ! Centre of curvature (m)
  REAL(wp),                 INTENT(in)    :: roc      ! Radius of curvature (m)
  INTEGER,                  INTENT(in)    :: w_ls     ! Large-scale smoothing
  INTEGER,                  INTENT(in)    :: w_smooth ! Smoothing above 7km
  INTEGER,                  INTENT(in)    :: w_low    ! Smoothing below 7km
  REAL(wp),                 INTENT(in)    :: hmax     ! Maximum height for WO
  CHARACTER(len=*),         INTENT(in)    :: filter   ! Filter method
  LOGICAL,                  INTENT(in)    :: opt_DL2  ! Degraded L2
  INTEGER,                  INTENT(in)    :: cff      ! Complex filter flag
  REAL(wp),                 INTENT(in)    :: dsh      ! Shadow border width (m)
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: impact   ! Impact parameters (m)
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: bangle   ! Bending angles (rad)
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: ba_cov   ! Bending angle covariance
  TYPE(PPdiag), OPTIONAL,   INTENT(inout) :: diag     ! Additional diagnostics

  ! 1.2 Local parameters

  COMPLEX(wp), PARAMETER :: Ci = (0.0_wp, 1.0_wp)  ! I = Sqrt(-1)
  INTEGER,  PARAMETER :: nv  = 5           ! Polynomial degree for velocity
  INTEGER,  PARAMETER :: np  = 3           ! Polynomial degree for filtering
  INTEGER,  PARAMETER :: nd  = 11          ! No. of points for differentiation
  REAL(wp), PARAMETER :: PL  = 1700.0_wp   ! Lower estimate minimum ray height
  REAL(wp), PARAMETER :: PLT = 7000.0_wp   ! Lower troposphere region height
  REAL(wp), PARAMETER :: ARH = 500.0_wp    ! Aperture for r-holo analysis
  REAL(wp), PARAMETER :: PSS0 = 10000.0_wp ! Upper height for lower-trop errors
  LOGICAL,  PARAMETER :: opt_SH = .TRUE.   ! Automatic shadow determination
  LOGICAL,  PARAMETER :: opt_QC = .FALSE.  ! Automatic QC from CT amplitude
  REAL(wp), PARAMETER :: dem = 0.004_wp    ! Limit bending angle error (rad)
  INTEGER,  PARAMETER :: nur = 32768

  ! 1.3 Local scalars

  INTEGER  :: i,j      ! Data index
  INTEGER  :: ic       ! Channel index
  INTEGER  :: icw      ! Index of worst quality channel
  INTEGER  :: iss      ! Sliding window index
  INTEGER  :: n        ! Number of input data
  INTEGER  :: nc       ! Number of channels
  INTEGER  :: nu       ! High-resolution grid dimension
  INTEGER  :: nf       ! Lowered-resolution grid dimension
  INTEGER  :: nr       ! Resolution ratio
  INTEGER  :: nrh      ! Number points for radio-holographic analysis
  INTEGER  :: nss      ! Number of slididng-spectra
  REAL(wp) :: E_dd     ! d(E)/d(d)
  REAL(wp) :: Ymin     ! Lower limit of Y-grid
  REAL(wp) :: Ymax     ! Upper limit of Y-grid
  REAL(wp) :: YC       ! Integration constant for Y(T)
  REAL(wp) :: stepY    ! Step of Y-grid
  REAL(wp) :: stepP    ! Step of P-grid
  REAL(wp) :: dY       ! Additional Y-interval
  INTEGER  :: wf       ! Smoothing window on reduced grid
  REAL(wp) :: Pmin     ! Minimum ray height estimate
  REAL(wp) :: Pmax     ! Maximum ray height estimate
  REAL(wp) :: Algt     ! Light-zone amplitude
  REAL(wp) :: Ashd     ! Shadow-zone amplitude
  REAL(wp) :: Athr     ! Amplitude threshold
  REAL(wp) :: Ascl     ! Scaled amplitude
  REAL(wp) :: dfI      ! Interpolated model Doppler
  REAL(wp) :: P0I      ! Interpolated model impact parameter
  REAL(wp) :: dEC      ! Correction of covariance
  INTEGER  :: IRHmin   ! Lower limit of sliding window
  INTEGER  :: IRHmid   ! Middle point of sliding window
  INTEGER  :: IRHmax   ! Upper bound of sliding window
  REAL(wp) :: dpw      ! Filter width of impact parameter [m]
  REAL(wp) :: cfw      ! Filter width for complex field [points]
  INTEGER  :: m        ! Spatial dimension index
  REAL(wp) :: t_norm   ! Normalized time

  REAL(wp) :: rgnsI    ! Interpolated GPS radius from r_coc
  REAL(wp) :: rleoI    ! Interpolated LEO radius from r_coc
  REAL(wp) :: thetaI   ! Interpolated satellite-to-satellite angle
  INTEGER  :: imax     ! Upper index of redefined grid
  INTEGER  :: imin     ! Lower index of redefined grid
  REAL(wp) :: dp       ! Step of low-res impact parameter grid
  REAL(wp) :: dpH      ! Step of hi-res impact parameter grid
  INTEGER  :: wh       ! Hi-res filter width
  INTEGER  :: whl      ! Hi-res filter width for lower troposphere
  INTEGER  :: ifb      ! Shadow border index
  REAL(wp) :: pfb      ! Filter border width [m]
  INTEGER  :: wfb      ! Window width near filter border [points]
  REAL(wp) :: dE       ! Mean channel-to-channel difference of bangle(p)
  INTEGER  :: ndE      ! Number of samples for estimation of dE
  CHARACTER(len = 1)   :: istr
  CHARACTER(len = 7)   :: nstr
  CHARACTER(len = 256) :: routine
  CHARACTER(len = 256) :: outstr

  ! 1.4 Local arrays

  REAL(wp), DIMENSION(3)      :: U0              ! GPS-LEO straight line dir
  REAL(wp), DIMENSION(3)      :: xgnsI           ! Interpolated GPS positions
  REAL(wp), DIMENSION(3)      :: xleoI           ! Interpolated LEO positions
  REAL(wp), DIMENSION(3)      :: vgnsI           ! Interpolated GPS velocities
  REAL(wp), DIMENSION(3)      :: vleoI           ! Interpolated LEO velocities
  REAL(wp), DIMENSION(6)      :: E_dr            ! d(E)/d(r_leo,r_gns)
  REAL(wp), DIMENSION(6)      :: P_dr            ! d(P)/d(r_leo,r_gns)
  REAL(wp), DIMENSION(0:nv,3) :: coeff_vleo      ! Regression coeffs for vleo
  REAL(wp), DIMENSION(0:nv,3) :: coeff_vgns      ! Regression coeffs for vgns

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo  ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo  ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns  ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns  ! GPS velocity by regression

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: k     ! Wave number [channel]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: d0    ! Vacuum relative Doppler shift
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: df    ! Filtered Doppler shift [ch,t]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: s0    ! Vacuum excess phase(m)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sf    ! Filtered excess phase [ch,t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dsf   ! Filtered derivative ex phase
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sm0   ! Model vacuum phase path (m)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sm    ! Model phase path (m)

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: P0    ! Impact parameter model
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: E0    ! Bending angle model (rad)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: P_dd  ! d(P)/d(d)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ft    ! P0(t)-d0(t)*P_dd

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Y     ! New coordinate
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fti   ! integral of ft over Y

  INTEGER,  DIMENSION(:),   ALLOCATABLE :: nh    ! Hi-res grid dimensions [ch]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: YH    ! Hi-res grids [ch, point]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: YHL   ! Hi-res grids lower trop
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: AH    ! Interpolated SNR [ch, point]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: Apf   ! Prefiltered amplitude
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: Af    ! Filtered amplitude
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: SMH   ! Interpolated eikonal [ch, pt]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: SMF   ! Filtered eikonal [ch, point]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PH    ! Impact parameters [ch, point]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: EH    ! Bending angles [ch, point]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: APY   ! Amplitude function of FIO
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: whv   ! Vector filter width

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: P_dh  ! Interpolated P_dd
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: th    ! Time as function of y [ch,pt]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: dh    ! Doppler function of p [ch,pt]
  COMPLEX(wp), DIMENSION(:),ALLOCATABLE :: UH    ! Interpolated complex field
  COMPLEX(wp), DIMENSION(:),ALLOCATABLE :: UR    ! Complex R.Hol with ref signal

  COMPLEX(wp), DIMENSION(:),ALLOCATABLE :: URS   ! Radio-hologram spectrum
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: YRH   ! Y-grid for radio-hologram
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PSS   ! P-grid for sliding windows
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ECSS  ! Estimate bangle covariance

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: AP    ! SNR transform field
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: A0    ! Vacuum SNR [ch]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PminC ! Estimate of shadow border
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: EW    ! Bangle for re-interpolation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ECW   ! Ba cov for re-interpolation
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: AW    ! Amplitude for re-interpol
  LOGICAL,  DIMENSION(:),   ALLOCATABLE :: MW    ! Mask for computation of dE
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: i0    ! Lower index of work area
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: i1    ! Upper index of work area
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: ib    ! Index of border point

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_DCT')

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

  ! 2.1 Determination of data sizes

  n  = SIZE(time)
  nc = SIZE(snr,1)

  ! 2.2 Array allocation

  ALLOCATE(df(nc,n))
  ALLOCATE(sm(nc,n))
  ALLOCATE(P0(nc,n))
  ALLOCATE(P_dd(nc,n))
  ALLOCATE(Y(nc,n))
  ALLOCATE(AP(nc,n))
  ALLOCATE(A0(nc))

!===============================================================================
! I. GEOMETRIC OPTICS PROCESSING [Section 3.4 Gorbunov2004]
!===============================================================================

!-------------------------------------------------------------------------------
! 3. Determination of vacuum model
!-------------------------------------------------------------------------------

  ! 3.1 Array allocation

  ALLOCATE(xleo(n,3))
  ALLOCATE(xgns(n,3))
  ALLOCATE(vleo(n,3))
  ALLOCATE(vgns(n,3))

  ALLOCATE(d0(n))
  ALLOCATE(s0(n))
  ALLOCATE(E0(nc,n))
  ALLOCATE(sf(nc,n))
  ALLOCATE(dsf(nc,n))
  ALLOCATE(sm0(nC,n))
  ALLOCATE(ft(nc,n))
  ALLOCATE(fti(nc,n))

  ! 3.2 Determination of satellite coordinates and velocities

  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,     &
                                    xgns, vgns, coeff_vleo, coeff_vgns)

  ! 3.3 Vacuum Doppler shift

  DO i=1,n
     U0    = (xleo(i,:) - xgns(i,:))/SQRT(SUM((xleo(i,:) - xgns(i,:))**2))
     d0(i) = (C_Light - DOT_PRODUCT(vleo(i,:),U0)) /     &
                 (C_Light - DOT_PRODUCT(vgns(i,:),U0)) - 1.0_wp
  ENDDO

  ! 3.4 Vacuum phase path

  s0(1) = 0.0_wp
  DO i=2,n
     s0(i) = s0(i-1) - C_Light*(d0(i-1) + d0(i))*(time(i) - time(i-1))/2.0_wp
  ENDDO

!-------------------------------------------------------------------------------
! 4. Determination of smooth model
!-------------------------------------------------------------------------------

  ! 4.1 Large-scale filtering of excess phase

  nr = CEILING(REAL(n,wp)/8000.0_wp)
  wf = CEILING(2.0_wp*w_ls/REAL(nr,wp))

  CALL ropp_pp_sliding_polynomial(time(1::nr), phase(:,1::nr), wf, np,    &
                                  sf(:,1::nr), dsf(:,1::nr))

  IF (nr > 1) THEN
     DO ic=1,nc
        CALL ropp_pp_interpol(time(1::nr), time(:), dsf(ic,1::nr), dsf(ic,:))
     ENDDO
  ENDIF

  ! 4.2 Filtered doppler shift [sigma0(t)]

  DO ic=1,nc
     df(ic,:) = d0(:) - dsf(ic,:)/C_Light
  ENDDO

  ! 4.3 Determination of impact parameter model [p0(t), f(t)]

  DO ic=1,nc
     DO i=1,n
        CALL ropp_pp_geometric_optics_adj(xleo(i,:)-r_coc(:), vleo(i,:),      &
                                          xgns(i,:)-r_coc(:), vgns(i,:),      &
                                          df(ic,i), P0(ic,i), E0(ic,i),       &
                                          P_dd(ic,i), P_dr(:), E_dd, E_dr(:))
        ft(ic,i) = P0(ic,i) - df(ic,i)*P_dd(ic,i)
      ENDDO
  ENDDO

  Pmax = MIN(hmax, MAXVAL(impact(1,:)) - roc)

  ! 4.4 Determination of new coordinate grid Y  [dY = d(sigma)/d(p0) dt]

  Y(:,1) = 0.0_wp

  DO ic=1,nc
     DO i=2,n
        Y(ic,i) = Y(ic,i-1) - C_Light*(1.0_wp/P_dd(ic,i-1) +      &
                       1.0_wp/P_dd(ic,i))*(time(i) - time(i-1))/2.0_wp
     ENDDO
     YC      = MINVAL(Y(ic,:))
     Y(ic,:) = Y(ic,:) - YC
  ENDDO

  ! 4.5 Determine accumulated phase path

  ! 4.5.1 Integral of f(t) over Y

  fti(:,1) = 0.0_wp
  DO i=2,n
     fti(:,i) = fti(:,i-1) + (ft(:,i-1) + ft(:,i))*(Y(:,i) - Y(:,i-1))/2.0_wp
  ENDDO

  ! 4.5.2 Model vacuum phase path
  DO ic=1,nc
     DO i=1,n
        sm0(ic,i) = s0(i) - roc*Y(ic,i) + fti(ic,i)
     ENDDO
  ENDDO

  ! 4.5.3 Model accumulated phase path [psi(t) = excess phase + psi0(t)]

  DO ic=1,nc
     DO i=1,n
        sm(ic,i) = phase(ic,i) + sm0(ic,i)
     ENDDO
  ENDDO

  ! 4.5 Deallocate finished arrays

  DEALLOCATE(xleo)
  DEALLOCATE(xgns)
  DEALLOCATE(vleo)
  DEALLOCATE(vgns)

  DEALLOCATE(d0)
  DEALLOCATE(s0)
  DEALLOCATE(E0)
  DEALLOCATE(sf)
  DEALLOCATE(dsf)
  DEALLOCATE(sm0)
  DEALLOCATE(ft)
  DEALLOCATE(fti)


!===============================================================================
! II. WAVE OPTICS PROCESSING
!===============================================================================


!-------------------------------------------------------------------------------
! 5. Determination of high-resolution grid
!-------------------------------------------------------------------------------

  ALLOCATE(k(nc))
  ALLOCATE(nh(nc))
  ALLOCATE(PminC(nc))
  ALLOCATE(i0(nc))
  ALLOCATE(i1(nc))
  ALLOCATE(ib(nc))

  ! 5.1 Computation of wave vectors

  k(1)   = 2*pi*f_L1/C_Light
  k(2)   = 2*pi*f_L2/C_Light

  ! 5.2 Definition of high resolution grid

  DO ic=1,nc
     dY     = ABS(Y(ic,n) - Y(ic,1))/4.0_wp
     Ymin   = MIN(Y(ic,1), Y(ic,n))
     Ymax   = MAX(Y(ic,1), Y(ic,n))
     nh(ic) = CEILING(k(ic)*(MAXVAL(impact(1,:)) - roc)*(Ymax - Ymin + 2*dY)/pi)
  ENDDO

!-------------------------------------------------------------------------------
! 6. Fourier integral operator
!-------------------------------------------------------------------------------

  Channels: DO ic=1,nc

    ! 6.1 Interpolation of complex field to high resolution Y-grid

    dY   = ABS(Y(ic,n) - Y(ic,1))/4.0_wp
    Ymin = MIN(Y(ic,1), Y(ic,N))
    Ymax = MAX(Y(ic,1), Y(ic,n))

    nu = ropp_pp_nearest_power2(3*nh(ic)/2)

    nf = MIN(NUR, nu)
    nr = nu/nf

    WRITE(istr, '(i1)') ic
    WRITE(nstr, '(i7)') nu
    CALL message(msg_diag, 'Channel ' //istr// '. Hi-res grid size = ' //nstr)

    ALLOCATE(YH(nu))
    ALLOCATE(smh(nu))
    ALLOCATE(P_dh(nu))
    ALLOCATE(UH(nu))
    ALLOCATE(AH(nu))

    DO i=1,nu

      YH(i) = ((Ymin-dY)*(nu-i) + (Ymax+dY)*(i-1))/REAL(nu-1,wp)

    ENDDO

    CALL ropp_pp_interpol(Y(ic,:), YH, sm(ic,:),   smh)
    CALL ropp_pp_interpol(Y(ic,:), YH, snr(ic,:),  AH)
    CALL ropp_pp_interpol(Y(ic,:), YH, P_dd(ic,:), P_dh)

    DO i=1,nu

      IF (YH(i) > Ymin .AND. YH(i) < Ymax) THEN

        UH(i) = ((Ymax - Ymin + 2*dY)/(nu-1.0_wp))*                      &
                   SQRT(2000.0_wp*ABS(-P_dh(i)/C_Light)*k(ic)/(2*pi))*   &
                   AH(i) * EXP(Ci*(MODULO(k(ic)*smh(i), 2*pi)))

      ELSE

        UH(i) = 0.0

      ENDIF
    ENDDO

    ! 6.2 Definition of impact parameter grid

    ALLOCATE(PH(nu))
    ALLOCATE(EH(nu))

    stepY = ABS(YH(nu) - YH(1))/(nu-1.0_wp)

    DO i=1,nu
      PH(i) = (i-1.0_wp)*2.0_wp*pi/(k(ic)*nu*stepY)
    ENDDO

    i0(ic) = MIN(SUM(MINLOC(PH(:), hmax >  PH(:) .AND. PH(:) > PL)),  &
                 SUM(MAXLOC(PH(:), hmax >  PH(:) .AND. PH(:) > PL)))
    i1(ic) = MAX(SUM(MINLOC(PH(:), hmax >  PH(:) .AND. PH(:) > PL)),  &
                 SUM(MAXLOC(PH(:), hmax >  PH(:) .AND. PH(:) > PL)))

    dpH = ABS(PH(i1(ic)) - PH(i0(ic)))/(i1(ic) - i0(ic))

    ! 6.3 Fourier transform

    call ropp_pp_FFT(UH, -1)

    ! 6.4 Calculation of amplitude and accumulated phase

    AH(:) = ABS(UH(:))

    WHERE (UH(:) /= 0.0_wp)
      smh(:) = ATAN2(AIMAG(UH(:)),REAL(UH(:)))
    ELSEWHERE
      smh(:) = 0.0
    END WHERE

    CALL Accumulate_Phase(smh(:),-1)

    DEALLOCATE(UH)

!-------------------------------------------------------------------------------
! 7. Complex field filtering
!-------------------------------------------------------------------------------

    IF (CFF >= 2) THEN

      CALL message(msg_info,"Complex field filtering \n")

      ! 7.1 Memory allocation

      ALLOCATE(smf(nu))
      ALLOCATE(UR(nu))

      ! 7.2 Filter width determination

!      dpw = REAL(w_smooth, wp)*(MAXVAL(impact(1,:)) - roc - Pmin)/REAL(n-1, wp)
      DPW = 250.0_wp
      cfw = dpw/(nr*dpH)

      WRITE(outstr, '(A,F10.3,2X,A,F7.1)') 'DPW = ', DPW, 'CFW = ', CFW
      CALL message(msg_diag, outstr)

      ! 7.3 Computation of phase model

      SELECT CASE(filter)
      CASE('optest')
        CALL ropp_pp_filter(nr*dpH, smh(1::nr), CEILING(cfw), nd, smf(1::nr))
      CASE('slpoly')
        CALL ropp_pp_sliding_polynomial(PH(1::nr), smh(1::nr), CEILING(cfw),  &
                                        np, smf(1::nr))
      END SELECT

      IF (nr > 1) THEN
        CALL ropp_pp_interpol(PH(1::nr), PH, smf(1::nr), smf)
      ENDIF

      ! 7.4 Multiplication with reference signal

      UR(:) = AH(:)*EXP(Ci*(MODULO(smh(:)-smf(:), 2*pi)))

      ! 7.5 Complex field filtering

      CALL ropp_pp_fourier_filter(UR, cfw)

      AH(:) = ABS(UR(:))

      IF (CFF == 2) THEN    ! Compute filtered CT phase

        WHERE (UR(:) /= 0.0_wp)
          smh(:) = ATAN2(AIMAG(UR(:)),REAL(UR(:)))
        ELSEWHERE
          smh(:) = 0.0
        END WHERE

        CALL Accumulate_Phase(smh(:))    ! Array of (accumulated) phase

        ! 7.6 Phase restoration

        smh(:) = smh(:) + smf(:)

      ENDIF

      ! 7.7 Memory deallocation

      DEALLOCATE(smf)
      DEALLOCATE(UR)

    ENDIF

!-------------------------------------------------------------------------------
! 8. Determination of shadow border
!-------------------------------------------------------------------------------

    ALLOCATE(Af(nu))

    IF (opt_SH) THEN

      ! 8.1 Determination of light and shadow amplitudes

      Algt = SQRT(SUM(AH(:)**2, &
                  Mask = (PH(:) > Pmax-5000.0) .AND. (PH(:) < Pmax)) /    &
                  COUNT(Mask = (PH(:) > Pmax-5000.0) .AND. (PH(:) < Pmax)))

      Ashd = SQRT(SUM(AH(:)**2, &
                  Mask = (PL-1000.0 < PH(:) .AND. PH(:) < PL)) /          &
                  COUNT(Mask = (PL-1000.0 < PH(:) .AND. PH(:) < PL)))

      ! 8.2 Determination of threshold and scaling amplitude

      Athr = 0.5_wp*(Algt + Ashd)

      Ascl = MIN(Athr, AH(I1(ic)) - Ashd)

      ! 8.3 Computation of correlation with step function

      Af(i1(ic)) = Ascl

      DO i=i1(ic)-1,i0(ic),-1
        Ascl  = MIN(Athr, AH(i)) - Ashd
        Af(i) = Af(i+1) + Ascl
      ENDDO

      DO i=i1(ic),i0(ic),-1
        Af(i) = Af(i)/SQRT(REAL(i1(ic)+1-i))
      ENDDO

      ! 8.4 Determination of shadow zone border and
      !     shifting it to nearest point of reduced grid

      ib(ic)    = i0(ic) + SUM(MAXLOC(Af(i0(ic):i1(ic)))) - 1
      ib(ic)    = ib(ic) + NINT(dSh/dpH)
      ib(ic)    = MIN(i1(ic)-nr, MAX(i0(ic),ib(ic)))
      ib(ic)    = nr*CEILING(REAL(ib(ic)-1)/REAL(nr)) + 1
      PminC(ic) = PH(ib(ic))

      WRITE(outstr,'(F10.3)') PminC(ic)
      CALL message(msg_diag, 'DP0 = '// outstr)

      ! 8.5 Determination of CT amplitude scintillations

      IF (opt_QC) THEN
        ib(ic)    = SUM(MAXLOC(PH(:), &
                        Mask = (AH(:) < 0.3*Algt) .AND. (PH(:) < Pmax)))
        IF ((ib(ic) < 1) .OR. (ib(ic) > size(PH))) THEN
          ib(ic) = 1
        ENDIF

        ib(ic) = nr*CEILING(REAL(ib(ic)-1)/REAL(nr)) + 1
        PminC(ic) = PH(ib(ic))
      ENDIF

    ELSE

      ! 8.6 Setting min impact height if no automatic shadow zone determination

      PminC(ic) = PL

    ENDIF

    ! 8.7 Setting minimum impact height from current or reference channel

    IF (opt_DL2) THEN
      Pmin = PminC(ic)
    ELSE
      Pmin = PminC(1)
    ENDIF

    WRITE(outstr, '(2(A,F10.3))') 'Pmin = ', Pmin, ' Pmax= ', Pmax
    CALL message(msg_diag, outstr)

!-------------------------------------------------------------------------------
! 9. Phase filtering and differentiation
!-------------------------------------------------------------------------------

    ALLOCATE(smf(nu))

    ifb = SUM(MinLoc(PH(:), PH(:) >= PminC(1)))
    ifb = MIN(i1(ic)-nr, MAX(i0(ic),ifb))
    ifb = nr*CEILING(REAL(ifb-1)/REAL(nr)) + 1

    dp  = (MAXVAL(impact(ic,:)) - roc - Pmin)/(N - 1)
    wh  = MAX(3, NINT(REAL(w_smooth)*dp/(nr*dpH)))

    SELECT CASE(filter)
    CASE('optest')
      CALL ropp_pp_filter(nr*dpH, smh(1::nr), wh, nd, smf(1::nr), YH(1::nr))

      IF ( w_low > 0 .AND. w_low /= w_smooth) THEN

        whl = MAX(3, NINT(REAL(w_low)*dp/(nr*dpH)))
        ALLOCATE(YHL(nu))
        CALL ropp_pp_filter(nr*dpH, smh(1::nr),whl, nd, smf(1::nr), YHL(1::nr))
        WHERE (PH(:) < PLT)   ! Lower troposphere region height
          YH(:) = YHL(:)
        ENDWHERE
        DEALLOCATE(YHL)

      ENDIF

    CASE('slpoly')

      ALLOCATE(whv(nu))
      whv(:) = wh

      IF ( w_low > 0 ) THEN
        whl = MAX(np, NINT(REAL(w_low)*dp/(nr*dpH)))
        WHERE (PH(:) < PLT)
          whv(:) = whl
        ENDWHERE
        PFB = REAL(w_low)*dp
      ELSE
        whl = wh
        PFB = REAL(w_smooth)*dp
      ENDIF

      wfb = MIN(whl, MAX(np, NINT(1000.0/(nr*dpH))))

      WHERE( (PminC(1) < PH(:)) .AND. (PH(:) <= PminC(1) + PFB))
        whv(:)=NINT((wfb*(PminC(1) + PFB - PH(:)) + whl*(PH(:)-PminC(1)))/PFB)
      ENDWHERE
      WHERE (PH(:) <= PminC(1))
        whv(:) = wfb
      ENDWHERE

      CALL ropp_pp_sliding_polynomial(PH(1:ifb:nr), smh(1:ifb:nr),       &
                                      whv(1:ifb:nr), np, smf(1:ifb:nr),  &
                                      DS=YH(1:ifb:nr))
      CALL ropp_pp_sliding_polynomial(PH(ifb::nr), smh(ifb::nr),       &
                                      whv(ifb::nr), np, smf(ifb::nr),  &
                                      DS=YH(ifb::nr))

      DEALLOCATE(whv)

    CASE DEFAULT

      CALL message(msg_error,'Filtering method '// TRIM(filter) //   &
                             ' not recognised. Check config%filter_method. ')

    END SELECT

    IF (nr > 1) THEN
      CALL ropp_pp_interpol(PH(1::nr), PH, YH(1::nr), YH)
    ENDIF

    smf(1) = 0
    DO i=2,nu
      smf(i) = smf(i-1) + dpH*(YH(i) + YH(i-1))/2
    ENDDO

!-------------------------------------------------------------------------------
! 10. Radio-holographic analysis
!-------------------------------------------------------------------------------

    ! 10.1 Computation of radiohologram with reference signal

    ALLOCATE(UR(nu))
    UR(:) = AH(:)*EXP(Ci*(MODULO(smh(:)-smf(:), 2*pi)))

    DEALLOCATE(smh)
    DEALLOCATE(smf)

    ! 10.2 Radio-holographic analysis

    stepP  = (PH(nu) - PH(1))/(nu-1)
    nrh = ropp_pp_nearest_power2(CEILING(ARH/stepP))

    ALLOCATE(YRH(nrh))
    ALLOCATE(URS(nrh))

    YRH(1) = 0.0_wp
    DO i=2,nrh/2+1
      YRH(i)       = (i-1)*2.0_wp*pi / (k(ic)*nrh*stepP)
      YRH(nrh-i+2) = -YRH(i)
    ENDDO

    YRH(:) = CSHIFT(YRH(:), nrh/2)

    nss = 2*(nu/nrh) - 1

    ALLOCATE(PSS(nss))
    ALLOCATE(ECSS(nss))

    DO iss=1,nss

      IRHmin   = (iss-1)*nrh/2 + 1
      IRHmid   = iss*nrh/2 + 1
      IRHmax   = (iss+1)*nrh/2

      PSS(iss) = PH(IRHmid)
      URS(:)   = UR(IRHmin:IRHmax)

      DO i=1,nrh
        URS(i) = COS(pi*REAL(i-1-nrh/2)/REAL(nrh))*URS(i)
      ENDDO

      ! 10.2.1 Fourier transform

      call ropp_pp_FFT(URS, -1)

      URS(:) = URS(:)/SQRT(REAL(nrh))
      URS(:) = CSHIFT(URS(:), nrh/2)

      ! 10.2.2 Bending angle error estimation [Gorbunov2006 Eq7]

      ECSS(iss) = SUM(YRH(:)**2*ABS(URS(:))**2, Mask = ABS(YRH(:)) < dem) / &
                         SUM(ABS(URS(:))**2, Mask = ABS(YRH(:)) < dem)

    ENDDO

    ! 10.3 Correction for natural line width

    CALL ropp_pp_interpol(PSS(:), PSS0, ECSS(:), dEC)

    ECSS(:) = MAX(0.0_wp, ECSS(:) - dEC)

    DEALLOCATE(YRH)
    DEALLOCATE(URS)
    DEALLOCATE(UR)

!-------------------------------------------------------------------------------
! 11. Determine impact parameter and bending angle
!-------------------------------------------------------------------------------

    ! 11.1 Computation of Y-coordinate

    YH(:) = -YH(:)/k(ic) - dY

    ! 11.2 Transform from y to t

    ALLOCATE(th(nu))

    CALL ropp_pp_interpol(Y(ic,:),YH(i0(ic):i1(ic)),time(:),th(i0(ic):i1(ic)))

    DEALLOCATE(YH)

    ! 11.3 Transform from p~ to Doppler

    ALLOCATE(dh(nu))

    DO i=i0(ic),i1(ic)
      CALL ropp_pp_interpol(time(:), th(i), df(ic,:), dfI)
      CALL ropp_pp_interpol(time(:), th(i), P0(ic,:), P0I)
      dh(i) = dfI + (roc + PH(i) - P0I)/(P_dh(i))
    ENDDO

    DEALLOCATE(P_dh)

    ! 11.4 Compute bending angle, impact parameter and amplitude function

    ALLOCATE(APY(nu))

    DO i=i0(ic),i1(ic)
      t_norm = (th(i) - time(1))/(time(n) - time(1))
      DO m=1,3
        CALL ropp_pp_polynomial(coeff_vgns(:,m), t_norm, xgnsI(m), vgnsI(m))
        CALL ropp_pp_polynomial(coeff_vleo(:,m), t_norm, xleoI(m), vleoI(m))
      ENDDO
      vgnsI = vgnsI/(time(n) - time(1))
      vleoI = vleoI/(time(n) - time(1))

      CALL ropp_pp_geometric_optics(xleoI-r_coc, vleoI, xgnsI-r_coc, vgnsI,   &
                                    dh(i), PH(i), EH(i))

      rgnsI  = SQRT(SUM((xgnsI(:) - r_coc(:))**2))
      rleoI  = SQRT(SUM((xleoI(:) - r_coc(:))**2))
      thetaI = vector_angle(xgnsI(:) - r_coc(:), xleoI(:) - r_coc(:))

      APY(i) = SQRT(SQRT(rgnsI**2 - PH(i)**2)*SQRT(rleoI**2 - PH(i)**2))
      APY(i) = APY(i)*SQRT(rgnsI*rleoI*SIN(thetaI))    ! amplitude fn a2(p,Y)
    ENDDO

    DEALLOCATE(th)
    DEALLOCATE(dh)

    ! 11.5 Re-define Pmin

    IF (opt_SH) THEN
      PminC(ic) = PH(ib(ic)) - roc
    ELSE
      PminC(ic) = PL
    ENDIF

    IF(opt_DL2) THEN
      Pmin = PminC(ic)
    ELSE
      Pmin = PminC(1)
    ENDIF

    WRITE(outstr, '(2(A,F10.3))') 'Pmin = ', Pmin, ' Pmax= ', Pmax
    CALL message(msg_diag, outstr)

    ! 11.5 Redefenition of impact parameter grid

    imax = SUM(MAXLOC(impact(ic,:), impact(ic,:) < roc + Pmax))

    IF (imax > n .OR. imax < 1) THEN
      CALL message(msg_warn, 'No data for CT2 processing - will not process')
      RETURN
    ENDIF

    IF (impact(ic,1) < roc + Pmax) THEN
      imin = 1
    ELSE
      imin = n
    ENDIF

    DO i=imin,imax,SIGN(1,imax-imin)
      impact(ic,i) = roc + (Pmin*REAL(imax-i,wp) + Pmax*REAL(i-imin,wp)) /    &
                     REAL(imax-imin,wp)
    ENDDO

    ! 11.6 Interpolation of bending angle and covariances onto impact grid

    ba_cov(ic,:) = 0.0_wp
    DO i=1,n
      IF (impact(ic,i) - roc < Pmax) THEN
        CALL ropp_pp_interpol(PH(ifb:i1(ic)), impact(ic,i),       &
                              EH(ifb:i1(ic)), bangle(ic,i))
        CALL ropp_pp_interpol(PSS(:), impact(ic,i)-roc, ECSS(:), ba_cov(ic,i))
      ENDIF
    ENDDO
    DEALLOCATE(PSS)
    DEALLOCATE(ECSS)

!-------------------------------------------------------------------------------
! 12. Computation of amplitude
!-------------------------------------------------------------------------------

    ! 12.1 Normalizing amplitude

    AH(i0(ic):i1(ic)) = AH(i0(ic):i1(ic)) * APY(i0(ic):i1(ic)) / APY(i0(ic))

    DEALLOCATE(APY)

    ! 12.2 Computation of output amplitude

!    IF (PRESENT(AP) .AND. PRESENT(A0)) THEN

      ! 12.2.1 Filtering

      dp  = ABS(impact(ic,n) - impact(ic,1)) / (n - 1)
      dpH = nr*ABS(PH(I1(ic)) - PH(i0(ic))) / (i1(ic) - i0(ic))
      wh  = NINT(REAL(w_smooth,wp)*dp / dpH)

      ALLOCATE(Apf(i0(ic):i1(ic)))

      DO i=i0(ic),i1(ic),nr
        imin   = MAX(1,      i-nr/2)
        imax   = MIN(nh(ic), i+nr/2)
        Apf(i) = SUM(AH(imin:imax))/SIZE(AH(imin:imax))
      ENDDO

      CALL ropp_pp_sliding_polynomial(PH(i0(ic):i1(ic):nr),          &
                                      Apf(i0(ic):i1(ic):nr),         &
                                      wh, np, Af(i0(ic):i1(ic):nr))

      DEALLOCATE(Apf)

      DO i=i0(ic),i1(ic)
        CALL ropp_pp_interpol(PH(i0(ic):i1(ic):nr), PH(i),           &
                              Af(i0(ic):i1(ic):nr), Af(i))
      ENDDO

      ! 12.2.1 Output CT amplitude as additional diagnostic

      IF (PRESENT(diag) .AND. ic == 1) THEN

        ALLOCATE(diag%CTimpact((i1(ic)-i0(ic)+10)/10))
        ALLOCATE(diag%CTamplitude((i1(ic)-i0(ic)+10)/10))
        ALLOCATE(diag%CTamplitude_smt((i1(ic)-i0(ic)+10)/10))

        j = 1
        DO i=i0(ic),i1(ic),10
          diag%CTimpact(j) = PH(i)
          diag%CTamplitude(j) = AH(i)
          diag%CTamplitude_smt(j) = Af(i)
          j = j+1
        ENDDO

      ENDIF

      IF (PRESENT(diag) .AND. ic == 2) THEN

        ALLOCATE(diag%CTimpactL2((i1(ic)-i0(ic)+10)/10))
        ALLOCATE(diag%CTamplitudeL2((i1(ic)-i0(ic)+10)/10))
        ALLOCATE(diag%CTamplitudeL2_smt((i1(ic)-i0(ic)+10)/10))

        j = 1
        DO i=i0(ic),i1(ic),10
          diag%CTimpactL2(j) = PH(i)
          diag%CTamplitudeL2(j) = AH(i)
          diag%CTamplitudeL2_smt(j) = Af(i)
          j = j+1
        ENDDO

      ENDIF

      ! 12.2.2 Interpolation to standard grid

      A0(ic) = SUM(AH(i0(ic):i1(ic)), &
                   Mask = (PH(i0(ic):i1(ic)) - roc > Pmax - 2000.0))/  &
                   COUNT(Mask = (PH(i0(ic):i1(ic)) - roc > Pmax - 2000.0))
      DO i=1,n
        IF (impact(ic,i) - roc < Pmax) THEN
          CALL ropp_pp_interpol(PH(i0(ic):i1(ic):nr), impact(ic,i),    &
                                Af(i0(ic):i1(ic):nr), AP(ic,i))
        ELSE
          AP(ic,i) = A0(ic)
        ENDIF
      ENDDO

      !    ENDIF


    DEALLOCATE(Af)
    DEALLOCATE(EH)
    DEALLOCATE(PH)
    DEALLOCATE(AH)

!-------------------------------------------------------------------------------

  ENDDO Channels

!-------------------------------------------------------------------------------
! 13. Synchronizing grids for all channels
!-------------------------------------------------------------------------------

  IF (opt_DL2) THEN

    ! 13.1 Determination of Pmin from best channel

    icw  = 1
    Pmin = PminC(icw)

    WRITE(outstr, '(2X,A,I2,A,F10.3)') 'ICW = ', ICW, '  Pmin = ', Pmin
    CALL message(msg_diag, outstr)

    ! 13.2 Grid redefinition and interpolation

    ALLOCATE(EW(n))
    ALLOCATE(ECW(n))
    ALLOCATE(AW(n))
    ALLOCATE(MW(n))

    DO ic=1,nc

      CALL ropp_pp_interpol(impact(ic,:), impact(icw,:), bangle(ic,:), EW(:), &
                            Cext=.true.)
      CALL ropp_pp_interpol(impact(ic,:), impact(icw,:), ba_cov(ic,:), ECW(:), &
                            Cext=.true.)

      MW(:) = (impact(icw,:) >= roc + PminC(ic) + 1000.0 .AND. &
               impact(icw,:) <= roc + PminC(ic) + 6000.0)

      ndE = COUNT(Mask = MW(:))

      IF (ndE > 20) THEN
        dE = SUM(EW(:) - bangle(icw,:), Mask = MW(:))/ndE
      ELSE
        dE = 0.0_wp
      ENDIF

      WHERE (impact(icw,:) < roc+PminC(ic))
        bangle(ic,:)  = bangle(icw,:) + dE
        ba_cov(ic,:) = ba_cov(icw,:)
      ELSEWHERE
        bangle(ic,:)  = EW(:)
        ba_cov(ic,:) = ECW(:)
      END WHERE

!      IF (PRESENT(diag)) THEN  !AP) .AND. PRESENT(A0)) THEN

        CALL ropp_pp_interpol(impact(ic,:), impact(icw,:), AP(ic,:), AW(:), &
                              Cext=.true.)

        WHERE (impact(icw,:) < PminC(ic))
          AP(ic,:) = 1e-6_wp
        ELSEWHERE
          AP(ic,:) = AW(:)
        END WHERE

!      ENDIF

      impact(ic,:) = impact(icw,:)

    ENDDO

    DEALLOCATE(EW)
    DEALLOCATE(ECW)
    DEALLOCATE(AW)
    DEALLOCATE(MW)

  ENDIF

!-------------------------------------------------------------------------------
! 14. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(df)
  DEALLOCATE(sm)
  DEALLOCATE(P0)
  DEALLOCATE(P_dd)

  DEALLOCATE(Y)
  DEALLOCATE(k)
  DEALLOCATE(nh)
  DEALLOCATE(PminC)
  DEALLOCATE(i0)
  DEALLOCATE(i1)
  DEALLOCATE(ib)
  DEALLOCATE(AP)
  DEALLOCATE(A0)

  CALL message_set_routine(routine)


CONTAINS

!-------------------------------------------------------------------------------
! 15. Transform phase to accumulated phase
!-------------------------------------------------------------------------------

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

END SUBROUTINE ropp_pp_DCT


