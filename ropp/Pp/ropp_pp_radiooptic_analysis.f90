! $Id: ropp_pp_radiooptic_analysis.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_radiooptic_analysis(time, r_leo, r_gns, r_coc, roc,    &
                                       phase_LM, phase, snr, PA, PD, OutRO, filnam)

!****s* Preprocessing/ropp_pp_radiooptic_analysis *
!
! NAME
!    ropp_pp_radiooptic_analysis - Radiooptical analysis of radio occultation
!                                  data - calculation of local spatial spectra
!                   
! SYNOPSIS
!    call ropp_pp_radiooptic_analysis(time, r_leo, r_gns, r_coc, roc, 
!                                         phase_LM, phase, snr, PA, PD, OutRO, filnam)
! 
! DESCRIPTION
!     Calculation of local spatial spectra
!
! INPUTS
!    real(wp), dimension(:)     :: time     ! relative time of samples (s)
!    real(wp), dimension(:,:)   :: r_leo    ! cartesian LEO coordinates (m)
!    real(wp), dimension(:,:)   :: r_gns    ! cartesian GPS coordinates (m)
!    real(wp), dimension(:)     :: r_coc    ! cartesian centre curvature (m)
!    real(wp)                   :: roc      ! radius of curvature (m)
!    real(wp), dimension(:)     :: phase_LM ! model excess phase (m)
!    real(wp), dimension(:,:)   :: phase    ! L1 and L2 excess phase (m)
!    real(wp), dimension(:,:)   :: snr      ! L1 and L2 amplitude
!    logical, optional          :: OutRO    ! Flag to output spectra results
!    character(len=*), optional :: filnam   ! Output file name root
!
! OUTPUT
!    real(wp), dimension(:,:)   :: PA       ! Average impact parameter (m)
!    real(wp), dimension(:,:)   :: PD       ! RMS deviation of impact parameter
!
!    Spectra as functions of (impact,bangle) written to output ASCII files
!    and netCDF files.
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

  USE typesizes,         ONLY: wp => EightByteReal
  USE ropp_utils,        ONLY: Get_IO_Unit
  USE ropp_io_types,     ONLY: ROprof
  USE ropp_io
  USE messages
! USE ropp_pp, not_this => ropp_pp_radiooptic_analysis
  USE ropp_pp
  USE ropp_utils,        ONLY: plane_coordinates

  IMPLICIT NONE

  TYPE(ROprof)                          :: ro_data      ! Hold spectra as RO extra data
  CHARACTER(LEN=256)                    :: text_filnam  ! Name of text file
  CHARACTER(LEN=256)                    :: nc_filnam    ! Name of netCDF file
  CHARACTER(LEN=256)                    :: cline        ! Dummy line in spectral file
  INTEGER, PARAMETER                    :: nsamp=10     ! Stride of output nc data
  INTEGER                               :: id_nx, ix, nx  ! 
  INTEGER                               :: id_ny, iy, ny  ! 
! INTEGER                               :: id_nz, iz, nz  ! ! Commented at 20 July, 2016
  REAL(wp), DIMENSION(:, :), ALLOCATABLE:: x,  y,  z      ! 2d fields from text file

  REAL(wp), DIMENSION(:),   INTENT(in)  :: time     ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo    ! Cartesian LEO coords (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns    ! Cartesian GPS coords (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc    ! Centre curvature coord (m)
  REAL(wp),                 INTENT(in)  :: roc      ! Radius of curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM ! Model excess phase (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: phase    ! L1 and L2 excess phase (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: snr      ! L1 and L2 amplitude
  LOGICAL,        OPTIONAL, INTENT(in)  :: OutRO    ! Flag to output RO spectra
  LOGICAL                               :: l_out    ! Output fields?
  CHARACTER(LEN=*), OPTIONAL,INTENT(in) :: filnam   ! Output file name root
  INTEGER                               :: lun      ! Unit number of output text files

  REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: PA  ! Avg impact parameter
  REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: PD  ! RMS deviation IP

  REAL(wp), DIMENSION(:), ALLOCATABLE   :: xleo     ! X coordinates of LEO
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: yleo     ! Y coordinates of LEO
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: xgns     ! X coordinates of GPS
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: ygns     ! Y coordinates of GPS

  INTEGER,  PARAMETER :: ncp = 2   ! Polynomial degree for calculating phase
  INTEGER,  PARAMETER :: nv  = 5   ! Polynomial degree for calculating velocity

  REAL(wp), PARAMETER :: DPM = 3000.0_wp ! Maximum impact parameter deviation
  REAL(wp), PARAMETER :: DE  = 0.0025_wp ! Spectral width of bending angle (rad)
  REAL(wp), PARAMETER :: DP  = 4000.0_wp ! Spectral width of impact param (m)
  REAL(wp)            :: DA              ! Sliding aperture step (m)
  REAL(wp)            :: zeta            ! Refractive parameter |1 - L de/dp|
  INTEGER             :: n               ! Number of data points
  INTEGER             :: nc              ! Number of channels
!  INTEGER             :: ny              ! Number of apertures
  INTEGER             :: nr              ! Number of points inside aperture
  INTEGER             :: np              ! Point counter
  INTEGER             :: ic              ! Channel number
  INTEGER             :: i, j            ! Array index

  REAL(wp), DIMENSION(3)                :: AX, AY    ! Occ plane XY-basis vector
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: t_norm     ! Normalised time
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KV         ! Regression matrix for v
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_vleo ! Regression coeffs vleo
  REAL(wp), DIMENSION(0:nv,3)           :: coeff_vgns ! Regression coeffs vgns
  REAL(wp)                              :: ymin, ymax ! Minimum, maximum Y
  REAL(wp), DIMENSION(3)                :: XG         ! GPS position
  REAL(wp), DIMENSION(3)                :: VG         ! GPS velocity
  REAL(wp), DIMENSION(3)                :: XL         ! LEO position
  REAL(wp), DIMENSION(3)                :: VL         ! LEO velocity
  REAL(wp), DIMENSION(3)                :: UGL        ! GPS-LEO direction

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: k    ! L1/L2 wave vectors
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: DY   ! Step of Y-grid in aperture
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: DT   ! Step of T-grid in aperture
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: FY   ! Y-aperture size for L1/L2    
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: FT   ! T-aperture size for L1/L2    
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: SMC  ! Model excess phase [ch,t]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: P    ! Accumulated phase [ch,t]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: TRH  ! Time grid for sliding aperture
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: YRH  ! Sliding aperture centers 
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: TR   ! T-grid for spectral analysis
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PR   ! Interpolated phase in aperture
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: SMR  ! Model excess phase in aperture
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: AR   ! Interpol amplitude in aperture
  COMPLEX(wp), DIMENSION(:), ALLOCATABLE :: UR  ! Complex field/spatial spectra
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: UA   ! Log(Abs(UR(:))
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: KY   ! Y-component of wave vector
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: PK   ! Impact parameter for KY
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: EK   ! Bending angle for KY
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PARH ! Average impact parameter
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: PDRH ! RMS deviation impact parameter
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: P2RH ! Spline coefficients
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: IP   ! Point numbers
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: BP   ! Matrix of basic polynomials
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: CP   ! Polynomial regression coeffs

  REAL(wp)           :: FYmax    ! Maximum aperture
  REAL(wp)           :: YRmin    ! Lower Y of aperture
  REAL(wp)           :: YRmax    ! Upper Y of aperture
  REAL(wp)           :: TRmin    ! Lower T of aperture
  REAL(wp)           :: TRmax    ! Upper T of aperture
  REAL(wp)           :: P0       ! Phase trend
  REAL(wp)           :: doppler  ! Relative Doppler frequency shift
  REAL(wp)           :: Umin     ! Minimum threshold of spectrum amplitude
  REAL(wp)           :: Umax     ! Spectrum maximum on grid
  LOGICAL            :: LMax     ! Indicator of local maximum
  LOGICAL            :: MMax     ! Indicator of main maximum
  REAL(wp)           :: KYext    ! Refined location of spectrum maximum
  REAL(wp)           :: UAext    ! Refined spectrum maximum
  REAL(wp)           :: Eext     ! Refraction angle at spectrum maximum
  REAL(wp)           :: Pext     ! Impact parameter at spectrum maximum

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: EARH ! Average impact parameter
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: EDRH ! RMS deviation impact parameter
  CHARACTER(len =  256) :: outstr
  CHARACTER(len = 256)                :: routine

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_radiooptic_analysis')

  l_out = .FALSE.
  IF ( PRESENT(OutRO) ) l_out = OutRO

!-------------------------------------------------------------------------------
! 3. Coordinate transforms
!-------------------------------------------------------------------------------

  ! 3.1. Array allocation
  
  n = SIZE(time)
  nc = SIZE(snr,1)

  IF (nc == 1) zeta = 300.0_wp
  IF (nc == 2) zeta = 1.0_wp

  IF (nc == 1) DA = 1000.0_wp
  IF (nc == 2) DA = 300.0_wp

  ALLOCATE(xgns(n))
  ALLOCATE(xleo(n))
  ALLOCATE(ygns(n))
  ALLOCATE(yleo(n))  
  ALLOCATE(CP(0:ncp))

  ! 3.2 Calculation of coordinates in occultation plane

  CALL plane_coordinates(r_leo,r_gns,r_coc,roc, xleo, yleo, xgns, ygns, ax, ay)
  yleo(:)  = yleo(:) - roc
  ygns(:)  = ygns(:) - roc

  ! 3.3. Calculation of regression coefficients
  
  ALLOCATE(t_norm(n))
  ALLOCATE(KV(n, 0:nv))

  t_norm(:) = (time(:) - time(1))/(time(n) - time(1))

  CALL ropp_pp_init_polynomial(t_norm, KV)
  
  DO i=1,3
     CALL ropp_pp_regression(KV, r_leo(:,i), coeff_vleo(:,i))
     CALL ropp_pp_regression(KV, r_gns(:,i), coeff_vgns(:,i))
  ENDDO

  DO i=1,3
    CALL ropp_pp_residual_regression(KV, t_norm, r_leo(:,i), coeff_vleo(:,i))
    CALL ropp_pp_residual_regression(KV, t_norm, r_gns(:,i), coeff_vgns(:,i))
  ENDDO

  DEALLOCATE(t_norm)
  DEALLOCATE(KV)

!-------------------------------------------------------------------------------
! 4. Calculation of complex field
!-------------------------------------------------------------------------------

  ! 4.1 Memory allocation
  
  ALLOCATE(P(nc,1:N))
  ALLOCATE(SMC(nc,1:N))
  ALLOCATE(k(nc))
  ALLOCATE(DY(nc))
  ALLOCATE(DT(nc))
  ALLOCATE(FY(nc))
  ALLOCATE(FT(nc))
  
  ! 4.2 Calculation of accumulated phase

  DO ic=1,nc
     IF (nc == 1) THEN
        k(ic) = 2.0_wp*Pi*f_L2/C_Light
     ELSE
        IF(ic == 1) k(ic) = 2.0_wp*Pi*f_L1/C_Light
        IF(ic == 2) k(ic) = 2.0_wp*Pi*f_L2/C_Light
     ENDIF
     P(ic,:) = k(ic)*phase(ic,:)
  ENDDO

  ! 4.3 Calculation of grid limits

  FYmax = MAXVAL(SQRT(zeta*xleo(n)*2.0_wp*Pi/k(:)))
  Ymin  = MINVAL(yleo) + FYmax/2.0_wp
  Ymax  = MAXVAL(yleo) - FYmax/2.0_wp
  ny    = 1 + CEILING((Ymax-Ymin)/DA)

  ! 4.4 Calculation of grid
  
  ALLOCATE(YRH(ny))
  ALLOCATE(TRH(ny))
  ALLOCATE(PARH(ny))
  ALLOCATE(PDRH(ny))
  ALLOCATE(P2RH(ny))
  ALLOCATE(EARH(ny))
  ALLOCATE(EDRH(ny))
  
  DO i=1,ny
     YRH(i)  = (Ymin*(ny-i) + Ymax*(i-1.0_wp))/(ny-1.0_wp)
  ENDDO
  
  CALL ropp_pp_interpol(yleo(1:n), yrh, time(1:n), trh)
  
  ! 4.5 Calculation of phase model   
    
  DO ic=1,nc
     IF (nc > 1) SMC(ic,:) = phase_LM(:)
     IF (nc == 1) SMC(ic,:) = phase(ic,:)
  ENDDO
  
!-------------------------------------------------------------------------------
! 5. Radiooptical analysis
!-------------------------------------------------------------------------------
  
  DO ic=1,nc

     ! 5.1 Calculation of internal and external scales

     DY(ic) = Pi/(2.0_wp*k(ic)*MAX(DE,DP/xleo(n)))
     FY(ic) = SQRT(zeta*xleo(n)*2.0_wp*Pi/k(ic))
     nr     = ropp_pp_nearest_power2(CEILING(FY(ic)/DY(ic)))
     Umin   = 0.0001_wp*nr*MAXVAL(snr(ic,:))

     WRITE (outstr,'(2X,A,I1,A,F10.3,2X,A,I5)')  &
          'L', ic, ' aperture: ', FY(ic),     &
          'FFT points: ', nr
     CALL message(msg_diag, outstr)

     IF ( l_out ) THEN

       lun = get_io_unit()

       IF (PRESENT(filnam)) THEN
         IF (TRIM(filnam) == "ropp_pp_spectra") THEN ! Default in ropp_pp_spectra_tool
            IF (ic == 1) text_filnam = "ROanalysis_ep_L1.dat"
            IF (ic == 2) text_filnam = "ROanalysis_ep_L2.dat"
         ELSE
            IF (ic == 1) text_filnam = TRIM(filnam)//"_ep_L1.dat"
            IF (ic == 2) text_filnam = TRIM(filnam)//"_ep_L2.dat"
         ENDIF
       ELSE ! Use original default output filenames
            IF (ic == 1) text_filnam = "ROanalysis_ep_L1.dat"
            IF (ic == 2) text_filnam = "ROanalysis_ep_L2.dat"
       ENDIF

       OPEN(UNIT=lun, FILE=TRIM(ADJUSTL(text_filnam)))

       WRITE(lun, '(A / A,I5,A,I5,A)')                                 &
                  'VARIABLES = "bangle (rad)", "impact (km)", "ln|U|"',&
                  'ZONE I=', nr, ' J=', ny, ' F=POINT'
     ENDIF

     ! 5.2 Array allocation
     
     ALLOCATE(TR(0:nr-1))
     ALLOCATE(PR(0:nr-1))
     ALLOCATE(SMR(0:nr-1))
     ALLOCATE(AR(0:nr-1))
     ALLOCATE(UR(0:nr-1))
     ALLOCATE(UA(ny,0:nr-1))
     ALLOCATE(KY(0:nr-1))
     ALLOCATE(PK(ny,0:nr-1))
     ALLOCATE(EK(ny,0:nr-1))
     ALLOCATE(IP(ny,0:nr-1))
     ALLOCATE(BP(0:nr-1,0:ncp))
     np = 0

     ! 5.3 Radiooptical analysis in apertures

     DO i=1,ny

        ! 5.3.1 Calculation of aperture limits
        
        YRmin  = YRH(i) - FY(ic)/2.0_wp
        YRmax  = YRH(i) + FY(ic)/2.0_wp

        CALL ropp_pp_interpol(yleo(1:n), YRmin, time(1:n), TRmin)
        CALL ropp_pp_interpol(yleo(1:n), YRmax, time(1:n), TRmax)

        DT(ic) = 2.0_wp*MIN(ABS(TRmax - TRH(i)), ABS(TRmin - TRH(i)))
        TRmax  = TRH(i) + DT(ic)/2.0_wp
        TRmin  = TRH(i) - DT(ic)/2.0_wp
                
        ! 5.3.2 Interpolation of phase and amplitude

        DO j=0,nr-1
           TR(j) = (TRmin*(nr-j) + TRmax*j)/nr - TRH(i)

           CALL ropp_pp_interpol(time, trh(i)+tr(j), P(ic,:), PR(j))
           CALL ropp_pp_interpol(time, trh(i)+tr(j), SMC(1,:), SMR(j))
           CALL ropp_pp_interpol(time, trh(i)+tr(j), snr(ic,:), AR(j))
           
        ENDDO
        
        ! 5.3.3 Phase regression

        CALL ropp_pp_init_polynomial(TR, BP)
        CALL ropp_pp_regression(BP, SMR, CP)
        CP(:) = k(ic)*CP(:)

        ! 5.3.4 Downconversion of frequency and calculation of complex signal

        DO j=0,nr-1
           
           IF (nc == 1) CALL ropp_pp_polynomial(CP, TR(j), P0)

           IF (nc == 2) P0 = k(ic)*SMR(j)

           PR(j) = PR(j) - P0
           UR(j) = AR(j)*EXP(CMPLX(0.0_wp, PR(j), wp))
           UR(j) = UR(j)*COS(Pi*TR(j)/(TRmax-TRmin))**2

        ENDDO
        
        ! 5.3.5 Calculation of frequencies

        KY(0) = CP(1)
        DO j=1,nr/2-1
           KY(j)    = CP(1) + 2.0_wp*Pi*j/(TRmax-TRmin)
           KY(nr-j) = CP(1) - 2.0_wp*Pi*j/(TRmax-TRmin)
        ENDDO
        KY(nr/2) = CP(1) - Pi*nr/(TRmax-TRmin)
        
        ! 5.3.6 Calculation of bending angles

        CALL ropp_pp_interpolate_trajectory(time, coeff_vleo, coeff_vgns,   &
                                            r_coc, TRH(i), XL, VL, XG, VG)

        UGL(:) = (XL(:) - XG(:))/SQRT(SUM((XL(:)-XG(:))**2))

        DO j=0,nr-1

           doppler = -KY(j)/(C_Light*k(ic)) +                             &
                         (DOT_PRODUCT(VG,UGL)-DOT_PRODUCT(VL,UGL)) /      &
                         (C_Light-DOT_PRODUCT(VG,UGL))
           
           CALL ropp_pp_geometric_optics(XL-r_coc, VL, XG-r_coc, VG,      &
                                         doppler, PK(i,j), EK(i,j))
           NP      = NP+1
           IP(i,j) = NP-1
           
        ENDDO

        ! 5.3.7 Fourier analysis

        CALL ropp_pp_FFT(UR, -1)

        WHERE(ABS(UR(:)) < Umin)
           UR(:) = Umin
        endwhere
        
        Umax    = MAXVAL(ABS(UR))
        UR(:)   = UR(:)/Umax
        UA(i,:) = LOG(ABS(UR(:)))
        
        UA(i,:) = CSHIFT(UA(i,:), -nr/2)
        UR(:)   = CSHIFT(UR,      -nr/2)
        KY(:)   = CSHIFT(KY,      -nr/2)
        EK(i,:) = CSHIFT(EK(i,:), -nr/2)
        PK(i,:) = CSHIFT(PK(i,:), -nr/2)
        
        ! 5.3.8 Computation of weighted values
      
        PARH(i) = SUM(PK(i,:)*ABS(UR(:))**2) / SUM(ABS(UR(:))**2)
        WHERE (ABS(PK(i,:)-PARH(i)) > DPM)
           UR(:) = 0.0
         endwhere
        PDRH(i) = SQRT(SUM((PK(i,:)-PARH(i))**2*ABS(UR(:))**2) /        &
                     SUM(ABS(UR(:))**2))

        EARH(i) = SUM(EK(i,:)*ABS(UR(:))**2) / SUM(ABS(UR(:))**2)
        EDRH(i) = SQRT(SUM((EK(i,:)-EARH(i))**2*ABS(UR(:))**2) /        &
                     SUM(ABS(UR(:))**2))

        IF ( l_out ) THEN
!          WRITE(lun, '(E13.5,1X,F10.4,1X,E13.5)') &
          WRITE(lun, '(E20.10,1X,E20.10,1X,E20.10)') &
                     (EK(i,j), (PK(i,j)-roc)/1000.0_wp, UA(i,j), j=0,nr-1)
        ENDIF

        ! 5.3.9 Writing spectral maxima

        DO j=1,nr-2

           MMax = (j == SUM(MAXLOC(UA(i,:))) - 1 + LBOUND(UA, Dim=2))
           LMax = (UA(i,j) > -2)             .AND. &
                (UA(i,j) > UA(i,j-1) + 0.05) .AND. &
                (UA(i,j) > UA(i,j+1) + 0.05)
           
           IF (MMax .OR. LMax) THEN

              CALL extremum(KY(j-1:j+1), UA(i,j-1:j+1), KYext, UAext)

              doppler = -KYext/(C_Light*k(ic)) +                        &
                         (DOT_PRODUCT(VG,UGL) - DOT_PRODUCT(VL,UGL)) /  &
                         (C_Light - DOT_PRODUCT(VG,UGL))
              
              CALL ropp_pp_geometric_optics(XL-r_coc, VL, XG-r_coc, VG, &
                                            doppler, Pext, Eext)

           END IF
           
         END DO
        
       END DO

   ! 5.4 Output argument computation
     
    IF (PRESENT(PA)) THEN

      CALL ropp_pp_init_spline(TRH(:), PARH(:), P2RH(:))

      DO i=1,n
        CALL ropp_pp_interpol_spline(TRH(:),PARH(:),P2RH(:),time(i),PA(ic,i))
      ENDDO

    ENDIF

    IF (PRESENT(PD)) THEN

      CALL ropp_pp_init_spline(TRH(:), PDRH(:), P2RH(:))
      
      DO i=1,n
        CALL ropp_pp_interpol_spline(TRH(:),PDRH(:),P2RH(:),time(i),PD(ic,i))
      ENDDO
      
    ENDIF

   ! 5.5 Array deallocation

    DEALLOCATE(TR)
    DEALLOCATE(PR)
    DEALLOCATE(SMR)
    DEALLOCATE(AR)
    DEALLOCATE(UR)
    DEALLOCATE(UA)
    DEALLOCATE(KY)
    DEALLOCATE(PK)
    DEALLOCATE(EK)
    DEALLOCATE(IP)
    DEALLOCATE(BP)

    ! 5.5 Output (thinned spectra) to netCDF files

    IF (l_out) THEN

      ! Simplest to read 2D fields from the test files just written

      REWIND (UNIT=lun)

      READ (lun, '(a80)') cline
      READ (lun, '(a80)') cline

      id_nx = INDEX( cline, '=')
      READ ( cline(id_nx+1:id_nx+6), '(i5)' ) nx

      id_ny = INDEX( cline(id_nx+1:), '=')
      READ ( cline(id_nx+id_ny+1:id_nx+id_ny+6), '(i5)' ) ny

      ALLOCATE (x(nx, ny), y(nx, ny), z(nx, ny))

      DO iy=1,ny
        DO ix=1,nx
          READ (lun, *) x(ix, iy), y(ix, iy), z(ix, iy)
        END DO
      END DO

      CLOSE (UNIT=lun)

      ! Store as (thinned) ROprof 'extra data'

      CALL ropp_io_init(ro_data%lev2c, 1)

      CALL ropp_io_addvar(ro_data,                                        &
                          name='bangle',                                  &
                          long_name='Bending angle',                      &
                          units='rad',                                    &
                          range=(/ -1000.0_wp, 1000.0_wp /),              &
                          data=x(1:nx, 1:ny:nsamp))

      CALL ropp_io_addvar(ro_data,                                        &
                          name='impact',                                  &
                          long_name='Impact parameter',                   &
                          units='km',                                     &
                          range=(/ -1000.0_wp, 1000.0_wp /),              &
                          data=y(1:nx, 1:ny:nsamp))

      CALL ropp_io_addvar(ro_data,                                        &
                          name='amp',                                     &
                          long_name='log of modulus of signal amplitude', &
                          units='',                                       &
                          range=(/ -1000.0_wp, 1000.0_wp /),              &
                          data=z(1:nx, 1:ny:nsamp))

      DEALLOCATE (x, y, z)

      ! Write out

      nc_filnam = TRIM(ADJUSTL(text_filnam(1:INDEX(text_filnam, '.', BACK=.TRUE.)))) // 'nc'

      CALL ropp_io_write(ro_data, nc_filnam, ranchk=.TRUE.)

      CALL ropp_io_free(ro_data)

    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(XLEO)
  DEALLOCATE(YLEO)
  DEALLOCATE(XGNS)
  DEALLOCATE(YGNS)

  DEALLOCATE(CP)
  DEALLOCATE(P)
  DEALLOCATE(SMC)
  DEALLOCATE(k)
  DEALLOCATE(DY)
  DEALLOCATE(DT)
  DEALLOCATE(FY)
  DEALLOCATE(FT)
  DEALLOCATE(YRH)
  DEALLOCATE(TRH)
  DEALLOCATE(PARH)
  DEALLOCATE(PDRH)
  DEALLOCATE(P2RH)
  DEALLOCATE(EARH)
  DEALLOCATE(EDRH)

  CALL message_set_routine(routine)

CONTAINS

!-------------------------------------------------------------------------------
! 7. Calculate extremum of square polynomial
!-------------------------------------------------------------------------------

  SUBROUTINE extremum(x, y, xext, yext)

    USE typesizes, ONLY: wp => EightByteReal
    IMPLICIT NONE

    ! 7.1 Declarations

    REAL(wp), DIMENSION(0:2), INTENT(in)  :: X     ! Arguments
    REAL(wp), DIMENSION(0:2), INTENT(in)  :: Y     ! Function values
    REAL(wp),                 INTENT(out) :: Xext  ! Extremum location
    REAL(wp),                 INTENT(out) :: Yext  ! Extremum value
    
    REAL(wp), DIMENSION(0:2)              :: A     ! Interpolation coefficients
    
    ! 7.2 Calculation of interpolation coefficients

    A(0) = Y(0)*X(1)*X(2)/((X(0)-X(1))*(X(0)-X(2))) + &
         Y(1)*X(0)*X(2)/((X(1)-X(0))*(X(1)-X(2))) + &
         Y(2)*X(0)*X(1)/((X(2)-X(0))*(X(2)-X(1)))
    
    A(1) = -Y(0)*(X(1)+X(2))/((X(0)-X(1))*(X(0)-X(2)))  &
         -Y(1)*(X(0)+X(2))/((X(1)-X(0))*(X(1)-X(2)))  &
         -Y(2)*(X(0)+X(1))/((X(2)-X(0))*(X(2)-X(1)))
    
    A(2) = Y(0)/((X(0)-X(1))*(X(0)-X(2))) + &
         Y(1)/((X(1)-X(0))*(X(1)-X(2))) + &
         Y(2)/((X(2)-X(0))*(X(2)-X(1)))

    ! 7.3 Calculation of extremum

    Xext = -A(1)/(2.0_wp*A(2))
    Yext = A(0) + Xext*(A(1) + Xext*A(2))
    
  END SUBROUTINE extremum

END SUBROUTINE ropp_pp_radiooptic_analysis

