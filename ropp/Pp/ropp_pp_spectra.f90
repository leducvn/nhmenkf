! $Id: ropp_pp_spectra.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_spectra(time, phase_L1, phase_L2, phase_LM,   &
                           impact_LM, config, OutRO, filnam)

!****s* Preprocessing/ropp_pp_spectra *
!
! NAME
!    ropp_pp_spectra - Calculation of local spatial spectra
!
! SYNOPSIS
!    call ropp_pp_spectra(time, phase_L1, phase_L2, phase_LM, impact_LM,
!                         config, OutRO, filnam)
!
! DESCRIPTION
!
! INPUTS
!    real(wp), dimension(:)     :: time      ! time of samples (s)
!    real(wp), dimension(:)     :: phase_L1  ! excess phase L1 (m)
!    real(wp), dimension(:)     :: phase_L2  ! excess phase L2 (m)
!    real(wp), dimension(:)     :: phase_LM  ! model excess phase (m)
!    real(wp), dimension(:)     :: impact_LM ! model impact parameter (m)
!    type(PPConfig)             :: config    ! Configuration options
!    logical, optional          :: OutRO     ! Flag to output spectra results
!    character(len=*), optional :: filnam    ! Output file name root
!
! OUTPUT
!    Spectra as functions of (doppler,time) written to output ASCII files
!    and netCDF files.
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
  USE ropp_utils,        ONLY: get_io_unit
  USE ropp_io_types,     ONLY: ROprof
  USE ropp_io
! USE ropp_pp, not_this => ropp_pp_spectra
  USE ropp_pp_types,     ONLY: PPConfig
  USE ropp_pp_constants, ONLY: pi, f_L1, f_L2, c_light

  IMPLICIT NONE

  TYPE(ROprof)                           :: ro_data      ! Hold spectra as RO extra data
  CHARACTER(LEN=256)                     :: text_filnam  ! Name of text file
  CHARACTER(LEN=256)                     :: nc_filnam    ! Name of netCDF file
  CHARACTER(LEN=256)                     :: cline        ! Dummy line in spectral file
  INTEGER, PARAMETER                     :: nsamp=10     ! Stride of output nc data
  INTEGER                                :: id_nx, ix, nx  ! 
  INTEGER                                :: id_ny, iy, ny  ! 
! INTEGER                                :: id_nz, iz, nz  ! ! Commented at 20 July, 2016
  REAL(wp), DIMENSION(:, :), ALLOCATABLE :: x,  y,  z      ! 2d fields from text file

  REAL(wp),     DIMENSION(:), INTENT(in) :: time      ! Time of samples (s)
  REAL(wp),     DIMENSION(:), INTENT(in) :: phase_L1  ! Excess phase L1 (m)
  REAL(wp),     DIMENSION(:), INTENT(in) :: phase_L2  ! Excess phase L2 (m)
  REAL(wp),     DIMENSION(:), INTENT(in) :: phase_LM  ! Model excess phase (m)
  REAL(wp),     DIMENSION(:), INTENT(in) :: impact_LM ! Model impact parameter (m)
  TYPE(PPconfig),             INTENT(in) :: config    ! Configuration options
  LOGICAL,          OPTIONAL, INTENT(in) :: OutRO     ! Flag to output RO spectra
  LOGICAL                                :: l_out     ! Output fields?
  INTEGER                                :: lun       ! Unit number of output text files
  CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: filnam    ! Output file name root

  COMPLEX(wp), DIMENSION(:), ALLOCATABLE :: U         ! Downconverted complx field
  COMPLEX(wp), DIMENSION(:), ALLOCATABLE :: US        ! Sliding spectrum
  REAL(wp),    DIMENSION(:), ALLOCATABLE :: dphase    ! Phase deviation
  REAL(wp),    DIMENSION(:), ALLOCATABLE :: doppler   ! Doppler frequency
  REAL(wp),    DIMENSION(:), ALLOCATABLE :: window    ! Window function
  REAL(wp)                               :: k         ! L1/L2 wave vectors
  INTEGER                                :: i, ic     ! Index counters
  INTEGER                                :: is        ! Sliding window index
  INTEGER                                :: is1       ! Start window position
  INTEGER                                :: isN       ! End window position
  INTEGER                                :: di        ! Index step
  INTEGER                                :: n         ! Number of data points
  INTEGER                                :: ns        ! Number of sliding windows
  INTEGER                                :: ocd       ! Occulatation direction

  COMPLEX(wp), PARAMETER                 :: Ci  = (0.0_wp, 1.0_wp)    ! Sqrt(-1)
  INTEGER,     PARAMETER                 :: nc  = 2                   ! Number of channels
  INTEGER,     PARAMETER                 :: nss = 64                  ! Sliding window width
  LOGICAL                                :: dummy

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

  n = SIZE(time)
  ALLOCATE(dphase(n))
  dummy = config%obs_ok  ! dummy use of otherwise unused argument

  l_out = .FALSE.
  IF ( PRESENT(OutRO) ) l_out = OutRO

!-------------------------------------------------------------------------------
! 3. Sliding spectral analysis
!-------------------------------------------------------------------------------

  ! 3.1 Occultation direction (-1 = setting, +1 = rising)

  ocd = NINT(SIGN(1.0_wp, impact_LM(n)-impact_LM(1)))

  SELECT CASE (ocd)
  CASE (-1)
    is1 = 1 + nss/2
    isN = n - (nss/2 - 1)
    di  = 1
  CASE (+1)
    is1 = n - nss/2
    isN = 1 + (nss/2 - 1)
    di  = -1
  END SELECT

  ns = di*isN - di*is1 + 1

  ALLOCATE(U(n))
  ALLOCATE(US(nss))
  ALLOCATE(doppler(nss))
  ALLOCATE(window(nss))

  ! 3.2 Doppler frequency

  doppler(1) = 0.0_wp
  DO i=2,nss/2+1
    doppler(i)       = (i-1.0_wp)/(nss*(time(2)-time(1)))
    doppler(nss-i+2) = -doppler(i)
  END DO
  doppler(:) = CSHIFT(doppler(:), nss/2)

  ! 3.3 Window function

  DO i=1,nss
    window(i) = (1.0 + COS(2.0_wp*Pi*REAL(i-1-nss/2,wp)/REAL(nss,wp)))/2.0_wp
  END DO

  ! 3.4 Compute spectrum

  DO ic=1,nc

    IF (l_out) THEN

      lun = get_io_unit()

      IF (PRESENT(filnam)) THEN
        IF (TRIM(filnam) == "ropp_pp_spectra") THEN ! Default in ropp_pp_spectra_tool
           IF (ic == 1) text_filnam = "ROanalysis_dt_L1.dat"
           IF (ic == 2) text_filnam = "ROanalysis_dt_L2.dat"
        ELSE
           IF (ic == 1) text_filnam = TRIM(filnam)//"_dt_L1.dat"
           IF (ic == 2) text_filnam = TRIM(filnam)//"_dt_L2.dat"
        ENDIF
      ELSE ! Use original default output filenames
           IF (ic == 1) text_filnam = "ROanalysis_dt_L1.dat"
           IF (ic == 2) text_filnam = "ROanalysis_dt_L2.dat"
      ENDIF

      OPEN(UNIT=lun, FILE=TRIM(ADJUSTL(text_filnam)))

      WRITE(lun, '(A / A,I5,A,I5,A,A,I2)')                         &
                 'VARIABLES = "doppler (Hz)", "time(s)", "ln|U|"', &
                 'ZONE I=', nss, ' J=', ns, ' F=POINT', ' OCD=', ocd

    ENDIF

    ! 3.4.1 Compute wave number

    IF (ic == 1) k = 2.0_wp * pi * f_L1 / C_Light
    IF (ic == 2) k = 2.0_wp * pi * f_L2 / C_Light

    ! 3.4.2 Compute phase deviation

    DO i=1,n
      IF (ic == 1) dphase(i) = k*(phase_L1(i) - phase_LM(i))
      IF (ic == 2) dphase(i) = k*(phase_L2(i) - phase_LM(i))
    END DO

    ! 3.4.3 Complex field

    U(:) = EXP(Ci*(dphase(:)))

    DO is = is1, isN, di

      US(:) = U(is - di*nss/2:is + di*nss/2 - di:di)*window(:)

      CALL ropp_pp_FFT(US, 1)

      US(:) = US(:)/SQRT(REAL(nss,wp))
      US(:) = US(:)/MAXVAL(ABS(US(:)))
      US(:) = CSHIFT(US(:), nss/2)

      IF (l_out) THEN
        DO i=1,nss
!          WRITE(lun,'(E13.5,1X,F10.4,1X,E13.5)')    &
          WRITE(lun, '(E20.10,1X,E20.10,1X,E20.10)') &
                     doppler(i), time(is), LOG(ABS(US(i)))
        END DO
      ENDIF

    END DO

    ! 3.4.4 Output (thinned spectra) to netCDF files

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
                          name='freq',                                    &
                          long_name='Frequency offset',                   &
                          units='Hz',                                     &
                          range=(/ -1000.0_wp, 1000.0_wp /),              &
                          data=x(1:nx, 1))

      CALL ropp_io_addvar(ro_data,                                        &
                          name='stime',                                   &
                          long_name='Time since start of occultation',    &
                          units='s',                                      &
                          range=(/ -1000.0_wp, 1000.0_wp /),              &
                          data=y(1, 1:ny:nsamp))

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

  END DO  ! loop over channels

  DEALLOCATE(dphase)
  DEALLOCATE(U)
  DEALLOCATE(US)
  DEALLOCATE(doppler)
  DEALLOCATE(window)


END SUBROUTINE ropp_pp_spectra
