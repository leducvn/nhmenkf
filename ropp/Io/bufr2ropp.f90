! $Id: bufr2ropp_mod.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE bufr2ropp

!****m* bufr2ropp/bufr2ropp *
!
! NAME
!   bufr2ropp    (bufr2ropp_mod.f90)
!
! SYNOPSIS
!   Module defining fixed values & common subroutines/functions for the
!   bufr2ropp main program (ECMWF or MetDB versions)
!
!   USE bufr2ropp
!
! USED BY
!   bufr2ropp  (bufr2ropp_ec or bufr2ropp_mo)
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

! Modules

  USE messages
  USE typesizes, dp => EightByteReal

! Fixed values

  INTEGER, PARAMETER :: RODescr = 310026    ! Table D descriptor for RO

  INTEGER, PARAMETER :: nb      =   500000  ! Max. length of BUFR message (bytes)
  INTEGER, PARAMETER :: nd      =   200000  ! Max. no. of descriptors
  INTEGER, PARAMETER :: no      =        1  ! Max. no. of observations
  INTEGER, PARAMETER :: nv      =    nd*no  ! Max. no. of data values

! Default (missing) data values

  INTEGER,  PARAMETER :: NMDFV = -9999999     ! Integer missing data flag value (MetDB)
  INTEGER,  PARAMETER :: NVIND = 2147483647   ! Integer missing data flag value (ECMWF)
  REAL,     PARAMETER :: RMDFV = -9999999.0   ! Real    missing data flag value (MetDB)
  REAL(dp), PARAMETER :: RVIND = 1.7E38_dp    ! Real    missing data flag value (ECMWF)

! MetDB I/O interface

  INTEGER, PARAMETER  :: Input  = 1        ! I/O input mode (r)

CONTAINS
!-------------------------------------------------------------------------

SUBROUTINE Usage()

!****s* bufr2ropp/Usage *
!
! NAME
!   Usage
!
! SYNOPIS
!   USE bufr2ropp
!   CALL Usage()
!
! INPUTS
!   None
!
! OUTPUTS
!   Summary usage text to stdout
!
! CALLED BY
!   GetOptions
!
! DESCRIPTION
!   Prints a summary of program usage (help) to stdout.
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

  PRINT *, 'Purpose:'
  PRINT *, '  Decode RO BUFR messages to ROPP (netCDF) format.'
  PRINT *, 'Usage:'
  PRINT *, ' > bufr2ropp bufr_file [bufr_file...] [-o ropp_file] [-m] [-a]'
  PRINT *, '                       [-f first] [-n number] [-d] [-h|?] [-v]'
  PRINT *, 'Input:'
  PRINT *, '  One or more files containing BUFR messages or WMO bulletins.'
  PRINT *, '  Any non-RO messages will be skipped.'
  PRINT *, 'Output:'
  PRINT *, '  One or more ROPP files.'
  PRINT *, 'Options:'
  PRINT *, '  -o ROPP netCDF output file name. Mandatory argument if used'
  PRINT *, '  -m write ROPP output as a multifile, i.e. if there'
  PRINT *, '     is more than one RO BUFR message in the input, they'
  PRINT *, '     are decoded into a single netCDF output file'
  PRINT *, '  -a append new profiles to file specified by -o. (-a implies -m)'
  PRINT *, '  -f specify the first RO message to decode'
  PRINT *, '  -n specify the number of RO messages to decode'
  PRINT *, '  -d outputs additional diagnostics to stdout'
  PRINT *, '  -h this help'
  PRINT *, '  -v version information'
  PRINT *, 'Defaults:'
  PRINT *, '  Input  file name : ropp.bufr'
  PRINT *, '  Output file name : from (first) occultation ID <occid>.nc'
  PRINT *, '  Output mode      : one netCDF output file per input RO BUFR message'
  PRINT *, '  Append mode      : exising file is overwritten'
  PRINT *, '  first   : 1'
  PRINT *, '  number  : 999999'
  PRINT *, 'See bufr2ropp(1) for details.'
  PRINT *, ''
END SUBROUTINE Usage
!-------------------------------------------------------------------------

SUBROUTINE GetOptions ( BUFRdsn,      & ! (out)
                        nfiles,       & ! (out)
                        ROPPdsn,      & ! (out)
                        multi,        & ! (out)
                        newfile,      & ! (out)
                        fMsgToDecode, & ! (out)
                        nMsgToDecode )  ! (out)

!****s* bufr2ropp/GetOptions *
!
! NAME
!   GetOptions
!
! SYNOPSIS
!   Get command line information & options or set defaults
!
!    USE bufr2ropp
!    CHARACTER (LEN=100), DIMENSION(:), ALLOCATABLE :: BUFRdsn
!    CHARACTER (LEN=100) :: roppdsn
!    INTEGER :: nfiles, fMsgToDecode, nMsgToDecode
!    LOGICAL :: multi, newfile
!    nfiles=IARGC()
!    ALLOCATE ( bufrdsn(nfiles) )
!    CALL getoptions ( bufrdsn, nfiles, roppdsn, &
!                      multi, newfile, fMsgToDecode, nMsgToDecode )
!   On the command line:
!    >  bufr2ropp bufr_file [...bufr_file...] [-o ropp_file]
!                           [-m] [-a]   [-f first] [-n number]
!                           [-d] [-h|?] [-v]
!
! INPUTS
!   None
!
! OUTPUTS
!   BUFRdsn       chr  BUFR input file name(s) (required)
!   nfiles        int  No. of BUFR input files
!   ROPPdsn       chr  ROPP output file name (optional)
!   multi         log  Multifile output (netCDF only) (default: single files)
!   newfile       log  .T. to start a new file, else append (default: new)
!   fMsgToDecode  int  First  RO message  to decode (default: 1)
!   nMsgToDecode  int  No. of RO messages to decode (default: all)
!
! CALLS
!   Usage
!   message
!   GETARG
!   IARGC
!
! CALLED BY
!   bufr2ropp
!
! DESCRIPTION
!   Provides a command line interface for the BUFR-to-ROPP
!   decoder application. See comments for main program bufr2ropp
!   for the command line details.
!
! SEE ALSO
!   bufr2ropp(1)
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN=*), PARAMETER :: Defdsn = "ropp.bufr" ! default input file name

! Argument list parameters

  CHARACTER (LEN=*), INTENT(OUT) :: BUFRdsn(:)    ! BUFR (input)  dataset name(s)
  CHARACTER (LEN=*), INTENT(OUT) :: ROPPdsn       ! ROPP (output) dataset name
  LOGICAL,           INTENT(OUT) :: multi         ! Single/Multifile output flag
  LOGICAL,           INTENT(OUT) :: newfile       ! New/Append output flag
  INTEGER,           INTENT(OUT) :: nfiles        ! No. of BUFR input files
  INTEGER,           INTENT(OUT) :: fMsgToDecode  ! First  ROPP message  to decode
  INTEGER,           INTENT(OUT) :: nMsgToDecode  ! No. of ROPP messages to decode

! Local parameters

  CHARACTER (LEN=256) :: carg    ! command line argument
  CHARACTER (LEN=256) :: Value   ! value extracted from argument
  INTEGER             :: ia, n   ! counters
  INTEGER             :: Narg    ! number of command line arguments
  INTEGER             :: ierr    ! error status code

  ! Some compilers may need the following declaration to be commented out
  INTEGER :: IARGC

!-------------------------------------------------------------
! 1. Initialise
!-------------------------------------------------------------

  BUFRdsn(:)   = Defdsn
  nfiles       = 0
  ROPPdsn      = " "
  multi        = .FALSE.
  newfile      = .TRUE.
  fMsgToDecode = 1
  nMsgToDecode = 999999

!-------------------------------------------------------------
! 2. Loop over all command line arguments.
!    If a switch has a trailing blank, then we need to get
!    the next string as it's argument.
!-------------------------------------------------------------

  ia   = 1
  narg = IARGC()

  DO WHILE ( ia <= Narg )

    CALL GETARG ( ia, carg )
    IF ( carg(1:1) == "?" .OR. &
         carg(1:6) == "--help"    ) carg = "-h"
    IF ( carg(1:9) == "--version" ) carg = "-v"

    IF ( carg(1:1) == "-" ) THEN   ! is this an option introducer?
                                   ! If so, which one?
      SELECT CASE (carg(2:2))

        CASE ("a","A")             ! Append (multifile) output requested
          newfile = .FALSE.
          multi   = .TRUE.

        CASE ("d","D")             ! debug/diagnostics wanted
          msg_MODE = VerboseMode

        CASE ("f","F" )            ! Start Message (first)
          Value = carg(3:)
          IF ( Value == " " ) THEN
            ia = ia + 1
            CALL GETARG ( ia, Value )
          END IF
          READ ( Value, *, IOSTAT=ierr ) n
          IF ( ierr == 0 .AND. &
                n    >  1 ) fMsgToDecode = n

        CASE ("h","H")             ! Help wanted
          CALL Usage()
          CALL EXIT(msg_exit_ok)

        CASE ("m","M")             ! Multifile output requested
          multi = .TRUE.

        CASE ("n","N" )            ! Number of Messages
          Value = carg(3:)
          IF ( Value == " " ) THEN
            ia = ia + 1
            CALL GETARG ( ia, Value )
          END IF
          READ ( Value, *, IOSTAT=ierr ) n
          IF ( ierr == 0 .AND. &
                n   >  0 ) nMsgToDecode = n

        CASE ("o","O" )            ! Output file
          ROPPdsn = carg(3:)
          IF ( ROPPdsn == " " ) THEN
            ia = ia + 1
            CALL GETARG ( ia, ROPPdsn )
          END IF

        CASE ("v","V" )            ! Only program version ID wanted
          CALL version_info()
          CALL EXIT(msg_exit_ok)

! Ignore anything else

        CASE DEFAULT               ! unknown option
      END SELECT

    ELSE
      nfiles = nfiles + 1
      BUFRdsn(nfiles) = carg       ! not an option - must be an input name
    END IF

    ia = ia + 1
  END DO         ! argument loop

  IF ( nfiles == 0 ) nfiles = 1    ! No input files - try a default name

END SUBROUTINE GetOptions
!----------------------------------------------------------------------------

SUBROUTINE ConvertBUFRtoROPP ( Values, & ! (in)
                               NFreq,  & ! (in)
                               ROdata, & ! (inout)
                               ierr )    ! (out)
!
!****s* bufr2ropp/ConvertBUFRtoROPP *
!
! NAME
!   ConvertBUFRtoROPP
!
! SYNOPSIS
!   Convert BUFR data to ROPP specification
!
!    USE ropp_io_types
!    USE bufr2ropp
!    TYPE (ROprof) ROdata
!    REAL*8  :: values(nv)
!    INTEGER :: nfreq, ierr
!    CALL ConvertBUFRtoROPP ( values, nfreq, ROdata, ierr )
!   where
!    ne is the number of elements (data items from BUFR)
!
! INPUTS
!   Values  dflt  Array of values from BUFR decoder
!   NFreq    int  No. of frequency sets in BUFR (1 or 3)
!   ROdata  dtyp  ROPP data - derived type
!
! OUTPUTS
!   ROdata  dtyp  ROPP data - derived type, updated
!   ierr     int  Return code: 0=OK
!                              1=sync error (mis-matched replic. factor)
!                              2=wrong frequency value
!                              3=wrong vertical significance code
!
! USES
!   ropp_io_types  - ROPP file I/O support
!   ropp_io        - i/o routines
!   geodesy        - geodesy routines
!   messages       - messages interface
!
! CALLS
!   ConvertCodes
!   ropp_io_occid
!   geometric2geopotential
!   messages
!   message_get_routine
!   message_set_routine
!
! CALLED BY
!   bufr2ropp
!
! DESCRIPTION
!   Converts decoded BUFR plain array data to ROPP derived type units etc.
!   This procedure is mostly scaling and/or range changing (e.g. Pa to hPa).
!   This routine also performs gross error checking, so that if data is not
!   valid (normally indicated by the BUFR 'missing data' flag value) that data
!   value is left as default "missing" in the ROPP structure.
!
! REFERENCES
!   1) ROPP User Guide - Part I
!      SAF/ROM/METO/UG/ROPP/002
!   2) WMO FM94 (BUFR) Specification for ROM SAF Processed Radio
!      Occultation Data. SAF/ROM/METO/FMT/BUFR/001
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

! Modules

  USE typesizes,     dp => EightByteReal
  USE ropp_io_types, ONLY: ROprof,          &
                           PCD_occultation
  USE ropp_io,       ONLY: ropp_io_occid
  USE geodesy,       ONLY: geometric2geopotential

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN=*), PARAMETER :: FmtNum = "(I20)"
  REAL(dp), PARAMETER :: MISSING = RVIND - 0.7E38_dp  ! Missing data check value (10^38)

! Argument list parameters

  REAL(dp),      INTENT(IN)    :: Values(:) ! Decoded values
  INTEGER,       INTENT(IN)    :: NFreq     ! No. of frequencies
  TYPE (ROprof), INTENT(INOUT) :: ROdata    ! RO data structure
  INTEGER,       INTENT(OUT)   :: ierr      ! Return code

! Local parameters

  CHARACTER (LEN=50) :: routine  ! saved routine name
  CHARACTER (LEN=80) :: outmsg   ! output message string
  CHARACTER (LEN=20) :: number   ! temporary string for numeric values
  CHARACTER (LEN=4)  :: Ccode    ! ICAO code associated with Ocode
  INTEGER            :: Gclass   ! GNSS satellite class value
  INTEGER            :: Gcode    ! GNSS PRN
  INTEGER            :: Lcode    ! LEO  code value
  INTEGER            :: Icode    ! Instrument code value
  INTEGER            :: Ocode    ! Origin. centre code value
  INTEGER            :: Scode    ! Sub-centre code value
  INTEGER            :: Bcode    ! B/G generator code value
! INTEGER            :: ProdType ! Product type code ! Commented at 20 July, 2016
  INTEGER            :: PCD      ! PCD bit flags (16-bit)
  INTEGER            :: in       ! loop counter for profile arrays
  INTEGER            :: IE       ! index to Values element
  INTEGER            :: RepFac   ! Replication Factor
! INTEGER            :: FOSsig   ! First-order statistics significance code ! Commented at 20 July, 2016
! INTEGER            :: TimeSig  ! Time signficance code ! Commented at 20 July, 2016
  INTEGER            :: ioerr    ! I/O error code
  REAL               :: SWver    ! Softwre version number
  REAL(dp)           :: lat      ! Nominal latitude of occultation
  REAL(dp)           :: ht       ! Sample geometric height
  REAL(dp)           :: MeanFreq ! Mean frequency (Hz)

!-------------------------------------------------
! 0. Initialise
!-------------------------------------------------

  CALL message_get_routine ( routine )
  CALL message_set_routine ( 'ConvertBUFRtoROPP' )

  ierr = 0

!-------------------------------------------------
! 1. ROPP Header block
!-------------------------------------------------

  Lcode  = NINT(Values(1))                      ! [001007] LEO ID
  Icode  = NINT(Values(2))                      ! [002019] RO Instrument
  Ocode  = NINT(Values(3))                      ! [001033] Proc.centre
  Scode  = 0
  Bcode  = Ocode
  Gclass = NINT(Values(21))                     ! [002020] GNSS class
  Gcode  = NINT(Values(22))                     ! [001050] GNSS PRN
  CALL ConvertCodes ( ROdata,        &
                      Gclass, Gcode, &
                      Lcode,  Icode, &
                      Ocode,  Scode, &
                      Ccode,  Bcode, &
                      -1 )
! ProdType = NINT(Values(4))                    ! [002172] Product type ! Commented at 20 July, 2016

  IF ( Values(5) < MISSING ) THEN               ! [025060] Software version
    SWver = Values(5) * 1E-3
    WRITE ( number,        &
            FMT="(F10.3)", &
            IOSTAT=ioerr ) SWver
    IF ( ioerr == 0 ) THEN
      IF ( SWver < 10.0 ) THEN
        ROdata%Software_Version = "V0" // ADJUSTL ( number )
      ELSE
        ROdata%Software_Version = "V"  // ADJUSTL ( number )
      END IF
    END IF
  END IF

! Date/time of start of occultation

! TimeSig             = NINT(Values(6))         ! [008021] Time.sig (17=start) ! Commented at 20 July, 2016
  ROdata%DTocc%Year   = NINT(Values(7))         ! [004001] Year
  ROdata%DTocc%Month  = NINT(Values(8))         ! [004002] Month
  ROdata%DTocc%Day    = NINT(Values(9))         ! [004003] Day
  ROdata%DTocc%Hour   = NINT(Values(10))        ! [004004] Hour
  ROdata%DTocc%Minute = NINT(Values(11))        ! [004005] Minute
  IF ( Values(12) < MISSING ) THEN              ! [004006] Seconds & MSecs
    ROdata%DTocc%Second = INT(Values(12))
    ROdata%DTocc%MSec   = NINT(MOD(Values(12),1.0_dp)*1E3)
  END IF

! Summary quality information. Only use 1st 16 bits in unswapped bit order

  IF ( Values(13) < MISSING ) THEN
    PCD = NINT(Values(13))
    ROdata%PCD = 0
    DO in = 0, 15
      IF ( BTEST(PCD, in) ) ROdata%PCD = IBSET(ROdata%PCD, 15-in)
    END DO
  END IF

  IF ( Values(14) < MISSING ) &                 ! [033007] Percent confidence
       ROdata%Overall_Qual = Values(14)

! We now have enough header info to generate the occultation ID

  CALL ropp_io_occid ( ROdata )

! Background information

  IF ( BTEST(ROdata%PCD,PCD_occultation) ) THEN
    ROdata%bg%Year   = ROdata%DTocc%Year
    ROdata%bg%Month  = ROdata%DTocc%Month
    ROdata%bg%Day    = ROdata%DTocc%Day
    ROdata%bg%Hour   = ROdata%DTocc%Hour
    ROdata%bg%Minute = ROdata%DTocc%Minute
  ELSE
    ROdata%bg%Source = "NONE"
  END IF

! Local Earth parameters

  IF ( Values(29) < MISSING ) &                 ! [004016] Time/start (s)
       ROdata%GeoRef%Time_Offset = Values(29)

  IF ( Values(30) < MISSING ) &                 ! [005001] Latitude (deg)
       ROdata%GeoRef%Lat = Values(30)

  IF ( Values(31) < MISSING ) &                 ! [006001] Longitude (deg)
       ROdata%GeoRef%Lon = Values(31)

  IF ( Values(32) < MISSING ) &                 ! [027031] CofC X (m)
       ROdata%GeoRef%r_CoC(1) = Values(32)

  IF ( Values(33) < MISSING ) &                 ! [028031] CofC Y (m)
       ROdata%GeoRef%r_CoC(2) = Values(33)

  IF ( Values(34) < MISSING ) &                 ! [010031] CofC Z (m)
       ROdata%GeoRef%r_CoC(3) = Values(34)

  IF ( Values(35) < MISSING ) &                 ! [010035] Radius value (m)
       ROdata%GeoRef%RoC = Values(35)

  IF ( Values(36) < MISSING ) &                 ! [005021] Line of sight bearing (degT)
       ROdata%GeoRef%Azimuth = Values(36)

  IF ( Values(37) < MISSING ) &                 ! [010036] Geoid undulation (m)
       ROdata%GeoRef%Undulation = Values(37)

!-------------------------------------------------
! 2. Level 1a data (SNR/Phase/POD)
!-------------------------------------------------

! LEO & GNSS POD - only one nominal value set in BUFR

  ROdata%Lev1a%Npoints  = 1
  ROdata%Lev1a%dtime(1) = 0.0                   ! [004016] Time offset from occ. start (s)

  IF ( Values(15) < MISSING ) &                 ! [027031] LEO X posn (m)
       ROdata%Lev1a%R_LEO(1,1) = Values(15)
  IF ( Values(16) < MISSING ) &                 ! [028031] LEO Y posn (m)
       ROdata%Lev1a%R_LEO(1,2) = Values(16)
  IF ( Values(17) < MISSING ) &                 ! [010031] LEO Z posn (m)
       ROdata%Lev1a%R_LEO(1,3) = Values(17)

  IF ( Values(18) < MISSING ) &                 ! [001041] LEO X vely (m/s)
       ROdata%Lev1a%V_LEO(1,1) = Values(18)
  IF ( Values(19) < MISSING ) &                 ! [001042] LEO Y vely (m/s)
       ROdata%Lev1a%V_LEO(1,2) = Values(19)
  IF ( Values(20) < MISSING ) &                 ! [001043] LEO Z vely (m/s)
       ROdata%Lev1a%V_LEO(1,3) = Values(20)

  IF ( Values(23) < MISSING ) &                 ! [027031] GNSS X posn (m)
       ROdata%Lev1a%R_GNS(1,1) = Values(23)
  IF ( Values(24) < MISSING ) &                 ! [028031] GNSS Y posn (m)
       ROdata%Lev1a%R_GNS(1,2) = Values(24)
  IF ( Values(25) < MISSING ) &                 ! [010031] GNSS Z posn (m)
       ROdata%Lev1a%R_GNS(1,3) = Values(25)

  IF ( Values(26) < MISSING ) &                 ! [001041] GNSS X vely (m/s)
       ROdata%Lev1a%V_GNS(1,1) = Values(26)
  IF ( Values(27) < MISSING ) &                 ! [001042] GNSS Y vely (m/s)
       ROdata%Lev1a%V_GNS(1,2) = Values(27)
  IF ( Values(28) < MISSING ) &                 ! [001043] GNSS Z vely (m/s)
       ROdata%Lev1a%V_GNS(1,3) = Values(28)

! If no valid POD, then skip Level 1a entirely

  IF ( ABS(ROdata%Lev1a%R_LEO(1,1)) < 0.1_dp .AND. &
       ABS(ROdata%Lev1a%V_LEO(1,1)) < 0.1_dp .AND. &
       ABS(ROdata%Lev1a%R_GNS(1,1)) < 0.1_dp .AND. &
       ABS(ROdata%Lev1a%V_GNS(1,1)) < 0.1_dp )     &
    ROdata%Lev1a%Npoints = 0

  IE = 37

!-------------------------------------------------
! 3. Level 1b data (bending angle profiles)
!-------------------------------------------------

  RepFac = NINT(Values(IE+1))                   ! [031002] Sample Replication factor
  IF ( RepFac /= ROdata%Lev1b%Npoints ) THEN
    WRITE ( outmsg, FMT="(A,2I10)" ) "Sync error: L1b RepFac /= NPoints", &
                                     repfac, ROdata%Lev1b%Npoints
    CALL message ( msg_error, TRIM(outmsg) )
    ierr = 1
  END IF

  DO in = 1, ROdata%Lev1b%Npoints

! Coordinates

    IF ( Values(IE+2) < MISSING ) &             ! [005001] Latitude (deg)
         ROdata%Lev1b%Lat_tp(in) = Values(IE+2)

    IF ( Values(IE+3) < MISSING ) &             ! [006001] Longitude (deg)
         ROdata%Lev1b%Lon_tp(in) = Values(IE+3)

    IF ( Values(IE+4) < MISSING ) &             ! [005021] Line of sight bearing (degT)
         ROdata%Lev1b%Azimuth_tp(in) = Values(IE+4)

! Extract L1+L2 if they were encoded

    RepFac = NINT(Values(IE+5))                 ! [031001] Frequency Replication factor
    IF ( RepFac /= nFreq ) THEN
      WRITE ( outmsg, FMT="(A,I12,3I6)" ) "Sync error: L1b RepFac /= nFreq", &
                                         repfac, nFreq,in,ie+5
      CALL message ( msg_error, TRIM(outmsg) )
      ierr = 1
    END IF

    IF ( NFreq == 3 ) THEN

! L1 data

      MeanFreq = Values(IE+6)                   ! [002121] Mean frequency (L1=1.5GHz)
      IF ( MeanFreq /= 1.5E9 ) THEN
        WRITE ( outmsg, FMT="(A,F5.1,A)") "Wrong L1 (1.5) frequency: ", &
                                          MeanFreq/1E9," GHz"
        CALL message ( msg_error, TRIM(outmsg) )
        ierr = 2
      END IF

      IF ( Values(IE+7) < MISSING ) &           ! [007040] Impact parameter (m)
           ROdata%Lev1b%Impact_L1(in) = Values(IE+7)

      IF ( Values(IE+8) < MISSING ) &           ! [015037] B/angle (rad)
           ROdata%Lev1b%BAngle_L1(in) = Values(IE+8)

!     FOSsig = NINT(Values(IE+9))               ! [008023] First-order stats signif. (13=RMS) ! Commented at 20 July, 2016

      IF ( Values(IE+10) < MISSING ) &          ! [015037] B/angle error (rad)
           ROdata%Lev1b%BAngle_L1_Sigma(in) = Values(IE+10)

!     FOSsig = NINT(Values(IE+11))              ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

! L2 data

      MeanFreq = Values(IE+12)                   ! [002121] Mean frequency (L2=1.2GHz)
      IF ( MeanFreq /= 1.2E9 ) THEN
        WRITE ( outmsg, FMT="(A,F5.1,A)" ) "Wrong L2 (1.2) frequency: ", &
                                           MeanFreq/1E9," GHz"
        CALL message ( msg_error, TRIM(outmsg) )
        ierr = 2
      END IF

      IF ( Values(IE+13) < MISSING ) &          ! [007040] Impact parameter (m)
           ROdata%Lev1b%Impact_L2(in) = Values(IE+13)

      IF ( Values(IE+14) < MISSING ) &          ! [015037] B/angle (rad)
           ROdata%Lev1b%BAngle_L2(in) = Values(IE+14)

!     FOSsig = NINT(Values(IE+15))              ! [008023] First-order stats signif. (13=RMS)

      IF ( Values(IE+16) < MISSING ) &          ! [015037] B/angle error (rad)
           ROdata%Lev1b%BAngle_L2_Sigma(in) = Values(IE+16)

!     FOSsig = NINT(Values(IE+17))              ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

    ELSE
      IE = IE - 12
    END IF

! Corrected bending angle (always encoded)

      MeanFreq = Values(IE+18)                   ! [002121] Mean frequency (Corrected=0)
      IF ( NINT(MeanFreq) /= 0 ) THEN
        WRITE ( outmsg, FMT="(A,F5.1,A)" ) "Wrong Corr (0.0) frequency: ", &
                                           MeanFreq/1E9," GHz"
        CALL message ( msg_error, TRIM(outmsg) )
        ierr = 2
      END IF

    IF ( Values(IE+19) < MISSING ) &            ! [007040] Impact parameter (m)
         ROdata%Lev1b%Impact(in) = Values(IE+19)

    IF ( Values(IE+20) < MISSING ) &            ! [015037] B/angle (rad)
         ROdata%Lev1b%BAngle(in) = Values(IE+20)

!   FOSsig = NINT(Values(IE+21))                ! [008023] First-order stats signif. (13=RMS) ! Commented at 20 July, 2016

    IF ( Values(IE+22) < MISSING ) &            ! [015037] B/angle error (rad)
         ROdata%Lev1b%BAngle_Sigma(in) = Values(IE+22)

!   FOSsig = NINT(Values(IE+23))                ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

    IF ( Values(IE+24) < MISSING ) &            ! [033007] Percent confidence
         ROdata%Lev1b%Bangle_Qual(in) = Values(IE+24)

    IE = IE + 23  ! L1b values per sample
  END DO
  IE = IE + 1     ! L1b RepFac

!-------------------------------------------------
! 4. Level 2a data (derived refractivity profile)
!-------------------------------------------------

  lat = 0.0
  IF ( ROdata%GEOref%lat >= -90.0_dp ) &
    lat = ROdata%GEOref%lat

  RepFac = NINT(Values(IE+1))                   ! [031002] Sample Replication factor
  IF ( RepFac /= ROdata%Lev2a%Npoints ) THEN
    WRITE ( outmsg, FMT="(A,2I10)" ) "Sync error: L2a RepFac /= NPoints", &
                                     repfac, ROdata%Lev2a%Npoints
    CALL message ( msg_error, TRIM(outmsg) )
    ierr = 1
  END IF

  DO in = 1, ROdata%Lev2a%Npoints

    IF ( Values(IE+2) < MISSING ) THEN          ! [007007] Height amsl (m)
         ROdata%Lev2a%Alt_Refrac(in) = Values(IE+2)

! Geopot ht is not in BUFR: re-create a value by applying
! latitude- & height-dependent gravity model

      ht = ROdata%Lev2a%Alt_Refrac(in)
      ROdata%Lev2a%Geop_Refrac(in) = geometric2geopotential(lat,ht) ! Geopot ht (gpm)
    END IF

    IF ( Values(IE+3) < MISSING ) &             ! [015036] Refractivity (N-units)
         ROdata%Lev2a%Refrac(in) = Values(IE+3)

!   FOSsig = NINT(Values(IE+4))                 ! [008023] First-order stats signif. (13=RMS) ! Commented at 20 July, 2016

    IF ( Values(IE+5) < MISSING ) &             ! [015036] Refractivity error (N-units)
         ROdata%Lev2a%Refrac_Sigma(in) = Values(IE+5)

!   FOSsig = NINT(Values(IE+6))                 ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

    IF ( Values(IE+7) < MISSING ) &             ! [033007] Percent confidence
         ROdata%Lev2a%Refrac_Qual(in) = Values(IE+7)

    IE = IE + 6   ! L2a values per sample
  END DO
  IE = IE + 1     ! L2a RepFac

!-------------------------------------------------
! 5. Level 2b data (retrieved P,T,q profile)
!-------------------------------------------------

  RepFac = NINT(Values(IE+1))                   ! [031002] Sample Replication factor
  IF ( RepFac /= ROdata%Lev2b%Npoints ) THEN
    WRITE ( outmsg, FMT="(A,2I10)" ) "Sync error: L2b RepFac /= NPoints", &
                                     repfac, ROdata%Lev2b%Npoints
    CALL message ( msg_error, TRIM(outmsg) )
    ierr = 1
  END IF

  DO in = 1, ROdata%Lev2b%Npoints

    IF ( Values(IE+2) < MISSING ) &             ! [007009] Geopot ht (gpm)
         ROdata%Lev2b%Geop(in) = Values(IE+2)

    IF ( Values(IE+3) < MISSING ) &             ! [010004] Pressure (Pa-->hPa)
         ROdata%Lev2b%Press(in) = Values(IE+3) * 1E-2

    IF ( Values(IE+4) < MISSING ) &             ! [012001] Temperature (K)
         ROdata%Lev2b%Temp(in) = Values(IE+4)

    IF ( Values(IE+5) < MISSING ) &             ! [013001] Spec/humidity (g/Kg)
         ROdata%Lev2b%SHum(in) = Values(IE+5) * 1E3

!   FOSsig = NINT(Values(IE+6))                 ! [008023] First-order stats signif. (13=RMS) ! Commented at 20 July, 2016

    IF ( Values(IE+7) < MISSING ) &             ! [010004] Pressure error (Pa-->hPa)
         ROdata%Lev2b%Press_Sigma(in) = Values(IE+7) * 1E-2

    IF ( Values(IE+8) < MISSING ) &             ! [012001] Temperature error (K)
         ROdata%Lev2b%Temp_Sigma(in) = Values(IE+8)

    IF ( Values(IE+9) < MISSING ) &             ! [013001] S/Hum error (g/Kg)
         ROdata%Lev2b%SHum_Sigma(in) = Values(IE+9) * 1E3

!   FOSsig = NINT(Values(IE+10))                ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

    IF ( Values(IE+11) < MISSING ) &            ! [033007] Percent confidence
         ROdata%Lev2b%Meteo_Qual(in) = Values(IE+11)

    IE = IE + 10   ! L2b values per sample
  END DO
  IE = IE + 1     ! L2b RepFac

!-------------------------------------------------
! 6. Level 2c data (retrieved surface params)
!-------------------------------------------------

  IF ( NINT(Values(IE+1)) == 0 ) THEN           ! [008003] Vertical signif. (0=surface)
    IF ( ROdata%Lev2c%Npoints == 1 ) THEN

      IF ( Values(IE+2) < MISSING ) &           ! [007009] Geoptot.Ht. (of surf) (gpm)
           ROdata%Lev2c%Geop_Sfc = Values(IE+2)

      IF ( Values(IE+3) < MISSING ) &           ! [010004] Surface pressure (Pa-->hPa)
           ROdata%Lev2c%Press_Sfc = Values(IE+3) * 1E-2

!     FOSsig = NINT(Values(IE+4))               ! [008023] First-order stats signif. (13=RMS) ! Commented at 20 July, 2016

      IF ( Values(IE+5) < MISSING ) &           ! [010004] S/press error (Pa-->hPa)
           ROdata%Lev2c%Press_Sfc_Sigma = Values(IE+5) * 1E-2

!     FOSsig = NINT(Values(IE+6))               ! [008023] First-order stats signif. (missing=off) ! Commented at 20 July, 2016

      IF ( Values(IE+7) < MISSING ) &           ! [033007] Percent confidence
           ROdata%Lev2c%Press_Sfc_Qual = Values(IE+7)

    END IF

  ELSE IF ( Values(IE+1) < MISSING ) THEN
    WRITE ( number, FMT=FmtNum ) NINT(Values(IE+1))
    CALL message ( msg_error, "Wrong vertical signficance, not surface (=0): "// &
                              TRIM(ADJUSTL(number)) )
    ierr = 3
  END IF

  CALL message_set_routine ( routine )

END SUBROUTINE ConvertBUFRtoROPP
!------------------------------------------------------------------------

END MODULE bufr2ropp
