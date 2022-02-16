! $Id: ropp2bufr_mod.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE ropp2bufr

!****m* ropp2bufr/ropp2bufr *
!
! NAME
!   ropp2bufr    (ropp2bufr_mod.f90)
!
! SYNOPSIS
!   Module defining fixed values & subroutines/functions for the
!   ropp2bufr main program (ECMWF or MetDB versions)
!
!   USE ropp2bufr
!
! USED BY
!   ropp2bufr   (ropp2bufr_ec or ropp2bufr_mo)
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

  USE typesizes, dp => EightByteReal
  USE messages

! Public fixed values for array sizes, etc

! GPSRO-specific GTS/RMDCN routing header elements

  CHARACTER (LEN=6), PARAMETER :: TTAAII  = "IUT?14"  ! Binary/UAir/Satellite/GPSRO

  INTEGER, PARAMETER :: NOhdrs = 0            ! No GTS routings headers
  INTEGER, PARAMETER :: ARhdrs = 1            ! GTS Abbreviated Routing headers only
  INTEGER, PARAMETER :: IPhdrs = 2            ! ARH plus IP support

! GPSRO-specific BUFR details

  INTEGER, PARAMETER :: RODescr     = 310026  ! Table D master descriptor for RO

  INTEGER, PARAMETER :: Edition     = 4       ! BUFR Edition (4)
  INTEGER, PARAMETER :: MasterTable = 0       ! BUFR Master Table (Meteorology)
  INTEGER, PARAMETER :: DataType    = 3       ! Table A (Sounding - satellite)
  INTEGER, PARAMETER :: IntlSubType = 50      ! International data sub-type (GPSRO)
  INTEGER, PARAMETER :: LoclSubType = 14      ! Local data sub-type
  INTEGER, PARAMETER :: VerMasTable = 12      ! Table version number (12)
  INTEGER, PARAMETER :: VerLocTable = 0       ! Local Table version  (not used)
  INTEGER, PARAMETER :: Sec3Type_mo = 1       ! Observed data, uncompressed (MetDB)
  INTEGER, PARAMETER :: Sec3Type_ec = 128     ! Observed data, uncompressed (ECMWF)

! Default (missing) data values

  INTEGER,  PARAMETER :: NMDFV = -9999999     ! Integer missing data flag value (MetDB)
  INTEGER,  PARAMETER :: NVIND = 2147483647   ! Integer missing data flag value (ECMWF)
  REAL,     PARAMETER :: RMDFV = -9999999.0   ! Real    missing data flag value (MetDB)
  REAL(dp), PARAMETER :: RVIND = 1.7E38_dp    ! Real    missing data flag value (ECMWF)

! MetDB I/O interface

  INTEGER, PARAMETER  :: Output = 2           ! I/O output mode (r+w, new)

! ECMWF SBYTE() packing & I/O interfaces

  INTEGER, PARAMETER  :: nbpc = 8                 ! bits per character (or byte)
  INTEGER, PARAMETER  :: nbpw = KIND(nbpc) * nbpc ! bits per word (default INTEGER)
  INTEGER, PARAMETER  :: nbsk = 0                 ! no. bits to skip when packing

CONTAINS
!--------------------------------------------------------------------

SUBROUTINE Usage()

!****s* ropp2bufr/Usage *
!
! NAME
!   Usage
!
! SYNOPIS
!   USE ropp2bufr
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
  PRINT *, '  Encode one or more ROPP-format netCDF files to WMO BUFR.'
  PRINT *, 'Usage:'
  PRINT *, '   > ropp2bufr ropp_file [ropp_file...] [-o bufr_file]'
  PRINT *, '                         [-g[i]] [-s csn_file]'
  PRINT *, '                         [-p thin_file|maxsamp] [-t time]'
  PRINT *, '                         [-u] [-l] [-d] [-m] [-h] [-v]'
  PRINT *, 'Input:'
  PRINT *, '  One or more files in ROPP netCDF format.'
  PRINT *, 'Output:'
  PRINT *, '  BUFR file, one message or bulletin per input RO profile'
  PRINT *, 'Options:'
  PRINT *, '  -o  BUFR output file name'
  PRINT *, '  -g  GTS routing headers/trailers required'
  PRINT *, '  -gi GTS headers preceded by 10-byte leading'
  PRINT *, '      size/type for GTS IP (FTP) transmission'
  PRINT *, '  -s  file containing last used channel sequence number'
  PRINT *, '      (updated on completion)'
  PRINT *, '  -p  thinning control file name or max. no. samples'
  PRINT *, '  -t  don''t encode data older than ''time'' ago (hh:mm)'
  PRINT *, '  -u  leave profiles unordered (i.e. in original height order)'
  PRINT *, '  -l  L1+L2 data (Level 1b) are not to be encoded,'
  PRINT *, '      only the ionospheric-corrected profile.'
  PRINT *, '  -m  met. data (Level 2c/d) are not to be encoded'
  PRINT *, '  -d  outputs additonal diagnostics to stdout'
  PRINT *, '  -h  this help'
  PRINT *, '  -v  version information'
  PRINT *, 'Defaults:'
  PRINT *, '  Input  file name        : none - at least one required'
  PRINT *, '  Output file name        : from (first) occultation ID <occid>.bufr'
  PRINT *, '  GTS routing headers     : not generated'
  PRINT *, '  Channel sequence nos.   : starts at 001'
  PRINT *, '  Reject time difference  : 00:00 (no rejection on time)'
  PRINT *, '   unless -g* option, when: 23:50 (assuming NRT on GTS)'
  PRINT *, '  Thinning                : sample to <= 375 levels'
  PRINT *, '  Re-ordering             : descending profiles re-ordered to ascending '
  PRINT *, 'See ropp2bufr(1) for details.'
  PRINT *, ''
END SUBROUTINE Usage
!--------------------------------------------------------------------

SUBROUTINE Usage_eum()

!****s* ropp2bufr/Usage_eum *
!
! NAME
!   Usage_eum
!
! SYNOPIS
!   USE ropp2bufr
!   CALL Usage_eum()
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
!   Prints a summary of eum2bufr program usage (help) to stdout.
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
  PRINT *, '  Encode one or more EUMETSAT netCDF files to WMO BUFR.'
  PRINT *, 'Usage:'
  PRINT *, '   > eum2bufr eum_file [eum_file...] [-o bufr_file]'
  PRINT *, '                       [-g[i]] [-s csn_file]'
  PRINT *, '                       [-p thin_file|maxsamp] [-t time] [-r resn]'
  PRINT *, '                       [-u] [-l] [-d] [-m] [-h] [-v]'
  PRINT *, 'Input:'
  PRINT *, '  One or more files in EUMETSAT netCDF-4 format.'
  PRINT *, 'Output:'
  PRINT *, '  BUFR file, one message or bulletin per input RO profile'
  PRINT *, 'Options:'
  PRINT *, '  -o  BUFR output file name'
  PRINT *, '  -g  GTS routing headers/trailers required'
  PRINT *, '  -gi GTS headers preceded by 10-byte leading'
  PRINT *, '      size/type for GTS IP (FTP) transmission'
  PRINT *, '  -s  file containing last used channel sequence number'
  PRINT *, '      (updated on completion)'
  PRINT *, '  -p  thinning control file name or max. no. samples'
  PRINT *, '  -t  don''t encode data older than ''time'' ago (hh:mm)'
  PRINT *, '  -u  leave profiles unordered (i.e. in original height order)'
  PRINT *, '  -l  L1+L2 data (Level 1b) are not to be encoded,'
  PRINT *, '      only the ionospheric-corrected profile.'
  PRINT *, '  -m  met. data (Level 2c/d) are not to be encoded'
  PRINT *, '  -r  resolution group of netCDF-4 EUM file'
  PRINT *, '  -d  outputs additonal diagnostics to stdout'
  PRINT *, '  -h  this help'
  PRINT *, '  -v  version information'
  PRINT *, 'Defaults:'
  PRINT *, '  Input  file name        : none - at least one required'
  PRINT *, '  Output file name        : from (first) occultation ID <occid>.bufr'
  PRINT *, '  GTS routing headers     : not generated'
  PRINT *, '  Channel sequence nos.   : starts at 001'
  PRINT *, '  Reject time difference  : 00:00 (no rejection on time)'
  PRINT *, '   unless -g* option, when: 23:50 (assuming NRT on GTS)'
  PRINT *, '  Thinning                : sample to <= 375 levels'
  PRINT *, '  Re-ordering             : descending profiles re-ordered to ascending '
  PRINT *, '  Resolution              : ''thinned'' '
  PRINT *, 'See eum2bufr(1) for details.'
  PRINT *, ''
END SUBROUTINE Usage_eum

!-------------------------------------------------------------------------

SUBROUTINE GetOptions ( centre,      & ! (in)
                        NCDFdsn,     & ! (out)
                        nfiles,      & ! (out)
                        BUFRdsn,     & ! (out)
                        CSNdsn,      & ! (out)
                        Thindsn,     & ! (out)
                        GTShdrType,  & ! (out)
                        RejTimeDiff, & ! (out)
                        CorrOnly,    & ! (out)
                        nomet,       & ! (out)
                        unordered,   & ! (out)
                        resolution   ) ! (out)

!****s* ropp2bufr/GetOptions *
!
! NAME
!   GetOptions
!
! SYNOPSIS
!   Get command line information & options or set defaults
!
!    USE ropp2bufr
!    CHARACTER (LEN=4) :: centre
!    CHARACTER (LEN=100) :: NCDFdsn(100), bufrdsn, csndsn
!    INTEGER :: nfiles, gtshdrtype, rejtimediff
!    LOGICAL :: corronly, nomet, unordered
!    CALL getOptions ( centre, NCDFdsn, nfiles, bufrdsn, csndsn, thindsn, &
!                      gtshdrtype, rejtimediff, &
!                      corronly, nomet, unordered, resolution )
!   On command line:
!   > ropp2bufr ropp_file [ropp_file...] [-o bufr_file]
!                         [-g[n]] [-s seq_file]
!                         [-p thin_file] [-t time]
!                         [-u] [-l] [-m] [-h|?] [-v] [-d]
!   > eum2bufr eum_file  [eum_file...] [-o bufr_file]
!                         [-g[n]] [-s seq_file]
!                         [-p thin_file] [-t time] [-r resol]
!                         [-u] [-l] [-m] [-h|?] [-v] [-d]
!
! INPUTS
!   centre       chr 'EUM' when called from eum2bufr to encode EUMETSAT netCDF
!                     files; anyhting else to encode ROPP netCDF files.
!
! OUTPUTS
!   NCDFdsn      chr  ROPP/EUM netCDF input file name(s)
!   nfiles       int  No. of ROPP input files
!   BUFRdsn      chr  BUFR output file name
!   CSNdsn       chr  Channel sequence number file name
!   Thindsn      chr  Thinning control file name
!   GTShdrType   int  GTS header type code
!   RejTimeDiff  int  Rejection time threshold (minutes)
!   CorrOnly     log  L1+L2 skip flag
!   nomet        log  Met data skip flag
!   unordered    log  Disable profile ordering flag
!   resolution   chr  resolution group to use in EUM netCDF-4 files
!
! CALLS
!   Usage
!   message
!   IARGC
!   GETARG
!
! CALLED BY
!   ropp2bufr
!   eum2bufr
!
! MODULES
!   DateTimeTypes   - Date & Time conversion definitions
!   GTShdrs         - GTS bulletin routing header support definitions
!
! DESCRIPTION
!   Provides a command line interface for the ROPP-to-BUFR and EUM-to-BUFR
!   encoder applications. See comments for main program ropp2bufr & eum2bufr
!   for the command line details.
!
! SEE ALSO
!   ropp2bufr(1), eum2bufr(1)
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

  USE DateTimeTypes, ONLY: nMinPerHour
  USE GTShdrs,       ONLY: GTSHDR_DEBUG

  IMPLICIT NONE

! Fixed parameters

  INTEGER,           PARAMETER :: DefRejTimeDiff = 1430   ! 23h50m in minutes

! Argument list parameters

  CHARACTER (LEN=*),   INTENT(IN)  :: centre      ! Centre (e.g. 'EUM')
  CHARACTER (LEN=*),   INTENT(OUT) :: NCDFdsn(:)  ! input  netCDF file name(s)
  CHARACTER (LEN=*),   INTENT(OUT) :: BUFRdsn     ! output BUFR file name
  CHARACTER (LEN=*),   INTENT(OUT) :: CSNdsn      ! Channel sequence number file name
  CHARACTER (LEN=*),   INTENT(OUT) :: Thindsn     ! thinning control file name
  INTEGER,             INTENT(OUT) :: nfiles      ! No. of ROPP input files
  INTEGER,             INTENT(OUT) :: GTShdrType  ! code for GTS header generation
  INTEGER,             INTENT(OUT) :: RejTimeDiff ! reject obs older than this
  LOGICAL,             INTENT(OUT) :: CorrOnly    ! .F. for L1+L2+C, .T. for C only
  LOGICAL,             INTENT(OUT) :: nomet       ! .F. for met data, .T. to skip
  LOGICAL,             INTENT(OUT) :: unordered   ! .T. to disable re-ordering of profiles
  CHARACTER(LEN=20),   INTENT(OUT) :: resolution  ! resolution group of EUM netCDF 4 files

! Local variables

  CHARACTER (LEN=256) :: carg                 ! command line argument
  INTEGER :: narg                             ! number of command line arguments
  INTEGER :: ia                               ! loop counter
  INTEGER :: ierr                             ! error status
  INTEGER :: hh, mm                           ! hours & minutes

! Some compilers may need the following declaration to be commented out
  INTEGER :: IARGC

!-------------------------------------------------------------
! 1. Initialise
!-------------------------------------------------------------

  NCDFdsn(:)  = " "
  nfiles      = 0
  BUFRdsn     = " "
  CSNdsn      = " "
  Thindsn     = "375"   ! to be interpreted as 'sample to no more than'
  GTShdrType  = NOhdrs
  RejTimeDiff = 0
  CorrOnly    = .FALSE.
  nomet       = .FALSE.
  unordered   = .FALSE.
  resolution  = 'thinned'

!-------------------------------------------------------------
! 2. Loop over all command line options.
!    If a switch has a trailing blank, then we need to get
!    the next string as it's argument.
!-------------------------------------------------------------

  ia   = 1
  narg = IARGC()

  DO WHILE ( ia <= narg )

    CALL GETARG ( ia, carg )
    IF ( carg(1:1) == "?" .OR. &
         carg(1:6) == "--help"    ) carg = "-h"
    IF ( carg(1:9) == "--version" ) carg = "-v"

    IF ( carg(1:1) == "-" ) THEN   ! is this an option introducer?
                                   ! If so, which one?
      SELECT CASE (carg(2:2))

        CASE ("d","D")             ! debug/diagnostics wanted
          msg_MODE = VerboseMode
          GTSHDR_DEBUG = .TRUE.

        CASE ("g","G")             ! GTS headers wanted - any extra IP?
          SELECT CASE (carg(3:3))
            CASE ("i","I")
              GTShdrType = IPhdrs   ! headers + IP
            CASE DEFAULT
              GTShdrType = ARhdrs   ! headers only
          END SELECT

        CASE ("h","H")             ! Help wanted
          narg = -1

        CASE ("l","L")             ! no L1/L2 (Corrected only)
          CorrOnly = .TRUE.

        CASE ("m","M")             ! no Met. (geophysical) data
          nomet = .TRUE.

        CASE ("o","O")             ! Output file name
          carg(1:2) = "  "
          IF ( carg(3:) == " " ) THEN
             ia = ia + 1
             CALL GETARG ( ia, carg )
          END IF
          BUFRdsn = ADJUSTL(carg)

        CASE ("p","P")             ! thinning control file name
          carg(1:2) = "  "
          IF ( carg(3:) == " " ) THEN
             ia = ia + 1
             CALL GETARG ( ia, carg )
          END IF
          Thindsn = ADJUSTL(carg)

        CASE ("r","R")             ! resolution group to use (EUM only)
          carg(1:2) = "  "
          IF ( carg(3:) == " " ) THEN
             ia = ia + 1
             CALL GETARG ( ia, carg )
          END IF
          IF ( centre == "EUM" ) resolution = ADJUSTL(carg)

        CASE ("s","S")             ! Channel sequence No. file name
          carg(1:2) = "  "
          IF ( carg(3:) == " " ) THEN
             ia = ia + 1
             CALL GETARG ( ia, carg )
          END IF
          CSNdsn = ADJUSTL(carg)

        CASE ("t","T")             ! Reject time difference (hh:mm)
          carg(1:2) = "  "
          IF ( carg(3:) == " " ) THEN
             ia = ia + 1
             CALL GETARG ( ia, carg )
          END IF
          carg = ADJUSTL(carg)
          READ ( carg, "(BN,I2,1X,I2)", IOSTAT=ierr ) hh, mm
          IF ( ierr == 0 ) RejTimeDiff = hh * nMinPerHour + mm

        CASE ("u","U")             ! Profile ordering
          unordered = .TRUE.

        CASE ("v","V")             ! Only program version ID wanted
          CALL version_info()
          CALL EXIT(msg_exit_ok)

        CASE DEFAULT               ! unknown option
      END SELECT

    ELSE                           ! not an option - must be an input name
      nfiles = nfiles + 1
      NCDFdsn(nfiles) = carg
    END IF

    ia = ia + 1
  END DO                               ! argument loop

  IF ( nfiles == 0 .AND. narg /= -1 ) THEN
    CALL message ( msg_error, "No input file(s) specified" )
    narg = 0
  END IF

  IF ( narg <= 0 ) THEN
    IF ( centre == "EUM" ) THEN
      CALL Usage_eum()
    ELSE
      CALL Usage()
    END IF
    CALL EXIT(msg_exit_status)
  END IF

!-------------------------------------------------------------
! 3. Set default time rejection if GTS routing headers to be
!    generated, on the assumption that the output is for NRT
!    GTS distribution.
!-------------------------------------------------------------

  IF ( GTShdrType /= NOhdrs .AND. &
       RejTimeDiff == 0 ) RejTimeDiff = DefRejTimeDiff

END SUBROUTINE GetOptions
!----------------------------------------------------------------------------

SUBROUTINE ConvertROPPtoBUFR ( ROdata,     & ! (in)
                               CorrOnly,   & ! (in)
                               OrigICAO,   & ! (out)
                               OrigCentre, & ! (out)
                               SubCentre,  & ! (out)
                               Values,     & ! (out)
                               nValues,    & ! (out)
                               RepFac,     & ! (out)
                               nRepFac )     ! (out)
!
!****s* ropp2bufr/ConvertROPPtoBUFR *
!
! NAME
!   ConvertROPPtoBUFR
!
! SYNOPSIS
!   Convert ROPP data to BUFR specification
!
!    USE ropp_io_types
!    USE ropp2bufr
!    TYPE (ROprof) rodata
!    CHARACTER (LEN=4) :: origicao
!    INTEGER :: origcentre, subcentre, nvalues, nrepfac, repfac(nr)
!    REAL(dp):: values(ne)
!    LOGICAL :: corronly
!    CALL convertropptobufr ( rodata, corronly, &
!                             origicao, origcentre, subcentre, &
!                             values, nvalues, repfac, nrepfac )
!   where
!    ne is the max. number of elements (data items for BUFR)
!    nr is the max. number of delayed replication factors
!
! INPUTS
!   ROdata     dtyp  RO data - derived type
!   CorrOnly    log  Flag for corrected Level 1b profile only
!
! OUTPUTS
!   ROdata     dtyp  RO data - derived type (potentially modified)
!   OrigICAO    chr  4-chr ICAO code associated with Orig.Centre
!   OrigCentre  int  Originating Centre code value
!   SubCentre   int  Originating subcentre (processing centre) code value
!   Values     dflt  Array(ne) of converted values for BUFR encoder
!   nValues     int  Total no. of values converted
!   RepFac      int  Array of Replication Factors
!   nRepFac     int  Total no. of Replication Factors
!
! MODULES
!   ropp_io_types  - ROPP file I/O support
!   ropp_utils     - ROPP utility functions & parameters
!
! CALLS
!   ConvertCodes
!   message
!   message_get_routine
!   message_set_routine
!
! CALLED BY
!   ropp2bufr
!
! DESCRIPTION
!   Converts RO data to BUFR units, etc, and returns converted data as a plain
!   1-D array. This procedure is mostly scaling and/or range changing (e.g
!   longitude from 0-360 to +/-180deg, hPa to Pa).
!   This routine also performs gross error checking, so that if data is not
!   valid (not within nominal range of BUFR bit width) that data value is set
!   "missing" in the output array.
!   The delayed replication factor counts are also returned for use with
!   the ECMWF BUFEN() encoder (not used with the MetDB ENBUFV4() encoder).
!   The processing (originating) centre's code & sub-centre code are
!   returned for insertion in BUFR Section 1, plus the ICAO Location
!   Indicator code for optonal use in a GTS routing header.
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

  USE ropp_io_types, ONLY: ROprof, &
                           PCD_occultation

  IMPLICIT NONE

! Fixed parameters

  REAL(dp), PARAMETER :: MISSING = RVIND     ! Missing data flag value

  INTEGER, PARAMETER :: ProdType = 2         ! Product type (limb sounding)
  INTEGER, PARAMETER :: TimeSig  = 17        ! Time significance (start)
  INTEGER, PARAMETER :: FOstats  = 13        ! First-order statistics (rms)

  REAL,    PARAMETER :: FreqL1   = 1.5E9     ! L1 frequency: 1.5GHz
  REAL,    PARAMETER :: FreqL2   = 1.2E9     ! L1 frequency: 1.2GHz
  REAL,    PARAMETER :: FreqLc   = 0.0       ! Corrected frequency (dummy)

  CHARACTER (LEN=*), PARAMETER :: numeric = "0123456789." ! valid for numerics

! Argument list parameters

  TYPE (ROprof),     INTENT(INOUT) :: ROdata
  LOGICAL,           INTENT(IN)    :: CorrOnly
  CHARACTER (LEN=4), INTENT(OUT)   :: OrigICAO
  INTEGER,           INTENT(OUT)   :: OrigCentre
  INTEGER,           INTENT(OUT)   :: SubCentre
  REAL(dp),          INTENT(OUT)   :: Values(:)
  INTEGER,           INTENT(OUT)   :: nValues
  INTEGER,           INTENT(OUT)   :: RepFac(:)
  INTEGER,           INTENT(OUT)   :: nRepFac

! Local parameters

  CHARACTER (LEN=10)  :: number  ! temporary strings for numeric values
  CHARACTER (LEN=256) :: routine ! temporary previously set message routine
  CHARACTER (LEN=4)   :: Ccode   ! ICAO code associated with Ocode
  INTEGER :: Gclass     ! GNSS class value
  INTEGER :: Gcode      ! GNSS PRN
  INTEGER :: Lcode      ! LEO  code value
  INTEGER :: Icode      ! Instrument code value
  INTEGER :: Ocode      ! Origin. centre code value
  INTEGER :: Scode      ! Sub-centre code value
  INTEGER :: Bcode      ! B/G generator code value
  INTEGER :: PCD        ! PCD bit flags (16-bit)
  INTEGER :: in         ! loop counter for profile arrays
  INTEGER :: IE         ! index offset to Values element
  INTEGER :: ierr       ! I/O error code
  REAL    :: SWver      ! Software version number

!-------------------------------------------------------------
! 1. Initialise
!------------------------------------------------------------

  CALL message_get_routine ( routine )
  CALL message_set_routine ( "ConvertROPPtoBUFR" )

  Values(:) = MISSING
  nValues   = 0

  RepFac(:) = 0
  nRepFac   = 0

!-------------------------------------------------------------
! 2. Convert ROPP character codes to BUFR numeric codes.
!-------------------------------------------------------------
! Possible Originating Centre & associated ICAO codes for GPSRO include:
!  007/KWBC - Washington (US) [UCAR/CDAAC]
!  074/EGRR - Exeter     (GB) [Met Office]
!  078/EDZW - Offenbach  (DE) [GFZ]
!  094/EKMI - Copenhagen (DK) [DMI/ROM SAF]
!  160/KNES - Washington (US) [NESDIS]
!  254/EUMS - Darmstadt  (DE) [EUMETSAT]
! See WMO BUFR code table 001033 (Common Code Table C-1 or C-11)
! for the full list of Originating Centres.
!-------------------------------------------------------------

  CALL ConvertCodes ( ROdata,        &
                      Gclass, Gcode, &
                      Lcode,  Icode, &
                      Ocode,  Scode, &
                      Ccode,  Bcode, &
                      1 )
  OrigICAO   = Ccode
  OrigCentre = Ocode
  SubCentre  = Scode

!-------------------------------------------------------------
! 3. Satellite data introducer
!-------------------------------------------------------------

  Values(1) = Lcode                                    ! [001007] LEO ID
  IF ( Values(1) <    0.0 .OR. &
       Values(1) > 1022.0 )    &
       Values(1) = MISSING

  Values(2) = Icode                                    ! [002019] RO Instrument
  IF ( Values(2) <    0.0 .OR. &
       Values(2) > 2046.0 )    &
       Values(2) = MISSING

  IF ( BTEST(ROdata%PCD,PCD_occultation) ) THEN
    Values(3) = Bcode                                  ! [001033] B/g gen.centre
  ELSE
    Values(3) = Ocode                                  ! [001033] Proc.centre
  END IF
  IF ( Values(3) <   0.0 .OR. &
       Values(3) > 254.0 )    &
       Values(3) = MISSING

  Values(4) = ProdType                                 ! [002172] Product type
                                                       !         (limb sounding)
  number = ROdata%Software_Version(1:10)
  DO in = 1, LEN_TRIM(number)
    IF ( INDEX ( numeric, number(in:in) ) == 0 ) number(in:in) = " "
  END DO
  READ ( number, FMT=*, IOSTAT=ierr) SWver
  IF ( ierr /= 0 ) SWver = -9.999
  Values(5) = SWver * 1E3                              ! [025060] Software version
  IF ( Values(5) <     0.0 .OR. &
       Values(5) > 16382.0 )    &
       Values(5) = MISSING

! Date/time of start of occultation (or background profile)

  Values(6)  = TimeSig                                 ! [008021] Time.sig (start)
  Values(7)  = ROdata%DTocc%Year                       ! [004001] Year
  Values(8)  = ROdata%DTocc%Month                      ! [004002] Month
  Values(9)  = ROdata%DTocc%Day                        ! [004003] Day
  Values(10) = ROdata%DTocc%Hour                       ! [004004] Hour
  Values(11) = ROdata%DTocc%Minute                     ! [004005] Minute
  Values(12) = ROdata%DTocc%Second &                   ! [004006] Seconds & MSecs
             + ROdata%DTocc%MSec * 1E-3
  IF ( Values(12) <  0.000 .OR. &
       Values(12) > 59.999 )    &
       Values(12) = MISSING

! Summary quality information

  PCD = 0
  DO in = 0, 15
    IF ( BTEST(ROdata%PCD, in) ) PCD = IBSET(PCD, 15-in) ! only use 1st 16 bits in swapped bit order
  END DO
  Values(13) = REAL(PCD)                               ! [033039] PCD
  IF ( Values(13) <     0.0 .OR. &
       Values(13) > 65534.0 )    &
       Values(13) = MISSING

  Values(14) = REAL(ROdata%Overall_Qual)               ! [033007] Percent confidence
  IF ( Values(14) <   0.0  .OR. &
       Values(14) > 100.0 )     &
       Values(14) = MISSING

! LEO & GNSS POD

  IF ( .NOT. ROdata%Lev1a%Missing ) THEN
    Values(15) = REAL(ROdata%Lev1a%R_LEO(1,1))         ! [027031] LEO X posn (m)
    IF ( ABS(Values(15)) > 10737418.23_dp ) &
             Values(15) = MISSING
    Values(16) = REAL(ROdata%Lev1a%R_LEO(1,2))         ! [028031] LEO Y posn (m)
    IF ( ABS(Values(16)) > 10737418.23_dp ) &
             Values(16) = MISSING
    Values(17) = REAL(ROdata%Lev1a%R_LEO(1,3))         ! [010031] LEO Z posn (m)
    IF ( ABS(Values(17)) > 10737418.23_dp ) &
             Values(17) = MISSING
    IF ( ABS(Values(15)) < 1.0 .AND. &
         ABS(Values(16)) < 1.0 .AND. &
         ABS(Values(17)) < 1.0 )     &
             Values(15:17) = MISSING
    Values(18) = REAL(ROdata%Lev1a%V_LEO(1,1))         ! [001041] LEO X vely (m/s)
    IF ( ABS(Values(18)) > 10737.41823_dp ) &
             Values(18) = MISSING
    Values(19) = REAL(ROdata%Lev1a%V_LEO(1,2))         ! [001042] LEO Y vely (m/s)
    IF ( ABS(Values(19)) > 10737.41823_dp ) &
             Values(19) = MISSING
    Values(20) = REAL(ROdata%Lev1a%V_LEO(1,3))         ! [001043] LEO Z vely (m/s)
    IF ( ABS(Values(20)) > 10737.41823_dp ) &
             Values(20) = MISSING
    IF ( ABS(Values(18)) < 1.0 .AND. &
         ABS(Values(19)) < 1.0 .AND. &
         ABS(Values(20)) < 1.0 )     &
             Values(18:20) = MISSING
  END IF

  Values(21) = Gclass                                  ! [002020] GNSS class
  IF ( Values(21) < 0 .OR. &
       Values(21) > 510 ) Values(21) = MISSING
  Values(22) = Gcode                                   ! [001050] GNSS PRN
  IF ( Values(22) < 0 .OR. &
       Values(22) > 131070 ) Values(22) = MISSING

  IF ( .NOT. ROdata%Lev1a%Missing ) THEN
    Values(23) = REAL(ROdata%Lev1a%R_GNS(1,1))         ! [027031] GNSS X posn (m)
    IF ( ABS(Values(23)) > 107374182.4_dp ) &
             Values(23) = MISSING
    Values(24) = REAL(ROdata%Lev1a%R_GNS(1,2))         ! [028031] GNSS Y posn (m)
    IF ( ABS(Values(24)) > 107374182.4_dp ) &
             Values(24) = MISSING
    Values(25) = REAL(ROdata%Lev1a%R_GNS(1,3))         ! [010031] GNSS Z posn (m)
    IF ( ABS(Values(25)) > 107374182.4_dp ) &
             Values(25) = MISSING
    IF ( ABS(Values(23)) < 1.0 .AND. &
         ABS(Values(24)) < 1.0 .AND. &
         ABS(Values(25)) < 1.0 )     &
             Values(23:25) = MISSING

    Values(26) = REAL(ROdata%Lev1a%V_GNS(1,1))         ! [001041] GNSS X vely (m/s)
    IF ( ABS(Values(26)) > 10737.41824_dp ) &
             Values(26) = MISSING
    Values(27) = REAL(ROdata%Lev1a%V_GNS(1,2))         ! [001042] GNSS Y vely (m/s)
    IF ( ABS(Values(27)) > 10737.41824_dp ) &
             Values(27) = MISSING
    Values(28) = REAL(ROdata%Lev1a%V_GNS(1,3))         ! [001043] GNSS Z vely (m/s)
    IF ( ABS(Values(28)) > 10737.41824_dp ) &
             Values(28) = MISSING
    IF ( ABS(Values(26)) < 1.0 .AND. &
         ABS(Values(27)) < 1.0 .AND. &
         ABS(Values(28)) < 1.0 )     &
             Values(26:28) = MISSING
  END IF

! Local Earth parameters

  Values(29) = REAL(ROdata%GeoRef%Time_Offset)         ! [004016] Time/start (s)
  IF ( Values(29) <    0.0 .OR. &
       Values(29) >  240.0 )    &
       Values(29) = MISSING

  Values(30) = REAL(ROdata%GeoRef%Lat)                 ! [005001] Latitude (deg)
  IF ( ABS(Values(30)) > 90.0 ) &
       Values(30) = MISSING

  Values(31) = REAL(ROdata%GeoRef%Lon)                 ! [006001] Longitude (deg)
  IF ( Values(31) > 180.0 ) &
       Values(31) = Values(31) - 360.0
  IF ( ABS(Values(31)) > 180.0 ) &
       Values(31) = MISSING

  Values(32) = REAL(ROdata%GeoRef%r_CoC(1))            ! [027031] CofC X (m)
  IF ( ABS(Values(32)) > 1000000.0_dp ) &
       Values(32) = MISSING

  Values(33) = REAL(ROdata%GeoRef%r_CoC(2))            ! [028031] CofC Y (m)
  IF ( ABS(Values(33)) > 1000000.0_dp ) &
       Values(33) = MISSING

  Values(34) = REAL(ROdata%GeoRef%r_CoC(3))            ! [010031] CofC Z (m)
  IF ( ABS(Values(34)) > 1000000.0_dp ) &
       Values(34) = MISSING

  Values(35) = REAL(ROdata%GeoRef%RoC)                 ! [010035] Radius value (m)
  IF ( Values(35) < 6200000.0_dp .OR. &
       Values(35) > 6600000.0_dp )    &
       Values(35) = MISSING

  Values(36) = REAL(ROdata%GeoRef%Azimuth)             ! [005021] Line of sight bearing (degT)
  IF ( Values(36) <    0.0 .OR. &
       Values(36) >= 360.0 )    &
       Values(36) = MISSING

  Values(37) = REAL(ROdata%GeoRef%Undulation)          ! [010036] Geoid undulation (m)
  IF ( ABS(Values(37)) > 163.82 ) &
       Values(37) = MISSING

  IE = 37

!-------------------------------------------------------------
! 4. Level 1b data (bending angle profile)
!-------------------------------------------------------------

! Interpolation thinning may generate a set of fixed impact height
! levels but no valid Level 1b profile data if no input L1b - reject such
! empty BA-profiles

  IF ( ROdata%Lev1b%Missing ) THEN
    CALL message(msg_diag, "Rejecting empty Level 1b (BA) profile")
    ROdata%Lev1b%Npoints = 0
  END IF

  Values(IE+1)    = ROdata%Lev1b%Npoints               ! [031002] Replication factor
  nRepFac         = nRepFac + 1
  RepFac(nRepFac) = NINT(Values(IE+1))

  DO in = 1, ROdata%Lev1b%Npoints

! Coordinates

    Values(IE+2) = REAL(ROdata%Lev1b%Lat_tp(in))       ! [005001] Latitude (deg)
    IF ( ABS(Values(IE+2)) > 90.0 ) &
         Values(IE+2) = MISSING

    Values(IE+3) = REAL(ROdata%Lev1b%Lon_tp(in))       ! [006001] Longitude (deg)
    IF ( Values(IE+3) > 180.0 ) &
         Values(IE+3) = Values(IE+3) - 360.0
    IF ( ABS(Values(IE+3)) > 180.0 ) &
         Values(IE+3) = MISSING

    Values(IE+4) = REAL(ROdata%Lev1b%Azimuth_tp(in))   ! [005021] Line of sight bearing (degT)
    IF ( Values(IE+4) <    0.0 .OR. &
         Values(IE+4) >= 360.0 )    &
         Values(IE+4) = MISSING

! Include L1+L2 or skip them?

    nRepFac = nRepFac + 1
    IF ( CorrOnly ) THEN
      Values(IE+5)    = 1                              ! [031001] Replication factor
      IE = IE - 12
      RepFac(nRepFac) = 1
    ELSE
      Values(IE+5)    = 3                              ! [031001] Replication factor
      RepFac(nRepFac) = 3

! L1 data

      Values(IE+6) = FreqL1                            ! [002121] L1=1.5Ghz

      Values(IE+7) = REAL(ROdata%Lev1b%Impact_L1(in))  ! [007040] Impact parameter (m)
      IF ( Values(IE+7) < 6200000.0_dp .OR. &
           Values(IE+7) > 6600000.0_dp )    &
           Values(IE+7) = MISSING

      Values(IE+8) = REAL(ROdata%Lev1b%BAngle_L1(in))  ! [015037] B/angle (rad)
      IF ( Values(IE+8) < -0.001 .OR. &
           Values(IE+8) >  0.08288 )  &
           Values(IE+8) = MISSING

      Values(IE+9) = FOstats                           ! [008023] 1st order stats (rms)

      Values(IE+10) = REAL(ROdata%Lev1b%BAngle_L1_Sigma(in)) ! [015037] B/angle error (rad)
      IF ( Values(IE+10) < 0.0 .OR. &
           Values(IE+10) > 0.009485 ) &                  ! 1/8 (-3 bits) from 015037
           Values(IE+10) = MISSING

      Values(IE+11) = MISSING                          ! [008023] 1st order stats (off)

! L2 data

      Values(IE+12) = FreqL2                           ! [002121] L2=1.2Ghz

      Values(IE+13) = REAL(ROdata%Lev1b%Impact_L2(in)) ! [007040] Impact parameter (m)
      IF ( Values(IE+13) < 6200000.0_dp .OR. &
           Values(IE+13) > 6600000.0_dp )    &
           Values(IE+13) = MISSING

      Values(IE+14) = REAL(ROdata%Lev1b%BAngle_L2(in)) ! [015037] B/angle (rad)
      IF ( Values(IE+14) < -0.001 .OR. &
           Values(IE+14) >  0.08288 )  &
           Values(IE+14) = MISSING

      Values(IE+15) = FOstats                          ! [008023] 1st order stats (rms)

      Values(IE+16) = REAL(ROdata%Lev1b%BAngle_L2_Sigma(in)) ! [015037] B/angle error (rad)
      IF ( Values(IE+16) < 0.0 .OR.   &
           Values(IE+16) > 0.009485 ) &                ! 1/8 (-3 bits) from 015037
           Values(IE+16) = MISSING

      Values(IE+17) = MISSING                          ! [008023] 1st order stats (off)
   END IF

! Corrected bending angle (always encoded)

    Values(IE+18) = FreqLc                             ! [002121] corrected

    Values(IE+19) = REAL(ROdata%Lev1b%Impact(in))      ! [007040] Impact parameter (m)
    IF ( Values(IE+19) < 6200000.0_dp .OR. &
         Values(IE+19) > 6600000.0_dp )    &
         Values(IE+19) = MISSING

    Values(IE+20) = REAL(ROdata%Lev1b%BAngle(in))      ! [015037] B/Ang (rad)
    IF ( Values(IE+20) < -0.001 .OR. &
         Values(IE+20) >  0.08288 )  &
         Values(IE+20) = MISSING

    Values(IE+21) = FOstats                            ! [008023] 1st order stats (rms)

    Values(IE+22) = REAL(ROdata%Lev1b%BAngle_Sigma(in)) ! [015037] Error in B/Ang (rad)
    IF ( Values(IE+22) < 0.0 .OR.   &
         Values(IE+22) > 0.009485 ) &                  ! 1/8 (-3 bits) from 015037
         Values(IE+22) = MISSING

    Values(IE+23) = MISSING                            ! [008023] 1st order stats (off)

    Values(IE+24) = REAL(ROdata%Lev1b%Bangle_Qual(in)) ! [033007] Percent confidence
    IF ( Values(IE+24) <   0.0  .OR. &
         Values(IE+24) > 100.0 )     &
         Values(IE+24) = MISSING

    IE = IE + 23
  END DO
  IE = IE + 1

!-------------------------------------------------------------
! 5. Level 2a data (derived refractivity profile)
!-------------------------------------------------------------

! Interpolation thinning may generate a set of fixed geometric height
! levels but no valid Level 2a profile data if no input L2a - reject such
! empty N-profiles

  IF ( ROdata%Lev2a%Missing ) THEN
    CALL message(msg_diag, "Rejecting empty Level 2a (N) profile")
    ROdata%Lev2a%Npoints = 0
  END IF

  Values(IE+1)    = ROdata%Lev2a%Npoints               ! [031002] Replication factor
  nRepFac         = nRepFac + 1
  RepFac(nRepFac) = NINT(Values(IE+1))

  DO in = 1, ROdata%Lev2a%Npoints

    Values(IE+2) = REAL(ROdata%Lev2a%Alt_Refrac(in))   ! [007007] Height amsl (m)
    IF ( Values(IE+2) <  -1000.0_dp .OR. &
         Values(IE+2) > 100000.0_dp )    &
         Values(IE+2) = MISSING

    Values(IE+3) = REAL(ROdata%Lev2a%Refrac(in))       ! [015036] Refrac (N-units)
    IF ( Values(IE+3) <   0.0 .OR. &
         Values(IE+3) > 524.28 )   &
         Values(IE+3) = MISSING

    Values(IE+4) = FOstats                             ! [008023] 1st order stats (rms)

    Values(IE+5) = REAL(ROdata%Lev2a%Refrac_Sigma(in)) ! [015036] Refrac error (N-units)
    IF ( Values(IE+5) <  0.0 .OR. &
         Values(IE+5) > 16.38 )   &                    ! 1/32 (-5 bits) from 015036
         Values(IE+5) = MISSING

    Values(IE+6) = MISSING                             ! [008023] 1st order stats (off)

    Values(IE+7) = REAL(ROdata%Lev2a%Refrac_Qual(in))  ! [033007] Percent confidence
    IF ( Values(IE+7) <   0.0  .OR. &
         Values(IE+7) > 100.0 )     &
         Values(IE+7) = MISSING

   IE = IE + 6
  END DO
  IE = IE + 1

!-------------------------------------------------------------
! 6. Level 2b data (retrieved P,T,q profile)
!-------------------------------------------------------------

! Interpolation thinning may generate a set of fixed geopotential height
! levels but no valid Level 2b profile data if no input L2b - reject such
! empty P,T,q-profiles


  IF ( ROdata%Lev2b%Missing ) THEN
    CALL message(msg_diag, "Rejecting empty Level 2b (T,q,P) profile")
    ROdata%Lev2b%Npoints = 0
  END IF

  Values(IE+1)    = ROdata%Lev2b%Npoints               ! [031002] Replication factor
  nRepFac         = nRepFac + 1
  RepFac(nRepFac) = NINT(Values(IE+1))

  DO in = 1, ROdata%Lev2b%Npoints

    Values(IE+2) = REAL(ROdata%Lev2b%Geop(in))         ! [007009] Geopot ht (gpm)
    IF ( Values(IE+2) <  -1000.0_dp .OR. &
         Values(IE+2) > 100000.0_dp )    &
         Values(IE+2) = MISSING

    Values(IE+3) = REAL(ROdata%Lev2b%Press(in)) * 1E2  ! [010004] Pressure (Pa)
    IF ( Values(IE+3) <=      0.0_dp .OR. &            ! Min. 0.1hPa
         Values(IE+3) >  150000.0_dp )    &
         Values(IE+3) = MISSING

    Values(IE+4) = REAL(ROdata%Lev2b%Temp(in))         ! [012001] Temperature (K)
    IF ( Values(IE+4) < 150.0 .OR. &
         Values(IE+4) > 350.0 )    &
         Values(IE+4) = MISSING

    Values(IE+5) = REAL(ROdata%Lev2b%SHum(in)) * 1E-3  ! [013001] Spec/humidity (Kg/Kg)
    IF ( Values(IE+5) <  0.0 .OR. &
         Values(IE+5) >  0.16 )   &
         Values(IE+5) = MISSING

    Values(IE+6) = FOstats                             ! [008023] 1st order stats (rms)

    Values(IE+7) = REAL(ROdata%Lev2b%Press_Sigma(in)) * 1E2 ! [010004] Pressure error (Pa)
    IF ( Values(IE+7) <   0.0 .OR. &
         Values(IE+7) > 620.0 )    &
         Values(IE+7) = MISSING

    Values(IE+8) = REAL(ROdata%Lev2b%Temp_Sigma(in))   ! [012001] Temperature error (K)
    IF ( Values(IE+8) < 0.0 .OR. &
         Values(IE+8) > 6.2 )    &
         Values(IE+8) = MISSING

    Values(IE+9) = REAL(ROdata%Lev2b%SHum_Sigma(in)) * 1E-3  ! [013001] S/Hum error (Kg/Kg)
    IF ( Values(IE+9) < 0.0 .OR. &
         Values(IE+9) > 0.0051 ) &
         Values(IE+9) = MISSING

    Values(IE+10) = MISSING                            ! [008023] 1st order stats (off)

    Values(IE+11) = REAL(ROdata%Lev2b%Meteo_Qual(in))  ! [033007] Percent confidence
    IF ( Values(IE+11) <   0.0  .OR. &
         Values(IE+11) > 100.0 )     &
         Values(IE+11) = MISSING

   IE = IE + 10
  END DO
  IE = IE + 1

!-------------------------------------------------------------
! 7. Level 2c data (retrieved surface params)
!-------------------------------------------------------------

  Values(IE+1) = 0                                     ! [008003] Vertical sig. (surf)

  VALUES(IE+2) = REAL(ROdata%Lev2c%Geop_Sfc)           ! [007009] Geoptot.Ht. (of surf) (gpm)
  IF ( Values(IE+2) < -1000.0_dp .OR. &
       Values(IE+2) > 10000.0_dp )    &
       Values(IE+2) = MISSING

  Values(IE+3) = REAL(ROdata%Lev2c%Press_Sfc) * 1E2    ! [010004] Surface pressure (Pa)
  IF ( Values(IE+3) <      0.0_dp .OR. &
       Values(IE+3) > 150000.0_dp )    &
       Values(IE+3) = MISSING

  Values(IE+4) = FOstats                               ! [008023] 1st order stats (rms)

  Values(IE+5) = REAL(ROdata%Lev2c%Press_Sfc_Sigma) * 1E2 ! [010004] S/press error (Pa)
  IF ( Values(IE+5) <   0.0 .OR. &
       Values(IE+5) > 620.0 )    &
       Values(IE+5) = MISSING

  Values(IE+6) = MISSING                               ! [008023] 1st order stats (off)

  Values(IE+7) = REAL(ROdata%Lev2c%Press_Sfc_Qual)     ! [033007] Percent confidence
  IF ( Values(IE+7) <   0.0  .OR. &
       Values(IE+7) > 100.0 )     &
       Values(IE+7) = MISSING

  nValues = IE + 7                                     ! Total no. of values

!-------------------------------------------------------------
! 8. Tidy up before return
!-------------------------------------------------------------

  CALL message_set_routine(routine)

END SUBROUTINE ConvertROPPtoBUFR
!--------------------------------------------------------------------

END MODULE ropp2bufr
