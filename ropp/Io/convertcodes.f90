! $Id: convertcodes.f90 4452 2015-01-29 14:42:02Z idculv $

SUBROUTINE ConvertCodes ( ROdata, & ! (inout)
                          Gclass, & ! (inout)
                          Gcode,  & ! (inout)
                          Lcode,  & ! (inout)
                          Icode,  & ! (inout)
                          Ocode,  & ! (inout)
                          Scode,  & ! (inout)
                          Ccode,  & ! (inout)
                          Bcode,  & ! (inout)
                          ind )     ! (in)

!****s* ropp2bufr/ConvertCodes *
!
! NAME
!   ConvertCodes
!
! SYNOPSIS
!  Convert header codes between ROPP and BUFR
!
!    USE ropp_io_types
!    TYPE (roprof) rodata
!    INTEGER :: gclass,gcode,lcode,icode,ocode,scode,bcode,ind
!    CHARACTER (LEN=4) :: ccode
!    ind = 1 ! to convert ROPP-->BUFR, ind = -1 for BUFR-->ROPP
!    CALL convertcodes(rodata,&
!                      gclass,gcode,lcode,icode,ocode,scode,ccode,bcode,&
!                      ind)
!
! INPUTS
!   ROdata  dtyp  RO data structure                     [ind>0]
!   Gclass   int  GNSS code (Satellite Class)           [ind<=0]
!   Gcode    int  GNSS PRN  (Platform Tx ID)            [ind<=0]
!   Lcode    int  LEO  code (Satellite ID)              [ind<=0]
!   Icode    int  Instrument code (Instrument ID)       [ind<=0]
!   Ocode    int  Originating (processing) Centre code  [ind<=0]
!   Scode    int  Originating Sub-centre code           [ind<=0]
!   Ccode    chr  Originating (GTS) centre ICAO code    [ind<=0]
!   Bcode    int  Background generating centre code     [ind<=0]
!   ind      int  ROPP-->BUFR if >0, else BUFR-->ROPP
!
! OUTPUTS
!   ROdata  dtyp  RO data structure                     [ind<=0]
!   Gclass   int  GNSS code (Satellite Class)           [ind>0]
!   Gcode    int  GNSS PRN  (Platform Tx ID)            [ind>0]
!   Lcode    int  LEO  code (Satellite ID)              [ind>0]
!   Icode    int  Instrument code (Instrument ID)       [ind>0]
!   Ocode    int  Originating (processing) Centre code  [ind>0]
!   Scode    int  Originating Sub-centre code           [ind>0]
!   Ccode    chr  Originating (GTS) centre ICAO code    [ind>0]
!   Bcode    int  Background generating centre code     [ind>0]
!
! MODULES
!   messages    - ROPP message library
!   ropp_utils  - ROPP utilities
!   ropp_io     - ROPP file I/O support
!
! CALLS
!   GETENV
!   Get_IO_Unit
!   message
!   message_get_routine
!   message_set_routine
!
! CALLED BY
!   ConvertBUFRtoROPP
!   ConvertROPPtoBUFR
!
! FILES
!   roppbufrcodes.nl  - in path BUFR_TABLES, BUFR_LIBRARY or PWD
!
! DESCRIPTION
!   Converts from character-based codes (as defined for ROPP) to numeric codes
!   suitable for BUFR encoding, if ind>0, else vice-versa.
!   The code conversion is driven by a set of look-up tables, which are read
!   from a NAMELIST file 'roppbufrcodes.nl' which is expected in the directory
!   path defined by at least one of the environment variables BUFR_TABLES
!   (ECMWF), BUFR_LIBRARY (MetDB) or PWD (searched in that order).
!   If this file cannot be found or opened, a warning is issued and an in-built
!   default set of tables is used instead.
!
! REFERENCES
!   1. Manual on Codes: International Codes, Part B & Part C.
!      WMO-No. 306, World Meteorological Organisation, Geneva.
!      http://www.wmo.int/pages/prog/www/WMOCodes/VolumeI2.html
!   2. Location Indicators. ICAO Document 7910/138
!      ISBN 978-92-9231-677-8, Ed.138, December 2010.
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
  USE ropp_utils,   ONLY: Get_IO_Unit
  USE ropp_io_types

  IMPLICIT NONE

! Fixed values

  INTEGER,  PARAMETER :: NVIND = 2147483647   ! Missing data flag value

!  NB: no. of elements given in NAMELIST file parameters must not
!  exceed these values - increase values below if necessary.

  INTEGER, PARAMETER :: ntx = 6               ! Max. no. of GNSS Tx types
  INTEGER, PARAMETER :: nrx = 25              ! Max. no. of LEO  Rx types
  INTEGER, PARAMETER :: noc = 10              ! Max. no. of orig. centre types
  INTEGER, PARAMETER :: nbg = 10              ! Max. no. of b/g centre types

  INTEGER, PARAMETER :: nep  = 3              ! No. of environment paths
  CHARACTER (LEN=*), PARAMETER :: NLenv(nep) = (/"BUFR_TABLES ", & ! N/L paths
                                                 "BUFR_LIBRARY", &
                                                 "PWD         " /)
  CHARACTER (LEN=*), PARAMETER :: NLdsn = "roppbufrcodes.nl" ! N/L file name

! Argument list parameters

  TYPE ( ROprof ),   INTENT(INOUT) :: Rodata  ! ROPP data structure
  INTEGER,           INTENT(INOUT) :: Gclass  ! GNSS class value
  INTEGER,           INTENT(INOUT) :: Gcode   ! GNSS PRN
  INTEGER,           INTENT(INOUT) :: Lcode   ! LEO  code value
  INTEGER,           INTENT(INOUT) :: Icode   ! Instrument code value
  INTEGER,           INTENT(INOUT) :: Ocode   ! Origin. centre code value
  INTEGER,           INTENT(INOUT) :: Scode   ! Sub-centre code value
  CHARACTER (LEN=*), INTENT(INOUT) :: Ccode   ! ICAO code
  INTEGER,           INTENT(INOUT) :: Bcode   ! B/G generator code value
  INTEGER,           INTENT(IN)    :: ind     ! RO->code if >1 else code->RO

! Define arrays for chararacter (ROPP) & numeric (BUFR code) lists.
! Set some defaults in case the NAMELISTs can't be read. NAMELIST
! values will overwrite these defaults. Include some dummy spares so
! that extra ones can be defined in the NAMELIST _without_ having to
! change the array sizes (up to the max. values) and rebuilding
! the program.

! Satellite Classification (GNSS Tx constellation) (Code Table 002020)
! NB: code 404 for Beidou is provisional

  CHARACTER (LEN=1),  DIMENSION(ntx) :: GNSlist = &
                      (/ "U",   "G", "R", "E", "B", "U"   /)
  INTEGER,            DIMENSION(ntx) :: GNScode = &
                      (/ NVIND, 401, 402, 403, 404, NVIND /)

! Satellite Identifier (LEO Rx mission) (Code Table 001007 or CCT C-5)
! and associated Instrument Type        (Code Table 002019 or CCT C-8)
! NB: instrument code 956 for GNOS (FY-3C) is provisional

  CHARACTER (LEN=4),  DIMENSION(nrx) :: LEOlist =        &
                      (/ "UNKN", "OERS", "CHMP",         &
                                 "SUNS", "SACC",         &
                                 "GRAA", "GRAB",         &
                                 "C001", "C002", "C003", &
                                 "C004", "C005", "C006", &
                                 "META", "METB", "METC", &
                                 "TSRX", "TDMX", "PAZE", &
                                 "OSAT", "CNOF", "FY3C", &
                                 "UNKN", "UNKN", "UNKN" /)
  INTEGER,            DIMENSION(nrx) :: LEOcode =        &
                        (/ NVIND,  040,    041,          &
                                   800,    820,          &
                                   722,    723,          &
                                   740,    741,    742,  &
                                   743,    744,    745,  &
                                   004,    003,    005,  &
                                   042,    043,    044,  &
                                   421,    786,    522,  &
                                  NVIND, NVIND, NVIND   /)
  INTEGER,            DIMENSION(nrx) :: Inscode =        &
                        (/ NVIND,  102,    102,          &
                                   102,    102,          &
                                   102,    102,          &
                                   102,    102,    102,  &
                                   102,    102,    102,  &
                                   202,    202,    202,  &
                                   103,    103,    103,  &
                                   287,    102,    956,  &
                                  NVIND, NVIND, NVIND   /)

! List of (BUFR) Originating Centre IDs & their BUFR codes
! (Code Table 001033, CCT C-1, or 001035, CCT C-11)
! The (Processing) Sub-centre code should be valid for the
! associated Originating Centre code (for which Sub-Centre is 0).
! (Code Table 001034, CCT C-12)
! plus associated ICAO Location Indicator codes (for GTS routing headers)
! (ICAO Document 7910: Location Indicators)
! NB: Origin code 'CMA' for FY-3C is provisional
! NB: ICAO code 'BABJ' for CMA is provisional
!
  CHARACTER (LEN=8), DIMENSION(noc) :: ORGlist = &
                      (/ "UNKNOWN ", "DMI     ", "GFZ     ", &
                                     "METO    ", "UCAR    ", &
                                     "NESDIS  ", "EUMETSAT", &
                                     "CMA     ",             &
                                     "UNKNOWN ", "UNKNOWN " /)
  INTEGER,            DIMENSION(noc) :: ORGcode =            &
                      (/ NVIND,      094,        078,        &
                                     074,        060,        &
                                     160,        254,        &
                                     038,                    &
                                     NVIND,      NVIND /)
  INTEGER,            DIMENSION(noc) :: Subcode =            &
                      (/ 000,        000,        173,        &
                                     000,        000,        &
                                     000,        000,        &
                                     000,                    &
                                     000,        000 /)
  CHARACTER (LEN=35), DIMENSION(noc) :: ORGname = &
                                    (/ "                                   ", &
                                       "(ROM SAF)                          ", &
                                       "Helmholtz Centre, Potsdam          ", &
                                       "Met Office, Exeter                 ", &
                                       "Boulder                            ", &
                                       "Washington                         ", &
                                       "Darmstatdt                         ", &
                                       "Beijing                            ", &
                                       "                                   ", &
                                       "                                   " /)
  CHARACTER (LEN=4), DIMENSION(noc)  :: ICAOcode = &
                      (/ "ZZZZ",     "EKMI",     "EDZW",  &
                                     "EGRR",     "KWBC",  &
                                     "KNES",     "EUMS",  &
                                     "BABJ",              &
                                     "ZZZZ",     "ZZZZ"  /)

! Orginating Centre (background profile)
! (Code Table 001033, CCT C-1, or 001035, CCT C-11)
! NB: Originating centre code 38 for CMA is provisional

  CHARACTER (LEN=8), DIMENSION(nbg) :: BGDlist = &
                      (/ "UNKNOWN ", "ECMWF   ", "DMI     ",             &
                                     "METO    ", "NCEP    ",             &
                                     "CMA     ", "NONE    ",             &
                                     "UNKNOWN ", "UNKNOWN ", "UNKNOWN " /)
  INTEGER,            DIMENSION(nbg) :: BGDcode =                        &
                       (/ NVIND,    98,        94,                       &
                                    74,         7,                       &
                                    38,     NVIND,                       &
                                 NVIND,     NVIND,     NVIND            /)

! Local variables

  CHARACTER (LEN=235) :: dir = " "       ! Translated BUFR directory (path)
  CHARACTER (LEN=255) :: FileSpec        ! Full sequence file name
  CHARACTER (LEN=50)  :: routine         ! Saved routine name
  INTEGER             :: NLunit          ! NAMELIST file unit no.
  INTEGER             :: i, j, l         ! loop counter/indices
  INTEGER             :: ierr            ! I/O error
  LOGICAL             :: exists          ! File exists flag
  LOGICAL             :: first = .TRUE.  ! First call flag

! Namelist parameters

  NAMELIST /GNScodes/ GNSlist, GNScode
  NAMELIST /LEOcodes/ LEOlist, LEOcode, Inscode
  NAMELIST /ORGcodes/ ORGlist, ORGcode, Subcode, ORGname, ICAOcode
  NAMELIST /BGDcodes/ BGDlist, BGDcode

  SAVE first

  CALL message_get_routine ( routine )
  CALL message_set_routine ( "ConvertCodes" )

!---------------------------------------------------
! 1. Find & open codes NAMELIST file, read lists
!---------------------------------------------------

  IF ( first ) THEN
    NLunit = Get_IO_Unit()
    DO i = 1, nep
      CALL GETENV ( TRIM(NLenv(i)), dir )
      l = LEN_TRIM(dir)
      IF ( l > 0 .AND. TRIM(dir) /= TRIM(NLenv(i)) ) THEN
        IF ( dir(l:l) /= "/" ) dir(l+1:l+1) = "/"
        FileSpec = ADJUSTL(TRIM(dir)//NLdsn)
        INQUIRE ( FILE=FileSpec, EXIST=exists )
        IF ( exists ) EXIT
      END IF
    END DO

    IF ( exists ) THEN
      OPEN ( UNIT=NLunit, FILE=FileSpec, ACTION="READ", IOSTAT=ierr )
      IF ( ierr == 0 ) THEN
        READ  ( UNIT=NLunit, NML=GNScodes, IOSTAT=ierr )
        IF ( ierr /= 0 ) &
          CALL message ( msg_warn, "Error loading NAMELIST GNScodes" )
        READ  ( UNIT=NLunit, NML=LEOcodes, IOSTAT=ierr )
        IF ( ierr /= 0 ) &
          CALL message ( msg_warn, "Error loading NAMELIST LEOcodes" )
        READ  ( UNIT=NLunit, NML=ORGcodes, IOSTAT=ierr )
        IF ( ierr /= 0 ) &
          CALL message ( msg_warn, "Error loading NAMELIST ORGcodes" )
        READ  ( UNIT=NLunit, NML=BGDcodes, IOSTAT=ierr )
        IF ( ierr /= 0 ) &
          CALL message ( msg_warn, "Error loading NAMELIST BGDcodes" )
        CLOSE ( UNIT=NLunit )
        IF ( ierr == 0 ) THEN
          CALL message ( msg_diag, "Loaded "//TRIM(FileSpec) )
        ELSE
          CALL message ( msg_warn, "Error loading ROPP-BUFR codes "// &
                                   "NAMELIST file" )
          CALL message ( msg_cont, "    ("//TRIM(FileSpec)//")" )
          CALL message ( msg_cont, "    Using default look-up tables for items "// &
                                   " not loaded" )
        END IF
      ELSE
        CALL message ( msg_warn, "ROPP-BUFR codes NAMELIST file"// &
                                 " could not be opened." )
        CALL message ( msg_cont, "    ("//TRIM(FileSpec)//")" )
        CALL message ( msg_cont, "    Using default look-up tables" )
      END IF
    ELSE
      CALL message ( msg_warn, "ROPP-BUFR codes NAMELIST file "// &
                               TRIM(NLdsn)//" could not be found." )
      CALL message ( msg_cont, "    Using default look-up tables" )
    END IF

    first = .FALSE.
  END IF

!---------------------------------------------------
! 2. Look up numeric (BUFR) code from character (ROPP)
!---------------------------------------------------

  IF ( ind >= 1 ) THEN

! Defaults

    Lcode  = NVIND
    Icode  = NVIND
    Gclass = NVIND
    Gcode  = NVIND
    Ocode  = NVIND
    Scode  = NVIND
    Ccode  = "ZZZZ"
    Bcode  = NVIND

! LEO Rx ID code (satellite & hence instrument)

    i = nrx
    DO WHILE ( i > 0 .AND. &
               LEOlist(i) /= ROdata%LEO_id )
      i = i - 1
    END DO
    IF ( i > 0 ) THEN
      Lcode = LEOcode(i)
      Icode = Inscode(i)
    END IF

! GNSS Tx ID code (satellite class) & separate PRN

    i = ntx
    DO WHILE ( i > 0 .AND. &
               GNSlist(i) /= ROdata%GNS_id(1:1) )
      i = i - 1
    END DO
    IF ( i > 0 ) Gclass = GNScode(i)
    READ ( ROdata%GNS_id(2:4), FMT=*, IOSTAT=ierr ) Gcode
    IF ( ierr /= 0 .OR. &
         Gcode < 0 .OR. &
         Gcode > 32 ) Gcode = NVIND

! Originating (encoding) centre code and associated
! sub-centre (processing) & ICAO (GTS node) codes

    i = noc
    DO WHILE ( i > 0 .AND. &
               ORGlist(i)(1:3) /= ROdata%Processing_Centre(1:3) )
      i = i - 1
    END DO
    Ocode = ORGcode(i)
    Scode = SUBcode(i)
    Ccode = ICAOcode(i)

! Look up background generator centre code

    i = nbg
    DO WHILE ( i > 0 .AND. &
               BGDlist(i)(1:3) /= ROdata%BG%Source(1:3) )
      i = i - 1
    END DO
    Bcode = BGDcode(i)

!---------------------------------------------------
! 3. Look up character (ROPP) code from numeric (BUFR)
!---------------------------------------------------

  ELSE

! Defaults

    ROdata%LEO_id            = "UNKN"
    ROdata%GNS_id            = "U999"
    ROdata%Processing_Centre = "UNKNOWN"
    ROdata%BG%Source         = "UNKNOWN"

! LEO Rx ID code (Satellite)

    i = nrx
    DO WHILE ( i > 0 .AND. &
               Lcode /= LEOcode(i) )
      i = i - 1
    END DO
    IF ( i > 0 ) ROdata%LEO_id = LEOlist(i)

! GNSS Tx ID code (from satellite class) & add PRN

    i = ntx
    DO WHILE ( i > 0 .AND. &
               Gclass /= GNScode(i) )
      i = i - 1
    END DO
    IF ( i > 0 ) ROdata%GNS_id(1:1) = GNSlist(i)
    IF ( Gcode < 0 .OR. Gcode > 999 ) Gcode = 999
    WRITE ( ROdata%GNS_id(2:4), &
            FMT="(I3.3)",       &
            IOSTAT=ierr ) Gcode

! Originating (RO processing) centre code

    i = noc
    DO WHILE ( i > 0 .AND. &
               Ocode /= ORGcode(i) )
      i = i - 1
    END DO
    IF ( i > 0 ) THEN
      j = MAX ( LEN_TRIM ( ORGlist(i) ), 4 )
      ROdata%Processing_Centre = ORGlist(i)(1:j) &
                       // " " // ORGname(i)
      Ccode = ICAOcode(i)
    END IF

! Background generating centre code

    i = nbg
    DO WHILE ( i > 0 .AND. &
               Bcode /= BGDcode(i) )
      i = i - 1
    END DO
    IF ( i > 0 ) ROdata%BG%Source = BGDlist(i)

  ENDIF

  CALL message_set_routine ( routine )

END SUBROUTINE ConvertCodes
