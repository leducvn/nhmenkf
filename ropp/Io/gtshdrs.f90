! $Id: gtshdrs.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE GTShdrs

!****m* BUFR/GTSHDRS/GTShdrs *
!
! NAME
!   GTShdrs    (gtshdrs.f90)
!
! SYNOPSIS
!   Module defining fixed values & subroutines for generating WMO GTS/RMDCN
!   routing headers, etc.
!
!   USE GTShdrs
!   INTEGER            :: CSN, DateTime(6), un
!   CHARACTER (LEN=50) :: CSNfile
!   CHARACTER (LEN=6)  :: TTAAII
!   CHARACTER (LEN=4)  :: CCCC
!   CHARACTER (LEN=10) :: cIPH
!   CHARACTER (LEN=31) :: cARH
!   CHARCATER (LEN=4)  :: cEOM
!   CHARACTER (LEN=1)  :: a2
!   INTEGER            :: lIPH, lARH, lEOM
!   INTEGER            :: Lat, Lon
!   REAL               :: Lats(n), Lons(n)
!   GTSHDR_DEBUG = .TRUE.  ! optionally enable debug mode
!   CALL GTShdrIPH ( cIPH, lIPH )
!   CALL GTShdrARH ( CSN, TTAAII, Lat, Lon, CCCC, DateTime, &
!                    cARH, lARH )
!   CALL GTShdrARH ( CSN, TTAAII, Lats, Lons, CCCC, DateTime, &
!                    cARH, lARH )
!   CALL GTShdrEOM ( cEOM, lEOM )
!   CALL GTShdrCSN ( CSNfile, CSN, 'R'|'W' )
!   un = GTShdrUNT ()
!   a2 = GTShdrGAD ( Lat, Lon )
!
! USED BY
!   BUFR or GRIB encoding applications requiring GTS headers
!
! DESCRIPTION
!   A 'naked' BUFR (or GRIB) "message" starts with the characters 'BUFR'
!   (or 'GRIB') and terminates with '7777' (Ref.1). However, such a message
!   cannot be transmitted over WMO GTS (or RMDCN) links without some form
!   of routing information. This is provided by prepending a standard
!   "Abbreviated Routing Header" (ARH) (Ref.2) and terminating with an
!   "end-of-message" (EOM) sequence. A message with ARH & EOM wrappers
!   is known as a "bulletin" and is the basic block of bytes transmitted
!   over the GTS/RMDCN.
!   Further, as protocols allow for more than one bulletin in a physical
!   file, each bulletin is usually (though not universally, depending on
!   local GTS node procedures) prefixed by a "Internet Protocol Header" (IPH)
!   or leader sequence which contains the total bulletin length in bytes
!   (including the IPH itself) and a data format identifier. The file is
!   terminated with a dummy IPH after the last bulletin (Ref.3).
!   This module provides routines to generate an IPH (actual or dummy), an
!   ARH and an EOM sequence. The routines return these sequences of bytes
!   as plain ASCII character strings (containing non-printing control
!   characters) which can be output 'as-is' using a write routine such as
!   METDB_CWRITE() from the MetDB BUFR library (Ref.4) or - with pre-packing
!   to integers - with PBWRITE() from the ECMWF BUFR library (Ref.5).
!
! EXAMPLES
!   1) Generate an ARH with CSN sequence cycling between runs (example
!      appropriate for E-GVAP data):
!
!     USE GTShdrs
!     INTEGER            :: CSN=0, DateTime(6) ! YR,MO,DY,HR,MN,SC [0,0,13,16,35,0]
!     CHARACTER (LEN=5)  :: TTAAII="ISX?14"    ! '?' is filled from lat/lon info
!     CHARACTER (LEN=4)  :: CCCC="EGRR"        ! UK Met Office
!     CHARACTER (LEN=31) :: cARH
!     CHARACTER (LEN=50) :: CSNfile = "~/data/gwv/csn.dat"
!     INTEGER            :: BUFRunit, lARH
!     INTEGER            :: Lat=50, Lon=10     ! Nominal European coverage
!     GTSHDR_DEBUG = .TRUE.                    ! optionally enable debug mode
!     CALL GTShdrCSN ( CSNfile, CSN, 'Read' )  ! read last used CSN (=8)
!     CALL GTShdrARH ( CSN, TTAAII, Lat, Lon, CCCC, DateTime, &
!                      cARH, lARH )            ! generate ARH
!     CALL GTShdrCSN ( CSNfile, CSN, 'Write' ) ! save last CSN (=9)
!     >> ARH = "....009...ISXD14.EGRR.131635..."
!
!   2) Write ARH to BUFR file using MetDB BUFR library:
!
!     INTEGER :: BUFRunit  ! from METDB_COPEN()
!     CALL METDB_CWRITE ( BUFRunit, cARH(1:lARH), lARH )
!
!   3) Write ARH to BUFR file using ECMWF BUFR library
!      (pack bytes from chr to int first):
!
!     INTEGER :: BUFRunit  ! from PBOPEN()
!     INTEGER :: uARH(31), pARH(8), ierr
!     DO i = 1, lARH
!       uARH(i) = IACHAR(cARH(i:i))
!     END DO
!     CALL SBYTES  ( pARH, uARH, 0, 8, 0, lARH )
!     CALL PBWRITE ( BUFRunit, pARH, lARH, ierr )
!
! REFERENCES
!   1) Manual on the Global Telecommunications System.
!      Vol.1, Part II. WMO, Geneva, 1986 (and updates)
!   2) WMO (2009). Operational Procedures for the GTS, Attachment II-5,
!      Data Designators T1T2A1A2ii in Abbreviated Headings, Table C3.
!      WMO, Geneva, 4 November 2009.
!   3) Nightingale, S. (2008). Exchange of data with the
!      Met Office FROST message switch.
!      FROST FTP Protocol Document, V2.1, 2 October 2008.
!   4) Met Office (2012). Decoding and Encoding BUFR messages.
!      MetDB Technote 1, Rev.1, 29/10/2012 [dmtn1.html].
!   5) Dragosavac, Milan (2009). BUFR User's Guide.
!      ECMWF Operations Department Technical Note, July 2009.
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!---------------------------------------------------------------------------------
!****

! Fixed local parameters

! INTEGER, PARAMETER, PRIVATE :: nbpc = 8                 ! bits per character (or byte)
! INTEGER, PARAMETER, PRIVATE :: nbpw = KIND(nbpc) * nbpc ! bits per word (default INTEGER)
! INTEGER, PARAMETER, PRIVATE :: nbsk = 0                 ! no. bits to skip when packing
  INTEGER, PARAMETER          :: GTSHDR_FLGVAL = 999      ! A2 flag value

  LOGICAL                     :: GTSHDR_DEBUG  = .FALSE.  ! local debug flag

! Explicit interfaces

  INTERFACE GTShdrIPH
    SUBROUTINE GTShdrIPH ( Lenm, cIPH, lIPH )
      INTEGER,           INTENT(IN)  :: Lenm
      CHARACTER (LEN=*), INTENT(OUT) :: cIPH
      INTEGER,           INTENT(OUT) :: lIPH
    END SUBROUTINE GTShdrIPH
  END INTERFACE GTShdrIPH

  INTERFACE GTShdrARH
    SUBROUTINE GTShdrARH_sca ( CSN, TTAAII, Lat, Lon, CCCC, DateTime, cARH, lARH )
      INTEGER,           INTENT(INOUT) :: CSN
      CHARACTER (LEN=*), INTENT(IN)    :: TTAAII
      INTEGER,           INTENT(IN)    :: Lat, Lon
      CHARACTER (LEN=*), INTENT(IN)    :: CCCC
      INTEGER,           INTENT(IN)    :: DateTime(:)
      CHARACTER (LEN=*), INTENT(OUT)   :: cARH
      INTEGER,           INTENT(OUT)   :: lARH
    END SUBROUTINE GTShdrARH_sca
    SUBROUTINE GTShdrARH_arr ( CSN, TTAAII, Lats, Lons, CCCC, DateTime, cARH, lARH )
      INTEGER,           INTENT(INOUT) :: CSN
      CHARACTER (LEN=*), INTENT(IN)    :: TTAAII
      REAL,              INTENT(IN)    :: Lats(:), Lons(:)
      CHARACTER (LEN=*), INTENT(IN)    :: CCCC
      INTEGER,           INTENT(IN)    :: DateTime(:)
      CHARACTER (LEN=*), INTENT(OUT)   :: cARH
      INTEGER,           INTENT(OUT)   :: lARH
    END SUBROUTINE GTShdrARH_arr
  END INTERFACE GTShdrARH

  INTERFACE GTShdrEOM
    SUBROUTINE GTShdrEOM ( cEOM, lEOM )
      CHARACTER (LEN=*), INTENT(OUT) :: cEOM
      INTEGER,           INTENT(OUT) :: lEOM
    END SUBROUTINE GTShdrEOM
  END INTERFACE GTShdrEOM

  INTERFACE GTShdrCSN
    SUBROUTINE GTShdrCSN ( CSNfile, CSN, RorW )
      CHARACTER (LEN=*), INTENT(IN)    :: CSNfile
      INTEGER,           INTENT(INOUT) :: CSN
      CHARACTER (LEN=*), INTENT(IN)    :: RorW
    END SUBROUTINE GTShdrCSN
  END INTERFACE GTShdrCSN

  INTERFACE GTShdrUNT
    FUNCTION GTShdrUNT() RESULT(unit)
      INTEGER :: unit
    END FUNCTION GTShdrUNT
  END INTERFACE GTShdrUNT

  INTERFACE GTShdrGAD
    FUNCTION GTShdrGAD ( Lat, Lon ) RESULT(A2)
      INTEGER, INTENT(IN) :: Lat, Lon
      CHARACTER (LEN=1)   :: A2
    END FUNCTION GTShdrGAD
  END INTERFACE GTShdrGAD

END MODULE GTShdrs
!---------------------------------------------------------------------------------

SUBROUTINE GTShdrIPH ( Lenm, & ! (in)
                       cIPH, & ! (out)
                       lIPH )  ! (out)

!---------------------------------------------------------------------------------
!****s* BUFR/GTSHDRS/GTShdrIPH *
!
! NAME
!   GTShdrIPH   (gtshdrs.f90)
!
! SYNOPIS
!   Generate IP start-of-message header (IPH) sequence
!
!   USE GTShdrsMod
!   CHARACTER (LEN=10) :: cIPH
!   INTEGER :: Lenm, lIPH
!   CALL GTShdrIPH ( Lenm, cIPH, lIPH )
!
! INPUT
!   Lenm  int  Length of GTS bulletin in bytes (including ARH but not IPH)
!              to set into IPH (0-99999999)
!
! OUTPUT
!   cIPH  chr  Chracter string of at least (10) for IPH sequence
!   lIPH  int  Length of IPH in bytes (normally 10)
!
! DESCRIPTION
!   Generates a 10-byte start-of-message header sequence for IP transmission
!   protocols according to WMO Format 00 (see Reference 1).
!   This header should be output immediately before the GTS bulletin
!   (BUFR message plus ARH & start/end wrappers generated by GTShdrARH &
!   GTShdrEOM).
!   The length of the bulletin, plus the IP header itself, is set into
!   the header. If Len is given as 0, this is assumed to be a dummy
!   IP header, and 10 '0's are returned in IPH. This dummy IPH is used to
!   terminate a file containing one or more bulletins intended for
!   transmission using IP (GTS/RMDCN etc). The IPH sequence is returned
!   as a plain ASCII string.
!
! REFERENCES
!   1) Nightingale, S. (2008). Exchange of data with the
!      Met Office FROST message switch.
!      FROST FTP Protocol Document, V2.1, 2 October 2008.
!
! SEE ALSO
!   GTShdrARH() - to generate a GTS Abbreviated Routing Header
!   GTShdrEOM() - to generate a GTS End-of-Message sequence
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!---------------------------------------------------------------------------------
!****
!
  IMPLICIT NONE

! Fixed values

  INTEGER, PARAMETER   :: niph = 10  ! Length of IPH

! Argument list parameters

  INTEGER,           INTENT(IN)  :: Lenm ! GTS bulletin length (incl. ARH+EOM)
  CHARACTER (LEN=*), INTENT(OUT) :: cIPH ! Plain character IPH
  INTEGER,           INTENT(OUT) :: lIPH ! Length of IPH (bytes)

! Check IPH is big enough - need at least 10 characters

  IF ( LEN(cIPH) < niph ) THEN
    cIPH = " "
    WRITE ( *, FMT="(A,I2)" ) "ERROR: [GTShdrIPH] String is too small for"//&
                              " IPH. LEN must be at least ",niph

! If Lenm<=0 generate dummy IPH with 8 zeros (no Format Identifier),
! else the message length is coded into bytes 1-8 with the Format Identifier
! (Bytes 9-10) always '00' (See Ref.1)

  ELSE
    IF ( Lenm >  0 .AND. &
         Lenm <= 999999999 ) THEN
      WRITE ( cIPH, FMT="(I8.8,A2)" ) Lenm, "00"
    ELSE
      cIPH = "00000000"
    END IF

  END IF

  lIPH = LEN_TRIM(cIPH)

END SUBROUTINE GTShdrIPH
!---------------------------------------------------------------------------------

SUBROUTINE GTShdrARH_sca ( CSN,      & ! (inout)
                           TTAAII,   & ! (in)
                           Lat,      & ! (in)
                           Lon,      & ! (in)
                           CCCC,     & ! (in)
                           DateTime, & ! (in)
                           cARH,     & ! (out)
                           lARH )      ! (out)

!---------------------------------------------------------------------------------
!****s* BUFR/GTSHDRS/GTShdrARH *
!
! NAME
!   GTShdrARH    (gtshdrs.f90)
!
! SYNOPSIS
!   Generate WMO GTS Abbreviated Routing Header (ARH) sequence
!
!   INTEGER            :: CSN, DateTime(6), lARH
!   CHARACTER (LEN=6)  :: TTAAII
!   CHARACTER (LEN=4)  :: CCCC
!   CHARACTER (LEN=31) :: cARH
!   INTEGER            :: Lat, Lon
!   REAL               :: Lats(n), Lons(n)
!   CALL GTShdrARH ( CSN, TTAAII, Lat,  Lon,  CCCC, DateTime, cARH, lARH )
!   CALL GTShdrARH ( CSN, TTAAII, Lats, Lons, CCCC, DateTime, cARH, lARH )
!
! INPUTS
!   CSN       int  Channel (GTS bulletin) sequence number (001-999)
!   TTAAII    chr  6-chr data ID (t1t2a1a2ii) Note that the a2 character
!                  will be set by this routine from (Lat,Lon) info.
!   EITHER:
!     Lat     int  Latitude  of data (-90  to +90 deg or 999)
!     Lon     int  Longitude of data (-180 to +180 or 0 to 360 deg or 999)
!   OR:
!     Lats    int  Array of Latitudes  of data (-90.0  to +90,0 deg)
!     Lons    int  Array of Longitudes of data (-180.0 to +180.0 or
!                  0.0 to 360.0 deg)
!   CCCC      chr  4-chr GTS source node ID (ICAO code) (cccc)
!   DateTime  int  6-element date/time (yr,mon,day,hr,min,sec). Only the
!                  day, hour and minute elements are used, which should
!                  be representative of the data timestamp(s) encoded in
!                  the message.
!
! OUTPUTS
!   CSN       int  Incremented Channel (GTS bulletin) sequence number (001-999)
!   cARH      chr  Character string of at least (31) for ARH sequence
!   lARH      int  Length of ARH in bytes (normally 31)
!
! USES
!   GTShdrGAD
!
! DESCRIPTION
!   Codes data routing header information into WMO/GTS standard 31-byte
!   Abbreviated Routing Header (ARH) (see references) of the form:
!     'sssrrriii'
!   where:
!     'sss' is a 10-byte message start sequence: 'SOH/CR/CR/LF/CSN/CR/CR/LF'
!     'rrr' is the 18-byte routing header:       't1t2a1a2ii cccc YYGGgg'
!     'iii' is a 3-byte data introducer iii:     'CR/CR/LF'.
!   The channel sequence number (CSN) should be a sequential 3-digit
!   integer cycling from 001 to 999, and is coded from the CSN value.
!   Data ID (t1t2a1 and ii) is coded as given.
!   A Geographical Area Designator code (a2) is derived from the given
!   nominal single lat/lon or array of lat/lon of the encoded data - normally
!   'A'-'L' for single observations or for local/regional coverage.
!   - For the scalar form, set |Lon|>=999 to force 'N' (Northern Hemisphere
!     or 'S' (Southern Hemisphere) depending on the sign of Lat;
!     set |Lat|>=999 to force 'T' if -45<=Lon<180 deg (sector 315-180 deg for
!     Europe/Asia) otherwise 'X' (Global); set both |Lat|>=999 and |Lat|>=999
!     to force 'X' (Global).
!   - For the array form, the lat/lon pair ranges are analysed to
!     automatically determine the appropriate nominal scalar lat/lon and/or
!     flag values.
!   - Tip: in place of a hard-coded flag value (999), use module parameter
!     GTSHDR_FLGVAL.
!   Source node (cccc) is the 4-chr ICAO code of the centre coding the
!   data (e.g. 'EGRR' for Met Office, Exeter).
!   Day of month (YY), hour (GG) and minute (gg) are coded from DateTime
!   values.
!   The ARH is returned as a plain ASCII string (containing non-printing
!   control characters). The message should then be terminated with an
!   end-of-message sequence, e.g. as generated with GTShdrEOM() to form a
!   complete GTS bulletin. An optional IPH may be needed for transmission of
!   files via TCP/IP protocol; this can be generated using GTShdrIPH().
!
! REFERENCES
!   1) WMO (2009). Operational Procedures for the GTS, Attachment II-5,
!      Data Designators T1T2A1A2ii in Abbreviated Headings, Table C3.
!      WMO, Geneva, 4 November 2009.
!   2) Nightingale, S. (2008). Exchange of data with the
!      Met Office FROST message switch.
!      FROST FTP Protocol Document, V2.1, 2 October 2008.
!
! SEE ALSO
!   GTShdrIPH() - to generate a GTS IP Header
!   GTShdrEOM() - to generate a GTS End-of-Message sequence
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!---------------------------------------------------------------------------------
!****
! Modules

  USE GTShdrs, ONLY: GTShdrGAD, &
                     GTSHDR_DEBUG

  IMPLICIT NONE

! Fixed values

  CHARACTER (LEN=1), PARAMETER :: SOH    = ACHAR(1)
  CHARACTER (LEN=3), PARAMETER :: CRCRLF = ACHAR(13)//ACHAR(13)//ACHAR(10)
  INTEGER,           PARAMETER :: narh   = 31  ! Length of ARH

! Argument list parameters

  INTEGER,           INTENT(INOUT) :: CSN         ! Channel sequence no.
  CHARACTER (LEN=*), INTENT(IN)    :: TTAAII      ! T1T2A1A2ii part of ARH
  INTEGER,           INTENT(IN)    :: Lat, Lon    ! Latitude & Longitude (deg)
  CHARACTER (LEN=*), INTENT(IN)    :: CCCC        ! Orininating centre ICAO code
  INTEGER,           INTENT(IN)    :: DateTime(:) ! Date/Time (Y,M,D,h,m,s)
  CHARACTER (LEN=*), INTENT(OUT)   :: cARH        ! Plain character ARH
  INTEGER,           INTENT(OUT)   :: lARH        ! ARH length (bytes)

! Local Variables

  CHARACTER (LEN=1)  :: A2
  CHARACTER (LEN=10) :: cv1, cv2

! Check ARH is big enough - need at least 31 characters

  IF ( LEN(cARH) < narh ) THEN
    cARH    = " "
    lARH    = 0
    WRITE ( *, FMT="(A,I2)" ) "ERROR: [GTShdrARH] String is too small for"//&
                              " ARH. LEN must be at least ",narh

  ELSE

! Start-of-message sequence (SOH/CR/CR/LF/CSN/CR/CR/LF) [10 bytes]
! Increment CSN (cycle in valid range 001-999) before use

    cARH(1:4)  = SOH // CRCRLF
    CSN = MOD ( MAX(CSN,0), 999 ) + 1
    WRITE ( cARH(5:7), FMT="(I3.3)" ) CSN
    cARH(8:10) = CRCRLF

! Data ID (T1T2A1A2ii) [6 bytes]
!  A2 Geographic Area Designator code: A-L for Longitude segments 0-90W,
!  90W-180, 180-90E, 90E-0 and for Latitude bands 90N-30N, 30N-30S, 30S-90S;
!  N,S for Northern/Southern Hemispheres; T for NH 45W-180E; X for Global
!  coverage.

    A2 = GTShdrGAD ( Lat, Lon )
    IF ( GTSHDR_DEBUG ) THEN
      WRITE ( cv1, "(I10)" ) Lat
      WRITE ( cv2, "(I10)" ) Lon
      WRITE ( *, "(A)" ) "DEBUG: [GTShdrARH]"          // &
                         " Lat=" // TRIM(ADJUSTL(cv1)) // &
                         " Lon=" // TRIM(ADJUSTL(cv2)) // &
                         " A2="  // A2
    END IF
    cARH(11:16) = TTAAII(1:3) // A2 // TTAAii(5:6)

! Originating centre (ICAO code CCCC) [6 bytes incl. spaces]

    cARH(17:22) = " " // CCCC // " "

! Day of month & time [6 bytes]

    WRITE ( cARH(23:28), FMT="(3I2.2)" ) DateTime(3:5)

! Data introducer sequence (CR/CR/LF) [3 bytes]

    cARH(29:31) = CRCRLF

    lARH = narh
  END IF

END SUBROUTINE GTShdrARH_sca
!---------------------------------------------------------------------------------

SUBROUTINE GTShdrARH_arr ( CSN,      & ! (inout)
                           TTAAII,   & ! (in)
                           Lats,     & ! (in)
                           Lons,     & ! (in)
                           CCCC,     & ! (in)
                           DateTime, & ! (in)
                           cARH,     & ! (out)
                           lARH )      ! (out)

! GTShdrARH_arr provides a lat/lon REAL array interface to GTShdrARH_sca
! by analysing the ranges of the input lat/lon pairs and determining an
! appropriate single lat/lon location or flag value(s) before calling
! GTShdrARH_sca

! Modules

  USE GTShdrs, ONLY: GTShdrARH_sca, &
                     GTShdrGAD,     &
                     GTSHDR_FLGVAL, &
                     GTSHDR_DEBUG

  IMPLICIT NONE

! Fixed values

  CHARACTER (LEN=1), PARAMETER :: SOH    = ACHAR(1)
  CHARACTER (LEN=3), PARAMETER :: CRCRLF = ACHAR(13)//ACHAR(13)//ACHAR(10)
  INTEGER,           PARAMETER :: narh   = 31  ! Length of ARH

! Argument list parameters

  INTEGER,           INTENT(INOUT) :: CSN         ! Channel sequence no.
  CHARACTER (LEN=*), INTENT(IN)    :: TTAAII      ! T1T2A1A2ii part of ARH
  REAL,              INTENT(IN)    :: Lats(:), Lons(:)  ! Latitudes & Longitudes (deg)
  CHARACTER (LEN=*), INTENT(IN)    :: CCCC        ! Orininating centre ICAO code
  INTEGER,           INTENT(IN)    :: DateTime(:) ! Date/Time (Y,M,D,h,m,s)
  CHARACTER (LEN=*), INTENT(OUT)   :: cARH        ! Plain character ARH
  INTEGER,           INTENT(OUT)   :: lARH        ! ARH length (bytes)

! Local variables

  REAL     :: MinLat, MaxLat
  REAL     :: MinLon, MaxLon
  INTEGER  :: B1, B2, Q1, Q2, N, Lat, Lon
!
  N = MIN ( SIZE(Lats), SIZE(Lons) ) ! just in case arrays not of equal size

  IF ( N > 1 ) THEN
    MinLat = MINVAL(MIN(MAX(Lats(1:N), -90.0), 90.0))  ! ensure all lats +/-90
    MaxLat = MAXVAL(MIN(MAX(Lats(1:N), -90.0), 90.0))
    MinLon = MINVAL(MOD(Lons(1:N)+360.0, 360.0))       ! ensure all lons 0-360
    MaxLon = MAXVAL(MOD(Lons(1:N)+360.0, 360.0))

    B1 = MAX(NINT(MINLAT+90.0) / 60, 0)   ! Bands (0=SH, 1=TR, 2=NH)
    B2 = MIN(NINT(MAXLAT+90.9) / 60, 2)
    Q1 = MAX(NINT(MINLON)      / 90, 0)   ! Quadrants (0-3 E of Greenwich)
    Q2 = MIN(NINT(MAXLON)      / 90, 3)
    IF ( GTSHDR_DEBUG ) THEN
      Lat = NINT(SUM(Lats)/N)
      Lon = NINT(SUM(Lons)/N)
      WRITE ( *,'(A,2F10.2,3I5,A,I3)' ) 'DEBUG: [GTShdrARH) Lat:',&
                                        MinLat,MaxLat,Lat,B1,B2,  &
                                        '   N =',N
      WRITE ( *,'(A,2F10.2,3I5,A,I3)' ) 'DEBUG: [GTShdrARH) Lon:',&
                                        MinLon,MaxLon,Lon,Q1,Q2,  &
                                        '   A2='//GTShdrGAD(Lat,Lon)
    END IF

    IF ( B1 == B2 .AND. Q1 == Q2 ) THEN   ! All in one Quadrant/Band box (A-L)
      Lat = NINT(MinLat)
      Lon = NINT(MinLon)

    ELSE IF ( MinLat >= 0.0 ) THEN        ! All in Northern Hemisphere (N,T)
      IF ( (MinLon <= 180.0 .OR. MinLon >= 315.0) .AND. &
           (MaxLon <= 180.0 .OR. MaxLon >= 315.0) ) THEN  ! 45W-180E (T)
        Lat = GTSHDR_FLGVAL
        Lon = NINT(MinLon)
      ELSE                                                ! else NH (N)
        Lat = NINT(MinLat)
        Lon = GTSHDR_FLGVAL
      END IF

    ELSE IF ( MaxLat < 0.0 ) THEN         ! All in Southern Hemisphere (S)
        Lat = NINT(MaxLat)
        Lon = GTSHDR_FLGVAL

    ELSE                                  ! Spans Equator (X)
      Lon = GTSHDR_FLGVAL
      Lat = GTSHDR_FLGVAL

    END IF

  ELSE IF ( N == 1 ) THEN                 ! Only one value (A-L)
    Lat = NINT(Lats(1))
    Lon = NINT(Lons(1))

  ELSE                                    ! Bad arrays (X) [should never happen!]
    WRITE ( *, '(A)' ) 'WARNING [GTShdrARH] Bad Lat/Lon arrays'
    Lat = GTSHDR_FLGVAL
    Lon = GTSHDR_FLGVAL

  END IF

  CALL GTShdrARH_sca ( CSN, TTAAII, Lat, Lon, CCCC, DateTime, cARH, lARH )

END SUBROUTINE GTShdrARH_arr
!---------------------------------------------------------------------------------

SUBROUTINE GTShdrEOM ( cEOM, & ! (out)
                       lEOM )  ! (out)

!---------------------------------------------------------------------------------
!****s* BUFR/GTSHDRS/GTShdrEOM *
!
! NAME
!   GTShdrEOM   (gtshdrs.f90)
!
! SYNOPIS
!   Generate WMO GTS end-of-message (EOM) sequence
!
!   USE GTShdrsMod
!   CHARACTER (LEN=4) :: cEOM
!   INTEGER :: lEOM
!   CALL GTShdrEOM ( cEOM, lEOM )
!
! INPUT
!   NONE
!
! OUTPUT
!   cEOM  chr  Character string of at least (4) for plain character EOM
!   lEOM  int  Length of EOM in bytes (normally 4)
!
! DESCRIPTION
!   Generates WMO GTS 4-byte end-of-message sequence 'CR/CR/LF/ETX' used to
!   terminate a WMO/GTS bulletin (such as BUFR or GRIB message with routing
!   headers) according to WMO Format 00 (see Reference 1).
!   The EOM is returned as a plain ASCII string (containing non-printing
!   control chracters).
!   The message should have been prefixed with an ARH (e.g. as generated
!   using GTShdrARH()) and optionally by an IP leader (GTShdrIPH()) if
!   intended for transmission over TCP/IP links.
!
! REFERENCES
!   1) Bridgeman, D. & Little, C. (2004). Exchange of data with the
!      Met Office FROST message switch.
!      FROST FTP Protocol Document, V1.0, 31 October 2004.
!
! SEE ALSO
!   GTShdrARH() - to generate a GTS Abbreviated Routing Header
!   GTShdrIPH() - to generate an IP Header for TCP/IP transmission
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!---------------------------------------------------------------------------------
!****
!
  IMPLICIT NONE

! Fixed values

  INTEGER,              PARAMETER :: neom = 4                       ! Length of EOM
  CHARACTER (LEN=neom), PARAMETER :: aEOM = ACHAR(13)//ACHAR(13)//& ! CR/CR/LF/ETX
                                            ACHAR(10)//ACHAR(3)

! Argument list parameters

  CHARACTER (LEN=*), INTENT(OUT) :: cEOM     ! Plain character EOM
  INTEGER,           INTENT(OUT) :: lEOM     ! Length of EOM (bytes)

! Check ARH is big enough - need at least 4 characters

  IF ( LEN(cEOM) < neom ) THEN
    cEOM    = " "
    lEOM    = 0
    WRITE ( *, FMT="(A,I2)" ) "ERROR: [GTShdrEOM] String is too small for"//&
                              " EOM. LEN must be at least ",neom

  ELSE
    cEOM = aEOM
    lEOM = neom
  END IF

END SUBROUTINE GTShdrEOM
!---------------------------------------------------------------------------

SUBROUTINE GTShdrCSN ( CSNfile, & !(in)
                       CSN,     & !(inout)
                       RorW )     !(in)

!****s* BUFR/GTSHDRS/GTShdrCSN *
!
! NAME
!   GTShdrCSN   (gtshdrs.f90)
!
! SYNOPSIS
!   Read or save a channel sequence number
!
!    CHARACTER (LEN=100) :: csnfile
!    INTEGER :: csn
!    CHARACTER (LEN=1) :: rw = ['R'|'W']
!    CALL GTShdrCSN ( csnfile, csn, rw )
!
! INPUTS
!   CSNfile  chr  Channel sequence number file name
!   CSN      int  Channel sequence number (001-999) (write)
!   RorW     chr  'R' or 'r' for input or 'W' or 'w' for output
!
! OUTPUTS
!   CSN      int  Channel sequence number (001-999) (read)
!
! ERRORS
!   If the CSN file cannot be opened for read, or the CSN is not a valid
!   integer, a default sequence initialisation CSN value will be returned as
!   zero. Otherwise, CSN values read from the file or passed to this routine
!   for write will always be returned in the range 001-999.
!
! CALLS
!   GTShdrUNT
!
! DESCRIPTION
!   Reads (if RorW='R') or writes (if RorW='W') a channel sequence number
!   (which should be in the range 001-999) from/to the given file. If the
!   file name is blank, or the first character of RorW is not 'R' or 'W',
!   nothing happens (RorW is case-insensitive).
!   Warning messages are written to stdout if any I/O error occurs (the value
!   of CSN is unchanged), but otherwise the action is silent unless the
!   debug mode is enabled.
!   The file name may include an absolute path, and should be accessible for
!   read & write. The file should also be in a stable location so as to be
!   available for subsequent program runs - don't use places such as /tmp or
!   /var/tmp, etc.
!   Note that GTShdrARH() will increment CSN (and keep it within the valid
!   range 1-999) before coding into the ARH; this routine merely reads the
!   value from the file and saves the given value back, unchanged (apart
!   from ensuring that the value is within valid range).
!
! SEE ALSO
!   GTShdrARH
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!---------------------------------------------------------------------------
!****
! Modules

  USE GTShdrs, ONLY: GTShdrUNT, &
                     GTSHDR_DEBUG

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN=*), PARAMETER :: fmt1 = "(A,I3.3,A)"
  CHARACTER (LEN=*), PARAMETER :: fmt2 = "(I3.3)"

! Argument list parameters

  CHARACTER (LEN=*), INTENT(IN)    :: CSNfile
  INTEGER,           INTENT(INOUT) :: CSN
  CHARACTER (LEN=*), INTENT(IN)    :: RorW

! Local variables

  CHARACTER (LEN=10) :: number
  INTEGER :: ierr, num
  INTEGER :: CSNunit
  LOGICAL :: exists

  IF ( CSNfile /= " " ) THEN
    CSNunit = GTShdrUNT()

!-------------------------------------------------------------
! 1. Read last used bulletin sequence number...
!-------------------------------------------------------------

    IF ( RorW(1:1) == "R" .OR. &
         RorW(1:1) == "r" ) THEN
      CSN = 0                              ! default return value

      INQUIRE ( FILE=CSNfile, &            ! does the file pre-exist?
                EXIST=exists )

      IF ( exists ) THEN                   ! yes, attempt to open it
        OPEN ( FILE=CSNfile, &
               UNIT=CSNunit, &
             ACTION="READ",  &
             IOSTAT=ierr )

        IF ( ierr == 0 ) THEN              ! if ok, attempt to read CSN
          READ ( UNIT=CSNunit, &
                  FMT="(A)",   &
               IOSTAT=ierr ) number
          READ ( number, &
                  FMT=*, &
               IOSTAT=ierr ) num

          IF ( ierr == 0 ) THEN            ! if ok, check range & return value
            CSN = MOD ( MAX(num,1)-1, 999 ) + 1
            IF ( GTSHDR_DEBUG ) THEN
              WRITE ( *, FMT=fmt1 ) &
                       "DEBUG: [GTShdrCSN] Read CSN ", CSN, " from "// &
                       TRIM(CSNfile)
            END IF

          ELSE
            WRITE ( *, FMT=fmt1 ) &
                       "WARNING: [GTShdrCSN] Invalid last value '"// &
                       TRIM(ADJUSTL(number))//"' - using CSN ", CSN
          END IF

        ELSE
          WRITE ( *, FMT=fmt1 ) "WARNING: [GTShdrCSN] Read file open error "// &
                                 TRIM(CSNfile)//" - using CSN ", CSN
        END IF

      ELSE
        WRITE ( *, FMT=fmt1 ) "WARNING: [GTShdrCSN] File not found "// &
                               TRIM(CSNfile)//" - using CSN ", CSN

      END IF

!-------------------------------------------------------------
! 2. ...or save current bulletin sequence number
!-------------------------------------------------------------

    ELSE IF ( RorW(1:1) == "W" .OR. &
              RorW(1:1) == "w" ) THEN

      OPEN ( FILE=CSNfile, &               ! attempt to open file (create if new)
             UNIT=CSNunit, &
           ACTION="WRITE", &
           IOSTAT=ierr )

      IF ( ierr == 0 ) THEN                ! if ok, check range & save value
        CSN = MOD ( MAX(CSN,1)-1, 999 ) + 1
        WRITE ( UNIT=CSNunit, &
                 FMT=fmt2,    &
              IOSTAT=ierr ) CSN

        IF ( ierr == 0 ) THEN
          IF ( GTSHDR_DEBUG ) THEN
            WRITE ( *, FMT=fmt1 ) &
                       "DEBUG: [GTShdrCSN] Saved CSN ", CSN, " to "// &
                       TRIM(CSNfile)
          ENDIF

        ELSE
          WRITE ( *, FMT=fmt1 ) &
                       "WARNING: [GTShdrCSN] Write error to "// &
                       TRIM(CSNfile)//" - CSN not saved"
        END IF

      ELSE
        WRITE ( *, FMT=fmt1 ) &
                       "WARNING: [GTShdrCSN] File open error for "// &
                       TRIM(CSNfile)//" - CSN not saved"
      END IF

    END IF

    CLOSE ( UNIT=CSNunit, &
          IOSTAT=ierr )

  END IF

END SUBROUTINE GTShdrCSN
!-----------------------------------------------------------------------

FUNCTION GTShdrUNT() RESULT(unit)

!****f* BUFR/GTSHDRS/GTShdrUNT *
!-----------------------------------------------------------------------
!
! NAME
!   GTShdrUNT   (gtshdrs.f90)
!
! SYNOPSIS
!   Obtain a free Fortran unit number
!
!   USE GTShdrs
!   INTEGER :: unit
!   unit = GTShdrUNT()
!
! INPUTS
!   None
!
! OUTPUTS
!   GTShdrUNT   (funcn)  Unit number (10-999 or 0 if no free numbers)
!
! ERRORS
!   If there are no free unit numbers in the range 10-999, then
!   a warning message is written to stdout and 'unit' is returned
!   with a value of zero.
!
! CALLED BY
!   GTShdrCSN
!
! DESCRIPTION
!   This function returns the first available unit number larger than 10, but
!   no larger than 999.
!   If no available unit number is found, a warning message will be written
!   to stdout, and the unit number will be returned as zero. This may be
!   valid for most compilers, so if free, the file open may still work.
!   Other compilers may give a file OPEN() 'invalid unit number' error.
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!-----------------------------------------------------------------------
!****

  IMPLICIT NONE

! Argument list parameters

  INTEGER :: unit

! Local variables

  LOGICAL :: unitopen

  unit     = 9
  unitopen = .TRUE.

  DO WHILE ( unitopen .AND. unit < 1000 )
    unit = unit + 1
    INQUIRE ( UNIT=unit, OPENED=unitopen )
  END DO

  IF ( unit > 999 ) THEN
    WRITE ( *, "(A)" ) "WARNING: [GTShdrUNT] No free units between 10 and 999"
    unit = 0
  ENDIF

END FUNCTION GTShdrUNT

!---------------------------------------------------------------------------------
FUNCTION GTShdrGAD ( Lat, Lon ) RESULT(A2)

!****f* BUFR/GTSHDRS/GTShdrA2code *
!-----------------------------------------------------------------------
!
! NAME
!   GTShdrGAD   (gtshdrs.f90)
!
! SYNOPSIS
!   Return the Geographic Area Designator (A2) code given a lat/lon pair.
!
!   USE GTShdrs
!   CHARACTER :: a2
!   REAL :: Lat, Lon
!   a2 = GTShdrGAD(Lat,Lon)
!
! INPUTS
!   Lat  int   Latitude  value in degrees (+/-90 deg or 999)
!   Lon  int   Longitude value in degrees (0-360 or +/-180 deg or 999)
!
! OUTPUTS
!   GTShdrGAD  (funcn) 'A2' Geographic Area Designator code (A-L,N,S,T,X)
!
! DESCRIPTION
!   This function returns the A2 'Geographical Area Designator' part of the
!   WMO GTS Abbreviated Routing Header (T1T2A1A2ii). The A2 code depends on
!   the given latitude and longitude pair of values, and will be one of:
!   the uppercase letters A,B,C...L denoting the longitude quadrant and
!   North, Tropic or South latitude bands; N for Northern Hemishere; S for
!   Southern Hemisphere; T for Europe/Asia or X for Global. See Reference 1.
!   When both Lat and Lon values are within their normal ranges, the
!   appropriate A-L code is determined; if |Lon|>=999 , N or S is
!   selected from the sign of Lat; if |Lat|>=999 and -45<=Lon<=180,
!   (sector 315-180 deg for Europe/Asia) then T is returned, otherwise X;
!   and if both |Lat|>=999 and |Lon|>=999, X is returned.
!   The N, S, T or X codes are appropriate when the encoded data has a large
!   geographic span, and a single A-L code is not representative of the
!   areal range of the data.
!
! REFERENCES
!   1) WMO (2009). Operational Procedures for the GTS, Attachment II-5,
!      Data Designators T1T2A1A2ii in Abbreviated Headings, Table C3.
!      WMO, Geneva, 4 November 2009.
!
! AUTHOR
!   D. Offiler, Met Office, Exeter, UK.
!
! COPYRIGHT
!   (c) Crown copyright 2013, Met Office. All rights reserved.
!
!   Use, duplication or disclosure of this code is subject to the restrictions
!   as set forth in the contract. If no contract has been raised with this
!   copy of the code, the use, duplication or disclosure of it is strictly
!   prohibited. Permission to do so must first be obtained in writing from
!   the Head of Satellite Applications at the following address:
!      Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!-----------------------------------------------------------------------
!****
! Modules

  USE GTShdrs, ONLY: GTSHDR_FLGVAL

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN=*), PARAMETER :: AtoL   = "ABCDEFGHIJKL"

! Argument list parameters

  INTEGER, INTENT(IN) :: Lat, Lon
  CHARACTER (LEN=1)   :: A2

! Local variables

  INTEGER :: LN, LL

  IF ( ABS(Lon) < GTSHDR_FLGVAL ) THEN
    LN = MOD(Lon+360, 360)                     ! ensure 0-360deg

    IF ( ABS(Lat) < GTSHDR_FLGVAL ) THEN
      LL = 4 - LN / 90                         ! Longitude quadrant
      IF ( Lat < -30 ) THEN                    ! Latitude band offset
        LL = LL + 8
      ELSE IF ( Lat < 30 ) THEN
        LL = LL + 4
      END IF
      A2 = AtoL(LL:LL)
    ELSE
      IF ( LN >= 315 .OR. LN <= 180 ) THEN       ! Europe/Asia
        A2 = "T"
      ELSE                                       ! Global
        A2 = "X"
      END IF
   END IF

  ELSE
    IF ( ABS(LAT) >= GTSHDR_FLGVAL ) THEN      ! Global
      A2 = "X"
    ELSE IF ( Lat >= 0 ) THEN                  ! Hemisphere
      A2 = "N"
    ELSE
      A2 = "S"
    END IF
  END IF

END FUNCTION GTShdrGAD
!---------------------------------------------------------------------------------
