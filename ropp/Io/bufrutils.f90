! $Id: bufrutils.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE BUFRutils

!****m* BUFR/BUFRUTILS/BUFRutils *
!
! NAME
!   BUFRutils    (BUFRutils.f90)
!
! SYNOPIS
!   Module providing BUFR utility routines
!
!   USE BUFRutils
!   CALL ConvertDescriptor ( Descr, FXXYYY,  F,X,Y, in )
!   CALL ConvertMOtoEC     ( nExpDescr, ExpDescr, nObs,
!                            Values1, Values2, pValues, cValues )
!   CALL ExtractMOinfo     ( cBUF, Supp, Sec0, Sec1, Sec2,
!                            Sec3, Sec4, Sec5, nDescr, Descr )
!   Cnt = ReplicationCount ( ExpDescr, Values, KeyDescr )
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
!****

CONTAINS
!---------------------------------------------------------------------

SUBROUTINE ConvertDescriptor ( Descr,  & ! (inout)
                               FXXYYY, & ! (inout)
                               F,X,Y,  & ! (inout)
                               in )      ! (in)

!****s* BUFR/BUFRUTILS/ConvertDescriptor *
!
! NAME
!   ConvertDescriptor    (BUFRutils.f90)
!
! SYNOPIS
!   Convert a BUFR descriptor between 16-bit integer, 6-digit readable
!   & individual element forms
!
!   INTEGER :: descr, fxxyyy, f,x,y, ind
!   CALL convertdescriptor( descr, fxxyyy, f,x,y, ind )
!
! INPUTS
!   Descr    int  descriptor in 16-bit integer (BUFR) form   (ind=1)
!   FXXYYY   int  descriptor in 6-digit readable form fxxyyy (ind=2)
!   F, X, Y  int  individual elements of descriptor          (ind=3)
!
! OUTPUTS
!   Descr    int  descriptor in 16-bit integer (BUFR) form   (ind!=1)
!   FXXYYY   int  descriptor in 6-digit readable form fxxyyy (ind!=2)
!   F, X, Y  int  individual elements of descriptor          (ind!=3
!
! DESCRIPTION
!   Converts a BUFR descriptor from one of:
!    1. a 16-bit integer (as encoded in Section 3 of a BUFR message),
!    2. a 6-digit descriptor suitable for printing (e.g. 301234) or
!    3. individual descriptor elements F,X,Y (e.g. 3,1,234)
!   to the other two forms. Which of these three forms to use as input
!   is specified by the 'ind' argument, 1, 2 or 3; any other value will
!   be silently ignored.
!
! REFERENCES
!   1) Manual on codes. International Codes, Vol 1.2, Part B - Binary Codes.
!      2010 Edition, WMO-306, WMO, Geneva.
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
!****

  IMPLICIT NONE

! Argument list parameters

  INTEGER, INTENT(INOUT) :: Descr     ! 1. 16-bit binary descriptor
  INTEGER, INTENT(INOUT) :: FXXYYY    ! 2. 6-digit readable descriptor
  INTEGER, INTENT(INOUT) :: F, X, Y   ! 3. Individual descriptor elements
  INTEGER, INTENT(IN)    :: in        ! Input indicator (1, 2, or 3)

! Local variables

  INTEGER :: FXX ! intermediate FXX value

  IF ( in == 1 ) THEN                ! Descr -> FXXYYY & F,X,Y
    FXX    = Descr / 256
    F      = MOD ( FXX / 64, 4 )
    X      = MOD ( FXX,     64 )
    Y      = MOD ( Descr,  256 )
    FXXYYY = F * 100000 + X * 1000 + Y

  ELSE IF ( in == 2 ) THEN           ! FXXYYY -> Descr & F,X,Y
    FXX   = FXXYYY / 1000
    F     = MOD ( FXX / 100, 4 )
    X     = MOD ( FXX,     100 )
    Y     = MOD ( FXXYYY, 1000 )
    Descr = F * 16384 + X * 256 + Y

  ELSE IF ( in == 3 ) THEN           ! F,X,Y -> Descr & FXXYYY
    F = MOD ( F,   4 ) ! 2 bits
    X = MOD ( X,  64 ) ! 6 bits
    Y = MOD ( Y, 256 ) ! 8 bits
    Descr  = F *  16384 + X *  256 + Y
    FXXYYY = F * 100000 + X * 1000 + Y

  END IF

END SUBROUTINE ConvertDescriptor
!-----------------------------------------------------------------------

SUBROUTINE ConvertMOtoEC ( nExpDescr, & ! (inout)
                           ExpDescr,  & ! (inout)
                           nObs,      & ! (in)
                           Values1,   & ! (in)
                           Values2,   & ! (out)
                           pValues,   & ! (in)
                           cValues )    ! (out)

!****s* BUFR/BUFRUTILS/ConvertMOtoEC *
!
! NAME
!   ConvertMOtoEC    (BUFRutils.f90)
!
! SYNOPSIS
!  Convert expanded descriptors and decoded values from MetDB to ECMWF format
!
!   CHARACTER (LEN=nn) :: pValues
!   CHRACATER (LEN=80) :: cvalues(nv)
!   INTEGER :: ExpDecr(nd), nExpDescr, nObs
!   REAL    :: Values1(nv)
!   REAL*8  :: Values2(nv)
!   CALL DEBUFR ( ExpDescr, Values1, pValues, nExpDescr, nObs, ... )
!   CALL DEBUFR ( ExpDescr, Values1, ... nExpDescr ...)
!   CALL ConvertMOtoEC ( nExpDescr, ExpDescr, nObs, &
!                        Values1, Values2, pValues, cValues )
!
! INPUTS
!   nExpDescr int  No. of expanded descriptors
!   ExpDescr  int  Array(nd) of expanded descriptors (16-bit format)
!   nObs      int  Number of observations
!   Values1   flt  Array(nv) of decoded numeric values
!   pValues   chr  String of packed dcoded character values
!
! OUTPUTS
!   nDescr    int  Updated no. of expanded descriptors
!   ExpDescr  int  Array(nd) of updated expanded descriptors (fxxyyy format)
!   Values2  dflt  Array(nv) of updated decoded numeric values
!   cvalues   chr  Array(nv) of updated decoded character values
!
! CALLS
!   ConvertDescriptor
!
! DESCRIPTION
!   Converts expanded descriptor list and decoded numeric and character values
!   from MetDB BUFR library DEBUFR() interface, to equivalents compatible
!   with those returned by ECMWF BUFREX() and BUSEL() routines (see Refs.)
!    1) convert MetDB raw 16-bit descriptor values to fxxyyy format;
!    For numeric data:
!    2) extract the replication count from any replication operator
!       descriptors (1xxyyy) and insert count into the Values array; replace
!       the replication operator descriptor with a replication factor
!       descriptor (03100y);
!    3) remove any change operator descriptors (2xxyyy)
!    4) Change MetDB missing data values (RMDFV) to ECMWF equivalent (RVIND);
!    For character data:
!    5) extract pointers from value array and unpack flat string to
!       character array
!   Hence for each replication, Values will lengthen by one and for each
!   change, ExpDescr will shorten by one until there is 1-to-1 correspondance
!   between ExpDescr and Values/cValues arrays.
!   Note that for correctly expanded descriptors, there will be no sequence
!   descriptors (3xxyyy), so these are not checked for.
!
! REFERENCES
!   1. Met Office (2012). Decoding and Encoding BUFR messages.
!      MetDB Technote 1, Rev.1, 29/10/2012 [dmtn1.html].
!   2. Dragosavac, Milan (2009). BUFR User's Guide.
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
!****

  IMPLICIT NONE

! Fixed parameters

  INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(P = 13, R = 307)

  INTEGER,  PARAMETER :: NMDFV = -9999999     ! Integer missing data flag value (MetDB)
  INTEGER,  PARAMETER :: NVIND = 2147483647   ! Integer missing data flag value (ECMWF)
  REAL,     PARAMETER :: RMDFV = -9999999.0   ! Real    missing data flag value (MetDB)
  REAL(dp), PARAMETER :: RVIND = 1.7E38_dp    ! Real    missing data flag value (ECMWF)

! Argument list parameters

  INTEGER,           INTENT(INOUT) :: nExpDescr   ! No. of expanded desriptors
  INTEGER,           INTENT(INOUT) :: ExpDescr(:) ! Expanded descriptor list
  INTEGER,           INTENT(IN)    :: nObs        ! No. of observations
  REAL,              INTENT(IN)    :: Values1(:)  ! Decoded   numeric values
  REAL(dp),          INTENT(OUT)   :: Values2(:)  ! Converted numeric values
  CHARACTER (LEN=*), INTENT(IN)    :: pValues     ! Decoded packed character data
  CHARACTER (LEN=*), INTENT(OUT)   :: cValues(:)  ! Unpacked character data

! Local variables

  INTEGER :: Descr                                ! decsriptor (16-bit format)
  INTEGER :: FXXYYY                               ! descriptor (fxxyy format)
  INTEGER :: F, X, Y                              ! descriptor in separate F,X,Y parts
  INTEGER :: RepFac                               ! replication factor
  INTEGER :: iVal                                 ! temporay values
  INTEGER :: pt, nc                               ! character index & no. of chrs
  INTEGER :: id, io, iv1, iv2                     ! loop counters & array indices
  INTEGER :: L                                    ! array length

!-------------------------------------------------------------
! 0. Initialize
!-------------------------------------------------------------

  Values2(:) = 0.0_dp
  cValues(:) = " "
  L = SIZE(ExpDescr)

  id  = 1
  iv1 = 0
  iv2 = 0

!-------------------------------------------------------------
! 1. Process descriptors
!-------------------------------------------------------------

  DO WHILE ( id <= nExpDescr )

!-------------------------------------------------------------
! 1.1 Convert descriptor from 16-bit to fxxyyy (& f,x,y) format
!-------------------------------------------------------------

    Descr = ExpDescr(id)
    CALL ConvertDescriptor ( Descr, FXXYYY, F,X,Y , 1 )
    ExpDescr(id) = FXXYYY

    IF ( Descr < 65536 ) THEN

!-------------------------------------------------------------
! 1.2 Normal Table B descriptor: copy values
!-------------------------------------------------------------

      IF ( F == 0 ) THEN
        DO io = 1, nObs
          iv1 = iv1 + 1
          iv2 = iv2 + 1
          Values2(iv2) = Values1(iv1)
        END DO

!-------------------------------------------------------------
! 1.2 If a replication operator descriptor (F=1); extract
!     count value (X*256+Y) from it and and insert this into
!     the Values array; replace operator descriptor  with an
!     apropriate replication factor descriptor (not foolproof
!     restoration of original)
!-------------------------------------------------------------

      ELSE IF ( F == 1 ) THEN
        RepFac = X * 256 + Y
        DO io = 1, nObs
          iv2 = iv2 + 1
          Values2(iv2) = RepFac
        END DO
        IF ( RepFac < 256 ) THEN
          ExpDescr(id) = 031001
        ELSE
          ExpDescr(id) = 031002
        END IF

!-------------------------------------------------------------
! 1.3 If a change operator descriptor (F=2), remove it
!-------------------------------------------------------------

      ELSE IF ( F == 2 ) THEN
        ExpDescr(id:L-1) = ExpDescr(id+1:L)
        ExpDescr(L)      = 0
        nExpDescr        = nExpDescr - 1
        id               = id        - 1

      END IF

!-------------------------------------------------------------
! 1.4 Unpack character values
!-------------------------------------------------------------

    ELSE
      DO io = 1, nObs
        iv1  = iv1 + 1
        iv2  = iv2 + 1
        iVal = NINT(Values1(iv1))
        Values2(iv2) = 0.0_dp
        pt   = MOD ( iVal, 65536 )
        nc   = iVal / 65536
        IF ( pt > 0 .AND. nc > 0 ) cValues(iv2) = pValues(pt:pt+nc-1)
      END DO

    ENDIF

    id = id + 1
  END DO

!-------------------------------------------------------------
! 2. Change MetDB missing data flag value to ECMWF equivalent
!-------------------------------------------------------------

  WHERE ( NINT(Values2) == NMDFV ) Values2 = RVIND

END SUBROUTINE ConvertMOtoEC
!-----------------------------------------------------------------

SUBROUTINE ExtractMOinfo ( cBUF,             & ! (in)
                           Supp,             & ! (out)
                           Sec0, Sec1, Sec2, & ! (out)
                           Sec3, Sec4, Sec5, & ! (out)
                           nDescr, Descr )     ! (out)

!****s* BUFR/BUFRUTILS/ExtractMOinfo *
!
! NAME
!   ExtractMOinfo   (BUFRutils.f90)
!
! SYNOPSIS
!   Extract information from MetDB BUFR Sections 0-5 (except encoded
!   Section 4 data)
!
!    CHARACTER (LEN=nnn) :: cBUF
!    INTEGER :: Supp(9), Sec10(3), Sec1(40),
!    INTEGER :: Sec2(4096), Sec3(4), Sec4(2), Sec5(2)
!    INTEGER :: nDescr, Descr(nd)
!    CALL ExtractMOinfo ( cBUF, Supp, Sec0, Sec1, Sec2, &
!                         Sec3, Sec4, Sec5, nDescr, Descr )
!
! INPUTS
!   cBUF    chr  BUFR message bit string (MetDB I/O interface)
!
! OUTPUTS
!   Supp    int  Array(9)    for Supplimentary info
!   Sec0    int  Array(3)    for Section 0 info
!   Sec1    int  Array(40)   for Section 1 info
!   Sec2    int  Array(4096) for Section 2 info
!   Sec3    int  Array(4)    for Section 3 info
!   Sec4    int  Array(2)    for Section 4 info
!   Sec5    int  Array(2)    for Section 5 info
!   nDecsr  int  No. of descriptors in Section 3
!   Descr   int  Array(nd)   of  descriptors (fxxyyy format)
!
! CALLS
!   ConvertDescriptor
!
! DESCRIPTION
!   Extracts BUFR Section 0-5 information (excep actual encoded Section 4
!   data) from a BUFR message using the MetDB character string I/O interface.
!   This routine, together with MetDB DEBUFR() provides the unexpanded and
!   expanded descriptor lists compatible with ECMWF BUSEL() and Section
!   header information returned by BUS123().
!
! REFERENCES
!   1. Manual on Codes: International Codes, Part B & Part C.
!      WMO-No. 306, World Meteorological Organisation, Geneva.
!      http://www.wmo.int/pages/prog/www/WMOCodes/VolumeI2.html
!   2. Met Office (2012). Decoding and Encoding BUFR messages.
!      MetDB Technote 1, Rev.1, 29/10/2012 [dmtn1.html].
!   3. Dragosavac, Milan (2009). BUFR User's Guide.
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
!****

  IMPLICIT NONE

! Argument list parameters

  CHARACTER (LEN=*), INTENT(IN)  :: cBUF       ! BUFR message bit string
  INTEGER,           INTENT(OUT) :: Supp(:)    ! Array for Supplimentary info
  INTEGER,           INTENT(OUT) :: Sec0(:)    ! Array for Section 0 info
  INTEGER,           INTENT(OUT) :: Sec1(:)    ! Array for Section 1 info
  INTEGER,           INTENT(OUT) :: Sec2(:)    ! Array for Section 2 info
  INTEGER,           INTENT(OUT) :: Sec3(:)    ! Array for Section 3 info
  INTEGER,           INTENT(OUT) :: Sec4(:)    ! Array for Section 4 info
  INTEGER,           INTENT(OUT) :: Sec5(:)    ! Array for Section 5 info
  INTEGER,           INTENT(OUT) :: nDescr     ! No. of descripors
  INTEGER,           INTENT(OUT) :: Descr(:)   ! Descriptor list

! Local parameters

  CHARACTER (LEN=4) :: Delim        ! message start/end delimiter

  INTEGER :: len0, len1, len2, len3 ! Lengths (bytes) of Sections 0, 1, 2 & 3
  INTEGER :: len4, len5, lenm       ! Lengths (bytes) of Sections 4, 5 & total
  INTEGER :: Edition                ! BUFR Edition in Section 0
  INTEGER :: Yr,Mo,Dy,Hr,Mn,Sc      ! Date & time in Section 1
  INTEGER :: DataCat                ! Data ategory in Section 1
  INTEGER :: IntnlSubCat            ! International sub-category in S1 (Ed.4)
  INTEGER :: LocalSubCat            ! Local sub-category in Section 1
  INTEGER :: is2                    ! Is Section 2 present? flag in Section 1
  INTEGER :: BufrMasterTable        ! Master table ID in Section 1
  INTEGER :: OCentre                ! Orig. centre code in Section 1
  INTEGER :: OSubCentre             ! Orig. sub-centre code in Section 1
  INTEGER :: UpdtSeqNo              ! Update sequence no. in Section 1
  INTEGER :: MasterTableVer         ! Master table version in Section 1
  INTEGER :: LocalTableVer          ! Local table version in Section 1
  INTEGER :: nObs                   ! No. of observations in Section 3
  INTEGER :: idata                  ! Observed data & compression flags in S3

  INTEGER :: Des                    ! Descriptor in raw 16-bit format
  INTEGER :: F, X, Y                ! Components of descriptor
  INTEGER :: K, L                   ! Pointers to current byte pos. in cBUF
! INTEGER :: Byte1                  ! Pointer to 'BUFR' in cBUF ! Commented at 20 July, 2016
  INTEGER :: i                      ! Loop counters

!-------------------------------------------------------------
! 0. Initialise
!-------------------------------------------------------------

  Supp(:)  = 0
  Sec0(:)  = 0
  Sec1(:)  = 0
  Sec2(:)  = 0
  Sec3(:)  = 0
  Sec4(:)  = 0
  Sec5(:)  = 0
  nDescr   = 0
  Descr(:) = 0

!-------------------------------------------------------------
! 1. Supplimentary information on Section array sizes
!-------------------------------------------------------------

  Supp(1) = SIZE(Sec1)
  Supp(2) = SIZE(Sec2)
  Supp(3) = SIZE(Sec3)
  Supp(4) = SIZE(Sec4)
  Supp(9) = SIZE(Sec0)

!---------------------------------------------------------------
! 2. Section 0 ('BUFR', Edition)
!---------------------------------------------------------------

  len0 = 8  ! fixed since Ed.2

! Check for start of message delimiter

  Delim = "BUFR"
  L = INDEX ( cBUF, Delim )
  IF ( L == 0 ) RETURN
! Byte1 = L ! Commented at 20 July, 2016

! Total BUFR message length

  lenm = ICHAR ( cBUF(L+4:L+4) ) * 65536 &
       + ICHAR ( cBUF(L+5:L+5) ) * 256   &
       + ICHAR ( cBUF(L+6:L+6) )

! BUFR Edition number

  Edition = ICHAR ( cBUF(L+7:L+7) )

  Supp(8) = lenm
  Sec0(1) = len0
  Sec0(2) = lenm
  Sec0(3) = Edition

  L = L + len0 - 1

!---------------------------------------------------------------
! 3. Section 1 (identification)
!---------------------------------------------------------------

  len1 = ICHAR ( cBUF(L+1:L+1) ) * 65536 &
       + ICHAR ( cBUF(L+2:L+2) ) * 256   &
       + ICHAR ( cBUF(L+3:L+3) )

  BufrMasterTable = ICHAR ( cBUF(L+4:L+4) )

  IF ( Edition == 3 ) THEN

! Orginating Centre & Sub-Centre and Update Sequence Number

    OSubCentre    = ICHAR ( cBUF(L+5:L+5) )
    OCentre       = ICHAR ( cBUF(L+6:L+6) )
    UpdtSeqNo     = ICHAR ( cBUF(L+7:L+7) )

! Data Type

    is2            = ICHAR ( cBUF(L+8:L+8) )
    DataCat        = ICHAR ( cBUF(L+9:L+9) )
    IntnlSubCat    = 255
    LocalSubCat    = ICHAR ( cBUF(L+10:L+10) )
    MasterTableVer = ICHAR ( cBUF(L+11:L+11) )
    LocalTableVer  = ICHAR ( cBUF(L+12:L+12) )

! Date & Time (pivot 2-digit century year to CCYY)

    Yr = ICHAR ( cBUF(L+13:L+13) )
    IF ( Yr < 50 ) THEN
      Yr = Yr + 2000
    ELSE IF ( Yr < 100 ) THEN
      Yr = Yr + 1900
    ELSE IF ( Yr == 100 ) THEN
      Yr = 2000
    END IF
    Mo = ICHAR ( cBUF(L+14:L+14) )
    IF ( Mo < 0 .OR. Mo > 12 ) Mo = 0
    Dy = ICHAR ( cBUF(L+15:L+15) )
    Hr = ICHAR ( cBUF(L+16:L+16) )
    Mn = ICHAR ( cBUF(L+17:L+17) )
    sc = 0

  ELSE   ! Ed.4 onwards

! Orginating Centre & Sub-Centre and Update Sequence Number

    OCentre       = ICHAR ( cBUF(L+5:L+5) ) * 256 &
                  + ICHAR ( cBUF(L+6:L+6) )
    OSubCentre    = ICHAR ( cBUF(L+7:L+7) ) * 256 &
                  + ICHAR ( cBUF(L+8:L+8) )
    UpdtSeqNo     = ICHAR ( cBUF(L+9:L+9) )


! Data Type

    is2            = ICHAR ( cBUF(L+10:L+10) )
    DataCat        = ICHAR ( cBUF(L+11:L+11) )
    IntnlSubCat    = ICHAR ( cBUF(L+12:L+12) )
    LocalSubCat    = ICHAR ( cBUF(L+13:L+13) )
    MasterTableVer = ICHAR ( cBUF(L+14:L+14) )
    LocalTableVer  = ICHAR ( cBUF(L+15:L+15) )

! Date & Time

    Yr = ICHAR ( cBUF(L+16:L+16) ) * 256 &
       + ICHAR ( cBUF(L+17:L+17) )
    Mo = ICHAR ( cBUF(L+18:L+18) )
    IF ( Mo < 0 .OR. Mo > 12 ) Mo = 0
    Dy = ICHAR ( cBUF(L+19:L+19) )
    Hr = ICHAR ( cBUF(L+20:L+20) )
    Mn = ICHAR ( cBUF(L+21:L+21) )
    Sc = ICHAR ( cBUF(L+22:L+22) )

  END IF

  Sec1(1)  = len1
  Sec1(2)  = Edition
  Sec1(3)  = OCentre
  Sec1(4)  = UpdtSeqNo
  Sec1(5)  = is2
  Sec1(6)  = DataCat
  Sec1(7)  = LocalSubCat
  Sec1(8)  = LocalTableVer
  Sec1(9)  = Yr
  Sec1(10) = Mo
  Sec1(11) = Dy
  Sec1(12) = Hr
  Sec1(13) = Mn
  Sec1(14) = BufrMasterTable
  Sec1(15) = MasterTableVer
  Sec1(16) = OSubCentre
  Sec1(17) = IntnlSubCat
  Sec1(18) = Sc

  L = L + len1

!---------------------------------------------------------------
! 4. Optional Section 2 (local data)
!---------------------------------------------------------------

  IF ( BTEST ( is2, 7 ) ) THEN
    len2 = ICHAR ( cBUF(L+1:L+1) ) * 65536 &
         + ICHAR ( cBUF(L+2:L+2) ) * 256   &
         + ICHAR ( cBUF(L+3:L+3) )
  ELSE
    len2 = 0
  END IF

  Sec2(1) = len2

  L = L + len2

!---------------------------------------------------------------
! 5. Section 3 (Descriptors)
!---------------------------------------------------------------

  len3 = ICHAR ( cBUF(L+1:L+1) ) * 65536 &
       + ICHAR ( cBUF(L+2:L+2) ) * 256   &
       + ICHAR ( cBUF(L+3:L+3) )

! Number of obs. in Section 4

  nObs = ICHAR ( cBUF(L+5:L+5) ) * 256  &
       + ICHAR ( cBUF(L+6:L+6) )

! Observed data & compression flags

  idata = ICHAR ( cBUF(L+7:L+7) )

! No. of descriptors & list

  nDescr = ( len3 - 7 ) / 2
  K = L + 7
  DO i = 1, nDescr
    Des = ICHAR ( cBUF(K+1:K+1) ) * 256 &
        + ICHAR ( cBUF(K+2:K+2) )
    CALL ConvertDescriptor ( Des, Descr(i), F,X,Y, 1 )
    K = K + 2
  END DO

  Sec3(1) = len3
  Sec3(3) = nObs
  Sec3(4) = idata

  L = L + len3

!---------------------------------------------------------------
! 6. Section 4 (data)
!---------------------------------------------------------------

  len4 = ICHAR ( cBUF(L+1:L+1) ) * 65536 &
       + ICHAR ( cBUF(L+2:L+2) ) * 256   &
       + ICHAR ( cBUF(L+3:L+3) )

  Sec4(1) = len4

  L = L + len4

!---------------------------------------------------------------
! 6. Section 5 ('7777')
!---------------------------------------------------------------

! Check for end-of-message delimiter

  Delim = "7777"
  len5  = 4
  IF ( Edition >= 3 ) THEN
    IF ( cBUF(L+1:L+4) /= Delim ) len5 = 0
  ELSE
    L = INDEX ( cBUF, Delim )
    IF ( L == 0 ) len5 = 0
  END IF

  Sec5(1) = len5
  L = L + len5
  Sec5(2) = L

END SUBROUTINE ExtractMOinfo
!-----------------------------------------------------------------------

INTEGER FUNCTION ReplicationCount ( ExpDescr,  & ! (in)
                                    Values,    & ! (in)
                                    KeyDescr )   ! (in)

!****f* BUFR/BUFRUTILS/ReplicationCount *
!
! NAME
!   ReplicationCount
!
! SYNOPSIS
!   Extract a replication count from list of expanded descriptors
!   and decoded values for a given key descriptor
!
!   INTEGER :: ExpDecr(nd), nDescr, KeyDescr
!   REAL*8  :: Values(nv)
!   CALL ReplicationCount ( ExpDescr, Vaues, KeyDescr )
!
! INPUTS
!   ExpDescr  int  Array(nDescr) of expanded descriptors (fxxyyy format)
!   Values   dflt  Array(nVals)  of decoded values corresponding 1-to-1
!                  with decriptors in ExpDescr
!   KeyDescr  int  Key descriptor (fxxyyy format)
!
! OUTPUTS
!   ReplicationCount  int  Replication count value (fn return)
!
! DESCRIPTION
!   Get a replication count by searching for a Replication Count Factor
!   descriptor in the expanded descriptor list which is followed by the key
!   descriptor, and extract the replication count from the associated element
!   of the decoded Values array.
!
! REFERENCES
!   1. Manual on Codes: International Codes, Part B & Part C.
!      WMO-No. 306, World Meteorological Organisation, Geneva.
!      http://www.wmo.int/pages/prog/www/WMOCodes/VolumeI2.html
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
!****

  IMPLICIT NONE

! Fixed parameters

  INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(P = 13, R = 307)

! Argument list parameters

  INTEGER,  INTENT(IN) :: ExpDescr(:)  ! expanded descriptor list (fxxyyy format)
  REAL(dp), INTENT(IN) :: Values(:)    ! decoded values
  INTEGER,  INTENT(IN) :: KeyDescr     ! key descriptor (fxxyyy format)

! Local variables

  INTEGER :: i                         ! loop counter

!-------------------------------------------------------------
! 0. Initialise
!-------------------------------------------------------------

  ReplicationCount = 0        ! default count

!-------------------------------------------------------------
! 1. Search for a replication factor descriptor (031001 or
!    031002) followed by the specified key descriptor.
!    If found, extract count from Values array.
!-------------------------------------------------------------

  DO i = 1, SIZE(ExpDescr)-1
    IF ( ExpDescr(i) == 0 ) THEN
      EXIT
    ELSE IF ( ( ExpDescr(i)   == 031001   .OR.  &
                ExpDescr(i)   == 031002 ) .AND. &
                ExpDescr(i+1) == KeyDescr ) THEN
      ReplicationCount = NINT(Values(i))
      EXIT
    END IF
  END DO

END FUNCTION ReplicationCount
!-----------------------------------------------------------------------

END MODULE BUFRutils
