!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CDA VERSION 2 COMPREHENSIVE UNPACK SUBROUTINES
!   INCLUDED MEMBERS
!    CDAUPGD
!    P4UNPK (FOR VER.2)
!    CDAUNPK
!    TWOUNPK200 (VER.2)
!      UNPSCL200 (VER.2)
!    TWOUNPK (VER.2.01)
!      UNPSCL (VER.2.01)
!    FOURUNPK
!    BETAUNPK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CDAUPGD (IC1,IC2,IC3,IC4,
     \                    NHT,HGT,NEL,MEL,
     \                    NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                    JGQC,IUSE,NOQC,JEQC,KEQC,
     \                    OBRP,OBQC,RELY,OBER,DVAL,
     \                    GUES,GSER,GSDH,GSDT,ANAL,
     \                    IVNM,IVUT,THIK,IVQC,RPHT,IRUT,
     \                    IQHE,
     \                    MISS,NXHT,NXEL,NXIH,NXIE,NXQC)

!---------------------------------------------------------------------
!     CDA VER.2 UNPACKING 'TUNNEL' SUBROUTINE
!     MOVE CDA VER.1 VARIABLES TO THOSE OF VER.2
!                                             2000.12.30 K.ONOGI
!---------------------------------------------------------------------

      PARAMETER (NXCOD=327)
      PARAMETER (NXTWC=62)

      INTEGER(2), DIMENSION (*) :: IC1,IC2,IC3,IC4

! --- CDAV2 ARRAY
      REAL   (4), DIMENSION (NXHT) :: HGT
      INTEGER(4), DIMENSION (NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION (NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION (NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION (NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION (NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION (NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION (NXQC,NXEL,NXHT) :: JEQC,KEQC

! --- CDAV1 ARRAY (EXCEPT COMMON ARRAYS WITH CDAV2)
      REAL   (4), DIMENSION (NXEL,NXHT) :: OBRP,OBQC,RELY,OBER,DVAL,
     \                                     GUES,GSER,GSDH,GSDT,ANAL
      REAL   (4), DIMENSION (NXHT) :: THIK,RPHT
      INTEGER(4), DIMENSION (NXHT) :: IVNM,IVUT,IVQC,IRUT


! === CLEAR NUMBERS
      NHT  = 0
      NEL  = 0
      NHCD = 0
      NECD = 0
      NOQC = 0


! === CDA VERSION 2 UNPACKING SUBROUTINE

      CALL CDAUNPK (IC1,IC2,IC3,IC4,
     \              NHT,HGT,NEL,MEL,
     \              NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \              JGQC,IUSE,NOQC,JEQC,KEQC,
     \              MISS,NXHT,NXEL,NXIH,NXIE,NXQC)


! === CLEAR CDA VERSION 1 VARIABLES

      DO NH=1,NHT
        IVNM(NH) = MISS
        IVUT(NH) = MISS
        THIK(NH) = MISS
        IVQC(NH) = MISS
        RPHT(NH) = MISS
        IRUT(NH) = MISS
        DO NE=1,NEL(NH)
          OBRP(NE,NH) = MISS
          OBQC(NE,NH) = MISS
          RELY(NE,NH) = MISS
          OBER(NE,NH) = MISS
          DVAL(NE,NH) = MISS
          GUES(NE,NH) = MISS
          GSER(NE,NH) = MISS
          GSDH(NE,NH) = MISS
          GSDT(NE,NH) = MISS
          ANAL(NE,NH) = MISS
        END DO
      END DO

! === TRANSFER CDA VER.2 VARIABLES TO CDA VER.1 VARIABLES

      IF (IQHE.EQ.0) THEN

!   --- HEIGHT INFORMATION

        DO NH=1,NHT

          DO NV=1,NHCD(NH)
            IF (KHCD(NV,NH).EQ.100) THEN
              IVNM(NH) = NINT(AHCD(NV,NH))
            ELSE IF (KHCD(NV,NH).EQ.200) THEN
              IVUT(NH) = NINT(AHCD(NV,NH))
            ELSE IF (KHCD(NV,NH).EQ.300) THEN
              THIK(NH) = AHCD(NV,NH)
            ELSE IF (KHCD(NV,NH).EQ.400) THEN
              IVQC(NH) = NINT(AHCD(NV,NH))
            ELSE IF (KHCD(NV,NH).EQ.500) THEN
              RPHT(NH) = AHCD(NV,NH)
            ELSE IF (KHCD(NV,NH).EQ.600) THEN
              IRUT(NH) = NINT(AHCD(NV,NH))
            END IF
          END DO

        END DO

!   --- ELEMENT INFORMATION

        DO NH=1,NHT
        DO NE=1,NEL(NH)

          DO NI=1,NECD(NE,NH)
            IF (KECD(NI,NE,NH).EQ.100) THEN
              OBRP(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.200) THEN
              OBQC(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.300) THEN
              RELY(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.400) THEN
              OBER(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.500) THEN
              DVAL(NE,NH) = AECD(NI,NE,NH)
!           KECD = 606 EXIST BEFORE 2004.07
            ELSE IF (KECD(NI,NE,NH)/100.EQ.6) THEN
              GUES(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.700) THEN
              GSER(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.810) THEN
              GSDH(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.910) THEN
              GSDT(NE,NH) = AECD(NI,NE,NH)
            ELSE IF (KECD(NI,NE,NH).EQ.1000) THEN
              ANAL(NE,NH) = AECD(NI,NE,NH)
            END IF
          END DO
!         GUES IS NOT SAVED IF OBQC & DVAL EXIST
          IF (NINT(OBQC(NE,NH)).NE.MISS .AND.
     \        NINT(DVAL(NE,NH)).NE.MISS .AND.
     \        NINT(GUES(NE,NH)).EQ.MISS)
     \      GUES(NE,NH) = OBQC(NE,NH) - DVAL(NE,NH)

        END DO
        END DO

      ELSE

!   --- HEIGHT INFORMATION

        DO NH=1,NHT

          DO NV=1,NHCD(NH)
            IF (KHCD(NV,NH).EQ.100) THEN
              IVNM(NH) = NINT(AHCD(NV,NH))
              KHCD(NV,NH) = MISS
            ELSE IF (KHCD(NV,NH).EQ.200) THEN
              IVUT(NH) = NINT(AHCD(NV,NH))
              KHCD(NV,NH) = MISS
            ELSE IF (KHCD(NV,NH).EQ.300) THEN
              THIK(NH) = AHCD(NV,NH)
              KHCD(NV,NH) = MISS
            ELSE IF (KHCD(NV,NH).EQ.400) THEN
              IVQC(NH) = NINT(AHCD(NV,NH))
              KHCD(NV,NH) = MISS
            ELSE IF (KHCD(NV,NH).EQ.500) THEN
              RPHT(NH) = AHCD(NV,NH)
              KHCD(NV,NH) = MISS
            ELSE IF (KHCD(NV,NH).EQ.600) THEN
              IRUT(NH) = NINT(AHCD(NV,NH))
              KHCD(NV,NH) = MISS
            END IF
          END DO

          NVZ = 0
          DO NV=1,NHCD(NH)
            IF (KHCD(NV,NH).NE.MISS) THEN
              NVZ = NVZ+1
              KHCD(NVZ,NH) = KHCD(NV,NH)
              AHCD(NVZ,NH) = AHCD(NV,NH)
            END IF
          END DO
          NHCD(NH) = NVZ

        END DO

!   --- ELEMENT INFORMATION

        DO NH=1,NHT
        DO NE=1,NEL(NH)

          DO NI=1,NECD(NE,NH)
            IF (KECD(NI,NE,NH).EQ.100) THEN
              OBRP(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.200) THEN
              OBQC(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.300) THEN
              RELY(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.400) THEN
              OBER(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.500) THEN
              DVAL(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
!           KECD = 606 EXIST BEFORE 2004.07
            ELSE IF (KECD(NI,NE,NH)/100.EQ.6) THEN
              GUES(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.700) THEN
              GSER(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.810) THEN
              GSDH(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.910) THEN
              GSDT(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            ELSE IF (KECD(NI,NE,NH).EQ.1000) THEN
              ANAL(NE,NH) = AECD(NI,NE,NH)
              KECD(NI,NE,NH) = MISS
            END IF
          END DO
!         GUES IS NOT SAVED IF OBQC & DVAL EXIST
          IF (NINT(OBQC(NE,NH)).NE.MISS .AND.
     \        NINT(DVAL(NE,NH)).NE.MISS .AND.
     \        NINT(GUES(NE,NH)).EQ.MISS)
     \      GUES(NE,NH) = OBQC(NE,NH) - DVAL(NE,NH)

          NIZ = 0
          DO NI=1,NECD(NE,NH)
            IF (KECD(NI,NE,NH).NE.MISS) THEN
              NIZ = NIZ+1
              KECD(NIZ,NE,NH) = KECD(NI,NE,NH)
              AECD(NIZ,NE,NH) = AECD(NI,NE,NH)
            END IF
          END DO
          NECD(NE,NH) = NIZ

        END DO
        END DO

      END IF

      RETURN
      END SUBROUTINE CDAUPGD
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE P4UNPK (NHT,HGT,NEL,MEL,
     \                   JGQC,IUSE,
     \                   OBRP,OBQC,RELY,OBER,DVAL,
     \                   GUES,GSER,GSDH,GSDT,ANAL,
     \                   NOQC,JEQC,KEQC,
     \                   IVNM,IVUT,IVSC,THIK,IVQC,RPHT,IRUT,IRSC,
     \                   IC4,IC1,MISS,
     \                   NXHT,NXEL,NXQC,MOPT)

!---------------------------------------------------------------------
!     CDA VER.2 UNPACKING WITH THE SAME SUBTRACTS OF P4UNPK
!                                             2000.12.30 K.ONOGI
!---------------------------------------------------------------------

      PARAMETER (NXIH  =  20)
      PARAMETER (NXIE  =  40)
      PARAMETER (NMXEL =  60)
      PARAMETER (NMXHT = 250)
      PARAMETER (LXD = 50)

      INTEGER(2), DIMENSION (*) :: IC1,IC4
      INTEGER(2), DIMENSION (LXD) :: JC2,JC3 ! DUMMY
      LOGICAL LX
      DATA LX/.TRUE./

! --- CDAV2 ARRAY
      REAL   (4), DIMENSION (NXHT) :: HGT
      INTEGER(4), DIMENSION (NXHT) :: NEL
      INTEGER(4), DIMENSION (NMXHT) :: NHCD
      INTEGER(4), DIMENSION (NXIH,NMXHT) :: KHCD
      REAL   (4), DIMENSION (NXIH,NMXHT) :: AHCD
      INTEGER(4), DIMENSION (NXEL,NXHT) :: MEL,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION (NXEL,NMXHT) :: NECD
      INTEGER(4), DIMENSION (NXIE,NMXEL,NMXHT) :: KECD
      REAL   (4), DIMENSION (NXIE,NMXEL,NMXHT) :: AECD
      INTEGER(4), DIMENSION (NXQC,NXEL,NXHT) :: JEQC,KEQC

! --- CDAV1 ARRAY (EXCEPT COMMON ARRAYS WITH CDAV2)
      REAL   (4), DIMENSION (NXEL,NXHT) :: OBRP,OBQC,RELY,OBER,DVAL,
     \                                     GUES,GSER,GSDH,GSDT,ANAL
      REAL   (4), DIMENSION (NXHT) :: THIK,RPHT
      INTEGER(4), DIMENSION (NXHT) :: IVNM,IVUT,IVSC,IVQC,IRUT,IRSC

      IF (LX) THEN
        JC2 = MISS  ! DUMMY
        JC3 = MISS  ! DUMMY

!   --- ARRAY SIZE CHECK
        IF (NMXHT.LT.NXHT .OR. NMXEL.NE.NXEL) THEN
          WRITE(*,*) 'DIMENSION SIZES INCONSISTENCY.'
          WRITE(*,*) 'SUBROUTINE P4PACK (FOR CDA VER.2)'
          WRITE(*,*) 'NMXHT,NXHT=',NMXHT,NXHT
          WRITE(*,*) 'NMXEL,NXEL=',NMXEL,NXEL
          STOP 180
        END IF

        LX = .FALSE.
      END IF


! === CDA VERSION CHECK
      IF (IC1(13)/100 == 2.AND.
     \   (IC1(14)==2.OR.IC1(14)==4.OR.IC1(14)==99)) THEN
        IVER = 2
      ELSE
        IVER = 1
        WRITE(*,*) 'CDA VER.1 IS NOT AVAILABLE'
        WRITE(*,*) (IC1(I),I=1,14)
        NHT = 0
        NEL = 0
        RETURN
      END IF


! === CDA VERSION 2 UNPACKING SUBROUTINE

      CALL CDAUNPK (IC1,JC2,JC3,IC4,
     \              NHT,HGT,NEL,MEL,
     \              NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \              JGQC,IUSE,NOQC,JEQC,KEQC,
     \              MISS,NXHT,NXEL,NXIH,NXIE,NXQC)

! --- IF MOPT=1 THEN ONLY NHT,HGT,NEL AND MEL ARE RETURNED.
      IF (MOPT.EQ.1) RETURN


!   === TRANSFER CDA VER.2 VARIABLES TO CDA VER.1 VARIABLES

!   --- CLEAR CDA VERSION 1 VARIABLES

      DO NH=1,NHT
        IVNM(NH) = MISS
        IVUT(NH) = MISS
        THIK(NH) = MISS
        IVQC(NH) = MISS
        RPHT(NH) = MISS
        IRUT(NH) = MISS
        DO NE=1,NEL(NH)
          OBRP(NE,NH) = MISS
          OBQC(NE,NH) = MISS
          RELY(NE,NH) = MISS
          OBER(NE,NH) = MISS
          DVAL(NE,NH) = MISS
          GUES(NE,NH) = MISS
          GSER(NE,NH) = MISS
          GSDH(NE,NH) = MISS
          GSDT(NE,NH) = MISS
          ANAL(NE,NH) = MISS
        END DO
      END DO

!   --- HEIGHT INFORMATION

      DO NH=1,NHT

        DO NV=1,NHCD(NH)
          IF (KHCD(NV,NH).EQ.100) THEN
            IVNM(NH) = NINT(AHCD(NV,NH))
          ELSE IF (KHCD(NV,NH).EQ.200) THEN
            IVUT(NH) = NINT(AHCD(NV,NH))
          ELSE IF (KHCD(NV,NH).EQ.300) THEN
            THIK(NH) = AHCD(NV,NH)
          ELSE IF (KHCD(NV,NH).EQ.400) THEN
            IVQC(NH) = NINT(AHCD(NV,NH))
          ELSE IF (KHCD(NV,NH).EQ.500) THEN
            RPHT(NH) = AHCD(NV,NH)
          ELSE IF (KHCD(NV,NH).EQ.600) THEN
            IRUT(NH) = NINT(AHCD(NV,NH))
          END IF
        END DO

      END DO

!   --- ELEMENT INFORMATION

      DO NH=1,NHT
      DO NE=1,NEL(NH)

        DO NI=1,NECD(NE,NH)
          IF (KECD(NI,NE,NH).EQ.100) THEN
            OBRP(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.200) THEN
            OBQC(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.300) THEN
            RELY(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.400) THEN
            OBER(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.500) THEN
            DVAL(NE,NH) = AECD(NI,NE,NH)
!         KECD = 606 EXIST BEFORE 2004.07
          ELSE IF (KECD(NI,NE,NH)/100.EQ.6) THEN
            GUES(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.700) THEN
            GSER(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.810) THEN
            GSDH(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.910) THEN
            GSDT(NE,NH) = AECD(NI,NE,NH)
          ELSE IF (KECD(NI,NE,NH).EQ.1000) THEN
            ANAL(NE,NH) = AECD(NI,NE,NH)
          END IF
        END DO
!       GUES IS NOT SAVED IF OBQC & DVAL EXIST
        IF (NINT(OBQC(NE,NH)).NE.MISS .AND.
     \      NINT(DVAL(NE,NH)).NE.MISS .AND.
     \      NINT(GUES(NE,NH)).EQ.MISS)
     \    GUES(NE,NH) = OBQC(NE,NH) - DVAL(NE,NH)

      END DO
      END DO

      RETURN
      END SUBROUTINE P4UNPK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CDAUNPK (IC1,IC2,IC3,IC4,
     \                    NHT,HGT,NEL,MEL,
     \                    NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                    JGQC,IUSE,NOQC,JEQC,KEQC,
     \                    MISS,NXHT,NXEL,NXIH,NXIE,NXQC)

!---------------------------------------------------------------------
!  CDA VER.2 COMPREHENSIVE UNPACKING SUBROUTINE
!                                             2000.12.29 K.ONOGI
!---------------------------------------------------------------------
!  IC1  (I/O;I2) : CDA PART 1
!  IC2  (I/O;I2) : CDA PART 2
!  IC3  (I/O;I2) : CDA PART 3
!  IC4  (IN ;I2) : CDA PART 4
!  NHT  (OUT;I4) : NUMBER OF HEIGHTS (LEVELS)
!  HGT  (OUT;R4) : HEIGHT (LEVEL)
!  NEL  (OUT;I4) : NUMBER OF PARAMETERS FOR A HEIGHT(LEVEL)
!  MEL  (OUT;I4) : ELEMENT ID NUMBER
!  NHCD (OUT;I4) : NUMBER OF HEIGHT INFORMATION
!  KHCD (OUT;I4) : HEIGHT INFORMATION CODE
!  AHCD (OUT:R4) : HEITHT INFORMATION
!  NECD (OUT;I4) : NUMBER OF ELEMENT INFORMATION
!  KECD (OUT;I4) : ELEMENT INFORMATION CODE
!  AECD (OUT:R4) : ELEMENT INFORMATION
!  JGQC (OUT;I4) : QC FINAL JUDGEMENT
!  IUSE (OUT;I4) : USED FLAG
!  NOQC (OUT;I4) : NUMBER OF QC INFORMATIONS
!  JEQC (OUT;I4) : QC JUSTICE OF EACH QC ITEM
!  KEQC (OUT;I4) : EACH QC ITEM
!  MISS (IN ;I4) : MISSING VALUE(=-2^15)
!  NXHT (IN ;I4) : ARRAY SIZE OF HEIGHT               250
!  NXEL (IN ;I4) : ARRAY SIZE OF ELEMENT               60
!  NXIH (IN ;I4) : ARRAY SIZE OF HEIGHT INFORMATION    20
!  NXIE (IN ;I4) : ARRAY SIZE OF ELEMENT INFORMATION   40
!  NXQC (IN ;I4) : ARRAY SIZE OF QC INFORMATION        20
!---------------------------------------------------------------------

      PARAMETER (NXPF4 =  65200)   ! CDA 4 BYTE WORK SPACES
      PARAMETER (NXPF2 =  NXPF4*2) ! CDA 2 BYTE WORK SPACES

! --- CDA ARRAY
      INTEGER(2), DIMENSION(*) :: IC1,IC2,IC3,IC4
      REAL   (4), DIMENSION(NXHT) :: HGT
      INTEGER(4), DIMENSION(NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION(NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION(NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION(NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION(NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION(NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION(NXQC,NXEL,NXHT) :: JEQC,KEQC

! --- CDA PART 4 WORK SPACE
      INTEGER(2), DIMENSION(NXPF2) :: IW4
      INTEGER(4), DIMENSION(NXPF4) :: JW4
      REAL   (4), DIMENSION(NXPF4) :: RW4
! BMT
!      EQUIVALENCE (IW4,JW4,RW4)
! BMT

! --- WORK
      CHARACTER(LEN=8) CSN

      LOGICAL LX/.TRUE./

      IF (LX) THEN

        IF (MISS.EQ.0) MISS = -32768

        LX = .FALSE.
      END IF

! --- IN CASE OF NO CDA DATA RECORD (IC1(7)<=200)
!     SKIP UNPACKING
      IF (IC1(7).LE.200) RETURN

      IF (IC3(1).GT.0.AND.
! BMT
     \    IC3(7).GE.0.AND.IC3(7).LE.9) THEN
        CALL cdatr_i2c(CSN, IC3(3:6), 4)
! BMT
      ELSE                   ! CALLED FROM P4UNPK
        CSN = '????????'
      END IF


! === EXPANDING IC4 AND PROVIDE VALUES TO EACH VARIABLES
      LEN4 = IC1(5)
      IF (IC1(5).LT.0) LEN4 = IC1(5)+32768*2

      DO L=1,LEN4
        IW4(L) = IC4(L)
      END DO

      IF (IC1(14).EQ.2) THEN

        IF (IC1(13) == 200) THEN
          CALL TWOUNPK200 (IW4,NHT,HGT,NEL,MEL,
     \                  NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                  JGQC,IUSE,NOQC,JEQC,KEQC,
     \                  IC1,LEN4,CSN,MISS,
     \                  NXHT,NXEL,NXIH,NXIE,NXQC,NXPF2)
        ELSE
          CALL TWOUNPK (IW4,NHT,HGT,NEL,MEL,
     \                  NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                  JGQC,IUSE,NOQC,JEQC,KEQC,
     \                  IC1,LEN4,CSN,MISS,
     \                  NXHT,NXEL,NXIH,NXIE,NXQC,NXPF2)
        END IF

      ELSE IF (IC1(14).EQ.4) THEN

        LEN4 = LEN4/2

        do i = 1, len4
           call cdatr_i2i4(jw4(i), iw4(2*i-1:2*i), 2)
           call cdatr_i2r4(rw4(i), iw4(2*i-1:2*i), 2)
        end do

        CALL FOURUNPK (JW4,RW4,NHT,HGT,NEL,MEL,
     \                 NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                 JGQC,IUSE,NOQC,JEQC,KEQC,
     \                 IC1,LEN4,CSN,MISS,
     \                 NXHT,NXEL,NXIH,NXIE,NXQC,NXPF4)

      ELSE IF (IC1(14).EQ.99) THEN

        LEN4 = LEN4/2
        do i = 1, len4
           call cdatr_i2i4(jw4(i), iw4(2*i-1:2*i), 2)
           call cdatr_i2r4(rw4(i), iw4(2*i-1:2*i), 2)
        end do

        CALL BETAUNPK (JW4,RW4,NHT,HGT,NEL,MEL,
     \                 NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                 JGQC,IUSE,NOQC,JEQC,KEQC,
     \                 IC1,LEN4,CSN,
     \                 NXHT,NXEL,NXIH,NXIE,NXQC,NXPF4)

      ELSE

        WRITE(*,*) 'SUBROUTINE CDAUNPK'
        WRITE(*,*) IC1(8),IC1(9),IC1(10),'   ',IC1(13),IC1(14)
        WRITE(*,*) 'IC1(14) IS NOT APPROPRIATE. =',IC1(14)
        STOP 201

      END IF


! --- MONITOR
!     IF (CSN.EQ.'03005   ' .AND.
!    \   (IC1(8).EQ.1200.OR.IC1(8).EQ.3000)) THEN
!       WRITE(*,*) 'IC1(14)',IC1(14)
!       WRITE(*,*) (IC1(I),I=1,14),'  ',CSN
!       IF (IC1(14).EQ.2) THEN
!         LW = IC1(5)
!       ELSE
!         LW = IC1(5)/2
!       END IF
!       DO L=1,LW,10
!         IF (MOD(L,100).EQ.1) WRITE(*,*)
!         LWE = MIN(L+9,LW)
!         IF (IC1(14).EQ.2) THEN
!           WRITE(*,8000) L,' :',(IW4(I),I=L,LWE)
!         ELSE
!           WRITE(*,8000) L,' :',(JW4(I),I=L,LWE)
!         END IF
!       END DO
!       WRITE(*,*)
!     END IF
!8000 FORMAT(1H ,I5,A2,10I7)


      RETURN
      END SUBROUTINE CDAUNPK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE TWOUNPK200 (IW4,NHT,HGT,NEL,MEL,
     \                    NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                    JGQC,IUSE,NOQC,JEQC,KEQC,
     \                    IC1,LEN4,CSN,MISS,
     \                    NXHT,NXEL,NXIH,NXIE,NXQC,NXPF2)

!---------------------------------------------------------------------
!     CDA VER.2
!     CDA PART 4  2 BYTE UNPACKING
!                                              2000.12.29  K.ONOGI
!---------------------------------------------------------------------

      PARAMETER (MZHGT = 2567) ! GEOPOTENTIAL HEIGHT ID (10-007)

      INTEGER(2), DIMENSION(*) :: IC1

      INTEGER(2), DIMENSION(NXPF2) :: IW4

      REAL   (4), DIMENSION(NXHT) :: HGT
      INTEGER(4), DIMENSION(NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION(NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION(NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION(NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION(NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION(NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION(NXQC,NXEL,NXHT) :: JEQC,KEQC

      INTEGER(4), DIMENSION(NXIH) :: IHSC
      REAL   (4), DIMENSION(NXIH) :: RHSC
      INTEGER(4), DIMENSION(NXIE) :: IESC
      REAL   (4), DIMENSION(NXIE) :: RESC

      CHARACTER(LEN=8) CSN

      DATA NBAD/0/


! === ALLOCATE DATA TO CDA 1-DIMENSIONAL ARRAY (2 BYTE)

      NHT = IW4(1)                           ! HEIGHT COUNT
      LW = 1

      DO NH=1,NHT

        IF (IW4(LW+1).NE.MISS) THEN
          IHGTSC = IW4(LW+2)
          CALL UNPSCL200 (RHGTSC, IHGTSC)
          HGT(NH) = IW4(LW+1)*RHGTSC
          IF (NINT(RHGTSC)==1.AND.
     \        HGT(NH).GE.-32000.AND.HGT(NH).LT.-2000.)
     \      HGT(NH) = 30000.-HGT(NH)
        ELSE
          HGT(NH) = MISS
        END IF
        NHCD(NH) = IW4(LW+3)                ! HEIGHT INFO. COUNT
        LW = LW+3

        DO NV=1,NHCD(NH)                    ! HEIGHT INFO. CODE
          KHCD(NV,NH) = IW4(LW+NV)/10
          IHSC(NV)    = MOD(IW4(LW+NV),10)
          CALL UNPSCL200 (RHSC(NV),IHSC(NV))
        END DO
        LW = LW+NHCD(NH)

        DO NV=1,NHCD(NH)                    ! HEIGHT INFO.
          AHCD(NV,NH) = IW4(LW+NV)*RHSC(NV)

          IF (KHCD(NV,NH).GE.300.AND.KHCD(NV,NH).LE.599) THEN
            KH = KHCD(NV,NH)/100
            IF (KH.EQ.3.OR.KH.EQ.5) THEN
              IF (AHCD(NV,NH).LT.-2000..AND.
     \            NINT(AHCD(NV,NH)).NE.MISS) THEN
                AHCD(NV,NH) = 30000.-AHCD(NV,NH)
              END IF
            END IF
          END IF

        END DO
        LW = LW+NHCD(NH)

        NEL(NH) = IW4(LW+1)                 ! ELEMENT COUNT
        LW = LW+1

        DO NE=1,NEL(NH)

          MEL(NE,NH)  = IW4(LW+1)           ! ELEMENT CODE
          IUSE(NE,NH) = IW4(LW+2)/1000      ! USED FLAG
          JGQC(NE,NH) = MOD(IW4(LW+2),1000) ! QC JUDGEMENT
          IF (JGQC(NE,NH) .EQ. 99) JGQC(NE,NH) = MISS
          NECD(NE,NH) = IW4(LW+3)           ! ELEMENT INFO. COUNT
          LW = LW+3

          DO NI=1,NECD(NE,NH)
            IF (IW4(LW+NI).GE.0) THEN
              KECD(NI,NE,NH) = IW4(LW+NI)/10  ! ELEMENT INFO. CODE
              IESC(NI) = MOD(IW4(LW+NI),10)
            ELSE
              KECD(NI,NE,NH) = (30000-IW4(LW+NI))/10
              IESC(NI) = MOD((30000-IW4(LW+NI)),10)
            END IF
            CALL UNPSCL200 (RESC(NI), IESC(NI))
          END DO
          LW = LW+NECD(NE,NH)

          DO NI=1,NECD(NE,NH)
            AECD(NI,NE,NH) = IW4(LW+NI)*RESC(NI)  ! ELEMENT INFO.

                             ! SPECIAL TREATMENT FOR Z OVER 32000M
            IF (MEL(NE,NH).EQ.MZHGT) THEN
              KE=KECD(NI,NE,NH)/100
                 ! EXCLUDE OBER DVAL GSER GSDH GSDT THTH
              IF (KE.NE.4 .AND. KE.NE.5 .AND.
     \            KE.NE.7 .AND. KE.NE.8 .AND. KE.NE.9 .AND.
     \            KE.NE.11) THEN
                IF (AECD(NI,NE,NH).LT.-2000..AND.
     \              NINT(AECD(NI,NE,NH)).NE.MISS) THEN
                  AECD(NI,NE,NH) = 30000.-AECD(NI,NE,NH)
                END IF
              END IF
            END IF

          END DO
          LW = LW+NECD(NE,NH)


!         IF (NH.EQ.6) THEN
!           WRITE(*,*) MEL(NE,NH)
!           WRITE(*,'(11I8)') (KECD(NI,NE,NH),NI=1,NECD(NE,NH))
!           WRITE(*,'(11F8.2)') (AECD(NI,NE,NH),NI=1,NECD(NE,NH))
!         END IF


          NOQC(NE,NH) = IW4(LW+1)                 ! QC INFO. COUNT
          LW = LW+1

          DO NQ=1,NOQC(NE,NH)
            JEQC(NQ,NE,NH) = IW4(LW+NQ)/1000      ! EACH QC JUDGE
            KEQC(NQ,NE,NH) = MOD(IW4(LW+NQ),1000) ! EACH QC ITEM
          END DO
          LW = LW+NOQC(NE,NH)

        END DO

      END DO

      IF (LW.NE.LEN4) THEN
        NBAD = NBAD+1
        WRITE(*,*) 'PART 4 LENGTH AND IC1(5) ARE NOT CONSISTENT',
     \             '  (SUBROUTINE TWOUNPK)'
        WRITE(*,*) LW,LEN4,'  ',IC1(7),IC1(8),IC1(9),IC1(10),'  ',CSN
        IF (NBAD.GT.100) STOP 197
      END IF

      RETURN
      END SUBROUTINE TWOUNPK200
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UNPSCL200(RSCL, ISCL)

!---------------------------------------------------------------------
!     UNPACKING SCALE FOR ELEMENT INFORMATION
!                                              2000.12.29 K.ONOGI
!---------------------------------------------------------------------

      IMPLICIT NONE

      REAL   (4) RSCL
      INTEGER(4) ISCL

      IF (ISCL.GE.4.AND.ISCL.LE.9) THEN
        ISCL = ISCL-10
      ELSE
        ISCL = ISCL
      END IF

      RSCL = 10.**ISCL

      RETURN
      END SUBROUTINE UNPSCL200
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE TWOUNPK (IW4,NHT,HGT,NEL,MEL,
     \                    NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                    JGQC,IUSE,NOQC,JEQC,KEQC,
     \                    IC1,LEN4,CSN,MISS,
     \                    NXHT,NXEL,NXIH,NXIE,NXQC,NXPF2)

!---------------------------------------------------------------------
!     CDA VER.2
!     CDA PART 4  2 BYTE UNPACKING
!                                              2000.12.29  K.ONOGI
!                                    Ver. 2.01 2004.06.04  Ohta Yukinari
!---------------------------------------------------------------------

      INTEGER(2), DIMENSION(*) :: IC1

      INTEGER(2), DIMENSION(NXPF2) :: IW4

      REAL   (4), DIMENSION(NXHT) :: HGT
      INTEGER(4), DIMENSION(NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION(NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION(NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION(NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION(NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION(NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION(NXQC,NXEL,NXHT) :: JEQC,KEQC

      INTEGER(4) :: IHGTSC
      REAL   (8) :: RHGTSC
      INTEGER(4), DIMENSION(NXIH) :: IHSC
      REAL   (8), DIMENSION(NXIH) :: RHSC
      INTEGER(4), DIMENSION(NXIE) :: IESC
      REAL   (8), DIMENSION(NXIE) :: RESC

      CHARACTER(LEN=8) CSN

      DATA NBAD/0/


! === ALLOCATE DATA TO CDA 1-DIMENSIONAL ARRAY (2 BYTE)

      NHT = IW4(1)                           ! HEIGHT COUNT
      LW = 1

      DO NH=1,NHT

        IF (IW4(LW+1).NE.MISS) THEN
          IHGTSC = IW4(LW+2)
          CALL UNPSCL (RHGTSC, IHGTSC)
          HGT(NH) = IW4(LW+1)*RHGTSC
        ELSE
          HGT(NH) = MISS
        END IF
        NHCD(NH) = IW4(LW+3)                ! HEIGHT INFO. COUNT
        LW = LW+3

        DO NV=1,NHCD(NH)                    ! HEIGHT INFO. CODE
          IF (IW4(LW+NV).GE.0) THEN
            KHCD(NV,NH) = IW4(LW+NV)/10
            IHSC(NV)    = MOD(IW4(LW+NV),10)
          ELSE
            KHCD(NV,NH) = (65536+IW4(LW+NV))/10
            IHSC(NV)    = MOD((65536+IW4(LW+NV)),10)
          END IF
          CALL UNPSCL (RHSC(NV),IHSC(NV))
        END DO
        LW = LW+NHCD(NH)

        DO NV=1,NHCD(NH)                    ! HEIGHT INFO.
          AHCD(NV,NH) = IW4(LW+NV)*RHSC(NV)
        END DO
        LW = LW+NHCD(NH)

        NEL(NH) = IW4(LW+1)                 ! ELEMENT COUNT
        LW = LW+1

        DO NE=1,NEL(NH)

          MEL(NE,NH)  = IW4(LW+1)           ! ELEMENT CODE
          IUSE(NE,NH) = IW4(LW+2)/1000      ! USED FLAG
          JGQC(NE,NH) = MOD(IW4(LW+2),1000) ! QC JUDGEMENT
          IF (JGQC(NE,NH) .EQ. 99) JGQC(NE,NH) = MISS
          NECD(NE,NH) = IW4(LW+3)           ! ELEMENT INFO. COUNT
          LW = LW+3

          DO NI=1,NECD(NE,NH)
            IF (IW4(LW+NI).GE.0) THEN
              KECD(NI,NE,NH) = IW4(LW+NI)/10  ! ELEMENT INFO. CODE
              IESC(NI) = MOD(IW4(LW+NI),10)
            ELSE
              KECD(NI,NE,NH) = (65536+IW4(LW+NI))/10
              IESC(NI) = MOD((65536+IW4(LW+NI)),10)
            END IF
            CALL UNPSCL (RESC(NI), IESC(NI))
          END DO
          LW = LW+NECD(NE,NH)

          DO NI=1,NECD(NE,NH)
            AECD(NI,NE,NH) = IW4(LW+NI)*RESC(NI)  ! ELEMENT INFO.
          END DO
          LW = LW+NECD(NE,NH)


!         IF (NH.EQ.6) THEN
!           WRITE(*,*) MEL(NE,NH)
!           WRITE(*,'(11I8)') (KECD(NI,NE,NH),NI=1,NECD(NE,NH))
!           WRITE(*,'(11F8.2)') (AECD(NI,NE,NH),NI=1,NECD(NE,NH))
!         END IF


          NOQC(NE,NH) = IW4(LW+1)                 ! QC INFO. COUNT
          LW = LW+1

          DO NQ=1,NOQC(NE,NH)
            JEQC(NQ,NE,NH) = IW4(LW+NQ)/1000      ! EACH QC JUDGE
            KEQC(NQ,NE,NH) = MOD(IW4(LW+NQ),1000) ! EACH QC ITEM
          END DO
          LW = LW+NOQC(NE,NH)

        END DO

      END DO

      IF (LW.NE.LEN4) THEN
        NBAD = NBAD+1
        WRITE(*,*) 'PART 4 LENGTH AND IC1(5) ARE NOT CONSISTENT',
     \             '  (SUBROUTINE TWOUNPK)'
        WRITE(*,*) LW,LEN4,'  ',IC1(7),IC1(8),IC1(9),IC1(10),'  ',CSN
        IF (NBAD.GT.100) STOP 197
      END IF

      RETURN
      END SUBROUTINE TWOUNPK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UNPSCL(RSCL, ISCL)

!---------------------------------------------------------------------
!     UNPACKING SCALE FOR ELEMENT INFORMATION
!                                              2000.12.29 K.ONOGI
!                                    Ver. 2.01 2004.06.04  Ohta Yukinari
!---------------------------------------------------------------------

      IMPLICIT NONE

      REAL   (8) RSCL
      INTEGER(4) ISCL

      REAL   (8), PARAMETER :: BASE = 10.D0


      IF (ISCL.GE.4.AND.ISCL.LE.9) THEN
        ISCL = ISCL-10
      ELSE
        ISCL = ISCL
      END IF

      RSCL = BASE**ISCL

      RETURN
      END SUBROUTINE UNPSCL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FOURUNPK (JW4,RW4,NHT,HGT,NEL,MEL,
     \                     NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                     JGQC,IUSE,NOQC,JEQC,KEQC,
     \                     IC1,LEN4,CSN,MISS,
     \                     NXHT,NXEL,NXIH,NXIE,NXQC,NXPF4)

!---------------------------------------------------------------------
!     CDA VER.2
!     CDA PART 4  4 BYTE UNPACKING
!                                              2000.12.29  K.ONOGI
!---------------------------------------------------------------------

      INTEGER(2), DIMENSION(*) :: IC1

      INTEGER(4), DIMENSION(NXPF4) :: JW4
      REAL   (4), DIMENSION(NXPF4) :: RW4

      REAL   (4), DIMENSION(NXHT) :: HGT
      INTEGER(4), DIMENSION(NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION(NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION(NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION(NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION(NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION(NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION(NXQC,NXEL,NXHT) :: JEQC,KEQC

      CHARACTER(LEN=8) CSN
      DATA NBAD/0/

! === ALLOCATE DATA TO CDA 1-DIMENSIONAL ARRAY (4 BYTE)

      NHT = JW4(1)                          ! HEIGHT COUNT
      LW = 1

      DO NH=1,NHT
        HGT(NH)  = RW4(LW+1)                ! HEIGHT
        NHCD(NH) = JW4(LW+2)                ! HEIGHT INFO. COUNT
        LW = LW+2

        DO NV=1,NHCD(NH)
          KHCD(NV,NH) = JW4(LW+NV)          ! HEIGHT INFO. CODE
        END DO
        LW = LW+NHCD(NH)

        DO NV=1,NHCD(NH)
          AHCD(NV,NH) = RW4(LW+NV)          ! HEIGHT INFORMATION
        END DO
        LW = LW+NHCD(NH)

        NEL(NH) = JW4(LW+1)                 ! ELEMENT COUNT
        LW = LW+1

        DO NE=1,NEL(NH)
          MEL(NE,NH)  = JW4(LW+1)           ! ELEMENT ID CODE
          IUSE(NE,NH) = JW4(LW+2)/1000      ! USED FLAG
          JGQC(NE,NH) = MOD(JW4(LW+2),1000) ! QC JUDGEMENT
          IF (JGQC(NE,NH) .EQ. 99) JGQC(NE,NH) = MISS
          NECD(NE,NH) = JW4(LW+3)           ! ELEMENT INFO. COUNT
          LW = LW+3
          DO NI=1,NECD(NE,NH)
            KECD(NI,NE,NH) = JW4(LW+NI)     ! ELEMENT INFO. CODE
          END DO
          LW = LW+NECD(NE,NH)

          DO NI=1,NECD(NE,NH)
            AECD(NI,NE,NH) = RW4(LW+NI)     ! ELEMENT INFORMATION
          END DO
          LW = LW+NECD(NE,NH)
          NOQC(NE,NH) = JW4(LW+1)           ! QC INFO. COUNT
          LW = LW+1

          DO NQ=1,NOQC(NE,NH)                     ! EACH QC INFO.
            JEQC(NQ,NE,NH) = JW4(LW+NQ)/1000      ! EACH QC JUDGE
            KEQC(NQ,NE,NH) = MOD(JW4(LW+NQ),1000) ! EACH QC ITEM
          END DO
          LW = LW+NOQC(NE,NH)

        END DO

      END DO

      IF (LW.NE.LEN4) THEN
        NBAD = NBAD+1
        WRITE(*,*) 'PART 4 LENGTH AND IC1(5) ARE NOT CONSISTENT',
     \             '  (SUBROUTINE FOURUNPK)'
        WRITE(*,*) LW,LEN4,'  ',IC1(7),IC1(8),IC1(9),IC1(10),'  ',CSN
        IF (NBAD.GT.100) STOP 197
      END IF

      RETURN
      END SUBROUTINE FOURUNPK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BETAUNPK (JW4,RW4,NHT,HGT,NEL,MEL,
     \                     NHCD,KHCD,AHCD,NECD,KECD,AECD,
     \                     JGQC,IUSE,NOQC,JEQC,KEQC,
     \                     IC1,LEN4,CSN,
     \                     NXHT,NXEL,NXIH,NXIE,NXQC,NXPF4)

!---------------------------------------------------------------------
!     CDA VER.2
!     CDA PART 4 BETA UNPACKING FOR TEMPORAL USE
!                                          2000.12.29  K.ONOGI
!---------------------------------------------------------------------

      INTEGER(2), DIMENSION(*) :: IC1

      INTEGER(4), DIMENSION(NXPF4) :: JW4
      REAL   (4), DIMENSION(NXPF4) :: RW4

      REAL   (4), DIMENSION(NXHT) :: HGT
      INTEGER(4), DIMENSION(NXHT) :: NEL,NHCD
      INTEGER(4), DIMENSION(NXIH,NXHT) :: KHCD
      REAL   (4), DIMENSION(NXIH,NXHT) :: AHCD
      INTEGER(4), DIMENSION(NXEL,NXHT) :: MEL,NECD,JGQC,IUSE,NOQC
      INTEGER(4), DIMENSION(NXIE,NXEL,NXHT) :: KECD
      REAL   (4), DIMENSION(NXIE,NXEL,NXHT) :: AECD
      INTEGER(4), DIMENSION(NXQC,NXEL,NXHT) :: JEQC,KEQC

      CHARACTER(LEN=8) CSN
      DATA NBAD/0/

! === ALLOCATE DATA TO CDA 1 DIMENSIONAL ARRAY

      NHT = JW4(1)
      LW = 1

      DO NH=1,NHT
        HGT(NH)  = RW4(LW+1)
        NEL(NH)  = JW4(LW+2)
        NHCD(NH) = JW4(LW+3)
        LW = LW+3
      END DO

      DO NH=1,NHT
        DO NV=1,NHCD(NH)
          KHCD(NV,NH) = JW4(LW+1)
          AHCD(NV,NH) = RW4(LW+2)
          LW = LW+2
        END DO
      END DO

      DO NH=1,NHT
        DO NE=1,NEL(NH)
          MEL(NE,NH)  = JW4(LW+1)
          NECD(NE,NH) = JW4(LW+2)
          JGQC(NE,NH) = JW4(LW+3)
          IUSE(NE,NH) = JW4(LW+4)
          NOQC(NE,NH) = JW4(LW+5)
          LW = LW+5
        END DO
      END DO

      DO NH=1,NHT
        DO NE=1,NEL(NH)
          DO NI=1,NECD(NE,NH)
            KECD(NI,NE,NH) = JW4(LW+1)
            AECD(NI,NE,NH) = RW4(LW+2)
            LW = LW+2
          END DO
        END DO
      END DO

      DO NH=1,NHT
        DO NE=1,NEL(NH)
          DO NQ=1,NOQC(NE,NH)
            JEQC(NQ,NE,NH) = JW4(LW+1)
            KEQC(NQ,NE,NH) = JW4(LW+2)
            LW = LW+2
          END DO
        END DO
      END DO

      IF (LW.NE.LEN4) THEN
        NBAD = NBAD+1
        WRITE(*,*) 'PART 4 LENGTH AND IC1(5) ARE NOT CONSISTENT',
     \             '  (SUBROUTINE BETAUNPK)'
        WRITE(*,*) IC1(7),IC1(8),IC1(9),IC1(10),'  ',CSN
        IF (NBAD.GT.100) STOP 197
      END IF

      RETURN
      END SUBROUTINE BETAUNPK
