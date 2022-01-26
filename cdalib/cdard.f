!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  INCLUDED MEMBERS
!     CDARD
!     CDAIN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CDARD (IC1,IC2,IC3,IC4, IEND,MCDA)

!---------------------------------------------------------------------
!     READ CDA 1 RECORD (SIMPLIFIED CDAIN SUBROUTINE)
!                                           2000.12.24   K.ONOGI
!---------------------------------------------------------------------

      INTEGER(2), DIMENSION(*) :: IC1,IC2,IC3,IC4

      IEND = 0     ! FILE END FLAG

! --- READ CDA 1 RECORD

      READ(MCDA,END=2000) IC1(1),IC1(2),
     \                   (IC1(I),I=3,IC1(2)),
     \                   (IC2(I),I=1,IC1(3)),
     \                   (IC3(I),I=1,IC1(4)),
     \                   (IC4(I),I=1,IC1(5))

* --- IN CASE THE CDA RECORD LENGTH IS OVER 32767
      LTOT= IC1(1)
      LC4 = IC1(5)
      IF(LTOT.LE.0) THEN
        LTOT=LTOT+32768*2
        IF(LC4.LE.0) THEN
          LC4=LC4+32768*2
        END IF !!TAB is replaced by 8 spaces. K. Ito 2012.09.09
        BACKSPACE MCDA     !! BACKSPACE 1 RECORD
        READ(MCDA,END=2000) IC1(1),IC1(2),
     \                     (IC1(I),I=3,IC1(2)),
     \                     (IC2(I),I=1,IC1(3)),
     \                     (IC3(I),I=1,IC1(4)),
     \                     (IC4(I),I=1,LC4)
      ENDIF

      RETURN

 2000 CONTINUE

      IEND = 1

      RETURN
      END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CDAIN (IC1,IC2,IC3,IC4,
     \                  LC1,LC2,LC3,LC4,IEND,MCDA)

!---------------------------------------------------------------------
!     CDA�t�@�C����1���R�[�h��ǂ�
!                                           1994.9.13   K.ONOGI
!---------------------------------------------------------------------

      INTEGER*2 IC1(*), IC2(*), IC3(*), IC4(*)
      LOGICAL   LX/.TRUE./

      IF (LX) THEN
        NRCMAX = 0
        LX = .FALSE.
      END IF

      IEND = 0     ! �t�@�C���I���t���O !

! --- CDA��1���R�[�h��ǂ�
      READ(MCDA,END=2000) IC1(1),IC1(2),
     \                   (IC1(I),I=3,IC1(2)),
     \                   (IC2(I),I=1,IC1(3)),
     \                   (IC3(I),I=1,IC1(4)),
     \                   (IC4(I),I=1,IC1(5))

! --- CDA�̊e���̒���
      LTOT= IC1(1)
      LC1 = IC1(2)
      LC2 = IC1(3)
      LC3 = IC1(4)
      LC4 = IC1(5)

* --- IN CASE THE CDA RECORD LENGTH IS OVER 32767
      IF(LTOT.LE.0) THEN
        LTOT=LTOT+32768*2
        IF(LC4.LE.0) THEN
          LC4=LC4+32768*2
        END IF  !!TAB is replaced by 8 spaces. K. Ito 2012.09.09
        BACKSPACE MCDA     !! BACKSPACE 1 RECORD
        READ(MCDA,END=2000) (IC1(I),I=1,LC1),
     \                      (IC2(I),I=1,LC2),
     \                      (IC3(I),I=1,LC3),
     \                      (IC4(I),I=1,LC4)
      ENDIF

! --- CDA�̍ő僌�R�[�h��
      IF (IC1(1).GT.NRCMAX) NRCMAX = IC1(1)

      RETURN

 2000 CONTINUE

      IEND = 1

      RETURN
      END
