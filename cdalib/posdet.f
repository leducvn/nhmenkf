      SUBROUTINE POSDET
     I  (OLAT,OLON,GLAT,IM,JM,
     O   II,JJ,X,Y)
C***************************************************************
C (3) �ϑ��_�̐������W���Βl�v�Z(OI,OJ)
C***************************************************************
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z,\)
      DIMENSION GLAT(IM,JM)
C*************** PROCEDURE *************************************
* === �ϑ��_�ʒu�̓���

* --- �ϑ��_�̈ܓx�o�x
C     OLAT = IC1(9)/100.  !�Q�o�C�g�f�[�^
C     OLON = IC1(10)/100. !�Q�o�C�g�f�[�^
      IF (OLON.LT.0.) OLON = OLON+360.

* --- �ɓ_�̈ܓx
*     (JUST90�x�ɂ���ƕʂ̈������K�v�Ȃ̂ł킸���ɂ��炷)
C     IF (IC1(9).EQ. 9000) OLAT =  89.99
C     IF (IC1(9).EQ.-9000) OLAT = -89.99
      IF(OLAT.GE. 90.) OLAT = 89.99
      IF(OLAT.LE.-90.) OLAT =-89.99

*   --- �ϑ��_��OI,OJ(�i�q�ł̐������W�ʒu) �����߂� (GAUSS

*   --- I����
        OI = OLON * IM / 360. + 1.
*   --- J����(�ܓx�����͓��Ԋu�łȂ��̂Ő��`���}�͌����ɂ͕s���m
        IF (OLAT.LE.GLAT(1,1) .AND. OLAT.GE.GLAT(1,JM)) THEN
          DO 120 J=2,JM
            IF (OLAT .GT. GLAT(1,J)) THEN
              OJ = (J-1) + (GLAT(1,J-1)-OLAT)/(GLAT(1,J-1)-GLAT(1,J))
              GOTO 130
            END IF
  120     CONTINUE
  130     CONTINUE
        ELSE IF (OLAT .GT. GLAT(1,1)) THEN       ! �k�ɗ̈�
          OJ =  (90.-OLAT)/(90.-GLAT(1,1))
        ELSE IF (OLAT .LT. GLAT(1,JM)) THEN      ! ��ɗ̈�
          OJ = JM + (GLAT(1,JM)-OLAT)/(GLAT(1,JM)-(-90.))
        END IF
C====================
C     �ϑ��_���}�ʒu
C====================
        II = OI
        X  = OI-II
        JJ = OJ
        Y  = OJ-JJ
      RETURN
      END
