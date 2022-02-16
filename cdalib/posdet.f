      SUBROUTINE POSDET
     I  (OLAT,OLON,GLAT,IM,JM,
     O   II,JJ,X,Y)
C***************************************************************
C (3) 観測点の水平座標相対値計算(OI,OJ)
C***************************************************************
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z,\)
      DIMENSION GLAT(IM,JM)
C*************** PROCEDURE *************************************
* === 観測点位置の特定

* --- 観測点の緯度経度
C     OLAT = IC1(9)/100.  !２バイトデータ
C     OLON = IC1(10)/100. !２バイトデータ
      IF (OLON.LT.0.) OLON = OLON+360.

* --- 極点の緯度
*     (JUST90度にすると別の扱いが必要なのでわずかにずらす)
C     IF (IC1(9).EQ. 9000) OLAT =  89.99
C     IF (IC1(9).EQ.-9000) OLAT = -89.99
      IF(OLAT.GE. 90.) OLAT = 89.99
      IF(OLAT.LE.-90.) OLAT =-89.99

*   --- 観測点のOI,OJ(格子での水平座標位置) を求める (GAUSS

*   --- I方向
        OI = OLON * IM / 360. + 1.
*   --- J方向(緯度方向は等間隔でないので線形内挿は厳密には不正確
        IF (OLAT.LE.GLAT(1,1) .AND. OLAT.GE.GLAT(1,JM)) THEN
          DO 120 J=2,JM
            IF (OLAT .GT. GLAT(1,J)) THEN
              OJ = (J-1) + (GLAT(1,J-1)-OLAT)/(GLAT(1,J-1)-GLAT(1,J))
              GOTO 130
            END IF
  120     CONTINUE
  130     CONTINUE
        ELSE IF (OLAT .GT. GLAT(1,1)) THEN       ! 北極領域
          OJ =  (90.-OLAT)/(90.-GLAT(1,1))
        ELSE IF (OLAT .LT. GLAT(1,JM)) THEN      ! 南極領域
          OJ = JM + (GLAT(1,JM)-OLAT)/(GLAT(1,JM)-(-90.))
        END IF
C====================
C     観測点内挿位置
C====================
        II = OI
        X  = OI-II
        JJ = OJ
        Y  = OJ-JJ
      RETURN
      END
