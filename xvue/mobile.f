      SUBROUTINE MOBILE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER LE MOBILE DE MEFISTO
C -----     UN PARTERRE DE TRIANGLES COLORIES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS   DECEMBRE 1997
C2345X7..............................................................012
      PARAMETER  (NBPAS=32, NBTIME=1000, A=1.7, R=2.0, Z=-0.75)
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/mecoit.inc"
      REAL      XYZ(3,3),XYZT(3,3)
C
C     LE TRIANGLE EQUILATERAL DE BARYCENTRE EN (0,0)
      RAC3 = SQRT( 3.0 )
C     SOMMET 1
      XYZT(1,1) = -A / 2.0
      XYZT(2,1) = -A * RAC3 / 6.0
      XYZT(3,1) =  Z
C     SOMMET 2
      XYZT(1,2) =  A / 2.0
      XYZT(2,2) =  XYZT(2,1)
      XYZT(3,2) =  Z
C     SOMMET 3
      XYZT(1,3) =  0.0
      XYZT(2,3) =  A / RAC3
      XYZT(3,3) =  Z
C
C     TRIANGLE
      ANGLE = 0.0
      PIS2  = ATAN(1.0) * 2
      PI    = PIS2 * 2
      PIP   = PI * 1.2
      PASA  = PIS2 * 4.0 / NBPAS
C
      XYZ(3,1) = Z
      XYZ(3,2) = Z
      XYZ(3,3) = Z
C
      DO 100 I=0,NBPAS
         ANGLE = PIP + I * PASA
         COSA = COS(ANGLE)
         SINA = SIN(ANGLE)
         XD = R * COSA
         YD = R * SINA
C        LA ROTATION + DEPLACEMENT
         DO 10 K=1,3
            XYZ(1,K) = XYZT(1,K)*COSA - XYZT(2,K)*SINA + XD
            XYZ(2,K) = XYZT(1,K)*SINA + XYZT(2,K)*COSA + YD
 10      CONTINUE
         NC = N1COUL + MOD(I,NDCOUL-N1COUL)
         CALL FACE3D( NC, NCNOIR, 3, XYZ )
 100  CONTINUE
      RETURN
      END
