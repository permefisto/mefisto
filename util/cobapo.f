      SUBROUTINE COBAPO( NBS, XYZP, XYZB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DU BARYCENTRE DU POLYGONE
C -----    DE NBS SOMMETS DE COORDONNEES XYZP
C
C ENTREES:
C --------
C NBS    : NOMBRE DE SOMMETS DU POLYGONE
C XYZP   : LES 3 COORDONNEES XYZ DES NBS SOMMETS DU POLYGONE
C
C SORTIES:
C --------
C XYZP   : LES 3 COORDONNEES XYZ DU BARYCENTRE DES
C          NBS SOMMETS DU POLYGONE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1994
C2345X7..............................................................012
      REAL   XYZP( 1:3, 1:NBS ), XYZB( 1:3 )
C
      DO 20 I=1,3
         R = 0.0
         DO 10 J=1,NBS
            R = R + XYZP(I,J)
 10      CONTINUE
         XYZB(I) = R / NBS
 20   CONTINUE
      END
