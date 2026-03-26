      SUBROUTINE COBAST( NBS, X, Y, Z, XYZB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DU BARYCENTRE DU POLYGONE
C -----    DE NBS SOMMETS DE COORDONNEES X, Y, Z

C ENTREES:
C --------
C NBS    : NOMBRE DE SOMMETS DU POLYGONE
C X, Y, Z: LES 3 COORDONNEES XYZ DES NBS SOMMETS DU POLYGONE

C SORTIES:
C --------
C XYZB   : LES 3 COORDONNEES XYZ DU BARYCENTRE DES
C          NBS SOMMETS DU POLYGONE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       JUIN 1996
C2345X7..............................................................012
      REAL   X(NBS), Y(NBS), Z(NBS), XYZB( 1:3 )

C     CALCUL DES COORDONNEES XYZB DU BARYCENTRE
      XYZB(1) = 0
      XYZB(2) = 0
      XYZB(3) = 0
      DO J=1,NBS
         XYZB(1) = XYZB(1) + X(J)
         XYZB(2) = XYZB(2) + Y(J)
         XYZB(3) = XYZB(3) + Z(J)
      ENDDO
      XYZB(1) = XYZB(1) / NBS
      XYZB(2) = XYZB(2) / NBS
      XYZB(3) = XYZB(3) / NBS

      RETURN
      END
