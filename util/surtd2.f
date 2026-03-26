      DOUBLE PRECISION FUNCTION SURTD2( P1 , P2 , P3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE LA SURFACE D'UN TRIANGLE DEFINI PAR 3 POINTS DE R**2
C -----
C PARAMETRES D ENTREE :
C ---------------------
C P1 P2 P3 : LES 3 FOIS 2 COORDONNEES DES SOMMETS DU TRIANGLE
C
C PARAMETRE RESULTAT :
C --------------------
C SURTD2 : SURFACE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      DOUBLE PRECISION  P1(2), P2(2), P3(2)
C
C     LA SURFACE DU TRIANGLE
      SURTD2 = ( ( P2(1)-P1(1) ) * ( P3(2)-P1(2) )
     %         - ( P2(2)-P1(2) ) * ( P3(1)-P1(1) ) ) * 0.5D0
      END
