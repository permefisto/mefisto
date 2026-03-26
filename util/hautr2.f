      REAL FUNCTION HAUTR2( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA HAUTEUR AU SOMMET P1 DU TRIANGLE P1 P2 P3
C -----
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C               SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C HAUTR2  : HAUTEUR ENTRE LE SOMMET P1 ET LE COTE P2 P3
C           0 SI P1=P2 ou P1=P3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
C     LA SURFACE DU TRIANGLE
      S = SURTR2( P1, P2, P3 )
C
C     LA LONGUEUR DE P2 P3
      HAUTR2 = SQRT( (P3(1)-P2(1))**2 + (P3(2)-P2(2))**2 )
      IF( HAUTR2 .EQ. 0 ) RETURN
C
C     LA HAUTEUR ISSUE DE P1
      HAUTR2 = 2 * ABS( S ) / HAUTR2
      END
