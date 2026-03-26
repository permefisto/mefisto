      REAL FUNCTION HAMNT2( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA HAUTEUR MINIMALE DU TRIANGLE P1 P2 P3
C -----
C
C ENTREES :
C ---------
C P1,P2,P3: LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C           SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C HAMNT2  : HAUTEUR MINIMALE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
      HAMNT2 = HAUTR2( P1, P2, P3 )
C
      H      = HAUTR2( P2, P3, P1 )
      IF( H .LT. HAMNT2 ) HAMNT2 = H
C
      H      = HAUTR2( P3, P1, P2 )
      IF( H .LT. HAMNT2 ) HAMNT2 = H
      END
