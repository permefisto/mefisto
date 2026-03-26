      REAL FUNCTION T2MEMN( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA MEDIANE MINIMALE DU TRIANGLE P1 P2 P3
C -----
C
C ENTREES :
C ---------
C P1,P2,P3: LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C           SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C T2MEMN  : MEDIANE MINIMALE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
      T2MEMN = TR2MED( P1, P2, P3 )
C
      H      = TR2MED( P2, P3, P1 )
      IF( H .LT. T2MEMN ) T2MEMN = H
C
      H      = TR2MED( P3, P1, P2 )
      IF( H .LT. T2MEMN ) T2MEMN = H
      END
