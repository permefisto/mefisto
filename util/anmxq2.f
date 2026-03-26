      SUBROUTINE ANMXQ2( P1, P2, P3, P4, ANGMAX, NSTMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LE PLUS GRAND DES 4 ANGLES AU SOMMET D'UN QUADRANGLE
C -----   CONVEXE
C
C ENTREES :
C ---------
C P1,P2,P3,P4 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C               SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C ANGMAX  : ANGLE MAXIMAL EN RADIANS DES 4 ANGLES DU QUADRANGLE CONVEXE
C NSTMAX  : NUMERO DE 1 A 4 DU SOMMET CORRESPONDANT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2),P4(2)
C
C     LES 4 ANGLES
      ANGMAX = ANGLE2( P1, P2, P4 )
      A      = ANGLE2( P2, P3, P1 )
      IF( A .GT. ANGMAX ) THEN
         ANGMAX = A
         NSTMAX = 2
      ELSE
         NSTMAX = 1
      ENDIF
C
      A = ANGLE2( P3, P4, P2 )
      IF( A .GT. ANGMAX ) THEN
         ANGMAX = A
         NSTMAX = 3
      ENDIF
C
      A = ANGLE2( P4, P1, P3 )
      IF( A .GT. ANGMAX ) THEN
         ANGMAX = A
         NSTMAX = 4
      ENDIF
      END
