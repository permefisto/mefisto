      REAL FUNCTION ANGMIT( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER L'ANGLE MINIMAL DU TRIANGLE (P1 P2 P3) EN RADIANS
C -----   ATTENTION LE TRIANGLE EST SUPPOSE ORIENTE DANS LE SENS DIRECT
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 2 COORDONNEES DES 3 SOMMETS DE L'ANGLE
C            SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C ANGMIT :  ANGLE MINIMAL DU TRIANGLE (P1P2,P1P3) EN RADIANS
C           0 SI P1=P2 ou P1=P3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
C     LES COTES
      X21 = P2(1) - P1(1)
      Y21 = P2(2) - P1(2)
      X31 = P3(1) - P1(1)
      Y31 = P3(2) - P1(2)
      X32 = P3(1) - P2(1)
      Y32 = P3(2) - P2(2)
C
C     LONGUEUR DES COTES
      A1 = X21 * X21 + Y21 * Y21
      A2 = X32 * X32 + Y32 * Y32
      A3 = X31 * X31 + Y31 * Y31
      IF( A1 * A2 * A3 .EQ. 0 ) THEN
         ANGMIT = 0
         RETURN
      ENDIF
C
C     COSINUS DE L'ANGLE
      C1 = (  X21 * X31 + Y21 * Y31 ) / SQRT( A1 * A3 )
      C2 = ( -X32 * X21 - Y32 * Y21 ) / SQRT( A2 * A1 )
      C3 = (  X31 * X32 + Y31 * Y32 ) / SQRT( A3 * A2 )
      C1 = MAX( C1, C2, C3 )
C
      IF( C1 .LE. -1.0 ) THEN
C        TILT SUR APOLLO SI ACOS( -1 -EPS )
         ANGMIT = ATAN( 1. ) * 4.
         RETURN
      ELSE IF( C1 .GE. 1.0 ) THEN
C        TILT SUR APOLLO SI ACOS( 1 + EPS )
         ANGMIT = 0
         RETURN
      ENDIF
C
C     L'ARC COSINUS DE CE COSINUS
      ANGMIT = ACOS( C1 )
      END
