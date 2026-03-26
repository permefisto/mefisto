      DOUBLE PRECISION FUNCTION ANGLED( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER L'ANGLE (P1P2,P1P3) EN RADIANS
C -----
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 2 COORDONNEES DES 3 SOMMETS DE L'ANGLE
C               SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C ANGLED :  ANGLE (P1P2,P1P3) EN RADIANS ENTRE [0 ET 2PI]
C           0 SI P1=P2 ou P1=P3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      DOUBLE PRECISION  P1(2),P2(2),P3(2),X21,Y21,X31,Y31,A1,A2,D,C
C
C     LES COTES
      X21 = P2(1) - P1(1)
      Y21 = P2(2) - P1(2)
      X31 = P3(1) - P1(1)
      Y31 = P3(2) - P1(2)
C
C     LONGUEUR DES COTES
      A1 = X21 * X21 + Y21 * Y21
      A2 = X31 * X31 + Y31 * Y31
      D  = SQRT( A1 * A2 )
      IF( D .EQ. 0 ) THEN
         ANGLED = 0
         RETURN
      ENDIF
C
C     COSINUS DE L'ANGLE
      C  = ( X21 * X31 + Y21 * Y31 ) / D
      IF( C .LE. -1.D0 ) THEN
C        TILT SUR APOLLO SI ACOS( -1 -EPS )
         ANGLED = ATAN( 1.D0 ) * 4.D0
         RETURN
      ELSE IF( C .GE. 1.D0 ) THEN
C        TILT SUR APOLLO SI ACOS( 1 + EPS )
         ANGLED = 0
         RETURN
      ENDIF
C
      ANGLED = ACOS( C )
      IF( X21 * Y31 - X31 * Y21 .LT. 0 ) THEN
C        DEMI PLAN INFERIEUR
         ANGLED = 8.D0 * ATAN( 1.D0 ) - ANGLED
      ENDIF
      END
