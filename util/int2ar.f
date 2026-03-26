      SUBROUTINE INT2AR( P1, P2, P3, P4, OUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LES 2 ARETES DE R**2 P1-P2  P3-P4 S'INTERSECTENT ELLES
C -----    ENTRE LEURS SOMMETS? (SOMMETS COMPRIS)
C
C ENTREES:
C --------
C P1,P2,P3,P4 : LES 2 COORDONNEES REELLES DES SOMMETS DES 2 ARETES
C
C SORTIE :
C --------
C OUI    : .TRUE. SI INTERSECTION, .FALSE. SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    OCTOBRE 1991
C2345X7..............................................................012
      DOUBLE PRECISION  P1(2),P2(2),P3(2),P4(2)
      DOUBLE PRECISION  X21,Y21,D21,X43,Y43,D43,D,X,Y,XX
      LOGICAL  OUI
C
C     LONGUEUR DES ARETES
      X21 = P2(1)-P1(1)
      Y21 = P2(2)-P1(2)
      D21 = X21**2 + Y21**2
C
      X43 = P4(1)-P3(1)
      Y43 = P4(2)-P3(2)
      D43 = X43**2 + Y43**2
C
C     LES 2 ARETES SONT-ELLES JUGEES PARALLELES ?
      D = X43 * Y21 - Y43 * X21
      IF( ABS(D) .LE. 0.001 * SQRT(D21 * D43) ) THEN
C        ARETES PARALLELES . PAS D'INTERSECTION
         OUI = .FALSE.
         RETURN
      ENDIF
C
C     LES 2 COORDONNEES DU POINT D'INTERSECTION
      X = ( P1(1)*X43*Y21 - P3(1)*X21*Y43 - (P1(2)-P3(2))*X21*X43 ) / D
      Y =-( P1(2)*Y43*X21 - P3(2)*Y21*X43 - (P1(1)-P3(1))*Y21*Y43 ) / D
C
C     COORDONNEES DE X,Y DANS LE REPERE NS1-NS2
      XX  = ( X - P1(1) ) * X21 + ( Y - P1(2) ) * Y21
C     LE POINT EST IL ENTRE P1 ET P2 ?
      OUI = -0.00001D0*D21 .LE. XX .AND. XX .LE. 1.00001D0*D21
C
C     COORDONNEES DE X,Y DANS LE REPERE NS3-NS4
      XX  = ( X - P3(1) ) * X43 + ( Y - P3(2) ) * Y43
C     LE POINT EST IL ENTRE P3 ET P4 ?
      OUI = OUI .AND. -0.00001D0*D43 .LE. XX .AND. XX .LE. 1.00001D0*D43
      END
