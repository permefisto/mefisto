      SUBROUTINE PTPLPL( P0, D, A, B, P, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DU POINT P A DISTANCE D
C -----    SUR LA DROITE INTERSECTION DES 2 PLANS DEFINIS PAR A ET B
C          ET CONTENANT LE POINT P0
C
C ENTREES:
C --------
C P0     : POINT APPARTENANT AUX 2 PLANS ET A LA DROITE INTERSECTION
C D      : DISTANCE P 0P SUR LA DROITE INTERSECTION
C A      : COEFFICIENTS DE L EQUATION DU PLAN 1
C          A(1) * X + A(2) * Y + A(3) * Z + A(4) = 0
C B      : COEFFICIENTS DE L EQUATION DU PLAN 2
C          B(1) * X + B(2) * Y + B(3) * Z + B(4) = 0
C
C SORTIE :
C --------
C P      : LES 3 COORDONNEES DU POINT SUR LA DROITE D'INTERSECTION
C          DES 2 PLANS (CHOIX DE LA RACINE POSITIVE)
C IERR   : 0 SI PAS D'ERREUR, 1 SI LES 2 PLANS SONT PARALLELES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1998
C2345X7..............................................................012
      DOUBLE PRECISION  P0(3), D, A(4), B(4), P(3), R1, R2, R3, DEN
C
C     LE VECTEUR ORTHOGONAL AUX VECTEURS NORMAUX DES 2 PLANS
      R1  = A(2) * B(3) - A(3) * B(2)
      R2  = A(3) * B(1) - A(1) * B(3)
      R3  = A(1) * B(2) - A(2) * B(1)
      DEN = R1**2 + R2**2 + R3**2
C
      IF( DEN .LT. 1D-12*D ) THEN
C
C        VECTEUR NORMAL NUL (DE LA DROITE INTERSECTION)
C        <=> 2 PLANS PARALLELES
         IERR = 1
C
      ELSE
C
C        CAS GENERAL P = P0 + D * (VECTEUR UNITAIRE DE LA DROITE INTERSECTION)
         DEN = SQRT( DEN )
         P(1) = P0(1) + D * R1 / DEN
         P(2) = P0(2) + D * R2 / DEN
         P(3) = P0(3) + D * R3 / DEN
C
      ENDIF
      RETURN
      END
