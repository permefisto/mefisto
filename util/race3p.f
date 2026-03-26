      SUBROUTINE RACE3P( P1, P2, P3, RAYON, CENTRE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE RAYON ET CENTRE DU CERCLE PASSANT PAR 3 POINTS DE R2
C -----
C
C ENTREES :
C ---------
C P1, P2, P3 : LES 2 COORDONNEES DES 3 POINTS
C
C SORTIES :
C ---------
C RAYON  : RAYON DU CERCLE
C CENTRE : LES 2 COORDONNEES DU CENTRE DU CERCLE
C IERR   : 0 SI CALCUL FAIT CORRECTEMENT
C          1 SI LES 3 POINTS SONT ALIGNES OU CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1993
C23456---------------------------------------------------------------012
      REAL    P1(2),P2(2),P3(2),CENTRE(2)
C
      X21 = P2(1) - P1(1)
      Y21 = P2(2) - P1(2)
      D21 = X21 * X21 + Y21 * Y21
      X31 = P3(1) - P1(1)
      Y31 = P3(2) - P1(2)
      D31 = X31 * X31 + Y31 * Y31
C
C     2 * AIRE DU TRIANGLE
      D   = X21 * Y31 - X31 * Y21
      IF( ABS(D) .LE. 1E-5 * (D21+D31) ) THEN
         IERR = 1
         RETURN
      ENDIF
C
C     LE CENTRE
      CENTRE(1) = ( D21 * Y31 - D31 * Y21 ) / D
      CENTRE(2) = ( X21 * D31 - X31 * D21 ) / D
C
      RAYON = SQRT( (P1(1)-CENTRE(1))**2  + (P1(2)-CENTRE(2))**2 )
      IERR  = 0
      END
