      SUBROUTINE TPPD2P1D( X, tPP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE INTEGRALE t[P1][P1] dx
C ----- POUR UN TRIANGLE 2P1D LAGRANGE DE DEGRE 1
C       AVEC INTEGRATION NUMERIQUE de tPP AUX 3 SOMMETS DU TRIANGLE
C       POUR OBTENIR UNE MATRICE DIAGONALE
C
C ENTREES:
C --------
C X      : COORDONNEES X ET Y DES 3 SOMMETS DU TRIANGLE
C
C SORTIES:
C --------
C tPP    : Integrale t[P][P] dX = TPP  est une MATRICE DIAGONALE
C tDPDP  : dt ALFA Integrale t[DP][DP] dX
C tPU2P  : dt BETA Integrale ( (P vn)**2 + (P wn)**2 )  t[P][P] dX
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      REAL              X(3,2)
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, S
      DOUBLE PRECISION  tPP(3)
      INTRINSIC         ABS
C
C     DETERMINANT DE LA TRANSFORMATION F: e ref --> e
C     -----------------------------------------------
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
C
      DELTA = ABS( X21 * Y31 - X31 * Y21 )
C
C     Integrale t[P][P] dX = TPP MATRICE DIAGONALE
C     --------------------------------------------
      S = DELTA / 6D0
      tPP(1) = S
      tPP(2) = S
      tPP(3) = S
C
      RETURN
      END
