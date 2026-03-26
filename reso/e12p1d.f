      SUBROUTINE E12P1D( X, F1, F2, DP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER NBPOLY, NPI ET
C -----    F1,F2  = COORDONNEES DES MILIEUX DES ARETES DU TRIANGLE
C
C ENTREE :
C --------
C X      : COORDONNEES DES 3 SOMMETS DU TRIANGLE COURANT
C
C SORTIES:
C --------
C F1, F2 : COORDONNEES X ET Y DU BARYCENTRE DU TRIANGLE
C DP     : GRADIENT (CONSTANT) DES POLYNOMES DE BASE SUE LE TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA
      DOUBLE PRECISION  F1(1), F2(1), DP(2,3)
      REAL              X(3,2)
C
      F1(1) = ( X(1,1) + X(2,1) + X(3,1) ) / 3D0
      F2(1) = ( X(1,2) + X(2,2) + X(3,2) ) / 3D0
C
C     LE GRADIENT
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
C
      DELTA = X21 * Y31 - X31 * Y21
      DELTA = 1D0 / DELTA
C
      DP(1,1) = -Y32 * DELTA
      DP(2,1) =  X32 * DELTA
C
      DP(1,2) =  Y31 * DELTA
      DP(2,2) = -X31 * DELTA
C
      DP(1,3) = -Y21 * DELTA
      DP(2,3) =  X21 * DELTA
      END
