      SUBROUTINE PE4INV( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INVERSER L'ORDRE 1 2 3 4 EN 1 4 3 2 DES 4 VALEURS DU TABLEAU
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL   X(4)
C
      U    = X(2)
      X(2) = X(4)
      X(4) = U
      END
