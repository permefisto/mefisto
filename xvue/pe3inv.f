      SUBROUTINE PE3INV( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INVERSER L'ORDRE 1 2 3 EN 1 3 2 DES 3 VALEURS DU TABLEAU
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(3)
C
      U    = X(2)
      X(2) = X(3)
      X(3) = U
      END
