      SUBROUTINE PECI3R( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PERMUTER CIRCULAIREMENT 3 VALEURS D'UN TABLEAU
C ----- 2 DEVIENT 1, 3 DEVIENT 2, 1 DEVIENT 3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(3)
C
      U    = X(1)
      X(1) = X(2)
      X(2) = X(3)
      X(3) = U
      END
