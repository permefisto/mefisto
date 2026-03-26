      SUBROUTINE PECM4R( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PERMUTER CIRCULAIREMENT 4 VALEURS D'UN TABLEAU
C ----- 4 DEVIENT 1, 3 DEVIENT 4, 2 DEVIENT 3, 1 DEVIENT 2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(4)
C
      U    = X(1)
      X(1) = X(4)
      X(4) = X(3)
      X(3) = X(2)
      X(2) = U
      END
