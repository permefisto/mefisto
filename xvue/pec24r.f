      SUBROUTINE PEC24R( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PERMUTER CIRCULAIREMENT 4 VALEURS D'UN TABLEAU
C ----- 2 DEVIENT 4, 3 DEVIENT 1, 4 DEVIENT 2, 1 DEVIENT 3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(4)
C
      U    = X(1)
      X(1) = X(3)
      X(3) = U
C
      U    = X(2)
      X(2) = X(4)
      X(4) = U
      END
