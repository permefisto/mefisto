      SUBROUTINE PE33IN( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INVERSION DU SENS DES SOMMETS ET TANGENTES DU TRIANGLE HCT
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(9)
C
      U    = X(2)
      X(2) = X(3)
      X(3) = U
C
      U    = X(4)
      X(4) = X(5)
      X(5) = U
C
      U    = X(6)
      X(6) = X(9)
      X(9) = U
C
      U    = X(7)
      X(7) = X(8)
      X(8) = U
      END
