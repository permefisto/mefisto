      SUBROUTINE PECM33( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PERMUTATION CIRCULAIRE DES 3 SOMMETS ET 6 TANGENTES CROISEES
C ----- DU TRIANGLE HCT LORSQUE LE SOMMET 3 DEVIENT LE SOMMET 1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL X(9)
C
      U    = X(1)
      X(1) = X(3)
      X(3) = X(2)
      X(2) = U
C
      U    = X(4)
      X(4) = X(8)
      X(8) = X(6)
      X(6) = U
C
      U    = X(7)
      X(7) = X(5)
      X(5) = X(9)
      X(9) = U
C
      END
