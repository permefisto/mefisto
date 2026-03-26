      SUBROUTINE PECI33( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  PERMUTATION CIRCULAIRE DES No DES 3 SOMMETS ET 6 TANGENTES CROISEES
C -----  DU TRIANGLE HCT LORSQUE LE SOMMET 2 DEVIENT LE SOMMET 1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      INTEGER  U, X(9)

      U    = X(1)
      X(1) = X(2)
      X(2) = X(3)
      X(3) = U

      U    = X(4)
      X(4) = X(6)
      X(6) = X(8)
      X(8) = U

      U    = X(5)
      X(5) = X(7)
      X(7) = X(9)
      X(9) = U

      RETURN
      END
