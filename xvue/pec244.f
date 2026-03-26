      SUBROUTINE PEC244( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  LE SOMMET 3 DU QUADRANGLE DVS DEVIENT LE SOMMET 1
C -----  ET DANS LE SENS DIRECT LES PERMUTATIONS DES SOMMETS
C        ET TANGENTES CROISEES SUIVENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL   X(12)
C
      U    = X(1)
      X(1) = X(3)
      X(3) = U
      U    = X(2)
      X(2) = X(4)
      X(4) = U
C
      U     = X( 5)
      X( 5) = X( 9)
      X( 9) = U
      U     = X( 6)
      X( 6) = X(10)
      X(10) = U
C
      U     = X( 7)
      X( 7) = X(11)
      X(11) = U
      U     = X( 8)
      X( 8) = X(12)
      X(12) = U
C
      END
