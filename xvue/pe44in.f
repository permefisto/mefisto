      SUBROUTINE PE44IN( X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INVERSER L'ORDRE 1 2 3 4 EN 1 4 3 2 DES 4 SOMMETS ET TANGENTES
C ----- CROISEES DU QUADRANGLE DVS REDUIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL   X(12)
C
      U    = X(2)
      X(2) = X(4)
      X(4) = U
C
      U    = X(5)
      X(5) = X(6)
      X(6) = U
C
      U     = X( 7)
      X( 7) = X(12)
      X(12) = U
C
      U     = X( 8)
      X( 8) = X(11)
      X(11) = U
C
      U     = X( 9)
      X( 9) = X(10)
      X(10) = U
      END
