       DOUBLE PRECISION FUNCTION VOLTET( P1 , P2 , P3 , P4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE VOLUME DU TETRAEDRE P1 P2 P3 P4
C -----    SOIT ENCORE LE DETERMINANT( P2-P1 , P3-P1 , P4-P1 ) / 6D0
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C VOLTET : LE VOLUME DU TETRAEDRE P1 P2 P3 P4 ou
C          LE DETERMINANT DE LA MATRICE ( P2-P1 , P3-P1 , P4-P1 ) / 6D0
C          LE RESULTAT PEUT ETRE >0. =0. OU <0.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS  MAI 1991
C2345X7..............................................................012
      DOUBLE PRECISION  P1(3), P2(3), P3(3), P4(3)
      DOUBLE PRECISION  X1,Y1,Z1, X12,Y12,Z12, X13,Y13,Z13, X14,Y14,Z14,
     %                  D, E, F

      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)

      X12 = X1 - P2(1)
      Y12 = Y1 - P2(2)
      Z12 = Z1 - P2(3)

      X13 = X1 - P3(1)
      Y13 = Y1 - P3(2)
      Z13 = Z1 - P3(3)

      X14 = X1 - P4(1)
      Y14 = Y1 - P4(2)
      Z14 = Z1 - P4(3)

      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13

      VOLTET = ( X12 * D + X13 * E + X14 * F ) / 6D0

      RETURN
      END
