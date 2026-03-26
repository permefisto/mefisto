       DOUBLE PRECISION FUNCTION VOLTER( P1 , P2 , P3 , P4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE VOLUME DU TETRAEDRE P1 P2 P3 P4
C -----    OU DETERMINANT( P2-P1 , P3-P1 , P4-P1 ) / 6D0
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C VOLTER : LE DETERMINANT DE LA MATRICE ( P2-P1 , P3-P1 , P4-P1 ) / 6D0
C          LE RESULTAT PEUT ETRE >0. =0. OU <0.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS  MAI 1987
C2345X7..............................................................012
      REAL              P1(3), P2(3), P3(3), P4(3)
      DOUBLE PRECISION  X12, Y12, Z12, X13, Y13, Z13, X14, Y14, Z14
      DOUBLE PRECISION  D, E, F
C
      X12 = P1(1) - P2(1)
      Y12 = P1(2) - P2(2)
      Z12 = P1(3) - P2(3)
C
      X13 = P1(1) - P3(1)
      Y13 = P1(2) - P3(2)
      Z13 = P1(3) - P3(3)
C
      X14 = P1(1) - P4(1)
      Y14 = P1(2) - P4(2)
      Z14 = P1(3) - P4(3)
C
      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13
C
      VOLTER = ( X12 * D + X13 * E + X14 * F ) / 6D0
C
      RETURN
      END
