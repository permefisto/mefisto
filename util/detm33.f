      DOUBLE PRECISION FUNCTION DETM33( A11, A12, A13,
     %                                  A21, A22, A23,
     %                                  A31, A32, A33 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU DETERMINANT DE LA MATRICE 3x3 A
C -----       A11 A12 A13
C             A21 A22 A23
C             A31 A32 A33
C ENTREES:
C --------
C A11,A12,A13 : LIGNE 1 DE LA MATRICE A
C A21,A22,A23 : LIGNE 2 DE LA MATRICE A
C A31,A32,A33 : LIGNE 3 DE LA MATRICE A

C SORTIES:
C --------
C DETM33 : DETERMINANT DE LA MATRICE 3x3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    DECEMBRE 1992
C2345X7..............................................................012
      DOUBLE PRECISION A11, A12, A13,
     %                 A21, A22, A23,
     %                 A31, A32, A33

      DETM33 = A11 * ( A22 * A33 - A23 * A32 )
     %       + A21 * ( A13 * A32 - A33 * A12 )
     %       + A31 * ( A12 * A23 - A22 * A13 )

      RETURN
      END
