      SUBROUTINE DETMAT33( A11, A12, A13,
     %                     A21, A22, A23,
     %                     A31, A32, A33,
     %                     DETM33 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETM33 EST LE DETERMINANT DE LA MATRICE A 3x3
C -----

C ENTREE :
C --------
C AIJ    : COEFFICIENT EN POSITION I J DE LA MATRICE A

C SORTIE :
C --------
C DETM33  : DETERMINANT DE LA MATRICE A 3x3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY   Janvier 2017
C23456...............................................................012
      DOUBLE PRECISION  A11, A12, A13,
     %                  A21, A22, A23,
     %                  A31, A32, A33, DETM33

C                       (A22     A23)            (A21     A23)
C     DETM33 = A11 * DET(           ) - A12 * DET(           )  +
C                       (A32     A33)            (A31     A33)
C
C                       (A21     A22)
C            + A13 * DET(           )
C                       (A31     A32)

      DETM33 = A11 * ( A22 * A33 - A23 * A32 )
     %       - A12 * ( A21 * A33 - A31 * A23 )
     %       + A13 * ( A21 * A32 - A31 * A22 )

      RETURN
      END
