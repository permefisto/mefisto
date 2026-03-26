      SUBROUTINE DETMAT44( A11, A12, A13, A14,
     %                     A21, A22, A23, A24,
     %                     A31, A32, A33, A34,
     %                     A41, A42, A43, A44,
     %                     DETM44 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETM44 EST LE DETERMINANT DE LA MATRICE A 4x4
C -----

C ENTREE :
C --------
C AIJ    : COEFFICIENT EN POSITION I J DE LA MATRICE A

C SORTIE :
C --------
C DETM44 : DETERMINANT DE LA MATRICE A 4x4
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY   Janvier 2017
C23456...............................................................012
      DOUBLE PRECISION  A11, A12, A13, A14,
     %                  A21, A22, A23, A24,
     %                  A31, A32, A33, A34,
     %                  A41, A42, A43, A44, DETM44
      DOUBLE PRECISION  DETM331, DETM332, DETM333, DETM334

C                       (A22 A23 A24)            (A21 A23 A24)
C     DETM44 = A11 * DET(A32 A33 A34) - A12 * DET(A31 A33 A34)  +
C                       (A42 A43 A44)            (A41 A43 A44)
C
C                       (A21 A22 A24)            (A21 A22 A23)
C            + A13 * DET(A31 A32 A34) - A14 * DET(A31 A32 A33)
C                       (A41 A42 A44)            (A41 A42 A43)

      CALL DETMAT33( A22, A23, A24,
     %               A32, A33, A34,
     %               A42, A43, A44, DETM331 )

C     LE SIGNE - EST PRIS EN COMPTE PAR LA PERMUTATION DE 2 COLONNES
      CALL DETMAT33( A23, A21, A24,
     %               A33, A31, A34,
     %               A43, A41, A44, DETM332 )

      CALL DETMAT33( A21, A22, A24,
     %               A31, A32, A34,
     %               A41, A42, A44, DETM333 )

C     LE SIGNE - EST PRIS EN COMPTE PAR LA PERMUTATION DE 2 COLONNES
      CALL DETMAT33( A22, A21, A23,
     %               A32, A31, A33,
     %               A42, A41, A43, DETM334 )

C     LE DETERMINANT DE LA MATRICE A 4x4
      DETM44 = A11 * DETM331 + A12 * DETM332
     %       + A13 * DETM333 + A14 * DETM334

      RETURN
      END
