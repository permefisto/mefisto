      DOUBLE PRECISION FUNCTION DETM44( A11, A12, A13, A14,
     %                                  A21, A22, A23, A24,
     %                                  A31, A32, A33, A34,
     %                                  A41, A42, A43, A44 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU DETM44 EST LE DETERMINANT DE LA MATRICE A 4x4
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
     %                  A41, A42, A43, A44

      CALL DETMAT44( A11, A12, A13, A14,
     %               A21, A22, A23, A24,
     %               A31, A32, A33, A34,
     %               A41, A42, A43, A44,
     %               DETM44 )

      RETURN
      END
