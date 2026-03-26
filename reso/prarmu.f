      SUBROUTINE PRARMU( Seuil, Coef, NBCOPG, PG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PG = PG * Coef AVEC MISE A ZERO SI VALEUR ABSOLUE TROP PETITE
C -----    TENTATIVE DE REDUIRE LES ERREURS D'ARRONDIS
C
C ENTREES:
C --------
C Seuil  : SEUIL AU DESSOUS DUQUEL LA VALEUR EST MISE A ZERO
C Coef   : VALEUR DU COEFFICIENT MULTIPLICATEUR
C NBCOPG : NOMBRE DE COEFFICIENTS DU TABLEAU PG
C
C MODIFIE:
C --------
C PG     : PG = PG * Coef AVEC MISE A ZERO SI VALEUR ABSOLUE TROP PETITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION   Seuil, Coef, PG(NBCOPG)
      INTEGER            NBCOPG, K
      INTRINSIC          ABS
C
      DO K = 1, NBCOPG
         IF( ABS( PG( K ) ) .LT. Seuil ) THEN
            PG(K) = 0D0
         ELSE
            PG(K) = PG( K ) * Coef
         ENDIF
      ENDDO
C
      RETURN
      END
