      SUBROUTINE COPTTE( COSOMM, NBPARE, NDIM, COPOIN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 2((N-1)**2+1) SOMMETS D UNE TETRAEDRISATION
C -----    REGULIERE DU SOUS-TETRAEDRE DE REFERENCE UNITE
C          OU SOUS-TETRAEDRES ISSU DU DECOUPAGE DU TETRAEDRE OU
C          PENTAEDRE OU HEXAEDRE DE REFERENCE UNITE
C
C ENTREES:
C --------
C COSOMM : COORDONNEES DES SOMMETS DU SOUS-TETRAEDRE PRIMAIRE
C NBPARE : NOMBRE DE POINTS SUR CHAQUE ARETE
C NDIM   : NOMBRE DE COORDONNEES DES POINTS + 1
C
C SORTIE :
C --------
C COPOIN : COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES INTERNES
C          DANS L'ELEMENT FINI DE REFERENCE GLOBAL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      REAL     COSOMM(3,4),
     %         COPOIN(NDIM,*)
      INTEGER  NOS1(6),
     %         NOS2(6)
      DATA     NOS1/ 1, 2, 3, 4, 4, 4 /
      DATA     NOS2/ 2, 3, 1, 1, 2, 3 /
C
      DO 1 J = 1, 4
         COPOIN( 1, J ) = COSOMM( 1, J )
         COPOIN( 2, J ) = COSOMM( 2, J )
         COPOIN( 3, J ) = COSOMM( 3, J )
 1    CONTINUE
C
      IF( NBPARE .GE. 3 ) THEN
C
C        CALCUL DU POINT MILIEU DES 6 ARETES
         DO 2 J =1, 6
            JJ = 4 + J
C           LE NUMERO DES 2 SOMMETS DE L'ARETE J
            J1 = NOS1( J )
            J2 = NOS2( J )
            COPOIN(1,JJ) = 0.5 * ( COSOMM(1,J1) + COSOMM(1,J2) )
            COPOIN(2,JJ) = 0.5 * ( COSOMM(2,J1) + COSOMM(2,J2) )
            COPOIN(3,JJ) = 0.5 * ( COSOMM(3,J1) + COSOMM(3,J2) )
 2       CONTINUE
C
      ENDIF
      END
