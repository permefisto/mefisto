      DOUBLE PRECISION FUNCTION DETMAT( N , A )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU DETERMINANT D'UNE MATRICE
C -----    PAR LA METHODE DE GAUSS PAR PIVOT TOTAL
C
C PARAMETRES D ENTREE :
C ---------------------
C N      : ORDRE DE LA MATRICE
C
C PARAMETRE D ENTREE ET RESULTAT :
C --------------------------------
C A      : MATRICE A REELLE DOUBLE PRECISION D ORDRE  N * N
C          RANGEE SOUS FORME D UN VECTEUR DE LONGUEUR N * N
C          COLONNE PAR COLONNE
C          EN SORTIE LA MATRICE EST FACTORISEE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  FEVRIER 1987
C ......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  A(N,N),P,S
C
C     LA SIGNATURE DU DETERMINANT
      DETMAT = 1.
C
C     LES N-1 ETAPES DE LA FACTORISATION DE GAUSS PIVOT TOTAL
C     =======================================================
      DO 100 K=1,N-1
C
C        RECHERCHE DU PIVOT TOTAL A(I0,J0)
         I0 = 0
         J0 = 0
         P  = 0.D0
         DO 20 I=K,N
            DO 10 J=K,N
               S = ABS( A(I,J) )
               IF( S .GT. P ) THEN
                  P  = S
                  I0 = I
                  J0 = J
               ENDIF
 10         CONTINUE
 20      CONTINUE
C
         IF( P .LE. 0D0 ) THEN
             DETMAT = 0D0
             RETURN
         ENDIF
C
C        LE CALCUL DU DETERMINANT
         DETMAT = DETMAT * A(I0,J0)
C
C        PERMUTATION DES COLONNES K ET J0
         IF( J0 .NE. K ) DETMAT = -DETMAT
         DO 40 I=K,N
            S       = A(I,J0)
            A(I,J0) = A(I,K)
            A(I,K)  = S
 40      CONTINUE
C
C        PERMUTATION DES LIGNES K ET I0
         P = 1D0 / A(I0,K)
         IF( I0 .NE. K ) DETMAT = -DETMAT
         DO 50 J=K,N
            S       = A(I0,J)
            A(I0,J) = A(K,J)
            A(K,J ) = S * P
 50      CONTINUE
C
C        TRAITEMENT DE L'ETAPE K
         K1 = K + 1
         DO 70 I=K1,N
            DO 60 J=K1,N
               A(I,J) = A(I,J) - A(I,K) * A(K,J)
 60         CONTINUE
 70      CONTINUE
100   CONTINUE
C
C     LE DETERMINANT EST DEFINITIVEMENT CALCULE
      DETMAT = DETMAT * A(N,N)
      END
