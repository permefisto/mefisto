      SUBROUTINE GAUSPT( N , M , A , IER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RESOUDRE UN SYSTEME LINEAIRE A * X  =  B DE M SECONDS MEMBRES
C ----- PAR LA METHODE DE GAUSS A PIVOT TOTAL
C
C ENTREES :
C ---------
C N   : ORDRE  DE LA MATRICE
C M   : NOMBRE DE SECONDS MEMBRES
C
C ENTREE ET RESULTAT :
C --------------------
C A   : TABLEAU RECTANGULAIRE DE N LIGNES ET N(+M) COLONNES
C       LA MATRICE A EST DEFINIE DANS LES N PREMIERES COLONNES
C       LES M SECONDS MEMBRES B   SONT DANS LES COLONNES N+1 A N+M DE A
C       EN SORTIE LES M SOLUTIONS SONT DANS LES COLONNES N+1 A N+M DE A
C                 LES N PREMIERES COLONNES CONTIENNENT LA FACTORISEE LU
C
C SORTIE :
C --------
C IER : 0 MATRICE INVERSE CORRECTEMENT CALCULEE
C       NON 0 => MATRICE NON INVERSIBLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  FEVRIER 1987
C ......................................................................
      DOUBLE PRECISION  A(N,*),P,S
C
C     LES N-1 ETAPES DE LA FACTORISATION DE GAUSS PIVOT TOTAL
C     =======================================================
      N1 = N + 1
      NM = N + M
C
C     K=N POUR EFFECTUER LA DIVISION PAR LE PIVOT A(N,N)
C         PROPRIETE EXPLOITEE LORS DE LA REMONTEE
C
      DO 100 K=1,N
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
             IER = 1
             RETURN
         ENDIF
C
C        PERMUTATION DES COLONNES K ET J0
         IF( K .NE. J0 ) THEN
            DO 40 I=K,N
               S       = A(I,J0)
               A(I,J0) = A(I,K)
               A(I,K)  = S
 40         CONTINUE
         ENDIF
C
C        PERMUTATION DES LIGNES K ET I0
         P = 1D0 / A(I0,K)
         DO 50 J=K,NM
            S       = A(I0,J)
            A(I0,J) = A(K,J)
            A(K,J ) = S * P
 50      CONTINUE
C
C        LE COEFFICIENT A(K,K)=1. .IL EST ECRASE PAR J0
C        NUMERO DE LA COLONNE NECESSAIRE LORS DES REMONTEES
         A(K,K) = J0
C
C        MISE A 0. DES A(I,K) POUR I=K+1 A N
         K1 = K + 1
         DO 70 I=K1,N
            DO 60 J=K1,NM
               A(I,J) = A(I,J) - A(I,K) * A(K,J)
 60         CONTINUE
 70      CONTINUE
100   CONTINUE
C
C     LES M REMONTEES DU SYSTEME
C     ==========================
      DO 200  I=N-1,1,-1
         I1 = I + 1
C
C        LE NUMERO J0 DE COLONNE PERMUTEE AVEC LA COLONNE I
         J0 = NINT( A(I,I) )
C
C        LA BOUCLE SUR LES M SECONDS MEMBRES
         DO 120 K=N1,NM,1
C
            S = A(I,K)
            DO 110 J=I1,N,1
               S = S - A(I,J) * A(J,K)
 110        CONTINUE
C
C           EN FAIT S=X(J0,K-N) DOIT ETRE PERMUTE AVEC X(I,K-N)
C           C-A-D A(J0,K) AVEC A(I,K)
            A(I ,K) = A(J0,K)
            A(J0,K) = S
 120     CONTINUE
 200  CONTINUE
C
C     RESOLUTION CORRECTE DU SYSTEME LINEAIRE
      IER = 0
      END
