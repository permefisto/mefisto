      SUBROUTINE AB1D(I,K,J,A,B,C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU PRODUIT MATRICIEL C = C + A * B
C -----
C ENTREES:
C --------
C I   : NBRE DE LIGNES DE A
C J   : NBRE DE COLONNES DE B
C K   : NBRE DE COLONNES DE A,ET DE LIGNES DE B
C A   : MATRICE I,K
C B   : MATRICE K,J
C
C SORTIES:
C --------
C C   : C = C + A * B
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1979
C23456---------------------------------------------------------------012
      DOUBLE PRECISION A(I,K), B(K,J), C(I,J), S
C
      DO 1 L=1,I
         DO 2 M=1,J
            S = C(L,M)
            DO 3 N=1,K
               S = S + A(L,N) * B(N,M)
    3       CONTINUE
            C(L,M) = S
    2    CONTINUE
    1 CONTINUE
      END
