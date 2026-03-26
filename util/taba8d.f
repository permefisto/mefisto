       SUBROUTINE TABA8D(I1,I2,A,B,C,AUX)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CALCUL DU PRODUIT MATRICIEL TRANSPOSE(A)*B*A DANS LE CAS OU
C ---  B EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE
C      C=TRANSPOSEE(A)*B*A
C
C ENTREES:
C --------
C I2     : NOMBRE DE LIGNES  DE A , DE LIGNES ET COLONNES DE B
C I1     : NOMBRE DE COLONNES DE A
C A      : MATRICE I2,I1
C B      : MATRICE (I2,I2) SYMETRIQUE STOCKEE DE HAUT EN BAS ET
C          DE LA GAUCHE VERS LA DROITE.LA PARTIE TRIANGULAIRE
C          SUPERIEURE EST STOCKEE
C AUX    : MATRICE AUXILIAIRE I1 I2 .CONTIENT LE PRODUIT
C          TRANSPOSEE(A)*B.PERMET DE CALCULER LE PRODUIT FINAL
C          EN 2*O(N**3) OPERATIONS
C
C SORTIE :
C --------
C C      : PRODUIT MATRICIEL TRANSPOSE(A)*B*A  C(I1,I1)
C          C EST SYMETRIQUE STOCKEE DE HAUT EN BAS ET DE LA GAUCHE
C          VERS LA DROITE . SEULE LA PARTIE TRIANGULAIRE
C          SUPERIEURE EST STOCKEE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
       DOUBLE PRECISION  S, B(*), A(I2,I1), C(*), AUX(I1,I2)
C
C     AUX=T(A)*B
C     ----------
      CALL TABSZD( I1, I2, A, B, AUX )
C
C     C=AUX*A
C     -------
      N = 0
      DO 10 J = 1 , I1
         DO 5 I = 1 , J
            S = 0.D0
            DO 2 K = 1 , I2
                S = S + AUX(I,K) * A(K,J)
 2          CONTINUE
            N    = N + 1
            C(N) = S
 5       CONTINUE
 10   CONTINUE
      RETURN
      END
