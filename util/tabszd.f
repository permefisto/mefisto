      SUBROUTINE TABSZD(I1,I2,A,B,C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:   CALCUL DU PRODUIT MATRICIEL TRANSPOSEE(A)*B DANS LE CAS OU
C ---    B EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE
C        C=TRANSPOSEE(A)*B
C
C PARAMETRES D ENTREE
C --------------------
C  I1   : NOMBRE DE COLONNES  DE A
C  I2   : NOMBRE DE LIGNES DE A ,DE LIGNES ET COLONNES DE B
C  A    : MATRICE I1,I2
C  B    : MATRICE SYMETRIQUE STOCKEE DE HAUT EN BAS ET
C        DE LA GAUCHE VERS LA DROITE.LA PARTIE TRIANGULAIRE
C        SUPERIEURE EST STOCKEE
C
C PARAMETRES DE SORTIE
C --------------------
C  C    : PRODUIT MATRICIEL A*B  C(I1,I2)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
       DOUBLE PRECISION A,B,C,S
       DIMENSION B(1),A(I2,I1),C(I1,I2)
C
       JC = 0
       DO 1 J=1,I2
          JJJ = JC + J
          J1  = J  + 1
          DO 2 I=1,I1
             S = 0D0
             DO 3 K=1,J
                S = S + B(JC + K) * A(K,I)
    3        CONTINUE
C
             IF(J .EQ. I2) GOTO 5
C
             JK = JJJ
             DO 4 K=J1,I2
                JK = JK + K - 1
                S  = S + B(JK) * A(K,I)
    4        CONTINUE
    5        C(I,J) = S
    2     CONTINUE
          JC = JJJ
    1   CONTINUE
        END
