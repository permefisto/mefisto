      SUBROUTINE ATB0D( I1 , I2 , I3 , A , B , C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                   S.P. ATB0D
C                   ----------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU PRODUIT MATRICIEL C = A * TRANSPOSEE( B )
C -----
C PARAMETRES D ENTREE :
C ----------------------
C I1  : NBRE DE LIGNES DE A
C I2  : NBRE DE COLONNES DE A ET DE B
C I3  : NBRE DE LIGNES DE B
C A   : MATRICE I1,I2
C B   : MATRICE I3,I2
C
C PARAMETRE RESULTAT :
C ----------------------
C C   : C = A * TRANSPOSEE ( B ) MATRICE I1,I3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET LAN189 PARIS ET INRIA     AVRIL 1982
C ......................................................................
      DOUBLE PRECISION S,A(I1,I2),B(I3,I2),C(I1,I3)
C
      DO 30 J=1,I3
         DO 20 I=1,I1
            S = 0.D0
            DO 10 K=1,I2
               S = S + A( I , K ) * B( J , K )
   10       CONTINUE
            C( I , J ) = S
   20    CONTINUE
   30 CONTINUE
      END
