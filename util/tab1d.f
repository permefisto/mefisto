      SUBROUTINE TAB1D(I1,I2,I3,A,B,C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : C = C + TA * B OU  C(I1,I3),A(I2,I1),B(I2,I3)
C ----
C PARAMETRES D ENTREE :
C ---------------------
C I1  : NOMBRE DE LIGNES DE C ET DE COLONNES DE A
C I2  : NOMBRE DE LIGNES DE A ET DE B
C I3  : NOMBRE DE COLONNES DE C ET DE B
C A   : MATRICE A
C B   : MATRICE B
C
C PARAMETRE RESULTAT :
C --------------------
C C   : MATRICE C + TA * B
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 ET INRIA PARIS           FEVRIER 1980
C ......................................................................
      DOUBLE PRECISION A,B,C,S
      DIMENSION A(I2,I1),B(I2,I3),C(I1,I3)
C
            DO 1 J = 1,I3
                 DO 2 I = 1 , I1
                 S = C(I,J)
                      DO 3 K = 1 , I2
                      S = S + A(K,I) * B(K,J)
    3                 CONTINUE
                 C(I,J) = S
    2            CONTINUE
    1       CONTINUE
      RETURN
      END
