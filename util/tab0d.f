      SUBROUTINE TAB0D(I1,I2,I3,A,B,C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : C = TA * B OU  C(I1,I3), A(I2,I1), B(I2,I3)
C ----
C ENTREES:
C --------
C I1  : NOMBRE DE LIGNES DE C ET DE COLONNES DEA
C I2  : NOMBRE DE LIGNES DE A ET DE B
C I3  : NOMBRE DE COLONNES DE C ET DE B
C A   : MATRICE A
C B   : MATRICE B
C
C SORTIE :
C --------
C C   : MATRICE TA * B
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 ET INRIA PARIS  FEVRIER 1980
C ......................................................................
      DOUBLE PRECISION  A(I2,I1), B(I2,I3), C(I1,I3), S
C
      DO 1 J = 1,I3
         DO 2 I = 1 , I1
            S = 0.D0
               DO 3 K = 1 , I2
                  S = S + A(K,I) * B(K,J)
    3          CONTINUE
            C(I,J) = S
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
