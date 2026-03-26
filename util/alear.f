      SUBROUTINE ALEAR( IC, N, A )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CALCUL DE N NOMBRES ALEATOIRES DE LA FORME : X=MOD(A*X+C,M)
C ----
C          M=2**IE+P
C          IP=-1,0,+1
C          IC PREMIER AVEC M
C          IA-1 MULTIPLE DES NOMBRES PREMIERS DIVISANT M
C
C ENTREES:
C --------
C IC     : NOMBRE PREMIER A CHOISIR
C N      : NOMBRE DE NOMBRES ALEATOIRES A CALCULER
C
C SORTIE :
C --------
C A      : LES N NOMBRES CALCULES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++12
C AUTEUR : DOMINIQUE BEGIS INRIA  VERSAILLES                JANVIER 1979
C ......................................................................
      REAL  A(N)
C
      IE = 15
      M  = 2 ** IE
      IA = 2 ** 8 + 1
      IX = 2 ** 7 - 29
C
      DO 1 I=1,N
         IX   = IX * IA + IC
         IX   = MOD(IX,M)
         A(I) = IX * 3.E-5
    1 CONTINUE
      END
