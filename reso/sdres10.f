       SUBROUTINE SDRES10(Q,LU,N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FACTORISATION LU DE LA MATRICE Q TRIDIAGONALE
C -----
C        Q(1,I) = SOUS-DIAGONALE PRINCIPALE
C        Q(2,I) = DIAGONALE
C        Q(3,I) = SUR-DIAGONALE PRINCIPALE
C
C ENTREE :
C ---------
C Q  : LA MATRICE A FACTORISER
C
C SORTIE :
C --------
C LU : LES FACTEURS DE GAUSS
C
C REMARQUES :
C      1) LA DIAGONALE DE L N'EST PAS STOCKEE
C      2) ON PEUT ECRASER Q PAR LU
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
C
      DOUBLE PRECISION  Q(3,N),LU(3,N)
C
C     LA PREMIERE LIGNE
      LU(1,1) = 0.D0
      LU(2,1) = Q(2,1)
      LU(3,1) = Q(3,1)
      DO 1 I  = 2 , N
         LU(1,I) = Q(1,I) / LU(2,I-1)
         LU(2,I) = Q(2,I) - LU(1,I) * LU(3,I-1)
         LU(3,I) = Q(3,I)
 1    CONTINUE
C     LA DERNIERE LIGNE
      LU(3,N) = 0.D0
C
      RETURN
      END
