      SUBROUTINE DIASYM( NB, DIAG, SYM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRANSFORMER UNE MATRICE DIAGONALE EN MATRICE SYMETRIQUE
C -----    A L'APPEL SYM PEUT ETRE IDENTIQUE A DIAG
C
C ENTREES:
C --------
C NB     : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE DIAGONALE
C DIAG   : LA MATRICE DIAGONALE STOCKEE SOUS LA FORME D'UN VECTEUR
C
C SORTIE :
C --------
C SYM    : LA MATRICE SYMETRIQUE STOCKEE SOUS LA FORME D'UN VECTEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     FEVRIER 1998
C2345X7..............................................................012
      DOUBLE PRECISION  DIAG(NB), SYM(NB*(NB+1)/2)
C
C     ADRESSE DU DERNIER COEFFICIENT DIAGONAL
      N = NB * (NB+1) / 2 + 1
C
      DO I = NB, 1, -1
C
C        ADRESSE DANS SYM
         N = N - 1
         SYM(N) = DIAG(I)
C
C        COMPLETION DE ZEROS
         DO J = I-1, 1, -1
            N = N - 1
            SYM(N) = 0D0
         ENDDO
C
      ENDDO
C
      RETURN
      END
