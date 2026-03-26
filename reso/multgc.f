      SUBROUTINE MULTGC(NTDL,LPLIGC,LPCOLC,AGC,B,X)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: MULTIPLICATION D'UN VECTEUR PAR LA MATRICE L TRANSPOSEE
C ---- AVEC  A |=| L * TL  FACTORISATION DE CHOLESKY INCOMPLETE
C
C      VERSION REELLE DOUBLE PRECISION
C
C PARAMETRES D ENTREE:
C --------------------
C NTDL   : ORDRE DE LA MATRICE ET DU SECOND MEMBRE
C LPLIGC :
C LPCOLC : LA MATRICE MORSE AGC (SYMETRIQUE)
C AGC    :
C B      : TABLEAU B(NTDL) DU  SECOND  MEMBRE
C
C PARAMETRE DE SORTIE:
C --------------------
C X      : TABLEAU X(NTDL) RESULTAT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         OCTOBRE 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AGC,B,X
      DIMENSION LPLIGC(NTDL+1),LPCOLC(1),AGC(1),B(NTDL),X(NTDL)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C
      DO 1 I = 1 , NTDL
         X(I) = B(I) * AGC(LPLIGC(I+1))
 1    CONTINUE
      DO 2 I = NTDL,1,-1
C        I LE NO DE LA COLONNE TRAITEE
         DO 3 K = LPLIGC(I)+1,LPLIGC(I+1)-1
C           J LE NO DE LA LIGNE DU COEFFICIENT AGC(I,J)
            J = LPCOLC(K)
            X(J) = X(J) + AGC(K) * B(I)
 3       CONTINUE
 2    CONTINUE
C
      RETURN
      END
