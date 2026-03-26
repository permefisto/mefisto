      SUBROUTINE DECHGC ( NTDL, LPLIGC, LPCOLC, AGC, B, X )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: DESCENTE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ---- LA FORME  A * X = L * TL * X = B  DE CHOLESKY INCOMPLETE
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
C X      : TABLEAU X(NTDL) SOLUTION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AGC,B,X,S
      DIMENSION LPLIGC(NTDL+1),LPCOLC(1),AGC(1),B(NTDL),X(NTDL)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     LA  DESCENTE
C     -------------
C
      DO 1 I = 1 , NTDL
C        I LE NO DE LA LIGNE TRAITEE
         S  = B(I)
         DO 2 K = LPLIGC(I)+1,LPLIGC(I+1)-1
C           J LE NO DE LA LIGNE DU COEFFICIENT AGC(I,J)
            J = LPCOLC(K)
            S = S - AGC(K) * X(J)
 2       CONTINUE
         X(I) = S / AGC(LPLIGC(I+1))
 1    CONTINUE
C
      RETURN
      END
