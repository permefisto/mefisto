      SUBROUTINE DRCHGC(NTDL,LPLIGC,LPCOLC,AGC,B,X)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     DESCENTE ET REMONTEE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ----     LA FORME  A * X = L * TL * X = B  DE CHOLESKY INCOMPLETE
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE ET DU SECOND MEMBRE
C LPLIGC : LE POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C LPCOLC : LE NUMERO DE LA COLONNE DES COEFFICIENTS
C AGC    : LA MATRICE MORSE AGC (SYMETRIQUE)
C B      : TABLEAU B(NTDL) DU  SECOND  MEMBRE
C
C SORTIES:
C --------
C X      : TABLEAU X(NTDL) SOLUTION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         OCTOBRE 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AGC(*), B(NTDL), X(NTDL), S
      INTEGER           LPLIGC(NTDL+1), LPCOLC(1)
C
C     LA  DESCENTE
C     -------------
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
C     LA REMONTEE
C     -----------
      DO 3 I = NTDL,1,-1
C        I LE NO DE LA COLONNE TRAITEE
         S = X(I) / AGC(LPLIGC(I+1))
         X(I) = S
         DO 4 K = LPLIGC(I)+1,LPLIGC(I+1)-1
C           J LE NO DE LA LIGNE DU COEFFICIENT AGC(I,J)
            J = LPCOLC(K)
            X(J) = X(J) - AGC(K) * S
 4       CONTINUE
 3    CONTINUE
C
      RETURN
      END
