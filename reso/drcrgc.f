      SUBROUTINE DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, B, X )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     DESCENTE ET REMONTEE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ----     LA FORME  A * X = L * D * TL * X = B  DE CROUT INCOMPLETE
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
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris OCTOBRE 2007
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AGC(*), B(NTDL), X(NTDL), S
      INTEGER           LPLIGC(NTDL+1), LPCOLC(*)
C
C     LA  DESCENTE
C     ------------
      DO 2 I = 1, NTDL
C        I NO DE LA LIGNE TRAITEE
         S = B(I)
         DO 1 K = LPLIGC(I)+1, LPLIGC(I+1)-1
C           J NO DE LA COLONNE DU COEFFICIENT AGC(K)
            J = LPCOLC(K)
            S = S - AGC(K) * X(J)
 1       CONTINUE
         X(I) = S
 2    CONTINUE
c
ccc      print *,'fin descente'
ccc      do 10 i=1,ntdl
ccc         print *,'x',i,'=',x(i)
ccc 10   continue
C
C     LA DIVISION PAR LA DIAGONALE
C     ----------------------------
      DO 5 I = 1, NTDL
C        I NO DE LA LIGNE TRAITEE
         X(I) = X(I) / AGC( LPLIGC(I+1) )
 5    CONTINUE
c
ccc      print *,'fin division diagonale'
ccc      do 20 i=1,ntdl
ccc         print *,'x',i,'=',x(i)
ccc 20   continue
C
C     LA REMONTEE
C     -----------
      DO 7 I = NTDL, 1, -1
C        I NO DE LA COLONNE TRAITEE
         S = X(I)
         DO 6 K = LPLIGC(I)+1, LPLIGC(I+1)-1
C           J NO DE LA LIGNE DU COEFFICIENT AGC(K)
            J = LPCOLC(K)
            X(J) = X(J) - AGC(K) * S
 6       CONTINUE
 7    CONTINUE
c
ccc      print *,'fin remontee'
ccc      do 30 i=1,ntdl
ccc         print *,'x',i,'=',x(i)
ccc 30   continue
C
      RETURN
      END
