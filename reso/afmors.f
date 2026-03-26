      SUBROUTINE AFMORS( NTDL , LPLIGN , LPCOLO , A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: AFFICHER LES COEFFICIENTS DE LA MATRICE MORSE A REEL2
C ----
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PARAMETRES D'ENTREE :
C  -------------------
C NTDL   : ORDRE DE LA MATRICE AG
C LPLIGN : POINTEUR SUR LE DERNIER COEFFICIENT DES LIGNES DE LPCOLO
C LPCOLO : NUMERO DE COLONNE DES VALEURS NON NULLES DE LA MATRICE
C A      :
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS      MAI 1990
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           LPLIGN(0:*),LPCOLO(1:*)
      DOUBLE PRECISION  A(1:*)
C
10000 FORMAT(1X,9(1H-),' D.L ---',I7,2X,105(1H-))
10001 FORMAT(4(2X,'A(',I6,')=',D14.7))
10002 FORMAT(1X,130(1H-)//)
C
       J1 = 1
       DO 100 I=1,NTDL
          J2 = LPLIGN(I)
          WRITE(IMPRIM,10000) I
C         REEL DOUBLE PRECISION
          WRITE(IMPRIM,10001) (LPCOLO(J),A(J),J=J1,J2)
          J1 = J2 + 1
100    CONTINUE
       WRITE(IMPRIM,10002)
C
       RETURN
       END
