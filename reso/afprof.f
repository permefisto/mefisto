      SUBROUTINE AFPROF( NP , LPDIAG , A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: AFFICHER LES NP COLONNES DE LA MATRICE PROFIL A REEL2
C ---- LPDIAG TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE A
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   MAI 1989
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           LPDIAG(1:*)
      DOUBLE PRECISION  A(1:*)
C
10006 FORMAT(1X,9(1H-),' D.L ---',I7,2X,105(1H-))
10001 FORMAT(1X,130(1H-)//)
10003 FORMAT(1H0,9(1H-),' MATRICE DIAGONALE DU D.L. 1 AU D.L.',I10,
     &5X,70(1H-))
10015  FORMAT(4(2X,'A(',I6,')=',D14.7))
C
      K=0
      IF(LPDIAG(2).EQ.1) THEN
C
C        MATRICE NON DIAGONALE
         DO 4 I=1,NP
            K1=K+1
            K2=K+LPDIAG(I+1)-LPDIAG(I)
            WRITE(IMPRIM,10006) I
C           REEL DOUBLE PRECISION
            WRITE(IMPRIM,10015) (L,A(L),L=K1,K2)
            K=K2
    4    CONTINUE
         WRITE(IMPRIM,10001)
C
      ELSE
C
C        MATRICE DIAGONALE
         WRITE(IMPRIM,10003) NP
         WRITE(IMPRIM,10015) (L,A(L),L=1,NP)
      ENDIF
      END
