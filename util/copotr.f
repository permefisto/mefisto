      SUBROUTINE COPOTR ( COSOMM , N , NDIM , COPOIN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES N * ( N + 1 ) / 2 SOMMETS D UNE TRIANGULATION
C -----    REGULIERE DU TRIANGLE RECTANGLE UNITE
C          DU HAUT EN BAS (PAR LIGNES) ET DE LA GAUCHE VERS LA DROITE
C
C ENTREES :
C ---------
C COSOMM : COORDONNEES DES SOMMETS DU TRIANGLE RECTANGLE UNITE
C N      : NOMBRE DE POINTS SUR ( 0 , 1 )
C NDIM   : NOMBRE DE COORDONNEES DES POINTS
C
C SORTIES :
C ---------
C COPOIN : COORDONNEES DES SOMMETS DES TRIANGLES INTERNES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      REAL COSOMM(2,3),COPOIN(NDIM,*)
C
      DELTA1 = 1. / ( N - 1 )
      Y      = 1.
      K      = 0
C
C     I DESIGNE LA LIGNE
      DO 20 I=1,N
         IF( I .EQ. N ) Y = 0.
         X = 0.
C
C        J DESIGNE LA COLONNE
         DO 10 J=1,I
C
C           K DESIGNE LE NO DU SOMMET = (I-1) * I / 2 + J
            K = K + 1
            Z = 1.0 - X - Y
            COPOIN(1,K) = COSOMM(1,1)*Z+COSOMM(1,2)*X+COSOMM(1,3)*Y
            COPOIN(2,K) = COSOMM(2,1)*Z+COSOMM(2,2)*X+COSOMM(2,3)*Y
            X = X + DELTA1
   10    CONTINUE
         Y = Y - DELTA1
   20 CONTINUE
      END
