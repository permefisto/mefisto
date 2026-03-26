       SUBROUTINE SDRES9(V,A,X,B,Y,N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : COMBINAISON LINEAIRE V = A * X + B * Y
C -----
C
C ENTREES :
C ---------
C X , Y  : VECTEURS DE LA COMBINAISON LINEAIRE
C A , B  : COEFFICIENTS DE LA COMBINAISON LINEAIRE
C N      : LONGUEUR DES VECTEURS
C
C SORTIES :
C ---------
C V      : VECTEUR RESULTAT (PEUT ETRE CONFONDU AVEC X OU Y)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
C
      DOUBLE PRECISION  A,B,X(N),Y(N),V(N)
C
      DO 1 I  = 1 , N
         V(I) = A * X(I) + B * Y(I)
 1    CONTINUE
C
      RETURN
      END
