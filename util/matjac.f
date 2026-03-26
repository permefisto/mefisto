      SUBROUTINE MATJAC( NDIM, NBPOLY, X, DP, DF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : A PARTIR DES DERIVEES DES POLYNOMES EN UN POINT DE L APPLICATION
C ----- F:ELEMENT DE REFERENCE=>ELEMENT COURANT ET DES COORDONNEES DES
C       POINTS,     CALCUL DE LA MATRICE JACOBIENNE DF EN CE POINT
C
C ENTREES :
C ---------
C NDIM   : DIMENSION DE L ESPACE (2 OU 3)
C NBPOLY : NOMBRE DE POLYNOMES
C X      : X(I,J) J-EME COORDONNEES DU I-EME POINT
C DP     : DERIVEES DES POLYNOMES AU POINT
C
C SORTIE :
C --------
C DF     : MATRICE JACOBIENNE DE LA TRANSFORMATION
C          DF(I,J) = SOMME   DPK/DXI * XKJ
C                  K=1,...,NBPOLY
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DF(NDIM,NDIM),
     %                  DP(NDIM,NBPOLY),
     %                  S
      REAL              X(NBPOLY,NDIM)
C
      DO 30 J=1,NDIM
         DO 20 I=1,NDIM
            S = 0D0
            DO 10 K=1,NBPOLY
               S = S + DP(I,K) * X(K,J)
   10       CONTINUE
            DF(I,J) = S
   20    CONTINUE
   30 CONTINUE
      END
