        SUBROUTINE TGVHCT( U, V, X, Y, Z,  TGXHCT, TGYHCT, TGZHCT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DES 3 COMPOSANTES DE LA DERIVEE EN V
C ----- DE L'INTERPOLATION AU POINT (U,V)
C       DU TRIANGLE RECTANGLE UNITE HSIEH-CLOUGH-TOCHER REDUIT
C       (REDUIT => DERIVEE NORMALE P1 SUR LES 3 COTES)
C       ET POUR LES 9 DEGRES DE LIBERTE STOCKES DANS X Y Z
C       CF ARTICLE DE M. BERNARDOU et K.HASSAN (1981)
C              ^
C              | V
C              |
C             S3=(0,1)
C              X
C              | \ \
C              |  \   ^
C              -   \    \
C              |    \     \
C              |     \      \
C              | E2   \       \
C              |       \    E1  \
C              ^     /   \        ^
C              |   /         \      \
C              | /      E3        \   \
C              X----->------------<-----X ----> U
C             S1=(0,0)                S2=(1,0)
C
C ENTREES:
C --------
C U, V   : VALEURS DES 2 PARAMETRES SUR LE TRIANGLE RECTANGLE UNITE
C          ( 0=<U,V<=1  ET V<=1-U )
C X, Y, Z: VALEURS DES 9 DEGRES DE LIBERTE RANGEES SELON K POUR X, Y, Z
C K=1  F(S1),        K=2  F(S2),        K=3  F(S3),
C K=4 DF(S1)(S2-S1), K=5 DF(S1)(S3-S1),
C K=6 DF(S2)(S3-S2), K=7 DF(S2)(S1-S2),
C K=8 DF(S3)(S1-S3), K=9 DF(S3)(S2-S3)
C
C SORTIES:
C --------
C TGXHCT : ABSCISSE AU POINT (U,V) DE LA TG
C TGYHCT : ORDONNEE AU POINT (U,V) DE LA TG
C TGZHCT : COTE     AU POINT (U,V) DE LA TG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1998
C2345X7..............................................................012
      REAL    U, V, X(9), Y(9), Z(9), DERVFB(9)
C
C     VALEUR DE LA DERIVEE SELON V DES 9 FONCTIONS DE BASE HCT AU POINT (U,V)
      CALL DDVHCT( U, V, DERVFB )
C
C     LES COMPOSANTES DE LA TG SELON V AU POINT (U,V) DE L'INTERPOLATION HCT
      TGXHCT = 0.0
      TGYHCT = 0.0
      TGZHCT = 0.0
      DO 10 I=1,9
         TGXHCT = TGXHCT + DERVFB(I) * X(I)
         TGYHCT = TGYHCT + DERVFB(I) * Y(I)
         TGZHCT = TGZHCT + DERVFB(I) * Z(I)
 10   CONTINUE
      END
