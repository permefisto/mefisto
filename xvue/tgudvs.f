        SUBROUTINE TGUDVS( U, V, X, Y, Z,  XTGDVS, YTGDVS, ZTGDVS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DE LA DERIVEE EN U DE L'INTERPOLATION
C ----- AU POINT (U,V) DU CARRE UNITE DE Fraeijs de VEUBEKE SANDER REDUIT
C       (REDUIT => DERIVEE NORMALE P1 SUR LES 3 COTES)
C       ET POUR LES 12 DEGRES DE LIBERTE STOCKES DANS X Y Z
C
C              ^ V
C              |
C              |
C              |
C             S4=(0,1)               S3=(1,1)
C              X----->------------<---X
C              | \     p1+p2+p3     / |
C              |   \              /   |
C              -     \          /     -
C              |       \      /       |
C              |          \ /         |
C              | p1+p3   /   \  p1+p2 |
C              |       /       \      |
C              ^     /           \    ^
C              |   /               \  |
C              | /        p1         \|
C              X----->----------<-----X  ------------> U
C             S1=(0,0)               S2=(1,0)
C
C ENTREE :
C --------
C U,V    : VALEURS DES 2 PARAMETRES SUR LE TRIANGLE RECTANGLE UNITE
C          ( 0=<U,V<=1  ET V<=1-U )
C X, Y, Z: VALEURS DES 12 DEGRES DE LIBERTE RANGEES SELON K POUR X Y Z
C
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C
C SORTIES:
C --------
C XTGDVS : ABSCISSE DE LA TG SELON U AU POINT (U,V) DU CARRE UNITE
C YTGDVS : ORDONNEE DE LA TG SELON U AU POINT (U,V) DU CARRE UNITE
C ZTGDVS : COTE     DE LA TG SELON U AU POINT (U,V) DU CARRE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1998
C2345X7..............................................................012
      REAL    U, V, X(12), Y(12), Z(12), DERUFB(12)
C
C     VALEUR DES 12 FONCTIONS DE BASE DVS AU POINT (U,V)
      CALL DDUDVS( U, V, DERUFB )
C
C     LES COORDONNEES X Y Z POUR LE POINT (U,V) DE L'INTERPOLATION DVS
      XTGDVS = 0.0
      YTGDVS = 0.0
      ZTGDVS = 0.0
      DO 10 I=1,12
         XTGDVS = XTGDVS + DERUFB(I) * X(I)
         YTGDVS = YTGDVS + DERUFB(I) * Y(I)
         ZTGDVS = ZTGDVS + DERUFB(I) * Z(I)
 10   CONTINUE
      END
