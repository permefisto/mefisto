        SUBROUTINE XYP3H( U, X, Y,  XP3H, YP3H )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DE L'INTERPOLATION P3 HERMITE
C ----- AU POINT U DE [0,1]
C
C ENTREES:
C --------
C U      : VALEUR DU PARAMETRE SUR [0,1]
C X, Y   : VALEURS DES 4 DEGRES DE LIBERTE RANGEES SELON K POUR X Y
C          K=1  F(S1),         K=2  F(S2),
C          K=3 DF(S1)(S2-S1),  K=4  DF(S2)(S1-S2)
C SORTIES:
C --------
C XP3H : ABSCISSE P3H AU POINT U DE [0,1]
C YP3H : ORDONNEE P3H AU POINT U DE [0,1]
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C2345X7..............................................................04
      REAL    U, X(4), Y(4), FB(4)
C
C     VALEUR DES 4 FONCTIONS DE BASE P3H AU POINT U
      CALL VFBP3H( U, FB )
C
C     LES COORDONNEES X Y Z POUR LE POINT U DE L'INTERPOLATION P3H
      XP3H = 0.0
      YP3H = 0.0
      DO 10 I=1,4
         XP3H = XP3H + FB(I) * X(I)
         YP3H = YP3H + FB(I) * Y(I)
 10   CONTINUE
      END
