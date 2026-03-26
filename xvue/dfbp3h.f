        SUBROUTINE DFBP3H( U, DFB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA VALEUR AU POINT U DE L'INTERVALLE UNITE
C -----    DE LA DERIVEE DES 4 FONCTIONS DE BASE (POLYNOMES P3)
C          DE L'INTERPOLATION C1 P3 HERMITE CROISEE
C
C  FB (U=0)(1) = 1   FB (U=0)(2) = 0   FB (U=0)(3) = 0   FB (U=0)(4) = 0
C  FB (U=1)(1) = 0   FB (U=1)(2) = 1   FB (U=1)(3) = 0   FB (U=1)(4) = 0
C  FB'(U=0)(1) = 0   FB'(U=0)(2) = 0   FB'(U=0)(3) = 1   FB'(U=0)(4) = 0
C -FB'(U=1)(1) = 0  -FB'(U=1)(2) = 0  -FB'(U=1)(3) = 0  -FB'(U=1)(4) = 1
C
C ENTREE :
C --------
C U      : VALEUR DU PARAMETRE SUR [0,1]
C
C SORTIE :
C --------
C DFB    : VALEUR DE LA DERIVEE AU POINT U DES 4 FONCTIONS DE BASE
C          POLYNOMES P3 EN U SELON L'ORDRE K DES DEGRES DE LIBERTE
C          K=1  F(0),       K=2   F(1),
C          K=3 DF(0)(1-0),  K=4  DF(1)(0-1)=-DF(1)=-FB'(1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      REAL  DFB(4)
C
      UU = U  * U
C
C     LA DERIVEE DU POLYNOME P3 ASSOCIE AU SOMMET S1 = 0
      DFB(1) = 6 * ( UU - U )
C
C     LA DERIVEE DU POLYNOME P3 ASSOCIE AU SOMMET S2 = 1
      DFB(2) = -DFB(1)
C
C     LA DERIVEE DU POLYNOME P3 ASSOCIE A LA DERIVEE(0)(1-0)
      DFB(3) = 1 - 4 * U + 3 * UU
C
C     LA DERIVEE DU POLYNOME P3 ASSOCIE A LA DERIVEE(1)(0-1)=-DERIVEE(1)
      DFB(4) = 2 * U - 3 * UU
C
      RETURN
      END
