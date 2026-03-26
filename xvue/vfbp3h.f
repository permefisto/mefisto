        SUBROUTINE VFBP3H( U, FB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA VALEUR AU POINT U DE L'INTERVALLE UNITE DES
C -----    4 FONCTIONS DE BASE DE L'INTERPOLATION C1 P3 HERMITE CROISEE
C
C  FB (U=0)(1) = 1   FB (U=0)(2) = 0   FB (U=0)(3) = 0   FB (U=0)(4) = 0
C  FB (U=1)(1) = 0   FB (U=1)(2) = 1   FB (U=1)(3) = 0   FB (U=1)(4) = 0
C  FB'(U=0)(1) = 0   FB'(U=0)(2) = 0   FB'(U=0)(3) = 1   FB'(U=0)(4) = 0
C -FB'(U=1)(1) = 0  -FB'(U=1)(2) = 0  -FB'(U=1)(3) = 0  -FB'(U=1)(4) = 1
C
C ENTREES:
C --------
C U      : VALEUR DU PARAMETRE SUR [0,1]
C
C SORTIES:
C --------
C FB     : VALEUR AU POINT U DES 4 FONCTIONS DE BASE
C          POLYNOMES P3 EN U SELON L'ORDRE K DES DEGRES DE LIBERTE
C          K=1  F(0),       K=2   F(1),
C          K=3 DF(0)(1-0),  K=4  DF(1)(0-1)=-DF(1)=-FB'(U=1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C2345X7..............................................................012
      REAL  U, FB(4)
C
      UU  = U  * U
      UUU = UU * U
C
C     LE POLYNOME P3 ASSOCIE AU SOMMET S1 = 0
      FB(1) = 1 - 3 * UU + 2 * UUU
C
C     LE POLYNOME P3 ASSOCIE AU SOMMET S2 = 1
      FB(2) =     3 * UU - 2 * UUU
C
C     LE POLYNOME P3 ASSOCIE A LA DERIVEE(0)(1-0)
      FB(3) = U - 2 * UU +     UUU
C
C     LE POLYNOME P3 ASSOCIE A LA DERIVEE(1)(0-1)=-DERIVEE(1)
      FB(4) =         UU -     UUU
C
      RETURN
      END
