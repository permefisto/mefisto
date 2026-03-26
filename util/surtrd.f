      DOUBLE PRECISION FUNCTION SURTRD( P1 , P2 , P3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE LA SURFACE D'UN TRIANGLE DEFINI PAR 3 POINTS DE R**3
C -----
C PARAMETRES D ENTREE :
C ---------------------
C P1 P2 P3 : LES 3 FOIS 3 COORDONNEES DES SOMMETS DU TRIANGLE
C
C PARAMETRE RESULTAT :
C --------------------
C SURTRD : LA SURFACE >=0 DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS DECEMBRE 1987
C ......................................................................
      DOUBLE PRECISION  P1(3),P2(3),P3(3), V(3,2), VN(3), PROSCD, SQRT
C
C     LES VECTEURS V1 = P2 - P1 , V2 = P3 - P1
C     ----------------------------------------
      DO 10 I=1,3
         V(I,1) = P2(I) - P1(I)
         V(I,2) = P3(I) - P1(I)
 10   CONTINUE
C
C     LE PRODUIT VECTORIEL DES 2 DIRECTIONS V1 V2
C     -------------------------------------------
      CALL PROVEC( V(1,1) , V(1,2) , VN )
C
C     LE PRODUIT SCALAIRE DE LA NORMALE AVEC ELLE MEME
C     ------------------------------------------------
      SURTRD = SQRT( PROSCD( VN , VN , 3 ) ) * 0.5D0
C
      RETURN
      END
