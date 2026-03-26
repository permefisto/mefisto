      SUBROUTINE PROVEC( V1 , V2 , V3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    V3 VECTEUR = PRODUIT VECTORIEL DE 2 VECTEURS DE R ** 3
C -----
C ENTREES:
C --------
C V1, V2 : LES 2 VECTEURS DE 3 COMPOSANTES
C
C SORTIE :
C --------
C V3     : VECTEUR = V1  PRODUIT VECTORIEL V2
CC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        MARS 1987
C2345X7..............................................................012
      DOUBLE PRECISION    V1(3), V2(3), V3(3)
C
      V3( 1 ) = V1( 2 ) * V2( 3 ) - V1( 3 ) * V2( 2 )
      V3( 2 ) = V1( 3 ) * V2( 1 ) - V1( 1 ) * V2( 3 )
      V3( 3 ) = V1( 1 ) * V2( 2 ) - V1( 2 ) * V2( 1 )
C
      RETURN
      END
