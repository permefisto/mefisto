      REAL FUNCTION SURTRR( P1 , P2 , P3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE LA SURFACE D'UN TRIANGLE DEFINI PAR 3 POINTS DE R**3
C -----
C ENTREES :
C ---------
C P1 P2 P3 : LES 3 FOIS 3 COORDONNEES DES SOMMETS DU TRIANGLE
C
C SORTIES :
C ---------
C SURTRR : SURFACE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS         DECEMBRE 1987
C23456---------------------------------------------------------------012
      REAL     P1(3),  P2(3), P3(3)
      REAL     V(3,2), VN(3)
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
      CALL PROVER( V(1,1) , V(1,2) , VN )
C
C     LE PRODUIT SCALAIRE DE LA NORMALE AVEC ELLE MEME
C     ------------------------------------------------
      SURTRR = PROSCR( VN , VN , 3 )
      SURTRR = SQRT( SURTRR ) * 0.5
      END
