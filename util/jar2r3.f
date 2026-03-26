      SUBROUTINE JAR2R3( DGL , DELTA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NORME EUCLIDIENNE DU PRODUIT VECTORIEL DE 2 VECTEURS
C -----   DGL(1,.) ET DGL(2,.)
C         (JACOBIEN D UNE TRANSFORMATION DE R ** 2 -> R ** 3 )
C
C ENTREE :
C --------
C DGL    : LES 3 COORDONNEES DES 2 VECTEURS
C
C SORTIE :
C --------
C DELTA  : LA NORME EUCLIDIENNE DU PRODUIT VECTORIEL DES 2 VECTEURS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DGL(2,3), DELTA
      INTRINSIC         SQRT
C
      DELTA = 0D0
      DO 10 I=1,3
         J = I + 1
         IF( J .GT. 3 ) J = J - 3
         K = I + 2
         IF( K .GT. 3 ) K = K - 3
         DELTA = (DGL(1,J) * DGL(2,K) - DGL(2,J) * DGL(1,K)) ** 2 +
     &            DELTA
   10 CONTINUE
      DELTA = SQRT( DELTA )

      RETURN
      END
