      INTEGER FUNCTION INVTRIAR( S1, S2, S3, VN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INVERSER LA NUMEROTATION DU TRIANGLE SI SA NORMALE VNT
C -----    FAIT UN ANGLE DE PLUS DE 90 DEGRES AVEC LE VECTEUR VN

C ENTREES :
C ---------
C S1,S2,S3: LES 3 SOMMETS DU TRIANGLE
C VN      : LE VECTEUR DE DIRECTION A COMPARER A LA NORMALE DU TRIANGLE

C SORTIE :
C --------
C INVTRIAR: = 1 SI VN et NORMALE A S1S2S3 ONT UN ANGLE <=90 DEGRES
C           =-1 SI VN et NORMALE A S1S2S3 ONT UN ANGLE > 90 DEGRES
CC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC LJLL & St PIERRE DU PERRAY  Octobre 2011
C2345X7..............................................................012
      REAL              S1(3), S2(3), S3(3)
      DOUBLE PRECISION  VN(3), VNT(3), VVN, PROSCD

C     VECTEUR NORMAL AU TRIANGLE S1 S2 S3
      CALL VECNOR3( S1, S2, S3, VNT )

C     PRODUIT SCALAIRE VNT . VN
      VVN = PROSCD( VNT, VN, 3 )

      IF( VVN .LT. -1D-4 ) THEN
         INVTRIAR = -1
      ELSE
         INVTRIAR = 1
      ENDIF

      RETURN
      END
