      SUBROUTINE VECNOR( DGL, NORMVN, DVN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 3 COMPOSANTES DU VECTEUR NORMAL UNITAIRE
C -----    A PARTIR DES 2 VECTEURS TANGENTS DEFINIS DANS DGL
C
C ENTREES:
C --------
C DGL    : DGL(2,3) LES 2 DERIVEES DE G EN CE POINT
C          DGL(I,J) = D(GL)(J) / DXI
C
C SORTIE :
C --------
C NORMVN : LA NORME EUCLIDIENNE DU VECTEUR NORMAL DGL(1) TENSORIEL DGL(2)
C DVN    : LES 3 COMPOSANTES DU VECTEUR NORMAL UNITAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1999
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DGL(2,3), DVN(3), NORMVN
C
C     CALCUL DE LA NORMALE: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
C     ------------------------------------------------------
      DVN(1) = DGL(2,2) * DGL(1,3) - DGL(1,2) * DGL(2,3)
      DVN(2) = DGL(2,3) * DGL(1,1) - DGL(1,3) * DGL(2,1)
      DVN(3) = DGL(2,1) * DGL(1,2) - DGL(1,1) * DGL(2,2)
C
C     LES COMPOSANTES DU VECTEUR NORMAL DE NORME 1
C     --------------------------------------------
      NORMVN = SQRT( DVN(1) ** 2 + DVN(2) ** 2 + DVN(3) ** 2 )
      IF( NORMVN .GE. 1D-78 ) THEN
         DVN(1) = DVN(1) / NORMVN
         DVN(2) = DVN(2) / NORMVN
         DVN(3) = DVN(3) / NORMVN
      ENDIF
C
      RETURN
      END
