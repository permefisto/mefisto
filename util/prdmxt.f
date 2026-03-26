      REAL FUNCTION PRDMXT(V1,V2,V3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PRODUIT MIXTE DES 3 VECTEURS V1 V2 V3 = V1 . ( V2 x V3 )
C -----
C ENTREES:
C --------
C V1, V2, V3 : LES 3 VECTEURS DE 3 COMPOSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 2000
C2345X7..............................................................012
      REAL V1(3),V2(3),V3(3), V(3)
C
      CALL PROVER(V2,V3,V)
      PRDMXT=PROSCR(V1,V,3)
      RETURN
      END
