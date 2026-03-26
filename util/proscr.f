      REAL FUNCTION PROSCR(X,Y,NC)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PRODUIT SCALAIRE REEL SIMPLE PRECISION DE 2 VECTEURS X Y REELS
C ----- DE NC COMPOSANTES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1979
C23456---------------------------------------------------------------012
      REAL  X(NC), Y(NC)

      PROSCR = 0.
      DO I=1,NC
         PROSCR = PROSCR + X(I) * Y(I)
      ENDDO

      RETURN
      END
