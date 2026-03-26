      DOUBLE PRECISION FUNCTION PROSCD( X, Y, N )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PRODUIT SCALAIRE DE 2 VECTEURS X Y REELS DOUBLE PRECISION
C ----- DE N COMPOSANTES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  MARS 1987
C ......................................................................
      DOUBLE PRECISION  X(N), Y(N)
      INTEGER           K, N
C
      PROSCD = 0.D0
      DO K = 1, N
         PROSCD = PROSCD + X(K) * Y(K)
      ENDDO
C
      RETURN
      END
