      SUBROUTINE NORME1( N, V, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NORMALISATION EUCLIDIENNE A 1 D UN VECTEUR V DE N COMPOSANTES
C -----
C ENTREES :
C ---------
C N       : NOMBRE DE COMPOSANTES DU VECTEUR
C
C MODIFIE :
C ---------
C V       : LE VECTEUR A NORMALISER A 1
C
C SORTIE  :
C ---------
C IERR    : 1 SI LA NORME DE V EST EGALE A 0
C           0 SI PAS D'ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS             MARS 1987
C ......................................................................
      DOUBLE PRECISION  V( N ), S, SQRT
C
      S = 0.0D0
      DO 10 I=1,N
         S = S + V( I ) * V( I )
   10 CONTINUE
C
C     TEST DE NULLITE DE LA NORME DU VECTEUR
C     --------------------------------------
      IF( S .LE. 0.0D0 ) THEN
C        NORME NULLE DU VECTEUR NON NORMALISABLE A 1
         IERR = 1
         RETURN
      ENDIF
C
      S = 1.0D0 / SQRT( S )
      DO 20 I=1,N
         V( I ) = V ( I ) * S
   20 CONTINUE
C
      IERR = 0
      END
