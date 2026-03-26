      SUBROUTINE NORMER( N, V, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    NORMALISATION EUCLIDIENNE A 1 D UN VECTEUR V DE N COMPOSANTES
C -----
C
C ENTREE :
C --------
C N      : NOMBRE DE COMPOSANTES DU VECTEUR V
C
C MODIFIE:
C --------
C V      : LE VECTEUR AVANT ET APRES NORMALISATION A 1
C
C SORTIE :
C --------
C IERR   : 1 SI NORME DE V NULLE
C          0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  MARS 1987
C ......................................................................
      REAL              V( N )
      DOUBLE PRECISION  S
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTRINSIC         SQRT
C
      S = 0.D0
      DO 10 I=1,N
         S = S + V( I ) * V( I )
   10 CONTINUE
C
C     TEST DE NULLITE DE LA NORME DU VECTEUR
C     --------------------------------------
      IF( S .LE. 0D0 ) THEN
C        NORME NULLE DU VECTEUR NON NORMALISABLE A 1'
         IERR = 1
         RETURN
      ENDIF
C
      S = 1.D0 / SQRT( S )
      DO 20 I=1,N
         V( I ) = REAL( V ( I ) * S )
   20 CONTINUE
C
      IERR = 0
      RETURN
      END
