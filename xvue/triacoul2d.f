      SUBROUTINE TRIACOUL2D( XY, COUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE TRIANGLE DE SOMMETS XY ET DE COULEURS COUL
C -----    SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C          ATTENTION (X,Y) EN COORDONNEES OBJETS 2D
C
C ENTREES:
C --------
C XY     : 2 COORDONNEES DES 3 SOMMETS
C COUL   : NUMERO REEL DE LA COULEUR AUX 3 SOMMETS DU TRIANGLE
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    OCTOBRE 1994
C2345X7..............................................................012
      include"./incl/minint.inc"
      REAL           XY(2,3), COUL(3)
      INTEGER*2      XYPX(2,3)
C
C     TRANSFORMATION EN PIXELS DANS LA FENETRE XV
      DO 10 I=1,3
         NX = NUPXEX( XY(1,I) )
         NY = NUPXEY( XY(2,I) )
C        SI LE NUMERO PIXEL EST INCORRECT ABANDON DU TRACE DE LA FACE
         IF( NX .EQ. MININT .OR. NY .EQ. MININT ) RETURN
         XYPX(1,I) = INT2( NX )
         XYPX(2,I) = INT2( NY )
 10   CONTINUE
C
C     TRACE EFFECTIF DU REMPLISSAGE DU TRIANGLE SELON DES COULEURS PROGRESSIVES
      CALL TRIACOUL( XYPX, COUL )
C
      RETURN
      END
