      SUBROUTINE TRIA2D( XY, NOCOUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE TRIANGLE DE SOMMETS XY SELON LA COULEUR NOCOUL
C -----
C          ATTENTION (X,Y) EN COORDONNEES OBJETS 2D
C
C ENTREES:
C --------
C XY     : 2 COORDONNEES DES 3 SOMMETS
C NOCOUL : NUMERO REEL DE LA COULEUR DU TRIANGLE
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 2007
C2345X7..............................................................012
      include"./incl/minint.inc"
      REAL           XY(2,3)
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
C     TRACE EFFECTIF DU REMPLISSAGE DU TRIANGLE SELON LA COULEUR
      CALL XVCOULEUR( NOCOUL )
      CALL XVFACE( 3, XYPX )
C
      RETURN
      END
