      SUBROUTINE TRAIT2D( NC, X1,Y1, X2,Y2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE TRAIT 2D DE (X1,Y1) A (X2,Y2) AVEC LA COULEUR NC
C -----  ATTENTION (X,Y) EN COORDONNEES OBJET 2D
C        LA TRANSFORMATION EN PIXELS EST ASSUREE DANS CE SP
C        LES AUTRES CARACTERISTIQUES DU TRACE SONT CELLES ACTUELLES
C
C ENTREES:
C --------
C NC     : NUMERO DE LA COULEUR DU TRAIT A TRACER
C X1,Y1  : 1-ERE EXTREMITE DU TRAIT
C X2,Y2  : 2-EME EXTREMITE DU TRAIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/minint.inc"
C
C     SI COULEUR NEGATIVE PAS DE TRACE
      IF( NC .LT. 0 ) RETURN
C
C     TRANSFORMATION EN PIXELS DANS LA FENETRE XV
      NX1 = NUPXEX( X1 )
      NX2 = NUPXEX( X2 )
      NY1 = NUPXEY( Y1 )
      NY2 = NUPXEY( Y2 )
      IF( NX1 .EQ. MININT .OR. NY1 .EQ. MININT .OR.
     %    NX2 .EQ. MININT .OR. NY2 .EQ. MININT ) RETURN
C
C     TRACE EFFECTIF
      CALL XVCOULEUR( NC )
      CALL XVTRAIT( NX1, NY1, NX2, NY2 )

      RETURN
      END
