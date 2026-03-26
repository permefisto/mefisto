      SUBROUTINE QUADCOUL2D( XY, COUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE QUADRANGLE DE SOMMETS XY ET DE COULEURS COUL
C -----    SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C          ATTENTION (X,Y) EN COORDONNEES OBJETS 2D
C
C ENTREES:
C --------
C XY     : 2 COORDONNEES DES 4 SOMMETS
C COUL   : NUMERO REEL DE LA COULEUR AUX 4 SOMMETS DU QUADRANGLE
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C2345X7..............................................................012
      include "./incl/minint.inc"
      REAL        XY(2,4), COUL(4), COULTR(3)
      INTEGER*2   XYPX(2,5), XYTR(2,3)
      INTRINSIC   INT2
C
C     TRANSFORMATION EN PIXELS DANS LA FENETRE XV
      DO 10 I=1,4
         NX = NUPXEX( XY(1,I) )
         NY = NUPXEY( XY(2,I) )
C        SI LE NUMERO PIXEL EST INCORRECT ABANDON DU TRACE DE LA FACE
         IF( NX .EQ. MININT .OR. NY .EQ. MININT ) RETURN
         XYPX(1,I) = INT2( NX )
         XYPX(2,I) = INT2( NY )
 10   CONTINUE
C
C     LE BARYCENTRE
      XYPX(1,5) = INT2( ( XYPX(1,1)+XYPX(1,2)+XYPX(1,3)+XYPX(1,4) ) /4 )
      XYPX(2,5) = INT2( ( XYPX(2,1)+XYPX(2,2)+XYPX(2,3)+XYPX(2,4) ) /4 )
      COULTR(3) = ( COUL(1) + COUL(2) + COUL(3) + COUL(4) ) / 4
C
C     TRACE DU TRIANGLE 125
      XYTR(1,1) = XYPX(1,1)
      XYTR(2,1) = XYPX(2,1)
      COULTR(1) = COUL(1)
C
      XYTR(1,2) = XYPX(1,2)
      XYTR(2,2) = XYPX(2,2)
      COULTR(2) = COUL(2)
C
      XYTR(1,3) = XYPX(1,5)
      XYTR(2,3) = XYPX(2,5)
      CALL TRIACOUL( XYTR, COULTR )
C
C     TRACE DU TRIANGLE 235
      XYTR(1,1) = XYPX(1,2)
      XYTR(2,1) = XYPX(2,2)
      COULTR(1) = COUL(2)
C
      XYTR(1,2) = XYPX(1,3)
      XYTR(2,2) = XYPX(2,3)
      COULTR(2) = COUL(3)
      CALL TRIACOUL( XYTR, COULTR )
C
C     TRACE DU TRIANGLE 345
      XYTR(1,1) = XYPX(1,3)
      XYTR(2,1) = XYPX(2,3)
      COULTR(1) = COUL(3)
C
      XYTR(1,2) = XYPX(1,4)
      XYTR(2,2) = XYPX(2,4)
      COULTR(2) = COUL(4)
      CALL TRIACOUL( XYTR, COULTR )
C
C     TRACE DU TRIANGLE 415
      XYTR(1,1) = XYPX(1,4)
      XYTR(2,1) = XYPX(2,4)
      COULTR(1) = COUL(4)
C
      XYTR(1,2) = XYPX(1,1)
      XYTR(2,2) = XYPX(2,1)
      COULTR(2) = COUL(1)
      CALL TRIACOUL( XYTR, COULTR )
C
      END
