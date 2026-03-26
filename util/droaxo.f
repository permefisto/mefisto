      SUBROUTINE DROAXO( X, Y, Z )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER JUSQU'AU AU POINT X Y Z EN AXONOMETRIE (CF SAPT3D)
C -----
C
C ENTREES :
C ---------
C X,Y,Z  : LES 3 COORDONNEES DU POINT OU ALLER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1990
C23456---------------------------------------------------------------012
      include"./incl/cm4vue.inc"
C
C     PASSAGE AUX COORDONNEES DE L'AXONOMETRIE
      CALL AXOYZ( X,Y,Z, R,S )
C
C     PASSAGE EN COORDONNEES CM
      XX = ( ( YAXMAX-R) * AXX1 + (R-YAXMIN) * AXX2 ) / (YAXMAX-YAXMIN)
      YY = ( ( ZAXMAX-S) * AXY1 + (S-ZAXMIN) * AXY2 ) / (ZAXMAX-ZAXMIN)
C
C     DEPLACEMENT AU POINT ECRAN XX YY
C     ATTENTION PAS DE CLIPPAGE !
      CALL TRAIT2D( NC, X0, Y0, XX, YY )
      END
