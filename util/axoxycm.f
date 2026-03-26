      SUBROUTINE AXOXYCM( X, Y, Z, XX, YY  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 2 COORDONNEES AXONOMETRIQUES PUIS
C -----    LES CONVERTIR EN CM DANS LA FENETRE
C
C ENTREES :
C ---------
C X,Y,Z  : LES 3 COORDONNEES DU POINT OU ALLER
C
C SORTIES :
C ---------
C XX, YY  : LES COORDONNEES CM DU POINT DANS LA FENETRE xvue
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       JUIN 1994
C23456---------------------------------------------------------------012
      include"./incl/cm4vue.inc"
C
C     PASSAGE AUX COORDONNEES DE L'AXONOMETRIE
      CALL AXOYZ( X,Y,Z, R,S )
C
C     PASSAGE EN COORDONNEES CM
      XX = ( ( YAXMAX-R) * AXX1 + (R-YAXMIN) * AXX2 ) / (YAXMAX-YAXMIN)
      YY = ( ( ZAXMAX-S) * AXY1 + (S-ZAXMIN) * AXY2 ) / (ZAXMAX-ZAXMIN)
      END
