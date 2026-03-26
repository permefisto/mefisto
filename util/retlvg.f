      SUBROUTINE RETLVG( X1, Y1, Y2 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE D'UNE LIGNE VERTICALE A GAUCHE D'UN RECTANGLE
C -----
C
C ENTREE :
C --------
C X1     : ABSCISSE DE LA LIGNE VERTICALE
C Y1, Y2 : ORDONNEES INITIALE ET FINALE DE LA LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1992
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
C
C     LIGNE EPAISSIE 3 FOIS ET DECALEE
      CALL XVEPAISSEUR( 3 )
      X = X1 - 0.1
      IF( Y1 .LE. Y2 ) THEN
         YMINI = Y1
         YMAXI = Y2
      ELSE
         YMINI = Y2
         YMAXI = Y1
      ENDIF
      CALL TRAIT2D( NCBLAN, X, YMINI-0.1,  X, YMAXI+0.1 )
      END
