      SUBROUTINE RETLHB( X1, X2, Y1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE D'UNE LIGNE HORIZONTALE BASSE ENTRE 2 LIGNES DE TEXTE
C -----
C
C ENTREE :
C --------
C X1, X2 : ABSCISSES DU DEBUT ET FIN DE LIGNE HORIZONTALE
C Y1     : ORDONNEE DE LA LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1992
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
C
C     LIGNE EPAISSIE 3 FOIS ET DECALEE
      CALL XVEPAISSEUR( 3 )
      YY = Y1 + 0.05
      IF( X1 .LE. X2 ) THEN
         XMINI = X1
         XMAXI = X2
      ELSE
         XMINI = X2
         XMAXI = X1
      ENDIF
      CALL TRAIT2D( NCNOIR, XMINI, YY, XMAXI+0.07, YY )
      END
