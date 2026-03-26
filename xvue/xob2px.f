      REAL FUNCTION XOB2PX( NXPX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LA VALEUR DE L'ABSCISSE OBJET 2D D'UN PIXEL D'ABSCISSE
C -----  NXPX DE LA FENETRE ACTUELLE

C        FONCTION INVERSE DE NUPXEX    VERSION xvue

C ENTREE :
C --------
C NXPX   : ABSCISSE PIXEL DU POINT

C SORTIE :
C --------
C XOB2PX : ORDONNEE OBJET 2D DU PIXEL DE LA FENETRE ACTUELLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/trvari.inc"

C     TRANSFORMATION PX X  => X OB
      XOB2PX = ( NXPX - BXOBPX ) / AXOBPX

      RETURN
      END
