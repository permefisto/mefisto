      REAL FUNCTION YOB2PX( NYPX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LA VALEUR DE L'ORDONNEE OBJET 2D D'UN PIXEL D'ORDONNEE
C -----  NYPX DE LA FENETRE ACTUELLE

C        FONCTION INVERSE DE NUPXEY   VERSION xvue

C ENTREE :
C --------
C NYPX   : ORDONNEE PIXEL DU POINT

C SORTIE :
C --------
C YOB2PX : ORDONNEE OBJET 2D DU PIXEL DE LA FENETRE ACTUELLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/trvari.inc"

C     TRANSFORMATION PX Y  => Y OB
      YOB2PX = ( NYPX - BYOBPX ) / AYOBPX

      RETURN
      END
