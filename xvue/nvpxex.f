      INTEGER FUNCTION NVPXEX( XVOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LE NOMBRE DE PIXELS DANS LA FENETRE DE L'ABSCISSE
C -----  DU VECTEUR D'ABSCISSE XVOB
C
C ENTREE :
C --------
C XVOB   : ABSCISSE COMPOSANTE DU VECTEUR
C
C SORTIE :
C --------
C NVPXEX : NOMBRE DE PIXELS OU ABSCISSE PX DU VECTEUR DE COMPOSANTE XVOB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    OCTOBRE 1995
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     TRANSFORMATION XVOB => PX X ET CONVERSION FINALE EN ENTIER
      NVPXEX = NINT( AXOBPX * XVOB )
      END
