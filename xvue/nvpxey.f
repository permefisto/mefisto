      INTEGER FUNCTION NVPXEY( YVOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LE NOMBRE DE PIXELS DANS LA FENETRE DE L'ORDONNEE
C -----  DU VECTEUR D'ORDONNEE XVOB
C
C ENTREE :
C --------
C YVOB   : ORDONNEE COMPOSANTE DU VECTEUR
C
C SORTIE :
C --------
C NVPXEY : NOMBRE DE PIXELS OU ORDONNEE PX DU VECTEUR DE COMPOSANTE YVOB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    OCTOBRE 1995
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     TRANSFORMATION YOB => PX Y ET CONVERSION FINALE EN ENTIER
      NVPXEY = NINT( AYOBPX * YVOB )
      END
