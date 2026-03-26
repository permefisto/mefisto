      INTEGER FUNCTION NCOQUA( QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NUMERO DE LA COULEUR EN FONCTION DE LA QUALITE
C -----
C
C ENTREE :
C --------
C QUALIT : QUALITE DE L'EF COMPRISE ENTRE 0 ET 1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1997
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     LA COULEUR DE L'EF VISUALISE SA QUALITE
      IF( QUALIT .LE. 0.1 ) THEN
         NCOQUA = NCROUG
      ELSE
         NCOQUA = N1COUL + 9 - INT( 10.0 * ( 1.0 - QUALIT ) )
      ENDIF
      END
