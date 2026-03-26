      INTEGER FUNCTION NUPXEX( XOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LE NUMERO DU PIXEL ( ABSCISSE ) DE XOB DANS LA FENETRE
C -----
C
C ENTREE :
C --------
C XOB    : ABSCISSE OBJET DU POINT
C
C SORTIE :
C --------
C NUPXEX : NUMERO DU PIXEL OU ABSCISSE PX DU POINT D'ABSCISSE OBJET XOB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/minint.inc"
      include"./incl/trvari.inc"
C
C     TRANSFORMATION XOB => PX X
      X = AXOBPX * XOB + BXOBPX
      IF( X .GT. MAXINT .OR. X .LT. MININT ) THEN
         NUPXEX = MININT
         RETURN
      ENDIF
C
C     CONVERSION FINALE EN ENTIER
      NUPXEX = NINT( X )
      END
