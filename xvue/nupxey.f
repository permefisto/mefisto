      INTEGER FUNCTION NUPXEY( YOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNE LE NUMERO DU PIXEL ( ORDONNEE ) DE YOB DANS LA FENETRE
C -----
C
C ENTREE :
C --------
C YOB    : ORDONNEE OBJET DU POINT
C
C SORTIE :
C --------
C NUPXEY : NUMERO DU PIXEL OU ORDONNEE PX DU POINT D'ORDONNEE OBJET YOB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/minint.inc"
      include"./incl/trvari.inc"
C
C     TRANSFORMATION YOB => PX Y
      Y = AYOBPX * YOB + BYOBPX
      IF( Y .GT. MAXINT .OR. Y .LT. MININT ) THEN
         NUPXEY = MININT
         RETURN
      ENDIF
C
C     CONVERSION FINALE EN ENTIER
      NUPXEY = NINT( Y )
      END
