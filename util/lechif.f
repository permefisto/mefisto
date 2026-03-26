      INTEGER FUNCTION LECHIF( CAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER L'ENTIER REPRESENTE PAR LE CARACTERE CAR
C -----    OU -1 SI CE N'EST PAS UN ENTIER 0 1 2 3 4 5 6 7 8 9
C
C ENTREE :
C --------
C CAR    : LE CARACTERE A DECRIPTE
C
C SORTIE :
C --------
C LECHIF : -1 SI CAR N'EST PAS 0 1 2 3 4 5 6 7 8 9
C          LE CHIFFRE ENTIER SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      CHARACTER*1       CAR
C
      IF( CAR .EQ. '0' ) THEN
         LECHIF = 0
      ELSE IF( CAR .EQ. '1' ) THEN
         LECHIF = 1
      ELSE IF( CAR .EQ. '2' ) THEN
         LECHIF = 2
      ELSE IF( CAR .EQ. '3' ) THEN
         LECHIF = 3
      ELSE IF( CAR .EQ. '4' ) THEN
         LECHIF = 4
      ELSE IF( CAR .EQ. '5' ) THEN
         LECHIF = 5
      ELSE IF( CAR .EQ. '6' ) THEN
         LECHIF = 6
      ELSE IF( CAR .EQ. '7' ) THEN
         LECHIF = 7
      ELSE IF( CAR .EQ. '8' ) THEN
         LECHIF = 8
      ELSE IF( CAR .EQ. '9' ) THEN
         LECHIF = 9
      ELSE
         LECHIF = -1
      ENDIF

      RETURN
      END
