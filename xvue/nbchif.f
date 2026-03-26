       INTEGER FUNCTION NBCHIF( LENTIE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LE NOMBRE DE CHIFFRES DE L'ENTIER LENTIE
C -----  LE SIGNE - COMPTE POUR UN CHIFFRE
C
C ENTREE:
C --------
C LENTIE : LE NOMBRE ENTIER POSITIF OU NEGATIF OU NUL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1994
C2345X7..............................................................012
      IF( LENTIE .GE. 0 ) THEN
         LEQUOT = LENTIE
         NBCHIF = 0
      ELSE
         LEQUOT = -LENTIE
C        LE SIGNE - COMPTE POUR UN CHIFFRE
         NBCHIF = 1
      ENDIF
C
 10   LEQUOT = LEQUOT / 10
      NBCHIF = NBCHIF + 1
      IF( LEQUOT .GT. 0 ) GOTO 10
      END
