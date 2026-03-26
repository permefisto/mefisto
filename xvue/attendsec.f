      SUBROUTINE ATTENDSEC( SECONDES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     ATTENDRE SECONDES DANS LE TEMPS DE L'UTILISATEUR
C -----
C
C ENTREE :
C --------
C SECONDES: NOMBRE DE SECONDES A ATTENDRE DANS LE TEMPS DE L'UTILISATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C2345X7..............................................................012
      DOUBLE PRECISION SECONDES, TINITIAL, TEMPSU, TFINAL
C
C     TEMPS INITIAL
      CALL Secondes1969( TINITIAL )
C
C     LE TEMPS A DEPASSER
      TFINAL = TINITIAL + SECONDES
C
C     TEMPS ACTUEL
      N  = 0
 10   CALL Secondes1969( TEMPSU )
      N = N + 1
      IF( TEMPSU .LT. TFINAL ) GOTO 10
C
      RETURN
      END
