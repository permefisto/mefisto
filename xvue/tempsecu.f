      SUBROUTINE TEMPSECU( SECONDES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     ATTENDRE SECONDES DEPUIS L'APPEL DE INITSEC DANS LE TEMPS
C -----     DE L'UTILISATEUR AVANT DE QUITTER CE SOUS-PROGRAMME
C
C ENTREE :
C --------
C SECONDES: NOMBRE DE SECONDES A ATTENDRE DANS LE TEMPS DE L'UTILISATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C2345X7..............................................................012
      include"./incl/trvari.inc"
      DOUBLE PRECISION SECONDES, TEMPSU, TD
C
C     LE TEMPS A DEPASSER DEPUIS L'APPEL DE INITSEC
      TD = TEMPUSER + SECONDES
C
C     TEMPS ACTUEL
 10   CALL Secondes1969( TEMPSU )
      IF( TEMPSU .LT. TD ) GOTO 10
C
      RETURN
      END
