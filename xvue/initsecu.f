      SUBROUTINE INITSECU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  INITIALISE TEMPUSER DANS LE TEMPS DE L'UTILISATEUR
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     TEMPS INITIAL
      CALL Secondes1969( TEMPUSER )
C
      RETURN
      END
