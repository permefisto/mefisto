      SUBROUTINE LXTSDS( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETRUIRE LE NOM KNOM TABLEAU MEMOIRE SECONDAIRE DU LEXIQUE NTLX
C -----
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A RETROUVER DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1985
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CHARACTER*(*) KNOM
C
      CALL LXNMDS( NTLX , KNOM )
      END
