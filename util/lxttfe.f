      SUBROUTINE LXTTFE( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FERMER LE NOM KNOM TABLEAU DE TABLEAUX DU LEXIQUE NTLX
C -----
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A RETROUVER DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1985
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      CHARACTER*(*) KNOM
      CHARACTER*4   KTYPE
C
C     RECUPERATION DE L'ADRESSE DU TABLEAU DE TABLEAUX
      CALL LXNMOU( NTLX , KNOM , KTYPE , NTTATA , MNTATA )
C
C     FERMETURE DES TABLEAUX OUVERTS (NUMERO <0 DANS LE TABLEAU DES NO)
      DO 10 N=4,MCN(MNTATA+3)+3
         NT = MCN( MNTATA + N )
         IF( NT .GE. 0 ) GOTO 10
C        LE TABLEAU OUVERT EST ALORS FERME
         NT = - NT
         CALL TAMSFE( NT )
         MCN( MNTATA + N ) = NT
 10   CONTINUE
C
C     FERMETURE DU TABLEAU DE TABLEAUX
      CALL TAMSFE( NTTATA )
      END
