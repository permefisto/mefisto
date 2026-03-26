      SUBROUTINE FICACH( NFFICA , CAR4 , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE DANS LE BUFFER KFICA LES 4 CARACTERES
C -----
C
C ENTREE :
C --------
C NFFICA : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C
C SORTIE :
C --------
C CAR4   : LA VALEUR DES 4 CARACTERES LU EN CARACTERES DANS KFICA
C NRETOU : 0 SI PAS D'ERREUR , >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   MARS 1989
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*4       CAR4
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC DANS KFICA
      CALL CAF1NB( NFFICA , NC , NRETOU )
      IF( NRETOU .NE. 0 ) RETURN
C
C     LE NUMERO DU DERNIER CARACTERE DES 4 CARACTERES
      NBC = INDEX( KFICA(NC:NCFICA) , ' ' )
      IF( NBC .GT. 0 ) THEN
         NBC = NC + NBC - 2
      ELSE
         NBC = NCFICA
      ENDIF
C
C     LE 4 CARACTERES
      CAR4 = KFICA(NC:NBC)
C
C     MISE A JOUR DU POINTEUR
      LCFICA = NBC
      END
