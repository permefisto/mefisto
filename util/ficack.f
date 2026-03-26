      SUBROUTINE FICACK( NFFICA , KNOM , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE DANS LE BUFFER KFICA LE NOM KNOM
C -----
C
C ENTREE :
C --------
C NFFICA : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C
C SORTIE :
C --------
C KNOM   : NOM A PORTER DANS LE BUFFER KFICA
C          EN FAIT SEULEMENT LES CARACTERES JUSQU'AU PREMIER BLANC
C NRETOU : 0 SI PAS D'ERREUR , >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     KNOM
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC DANS KFICA
      CALL CAF1NB( NFFICA , NC , NRETOU )
      IF( NRETOU .NE. 0 ) RETURN
C
C     LE NUMERO DU DERNIER CARACTERE DU NOM
      NBC = INDEX( KFICA(NC:NCFICA) , ' ' )
      IF( NBC .GT. 0 ) THEN
         NBC = NC + NBC - 2
      ELSE
         NBC = NCFICA
      ENDIF
C
C     LE NOM
      KNOM = KFICA(NC:NBC)
C
C     MISE A JOUR DU POINTEUR
      LCFICA = NBC
      END
