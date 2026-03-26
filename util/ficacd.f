      SUBROUTINE FICACD( NFFICA , REEL2 , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE DANS LE BUFFER KFICA LE REEL2
C -----
C
C ENTREE :
C --------
C NFFICA : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C
C SORTIE :
C --------
C REEL2  : LA VALEUR DU REEL2 LU EN CARACTERES DANS KFICA
C NRETOU : 0 SI PAS D'ERREUR , >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  REEL2
      CHARACTER*2       KNB
      CHARACTER*8       KFORMA
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC DANS KFICA
      CALL CAF1NB( NFFICA , NC , NRETOU )
      IF( NRETOU .NE. 0 ) RETURN
C
C     LE NUMERO DU DERNIER CARACTERE DU REEL2
      NBC = INDEX( KFICA(NC:NCFICA) , ' ' )
      IF( NBC .GT. 0 ) THEN
         NBC = NC + NBC - 2
      ELSE
         NBC = NCFICA
      ENDIF
C     LE NOMBRE DE CARACTERES
      NB = NBC - NC + 1
C
C     CONVERSION DE NB EN CARACTERES
      WRITE( KNB , '(I2)' ) NB
C
C     GENERATION DU FORMAT DE LECTURE EN  D .17
      KFORMA = '(D' // KNB // '.17)'
C
C     LE REEL2
      READ( KFICA(NC:NBC) , KFORMA , IOSTAT=NRETOU ) REEL2
C
C     MISE A JOUR DU POINTEUR
      LCFICA = NBC
      END
