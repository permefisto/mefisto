      SUBROUTINE CAF1NB( NFFICA , NCA1NB , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TROUVER LE NUMERO DANS KFICA DU PREMIER CARACTERE NON BLANC
C -----
C ENTREE :
C --------
C NFFICA : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C
C SORTIE :
C --------
C NCA1NB : NUMERO DANS KFICA DU PREMIER CARACTERE NON BLANC
C NRETOU : 0 SI PAS D'ERREUR , >0 SINON
C          1001 SI LA FIN DE FICHIER DES CARACTERES EST ATTEINTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC
  5   DO 10 NCA1NB = LCFICA+1 , NCFICA
         IF( KFICA(NCA1NB:NCA1NB) .NE. ' ' ) RETURN
 10   CONTINUE
C
C     ICI LE BUFFER EST EPUISE
C     PASSAGE A L'ENREGISTREMENT SUIVANT DU FICHIER NFFICA
      READ( UNIT=NFFICA , FMT='(A)' , IOSTAT=NRETOU ) KFICA
      IF( NRETOU .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'CAF1NB: FICHIER CARACTERES EPUISE'
         CALL LEREUR
         NRETOU = 1001
         RETURN
      ENDIF
C
C     MISE A JOUR DU POINTEUR SUR LE BUFFER
      LCFICA = 0
      GOTO 5
      END
