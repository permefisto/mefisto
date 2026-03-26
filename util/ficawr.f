      SUBROUTINE FICAWR( NFFICA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     ECRIRE LE BUFFER KFICA SUR LE FICHIER NFFICA
C -----
C
C ENTREE :
C --------
C NFFICA : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     UN BUFFER VIDE N'EST PAS ECRIT
      IF( LCFICA .LE. 0 ) RETURN
C
C     MISE A BLANC AU DELA DU POINTEUR
      IF( LCFICA .LT. NCFICA ) THEN
         KFICA(LCFICA+1:NCFICA) = ' '
      ENDIF
C
C     MISE SUR LE FICHIER
      WRITE( NFFICA , '(A)' , IOSTAT=I ) KFICA
C
C     SI ERREUR IMPRESSION
      IF( I .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NFFICA
         KERR(1) =  'FICAWR: ERREUR EN ECRITURE FICHIER'//KERR(2)(1:4)
         CALL LEREUR
         RETURN
      ENDIF
C
C     MISE A JOUR DU POINTEUR
      LCFICA = 0
      END
