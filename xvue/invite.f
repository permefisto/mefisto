      SUBROUTINE INVITE( NOTEXT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INVITER A ENTRER UNE DONNEE A PARTIR DU NUMERO DE LA LIGNE
C -----    DANS LE TABLEAU KINVTX
C
C ENTREE :
C --------
C NOTEXT : NUMERO DE LA LIGNE DE TEXTE DANS LE TABLEAU KINVTX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1999
C23456---------------------------------------------------------------012
      include"./incl/iinvtx.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      IF( INTERA .GT. 0 .AND. LHLECT .EQ. 1 ) THEN
C
C        CONTROLE DU NUMERO
         IF( NOTEXT .LE. 0 .OR. NOTEXT .GT. NINVTX ) THEN
            NBLGRC( NRERR ) = 3
            KERR(1) = 'NO DE LIGNE D''INVITE INCORRECT'
            KERR(2) = 'CORRIGER l''APPEL AVANT de RELANCER'
            KERR(3) = 'ou AUGMENTER NINVTX dans incl/iinvtx.inc'
            CALL LEREUR
            NOTEXT = 0
         ENDIF
C
C        MODE INTERACTIF SUPPRESSION DES DOUBLES BLANCS
         L = NUDCNB( KINVTX(NOTEXT) )
C
         IF( INTERA .LT. 3 ) THEN
C
C           ENTREE CLAVIER
            WRITE(IMPRIM,10000) KINVTX(NOTEXT)(1:L)
10000       FORMAT(A,'?')
C
         ELSE
C
C           ENTREE MENU ET SOURIS : LA PREMIERE LIGNE DE L'INVITE
            NBLGRC(NRINVI) = 1
            KINVI(1) = KINVTX(NOTEXT)(1:L) // '? '
C           LE NUMERO DU CARACTERE AU DELA DU DERNIER NON BLANC
            L = MIN( NUDCNB( KINVI(1) ), NBCAIN )
            MDLGRC(NRINVI) = L
            MDRECT(NRINVI) = L
C
         ENDIF
      ENDIF
C
      RETURN
      END
