      SUBROUTINE INVITR( TEXTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   INVITER A ENTRER UNE DONNEE ET SAISIR ENSUITE DANS LA FENETRE
C -----   GRAPHIQUE ET NON PAS DANS UNE FENETRE AVEC MENU
C ENTREE :
C --------
C TEXTE : LE TEXTE INVITANT A ENTRER LA DONNEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*(*)     TEXTE
C
      IF( INTERA .GT. 0 .AND. LHLECT .EQ. 1 ) THEN
C
C        MODE INTERACTIF
         IF( INTERA .LT. 3 ) THEN
C
C           ENTREE CLAVIER
            WRITE(IMPRIM,10000) TEXTE
10000       FORMAT(A,'?')
C
         ELSE
C
C           L'INVITE PRECEDENTE EST EFFACEE
            CALL RECTEF( NRINVI )
C
C           ENTREE MENU ET SOURIS : LA PREMIERE LIGNE DE L'INVITE
            NBLGRC(NRINVI) = 1
            KINVI(1) = TEXTE
C
C           LE NUMERO DU CARACTERE AU DELA DU DERNIER NON BLANC
            L             = MIN( NUDCNB( KINVI(1) ) + 1 , NBCAIN )
            KINVI(1)(L:L) = '?'
            MDLGRC(NRINVI) = L
            MDRECT(NRINVI) = L
            CALL LESREC
C
         ENDIF
      ENDIF
      END
