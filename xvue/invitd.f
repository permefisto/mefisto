      SUBROUTINE INVITD( TEXTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   INVITER A ENTRER UNE DONNEE A PARTIR D'UN TEXTE DE 1 LIGNE
C -----
C ENTREE :
C --------
C TEXTE : LE TEXTE INVITANT A ENTRER LA DONNEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     TEXTE
C
      IF( INTERA .GT. 0 .AND. LHLECT .EQ. 1 ) THEN
C
C        MODE INTERACTIF SUPPRESSION DES DOUBLES BLANCS
         L = NUDCNB( TEXTE )
C
         IF( INTERA .LT. 3 ) THEN
C
C           ENTREE CLAVIER
            WRITE(IMPRIM,10000) TEXTE(1:L)
10000       FORMAT(A,'?')
C
         ELSE
C
C           ENTREE MENU ET SOURIS : LA PREMIERE LIGNE DE L'INVITE
            NBLGRC(NRINVI) = 1
            KINVI(1) = TEXTE(1:L) // '? '
C           LE NUMERO DU CARACTERE AU DELA DU DERNIER NON BLANC
            L = MIN( NUDCNB( KINVI(1) ), NBCAIN )
            MDLGRC(NRINVI) = L
            MDRECT(NRINVI) = L
C
         ENDIF
      ENDIF
      END
