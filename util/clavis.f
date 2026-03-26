      SUBROUTINE CLAVIS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFACER LE CLAVIER ET LA ZONE TEXTE ET LIGNES
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C23456---------------------------------------------------------------012
CCC      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
C     EFFACAGE DE L'INVITE
CCC      CALL RECTEF( NRINVI )
C
C     EFFACAGE DE LA LIGNE DE SAISIE
CCC      CALL RECTEF( NRLGSA )
C
C     EFFACAGE DES LIGNES LUES
CCC      CALL RECTEF( NRLGLU )
      WRITE(IMPRIM,*) 'ENTREE DANS CLAVIS: APPEL A SUPPRIMER'
      END
