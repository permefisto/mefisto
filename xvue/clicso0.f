      SUBROUTINE CLICSO0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ATTENDRE UN CLIC SUR UNE DES TOUCHES DE LA SOURIS
C------    OU BIEN UN CARACTERE AU CLAVIER  VERSION xvue
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      IF( INTERA .GE. 3 ) THEN
C
C        PREVENIR POUR LE CLIC SOURIS OU L'ENTREE D'UN CARACTERE AU CLAVIER
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CLIQUER pour CONTINUER ...'
         ELSE
            KERR(1) = 'CLICK to CONTINUE ...'
         ENDIF
         CALL LERESU
C
C        SAISIE D'UN POINT PAR CLIC DE LA SOURIS OU ENTREE D'UN CARACTERE
         CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
C
C        LE RECTANGLE DE L'INVITE EST EFFACE
         CALL RECTEF( NRERR )
      ENDIF
C
      RETURN
      END
