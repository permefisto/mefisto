      SUBROUTINE CLICSO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ATTENDRE UN CLIC SUR UNE DES TOUCHES DE LA SOURIS
C------    OU BIEN UN CARACTERE AU CLAVIER  VERSION xvue
C          LORBITE EST MIS A ZERO SI CARACTERE Abandon ou @ FRAPPE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      IF( INTERA .GE. 3 ) THEN
C
C        PREVENIR POUR LE CLIC SOURIS OU L'ENTREE D'UN CARACTERE AU CLAVIER
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CLIQUER pour CONTINUER'
            KERR(2) = 'ou FRAPPER Escape'
         ELSE
            KERR(1) = 'CLICK to CONTINUE ...'
            KERR(2) = 'or TYPE Escape'
         ENDIF
         CALL LERESU
C
C        POUR VIDER LE BUFFER DE X11
         CALL XVVOIR
C
C        SAISIE D'UN POINT PAR CLIC DE LA SOURIS OU ENTREE D'UN CARACTERE
         CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
C        NOTYEV: NUMERO DU TYPE DE L'EVENEMENT
C         0        => ABANDON DEMANDE
C         1 2 ou 3 => NUMERO DU BOUTON DE LA SOURIS ET
C                       NOPXX, NOPXY INITIALISES AU POINT CLIC DE LA SOURIS
C        -1        => CARACTERE TAPE AU CLAVIER AVEC NOCHAR
C        NOPXX, NOPXY : SI NOTYEV>0 LES COORDONNEES PIXELS DU POINT CLIQUE
C        NOCHAR       : SI NOTYEV<0 NUMERO ASCII DU CARACTERE TAPE AU CLAVIER
C        Le caractere Echappement EST en position 27 dans la table ASCII
C        Le caractere @           EST en position 64 dans la table ASCII

C        UN CLIC SOURIS ou LA FRAPPE d'UN CARACTERE a ete FAIT
C        => ARRET de LORBITE
         LORBITE = 0

C        LE RECTANGLE DE L'INVITE EST EFFACE
         CALL RECTEF( NRERR )

      ENDIF
C
      RETURN
      END
