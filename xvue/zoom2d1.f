      SUBROUTINE ZOOM2D1( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TRACE 2D AVEC ZOOM ET TRANSLATION
C -----    UN BOUTON ENFONCE et DEPLACE POUR REDEFINIR LA VISEE
C SORTIE :
C --------
C NOTYEV : NO DE L'EVENEMENT =0 SI ABANDON DEMANDE, NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X...............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
C
      REAL  DX, DY
C
C     INVITE A CLIQUER UN BOUTON
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Taper ECHAPPEMENT ou Cliquer le BOUTON 1 ou 3'
      ELSE
         KERR(1) = 'Type ESCAPE or Push the BUTTON 1 or 3'
      ENDIF
      CALL LERESU
C
C     EN ATTENTE D'EVENEMENT
 10   CALL XVSOURIS( NOTYEV, NBC, NOPX, NOPY )
C
      IF( NOTYEV .EQ. 0 .OR. NOTYEV .EQ. 2 ) THEN
C        ABANDON DEMANDE PAR FRAPPE AU CLAVIER OU BOUTON 2
         NOBOUTON = 2
         NOETATBOUTON = 0
         NOTYEV  = 0
         NORBITE = 0
C        NOPX0, NOPY0: NUMERO EN X ET Y DU PIXELS
C        VALEUR INDIQUANT UNE NON INITIALISATION
         NOPX0 = -1
         NOPY0 = -1
         RETURN
      ENDIF
C
      IF( NOETATBOUTON .EQ. 0 ) THEN
C        EN ATTENTE D'UN BOUTON ENFONCE
         IF( NOTYEV .EQ. -1 ) THEN
C           LE NO DE BOUTON ENFONCE ET L'ETAT CHANGE
            NOBOUTON = NBC
            NOETATBOUTON = 1
            NOPX0 = NOPX
            NOPY0 = NOPY
C           INVITE A CLIQUER UN BOUTON
            NBLGRC(NRERR) = 1
            WRITE(KERR(2)(1:1),'(I1)') NOBOUTON
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='RELACHER le BOUTON '//KERR(2)(1:1)// ' ou Echap'
            ELSE
               KERR(1)='RELAX the BUTTON '//KERR(2)(1:1) // ' or Escape'
            ENDIF
            CALL LERESU
         ENDIF
         GOTO 10
      ENDIF
C
C     LE BOUTON NOBOUTON EST ACTUELLEMENT ENFONCE
      IF( NOTYEV .EQ. 1 ) THEN
C        LE NO DE BOUTON EST RELACHE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 0
         NOPX0 = NOPX
         NOPY0 = NOPY
         GOTO 10
C
      ELSE IF( NOTYEV .EQ. -2 ) THEN
C
C        SOURIS BOUGEE BOUTON NOBOUTON ENFONCE
         IF( NOBOUTON .EQ. 1 ) THEN
C
C           BOUTON 1 ENFONCE: => TRANSLATION SELON CELLE EN PIXELS
C           TRANSLATION EN X
            DX = ( NOPX0 - NOPX ) * (XOBMAX-XOBMIN) / FLOAT(LAPXFE)
C           TRANSLATION EN Y
            DY = ( NOPY - NOPY0 ) * (YOBMAX-YOBMIN) / FLOAT(LHPXFE)
C           DECALAGE DE LA FENETRE
            XOBMIN = XOBMIN + DX
            XOBMAX = XOBMAX + DX
            YOBMIN = YOBMIN + DY
            YOBMAX = YOBMAX + DY
C
         ELSE IF( NOBOUTON .EQ. 3 ) THEN
C
C           BOUTON 3 ENFONCE: => ZOOM => CALCUL DE LA LARGEUR ET HAUTEUR
C           RAPPORT DU DEPLACEMENT DE LA SOURIS EN Y
            DY = ( NOPY0 - NOPY ) * (YOBMAX-YOBMIN) / FLOAT(LHPXFE)
            XOBMIN = XOBMIN - DY
            XOBMAX = XOBMAX + DY
            YOBMIN = YOBMIN - DY
            YOBMAX = YOBMAX + DY
C
         ENDIF
C
      ELSE
C        AUTRE VALEUR DE NOTYEV
         GOTO 10
      ENDIF
C
      NOPX0 = NOPX
      NOPY0 = NOPY
C
C     FENETRE MAXIMALE POUR VOIR L'OBJET
      CALL ISOFENETRE( XOBMIN, XOBMAX, YOBMIN, YOBMAX )
C
      AXOPTV(1) = ( XOBMIN + XOBMAX ) / 2
      AXOPTV(2) = ( YOBMIN + YOBMAX ) / 2
      AXOPTV(3) = 0
C
C     LA MEMOIRE PIXELS EST EFFACEE
      CALL EFFACEMEMPX
C
C     NOMBRE DE PASSAGES SUR ORBITE1
      NORBITE = NORBITE + 1
C
      RETURN
      END
