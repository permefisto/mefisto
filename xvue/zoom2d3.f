      SUBROUTINE ZOOM2D3( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TRACE 2D AVEC ZOOM ET TRANSLATION
C -----    UN BOUTON ENFONCE et DEPLACE et RELACHE POUR DEFINIR LA VISEE
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
         KERR(1) = 'CLIQUER un BOUTON 1 ou 3 ou Echap'
      ELSE
         KERR(1) = 'PUSH the BUTTON 1 or 3 or Escape'
      ENDIF
      CALL LERESU
C
C     EN ATTENTE D'EVENEMENT D'ABANDON OU DE BOUTON ENFONCE
C     -----------------------------------------------------
 10   CALL XVSOURIS( NOTYEV, NBC, NOPX0, NOPY0 )
C
      IF( NOTYEV .EQ. 0 .OR. NOTYEV .EQ. 2 ) THEN
C        ABANDON DEMANDE PAR FRAPPE AU CLAVIER
         NOBOUTON = 0
         NOETATBOUTON = 0
         NOTYEV  = 0
         NORBITE = 0
         NOPX0 = -1
         NOPY0 = -1
         RETURN
      ENDIF
C
C     EN ATTENTE D'UN BOUTON ENFONCE
      IF( NOTYEV .EQ. -1 ) THEN
         IF( NBC .EQ. 2 ) GOTO 10
C        LE NO DE BOUTON 1 ou 3 EST ENFONCE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 1
         GOTO 20
      ENDIF
      GOTO 10
C
C     EN ATTENTE D'EVENEMENT A PARTIR D'UN BOUTON ENFONCE
C     ET QUE CE BOUTON SE RELACHE
C     ---------------------------------------------------
C
C     INVITE A CLIQUER UN BOUTON
 20   NBLGRC(NRERR) = 1
      WRITE(KERR(2)(1:1),'(I1)') NOBOUTON
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'RELACHER le BOUTON ' // KERR(2)(1:1) // ' ou Echap'
      ELSE
         KERR(1) = 'RELAX the BUTTON ' // KERR(2)(1:1) // ' or Escape'
      ENDIF
      CALL LERESU
C
      CALL XVSOURIS( NOTYEV, NBC, NOPX, NOPY )
C
      IF( NOTYEV .EQ. 0 .OR. NOTYEV .EQ. 2 ) THEN
C
C        ABANDON DEMANDE PAR FRAPPE AU CLAVIER
         NOBOUTON = 0
         NOETATBOUTON = 0
         NOTYEV = 0
         NORBITE = 0
         RETURN
C
      ENDIF
C
      IF( NOTYEV .EQ. 1 ) THEN
C
C        LE NO DE BOUTON EST RELACHE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 0
C
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
C
C        AUTRE VALEUR DE NOTYEV
         GOTO 20
C
      ENDIF
C
C     SORTIE APRES INITIALISATION DE LA NOUVELLE VISEE
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
