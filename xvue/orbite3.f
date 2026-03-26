      SUBROUTINE ORBITE3( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TRACE 3D AVEC TRANSLATION ORBITE TRANSLATION
C -----    UN BOUTON ENFONCE et DEPLACE et RELACHE POUR REDEFINIR LA VISEE
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
      REAL  RLONGITUDE, RLATITUDE, DX, DY, XYZ(3)
C
cccC     INVITE A CLIQUER UN BOUTON
ccc      NBLGRC(NRERR) = 1
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         KERR(1) = 'CLIQUER un BOUTON 1 ou 2 ou 3 ou Echap'
ccc      ELSE
ccc         KERR(1) = 'PUSH the BUTTON 1 or 2 or 3 or Escape'
ccc      ENDIF
ccc      CALL LERESU
C
C     EN ATTENTE D'EVENEMENT D'ABANDON OU DE BOUTON ENFONCE
C     -----------------------------------------------------
 10   CALL XVSOURIS( NOTYEV, NBC, NOPX0, NOPY0 )
c     notyev: = 0 Si ABANDON par frappe de la touche Echappement ou @
c             = 1 Si CLIC ENFONCE et RELACHE D'UN BOUTON DE LA SOURIS => nopx no
c             =-1 Si CLIC SEULEMENT ENFONCE  D'UN BOUTON DE LA SOURIS => nopx no
c             =-2 Si le pointeur de la souris a bouge                 => nopx no
c             = 2 Si FRAPPE D'UN CARACTERE AU CLAVIER
c     nbc    : seulement actif si notyev est non nul
c              si notyev=-+1  nbc=numero du bouton 1 ou 2 ou 3
c              si notyev= -2  nbc=0 (pas de bouton designe)
c              si notyev= +2  nbc=numero du caractere dans la table ASCII
c   nopx,nopy: seulement actif si notyev=+-1 ou -2
c              coordonnees pixels du point clique par rapport au coin
c              superieur gauche de la fenetre
C
      IF( NOTYEV .EQ. 0 .OR. NOTYEV .EQ. 2 ) THEN
C
C        ABANDON DEMANDE PAR FRAPPE AU CLAVIER
         NOBOUTON = 0
         NOETATBOUTON = 0
         NOTYEV  = 0
         NORBITE = 0
         NOPX0 = -1
         NOPY0 = -1
         RETURN
C
      ENDIF
C
C     EN ATTENTE D'UN BOUTON ENFONCE
      IF( NOTYEV .EQ. -1 ) THEN
C        LE NO DE BOUTON ENFONCE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 1
         GOTO 20
      ENDIF
      GOTO 10
C
C     EN ATTENTE D'EVENEMENT A PARTIR D'UN BOUTON ENFONCE
C     ET QUE CE BOUTON SE RELACHE
C     ---------------------------------------------------
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
C     LE BOUTON NOBOUTON EST ACTUELLEMENT ENFONCE
      IF( NOTYEV .EQ. 1 ) THEN
C
C        LE BOUTON NBC EST RELACHE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 0
C
         IF( NOBOUTON .EQ. 1 ) THEN
C
C           BOUTON 1: => TRANSLATION SELON CELLE EN PIXELS
C           TRANSLATION EN X
ccc            DX = ( NOPX0 - NOPX ) * (5*AXOLAR) / FLOAT(LAPXFE)
            DX = ( NOPX0 - NOPX ) * (XOBMAX-XOBMIN) / FLOAT(LAPXFE)
ccc            DX = ( NOPX0 - NOPX ) * DIAMOB / FLOAT(LHPXFE)
C           TRANSLATION EN Y
ccc            DY = ( NOPY - NOPY0 ) * (5*AXOHAU) / FLOAT(LHPXFE)
            DY = ( NOPY - NOPY0 ) * (YOBMAX-YOBMIN) / FLOAT(LHPXFE)
ccc            DY = ( NOPY - NOPY0 ) * DIAMOB / FLOAT(LHPXFE)
C           DECALAGE DU POINT VU
            XYZ(1) = DX
            XYZ(2) = DY
            XYZ(3) = 0
            CALL AXOXYZ( XYZ, AXOPTV )
C           DECALAGE DE L'OEIL
            CALL XYZAXO( AXOEIL, XYZ )
            XYZ(1) = XYZ(1) + DX
            XYZ(2) = XYZ(2) + DY
            CALL AXOXYZ( XYZ, AXOEIL )
            CALL AXONOMETRIE( AXOPTV, AXOEIL,
     %                        AXOLAR, AXOHAU,
     %                        AXOARR, AXOAVA )
C           LONGITUDE ET LATITUDE ACTUELLE
            CALL LOLARA( AXOPTV(1), AXOPTV(2), AXOPTV(3),
     %                   AXOEIL(1), AXOEIL(2), AXOEIL(3),
     %                   RLONGITUDE0, RLATITUDE0, LEPOLE )
C
         ELSE IF( NOBOUTON .EQ. 2 ) THEN
C
C           BOUTON 2: => ORBITE => CALCUL DE LA LONGITUDE ET LATITUDE
C           LAPXFE REPRESENTE 2 PI EN LONGITUDE
C           LHPXFE REPRESENTE 1 PI EN LATITUDE
            PI = ATAN(1.0) * 4.0
            RLONGITUDE= RLONGITUDE0 + 2*PI *(NOPX0-NOPX)/FLOAT(LAPXFE)
            RLATITUDE = RLATITUDE0  -   PI *(NOPY0-NOPY)/FLOAT(LHPXFE)
            CALL PTVLONLAT(AXOPTV, RLONGITUDE*180/PI, RLATITUDE*180/PI)
            RLONGITUDE0 = RLONGITUDE
            RLATITUDE0  = RLATITUDE
C
         ELSE IF( NOBOUTON .EQ. 3 ) THEN
C
C           BOUTON 3: => ZOOM => CALCUL DE LA LARGEUR ET HAUTEUR
C           RAPPORT DU DEPLACEMENT DE LA SOURIS EN Y
            IF( NOPY .GT. NOPY0 ) THEN
C              ZOOM +: L'OBJET SE RAPPROCHE
ccc               N = 5*(NOPY-NOPY0)
ccc               AXOLAR = ( AXOLAR * (LHPXFE-N) + AXLMIN * N ) / LHPXFE
ccc               AXOHAU = ( AXOHAU * (LHPXFE-N) + AXHMIN * N ) / LHPXFE
ccc               DEPL = YOB2PX( NOPY - NOPY0 )
               DEPL = (NOPY-NOPY0) * (YOBMAX-YOBMIN) / LHPXFE
ccc               DEPL = (NOPY-NOPY0) * DIAMOB / LHPXFE
               AXOHAU = AXOHAU - DEPL
               AXOLAR = AXOLAR - DEPL * LAPXFE / LHPXFE
            ELSE IF( NOPY .LT. NOPY0 ) THEN
C              ZOOM -: L'OBJET S'ELOIGNE
ccc               N = 5*(NOPY0-NOPY)
ccc               AXOLAR = ( AXOLAR * (LHPXFE-N) + AXLMAX * N ) / LHPXFE
ccc               AXOHAU = ( AXOHAU * (LHPXFE-N) + AXHMAX * N ) / LHPXFE
ccc               DEPL = YOB2PX( NOPY0 - NOPY )
               DEPL = (NOPY0-NOPY) * (YOBMAX-YOBMIN) / LHPXFE
ccc               DEPL = (NOPY0-NOPY) * DIAMOB / LHPXFE
               AXOHAU = AXOHAU + DEPL
               AXOLAR = AXOLAR + DEPL * LAPXFE / LHPXFE
            ELSE
               GOTO 10
            ENDIF
         ENDIF
C
C        REDUCTION OU AGRANDISSEMENT DE LA FENETRE VUE
         CALL AXONOMETRIE( AXOPTV,  AXOEIL,
     %                     AXOLAR,  AXOHAU,
     %                     AXOARR,  AXOAVA )
C        LONGITUDE ET LATITUDE ACTUELLE EN RADIANS
         CALL LOLARA( AXOPTV(1), AXOPTV(2), AXOPTV(3),
     %                AXOEIL(1), AXOEIL(2), AXOEIL(3),
     %                RLONGITUDE0, RLATITUDE0,  LEPOLE )
C
      ELSE
C        AUTRE VALEUR DE NOTYEV
         GOTO 20
      ENDIF
C
C     SORTIE APRES INITIALISATION DE LA NOUVELLE VISEE
C     LA MEMOIRE PIXELS EST EFFACEE
      CALL EFFACEMEMPX
C
C     NOMBRE DE PASSAGES SUR ORBITE1
      NORBITE = NORBITE + 1
C
      RETURN
      END
