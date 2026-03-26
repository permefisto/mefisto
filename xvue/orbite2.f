      SUBROUTINE ORBITE2( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TRACE 3D AVEC TRANSLATION ORBITE TRANSLATION
C -----    + DEPLACEMENT DE LA SOURIS SANS CLIC D'UN BOUTON
C SORTIE :
C --------
C NOTYEV : NO DE L'EVENEMENT =0 SI ABANDON DEMANDE, NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X...............................................................012
      include"./incl/trvari.inc"
C
      REAL  RLONGITUDE, RLATITUDE, DX, DY, XYZ(3)
C
      PI = ATAN(1.0) * 4.0
C
C     EN ATTENTE D'EVENEMENT
 12   CALL XVSOURIS( NOTYEV, NBC, NOPX, NOPY )
C
      IF( (NOTYEV.EQ.0.AND.NBC.GT.3) .OR. NOTYEV .EQ. 2 ) THEN
C        ABANDON DEMANDE PAR FRAPPE AU CLAVIER OU BOUTON 2
         NOBOUTON = 2
         NOETATBOUTON = 0
         NOTYEV = 0
         NORBITE = 0
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
            GOTO 12
         ELSE IF( NOTYEV .LE. -2 ) THEN
            NOETATBOUTON = 0
            NOPX1 = NOPX
            NOPY1 = NOPY
            NOTYEV = -3
            RETURN
         ENDIF
      ENDIF
C
C     LE BOUTON NOBOUTON EST ACTUELLEMENT ENFONCE
      IF( NOTYEV .EQ. 1 ) THEN
C        LE NO DE BOUTON EST RELACHE ET L'ETAT CHANGE
         NOBOUTON = NBC
         NOETATBOUTON = 0
         GOTO 12
      ELSE IF( NOTYEV .EQ. -2 ) THEN
C
C        SOURIS BOUGEE BOUTON NOBOUTON ENFONCE
C
         IF( NOBOUTON .EQ. 1 ) THEN
C
C           BOUTON 1: => TRANSLATION SELON CELLE EN PIXELS
C           TRANSLATION EN X
            DX = ( NOPX0 - NOPX ) * (4*AXOLAR) / FLOAT(LAPXFE)
C           TRANSLATION EN Y
            DY = ( NOPY - NOPY0 ) * (4*AXOHAU) / FLOAT(LHPXFE)
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
C           LHPXFE REPRESENTE   PI EN LATITUDE
            RLONGITUDE= RLONGITUDE0 + 2*PI*(NOPX0-NOPX)/FLOAT(LAPXFE)
            RLATITUDE = RLATITUDE0  -   PI*(NOPY0-NOPY)/FLOAT(LHPXFE)
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
               DEPL = (NOPY-NOPY0) * DIAMOB / LHPXFE
               AXOHAU = AXOHAU - DEPL
               AXOLAR = AXOLAR - DEPL * LAPXFE / LHPXFE
            ELSE IF( NOPY .LT. NOPY0 ) THEN
C              ZOOM -: L'OBJET S'ELOIGNE
ccc               N = NOPY0-NOPY
ccc               AXOLAR = ( AXOLAR * (LHPXFE-N) + AXLMAX * N ) / LHPXFE
ccc               AXOHAU = ( AXOHAU * (LHPXFE-N) + AXHMAX * N ) / LHPXFE
               DEPL = (NOPY-NOPY0) * DIAMOB / LHPXFE
               AXOHAU = AXOHAU - DEPL
               AXOLAR = AXOLAR - DEPL * LAPXFE / LHPXFE
            ELSE
               GOTO 12
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
         GOTO 12
      ENDIF
C
      NOPX0 = NOPX
      NOPY0 = NOPY
C
C     LA MEMOIRE PIXELS EST EFFACEE
      CALL EFFACEMEMPX
C
C     NOMBRE DE PASSAGES SUR ORBITE1
      NORBITE = NORBITE + 1
C
C     LE BOUTON 2 ICI SERT A L'ORBITE ET NON A L'ABANDON
C     LE NO D'EVENEMENT EST DONC MODIFIE
      IF( NOBOUTON .EQ. 2 ) NOTYEV = 1
C
      RETURN
      END
