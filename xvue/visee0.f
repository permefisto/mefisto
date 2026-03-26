      SUBROUTINE VISEE0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   INITIALISER LA VISEE PAR DEFAUT
C -----   VERSION xvue
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C ......................................................................
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      REAL      HEXSEC(6,2)
C
C     PAR DEFAUT L'HEXAEDRE DE SECTION EST L'HEXAEDRE ENGLOBANT
C     =========================================================
      IF( COOEXT(1,1) .LT. RINFO( 'GRAND' ) ) THEN
C
C        COOEXT EST INITIALISE CORRECTEMENT
         CALL TRTATA( COOEXT, HEXSEC, 6*2 )
         INIEXT = 1
C
      ELSE
C
C        COOEXT INCORRECT POUR LA SUITE => C'EST UN RECTANGLE CENTRE
         HEXSEC(1,1) =-1.0
         HEXSEC(1,2) = 1.0
         HEXSEC(2,1) =-0.6
         HEXSEC(2,2) = 0.6
         HEXSEC(3,1) = 0.0
         HEXSEC(3,2) = 0.0
C
      ENDIF
C
C     LA DIMENSION DE L'ESPACE DE L'OBJET
      IF( HEXSEC(3,1) .NE. 0.0  .OR. HEXSEC(3,2) .NE. 0.0 ) THEN
C        DIMENSION 3
         NDIMLI = 3
      ELSE IF( HEXSEC(2,1) .NE. 0.0  .OR. HEXSEC(2,2) .NE. 0.0 ) THEN
C        DIMENSION 2
         NDIMLI = 2
      ELSE
C        DIMENSION 1
         NDIMLI = 1
      ENDIF
C
C     RECHERCHE DE LA DIMENSION MAXIMALE ECAMAX DE L'OBJET
      ECAMAX = 0.0
      ECAMOY = 0.0
      DO 10  I=1,NDIMLI
         R = HEXSEC(I,2) - HEXSEC(I,1)
         IF( R .LT. 0 ) THEN
            RAPECR      = HEXSEC(I,1)
            HEXSEC(I,1) = HEXSEC(I,2)
            HEXSEC(I,2) = RAPECR
            R           =-R
         ENDIF
         ECAMAX = MAX( ECAMAX , R )
         ECAMOY = ECAMOY + R
 10   CONTINUE
      IF( ECAMAX .LE. 0.0 ) ECAMAX = 1.0
      ECAMOY = ECAMOY / NDIMLI
C
C     RAPPORT LARGEUR FENETRE / HAUTEUR FENETRE
      R = FLOAT(LAPXFE) / FLOAT(LHPXFE) * CYMMPX / CXMMPX
C
C     DEFINITION DE L'ECRAN EN COORDONNEES OBJET
C     ==========================================
      IF( NDIMLI .EQ. 1 ) THEN
C
C        OBJET UNIDIMENSIONNEL
C        ---------------------
C        AJUSTAGE SUR Y DE L'ECART MAXIMAL
         AXOPTV(1) = ( HEXSEC(1,1) + HEXSEC(1,2) ) * 0.5
         AXOPTV(2) = ( HEXSEC(2,1) + HEXSEC(2,2) ) * 0.5
         ECARTX    = ( HEXSEC(1,2) - HEXSEC(1,1) ) / 20
         AXOLAR    = ( HEXSEC(1,2) - HEXSEC(1,1) ) * 0.5 + ECARTX
         ECARTY    = ( HEXSEC(2,2) - HEXSEC(2,1) ) / 20
         AXOHAU    = ( HEXSEC(2,2) - HEXSEC(2,1) ) * 0.5 + ECARTY
         IF( AXOLAR .LE. 0 ) AXOLAR = 1
         IF( AXOHAU .LE. 0 ) AXOHAU = AXOLAR / R
         AXOARR = 0
         AXOAVA = 0
C
C        L'ECRAN VU EN COORDONNEES OBJET 2D EN UNITES UTILISATEUR
         CALL FENETRE( AXOPTV(1)-AXOLAR, AXOPTV(1)+AXOLAR,
     %                 AXOPTV(2)-AXOHAU, AXOPTV(2)+AXOHAU )
C
C        LE TYPE DE LA VISEE DEVIENT
         NOTYVI = 1
C
      ELSE IF( NDIMLI .EQ. 2 ) THEN
C
C        OBJET BIDIMENSIONNEL
C        --------------------
C        AJUSTAGE SUR Y DE L'ECART MAXIMAL
         ECART  = ECAMAX * 0.017
         AXOPTV(1) = ( HEXSEC(1,1) + HEXSEC(1,2) ) * 0.5
         AXOPTV(2) = ( HEXSEC(2,1) + HEXSEC(2,2) ) * 0.5
         AXOLAR = ( HEXSEC(1,2) - HEXSEC(1,1) ) * 0.5 + ECART
         AXOHAU = ( HEXSEC(2,2) - HEXSEC(2,1) ) * 0.5 + ECART
C        22/3/1999
         IF( AXOLAR .LT. AXOHAU ) AXOLAR = AXOHAU
C        22/3/1999
         IF( AXOHAU .LE. 0 ) AXOHAU = AXOLAR / R
         AXOARR = 0
         AXOAVA = 0
C
C        L'ECRAN VU EN COORDONNEES OBJET 2D EN UNITES UTILISATEUR
         CALL ISOFENETRE( AXOPTV(1)-AXOLAR, AXOPTV(1)+AXOLAR,
     %                    AXOPTV(2)-AXOHAU, AXOPTV(2)+AXOHAU )
C
C        LE TYPE DE LA VISEE DEVIENT
         NOTYVI = 1
C
      ELSE
C
C        OBJET TRIDIMENSIONNEL
C        ---------------------
C        AXONOMETRIE PAR AXOPTV AXOEIL ET AXOARR AXOAVA
         NOTYVI = 11
C
C        LA MARGE
         AXOLAR = ECAMOY * 0.98
         AXOHAU = AXOLAR / R
C
C        A PRIORI PAS DE TRONCATURE ENTRE LES PLANS ARRIERE ET AVANT
         AXOAVA = 0
         AXOARR = 0
C
C        LE POINT VISE
         AXOPTV(1) = ( HEXSEC(1,1) + HEXSEC(1,2) ) / 2
         AXOPTV(2) = ( HEXSEC(2,1) + HEXSEC(2,2) ) / 2
         AXOPTV(3) = ( HEXSEC(3,1) + HEXSEC(3,2) ) / 2

C        LONGUEUR DE LA DIAGONALE DES SOMMETS EXTREMES DE L'HEXAEDRE
         AXODIS = SQRT( (HEXSEC(1,2) - HEXSEC(1,1)) ** 2
     %                + (HEXSEC(2,2) - HEXSEC(2,1)) ** 2
     %                + (HEXSEC(3,2) - HEXSEC(3,1)) ** 2 )
C
C        LONGITUDE ET LATITUDE EN RADIANS DE LA POSITION DE L'OEIL PAR DEFAUT
         R = ATAN(1.0) / 45.0

cccC        pour le rouleau dans une boite vue de l'ecoulement en X
ccc         AXOLON = -80

         AXOLON = 30.0
         AXOLAT = 10.0
         RADLON = AXOLON * R
         RADLAT = AXOLAT * R
         COSLON = COS( RADLON )
         SINLON = SIN( RADLON )
         COSLAT = COS( RADLAT )
         SINLAT = SIN( RADLAT )

C        RAYON DANS LE PLAN XY
         RXY       = AXODIS * COSLAT
         AXOEIL(1) = AXOPTV(1) + RXY    * COSLON
         AXOEIL(2) = AXOPTV(2) + RXY    * SINLON
         AXOEIL(3) = AXOPTV(3) + AXODIS * SINLAT * 3
C
C        LA MATRICE DE L'AXONOMETRIE
         CALL MATAXO
         CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )

      ENDIF
C
      RETURN
      END
