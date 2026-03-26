      SUBROUTINE TRVI3D( KNOMOB,   NBPOI,    XYZPOI,
     %                   SEMINVIT, SEMAXVIT, CMVITE,
     %                   NCAS0,    NCAS1,    vitx,   vity,   vitz,
     %                   VITMOY,   VITMIN,   VITMAX, TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES VITESSES EN TOUS LES POINTS DE LA TETRAEDRISATION
C -----
C ENTREES:
C --------
C KNOMOB :  NOM DE L'OBJET 3D
C NBPOI  :  NOMBRE DE POINTS = NOEUDS DES TETRAEDRES
C XYZPOI :  3 COORDONNEES DES NBPOI POINTS DES EF
C SEMINVIT: MIN DE LA NORME DE LA VITESSE AU DESSOUS DUQUEL LA
C           FLECHE NE DOIT PAS ETRE TRACEE
C SEMAXVIT: MAX DE LA NORME DE LA VITESSE AU DESSUS  DUQUEL LA
C           FLECHE NE DOIT PAS ETRE TRACEE
C CMVITE :  NOMBRE DE CM PAR UNITE DE VITESSE

C NCAS0:NCAS1 : NOMBRE DE VECTEURS VITESSE STOCKES  (=NCAS1-NCAS0+1)
C NCAS0  : NUMERO DU PREMIER CAS A TRACER PARMI LES NCAS0:NCAS1 VECTEURS
C NCAS1  : NUMERO DU DERNIER CAS A TRACER PARMI LES NCAS0:NCAS1 VECTEURS

C VITX   : COMPOSANTE X DE LA VITESSE EN CHAQUE POINT DES TETRAEDRES
C VITY   : COMPOSANTE Y DE LA VITESSE EN CHAQUE POINT DES TETRAEDRES
C VITZ   : COMPOSANTE Z DE LA VITESSE EN CHAQUE POINT DES TETRAEDRES

C VITMOY : NORME MOYENNE  DE LA VITESSE EN UN NOEUD (=POINT)
C VITMIN : NORME MINIMALE DE LA VITESSE EN UN NOEUD
C VITMAX : NORME MAXIMALE DE LA VITESSE EN UN NOEUD
C TIMES  : TEMPS DU CALCUL DES NBVPFILE VECTEURS STOKES SUR FICHIERS

C REMARQUE:
C          LA FONCTION UTILISATEUR Region(t,x,y,z) PEUT LIMITER LA
C          REGION OU TRACER LA FLECHE DE LA VITESSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris  Fevrier 2007
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY              Mars 2021
C23456...............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/ctemps.inc"
      include"./incl/xvfontes.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*24      NOMFGIF
      REAL              XYZPOI(3,NBPOI), TIMES(NCAS0:NCAS1)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity, vitz

      DOUBLE PRECISION  SEMINVIT, SEMAXVIT, VITMOY, VITMIN, VITMAX
      INTRINSIC         SQRT
C
      NBCOOR = 3
      MNXYZF = 0
      MNBAFR = 0
      MNNBFR = 0
      MNAREF = 0
      NBCOUL = NDCOUL - N1COUL + 1
C
C     NOM DU FICHIER VIDEO  SELON MODECO // 'Arow'
      CALL VIDEONM( 13, 'Arow', NOMFGIF )
C
C     LES COORDONNEES EXTREMES SONT CELLES DE CET OBJET
C     -------------------------------------------------
      CALL MIMXPT( NBCOOR, NBPOI, XYZPOI, COOEXT )
C     CADRE AVEC 20% EN PLUS
      DO I=1,NBCOOR
         EC = ( COOEXT(I,2) - COOEXT(I,1) ) * 1.1 / 2
         CM = ( COOEXT(I,1) + COOEXT(I,2) ) / 2
         COOEXT(I,1) = CM - EC
         COOEXT(I,2) = CM + EC
      ENDDO
C
cccC     DEFINITION DE LA VISEE AXONOMETRIQUE
cccC     ------------------------------------
cccC     POINT VU LE CENTRE DE L'HEXAEDRE ENGLOBANT
ccc      AXOPTV(1) = ( COOEXT(1,1) + COOEXT(1,2) ) * 0.5
ccc      AXOPTV(2) = ( COOEXT(2,1) + COOEXT(2,2) ) * 0.5
ccc      AXOPTV(3) = ( COOEXT(3,1) + COOEXT(3,2) ) * 0.5
cccC     POSITION DE L'OEIL
ccc      AXOEIL(1) = AXOPTV(1) + 4.0 * ( COOEXT(1,2) - COOEXT(1,1) )
ccc      AXOEIL(2) = AXOPTV(2) + 3.0 * ( COOEXT(2,2) - COOEXT(2,1) )
ccc      AXOEIL(3) = AXOPTV(3) + 3.5 * ( COOEXT(3,2) - COOEXT(3,1) )
cccC     LARGEUR/2 et HAUTEUR/2 de la FENETRE
ccc      AXOLAR = (COOEXT(1,2)-COOEXT(1,1) + COOEXT(2,2)-COOEXT(2,1))*0.45
ccc      AXOHAU = ( COOEXT(3,2) - COOEXT(3,1) ) * 0.85
cccC     PAS DE PLAN ARRIERE ET AVANT
ccc      AXOARR = 0
ccc      AXOAVA = 0
ccc      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )
C
C     ECART MAXIMAL MIN-MAX = DIAGONALE DE L'HEXAEDRE
      ECMX = SQRT( (COOEXT(1,2)-COOEXT(1,1))**2 +
     %             (COOEXT(2,2)-COOEXT(2,1))**2 +
     %             (COOEXT(3,2)-COOEXT(3,1))**2 )
C
C     CONSTRUCTION DES FACES FRONTALIERES ET DE LEURS BARYCENTRES
C     ===========================================================
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE
      CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     CREATION DU HACHAGE DES ARETES DES FACES FRONTALIERES DE L'OBJET
      CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %             NTAFOB, MNAFOB, I )
C
C     LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
      MOARFR = MCN( MNAFOB + WOARFR )
C     LA MAJORATION DU NOMBRE DES ARETES FRONTALIERES
      MXARFR = MCN( MNAFOB + WXARFR )
C     LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
      L1ARFR = MCN( MNAFOB + W1ARFR )
C     LE NOMBRE D'ARETES FRONTALIERES DANS LE CHAINAGE
      NBARFR = MCN( MNAFOB + WBARFR )
C
C     RESERVATION DES TABLEAUX NO XYZ POUR LE TRI
      CALL TNMCDC( 'ENTIER',  NBARFR, MNNBFR )
      CALL TNMCDC( 'REEL'  ,3*NBARFR, MNBAFR )
C
C     CALCULER LE BARYCENTRE DE CHAQUE FACE FRONTALIERE
C              LE NUMERO DE LA FACE POUR LE TRI PAR TAS ULTERIEUR
      CALL TRIARFR( MOARFR, L1ARFR, MCN(MNAFOB+WAREFR), XYZPOI,
     %              NBARFR, MCN(MNNBFR), RMCN(MNBAFR) )
C
C     -----------------------------------------------------------
C     OPTIONS DE LA VISEE POUR VOIR L'OBJET
C     -----------------------------------------------------------
 150  IF( INTERA .GT. 0 ) CALL CHOIXFONTE( NPHFCO )
      CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
C     -------------------------------------------
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 150
      ENDIF
C
C     TRACE EFFECTIF DES AXES, ARETES FRONTALIERES ET DES FLECHES
C     ===========================================================
 210  DO NCAS = NCAS0, NCAS1

C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX
C
C        TEMPS DE CALCUL DU VECTEUR NCAS
         TEMPS = TIMES( NCAS )
C
C        TRACE DES AXES 3D
         CALL TRAXE3
C
C        TRACE DES FLECHES VITESSES
         CALL TRVITE3D( NBPOI,    XYZPOI,
     %                  MOARFR,   MXARFR,   L1ARFR, MCN(MNAFOB+WAREFR),
     %                  SEMINVIT, SEMAXVIT, CMVITE,
     %                  NCAS,
     %                  vitx(NCAS)%dptab,
     %                  vity(NCAS)%dptab,
     %                  vitz(NCAS)%dptab,
     %                  VITMAX )
C
C        LE TRACE DU TITRE FINAL
         CALL LEGVIT( KNOMOB, NCAS, VITMOY, VITMIN, VITMAX, CMVITE )
C
C        MISE SUR FICHIER NomfgifBoImage.xwd puis NomfgifNoImage.jpg
C        DE LA PIXMAP de la FENETRE X11 ACTUELLE
         CALL VIDEO1( NOMFGIF, NCAS )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C        FIN DE LA BOUCLE SUR LES CAS
      ENDDO


C     CONSTRUCTION FINALE DU FICHIER.gif
      CALL VIDEOFIN( NOMFGIF )
C
C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ORBITE BOUTON ENFONCE et DEPLACE
            CALL ORBITE1( NOTYEV )
         ELSE
C           ORBITE BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ORBITE3( NOTYEV )
         ENDIF
         IF( NOTYEV .EQ. 0 ) GOTO 150
         GOTO 210
      ELSE
C        POUR LIRE LE TRACE AVANT D'AFFICHER UN MENU
         CALL CLICSO
      ENDIF
      GOTO 150
C
C     SORTIE DU TRACE
C     ===============
 9000 IF( MNBAFR .GT. 0 ) CALL TNMCDS( 'REEL',  NBARFR*3, MNBAFR )
      IF( MNNBFR .GT. 0 ) CALL TNMCDS( 'ENTIER',NBARFR,   MNNBFR )
C
      RETURN
      END
