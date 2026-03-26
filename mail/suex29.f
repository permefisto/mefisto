      SUBROUTINE SUEX29( NUSUFI, NTLXSU, LADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MODIFIER INTERACTIVEMENT une QUADRAN-TRIANGULATION devenue
C -----    une TRIANGULATION d'une SURFACE 2D ou 3D
C          a l'aide de MENUS et une SOURIS avec 3 BOUTONS pour POINTER
C          des SOMMETS, des ARETES, des TRIANGLES

C ENTREES:
C --------
C NUSUFI : NUMERO DE LA SURFACE FINALE DANS LE LEXIQUE DES SURFACES
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE FINALE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE FINALE
C          CF $MEFISTO/td/d/a_surface__definition

C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DE LA SURFACE FINALE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DE LA SURFACE FINALE
C          CF $MEFISTO/td/d/a___nsef
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE FINALE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE FINALE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          1 SI TRACE EN MODE BATCH I.E. NON INTERACTIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: LEVI-OVSIOUK DEA ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2000
C MODIFS : Alain PERRONNET Labo J-L. LIONS    UPMC PARIS  SEPTEMBRE 2007
C AJOUTS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C AJOUTS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY       MAI 2016
C AJOUTS : Alain PERRONNET Saint PIERRE du PERRAY             AVRIL 2017
C AJOUTS : Alain PERRONNET Saint PIERRE du PERRAY             AVRIL 2020
C AJOUTS : Alain PERRONNET Saint PIERRE du PERRAY           JANVIER 2021
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___trace.inc"
      include"./incl/traaxe.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pilect.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a_surface__definition.inc"

C     SEUILS DE PRECISION POUR IDENTIFIER 2 POINTS ou SOMMETS
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / TRACE1 / PTCOUR(3), MUETR1(7)

      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      (RMCN(1),MCN(1))

C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)

      CHARACTER*24      KNOMPO, KNMSUFI
      CHARACTER*64      KNOM
      CHARACTER*80      LETITRE
      REAL              CADRSA(2,3),CADRS0(1:2,1:3)
      INTEGER           LADEFI(0:*), NUMST0, NBARXF(3), NOSOTQ(4)
      REAL              X0, Y0, XYZPT0(3)

C     LES VARIABLES LOCALES
      REAL              X(12),  Y(12)
      SAVE              NBPOIN, KNOMPO
      DATA              NBPOIN / 0 /
      DATA              KNOMPO / 'P' /

      IF( INTERA .LT. 3 ) THEN
C        MODE INTERACTIF INTERDIT
         IERR = 1
         RETURN
      ENDIF
      IERR = 0

C     LE NOM DE LA SURFACE FINALE A MODIFIER
C     ======================================
      CALL NMOBNU( 'SURFACE', NUSUFI, KNMSUFI )

      MNARFA = 0
      L1ARFA = 0
      L2ARFA = 0

      MXAR1F = 0
      MXAR2F = 0
      MXA2SF = 0
      MXAR2FCOL = 0
      MXAR3F = 0

      MNAR1F = 0
      MNAR2F = 0
      MNA2SF = 0
      MNAR2FCOL=0
      MNAR3F = 0

      MXITDI = 0
      MNITDI = 0
      MXITAN = 0
      MNITAN = 0
      MNNEWS = 0
      MXNEWS = 0

      NOTYVI = 0
      LORBITE= 1
      NETAXE = 0
      NUMST0 = 0

C     NUMERO DU TQ DETRUIT NON INITIALISE
      NOSOTQ(1) = 0

C     LE NUMERO DE LA SURFACE INITIALE
      NUSUIN = LADEFI(WUSUIN)

C     CREE UNE TRIANGULATION NON STRUCTUREE DE LA SURFACE
C     AVEC 8 FOIS PLUS DE PLACES QU'AU DEPART
C     ---------------------------------------------------
      CALL AUGTRIANG( 8,      NTLXSU, NUSUIN,
     %                NTXYZS, MNXYZS, NTNSEF, MNNSEF,
     %                MOXYZS, MONSEF, IERR )
      IF( IERR .NE. 0 ) GOTO 9990


C     MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN
C     IDENTIFIANT LES SOMMETS PROCHES
C     RENUMEROTANT LES SOMMETS
C     ELIMINANT LES EF DESACTIVES ou AYANT 2 SOMMETS DE MEME NUMERO
C     -------------------------------------------------------------
 5    IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNEWS, MNNEWS )
      MXNEWS = MCN(MNXYZS+WNBSOM) + 1
      CALL TNMCDC( 'ENTIER', MXNEWS, MNNEWS )
      CALL MAJXYZNSE( MCN(MNXYZS+WNBSOM), RMCN(MNXYZS+WYZSOM),
     %                MCN(MNNEWS), MCN(MNNSEF+WBSOEF),
     %                MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF) )
      IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNEWS, MNNEWS )

C     MISE A JOUR DES COORDONNEES EXTREMES DE CETTE SURFACE
 10   CALL MAJEXT( MNXYZS )
      CALL DIMCOO( MCN(MNXYZS+WNBSOM), RMCN(MNXYZS+WYZSOM), NDIMLI )


C     INITIALISATION DE LA VISEE INITIALE POUR LES TRACES
C     ---------------------------------------------------
      IF( NOTYVI .EQ. 1 .OR. NOTYVI .EQ. 11 ) THEN
C        VISEE PAR DEFAUT DEJA FAITE
         GOTO 50
      ENDIF

      IF( NDIMLI .LE. 2 ) THEN

C        AJUSTAGE DE LA VISEE EN 2D
C        CADRE DE SAISIE CADRSA(X Y, MIN MAX PAS) EN 2D
         IF( COOEXT(1,1) .LT. RINFO('GRAND') )  THEN
            ECAMAX=MAX(COOEXT(1,2)-COOEXT(1,1), COOEXT(2,2)-COOEXT(2,1))
            IF( ECAMAX .LE. 0. ) ECAMAX = 1.
            ECAMAX = ECAMAX * 0.1
C           ABSCISSES
            CADRS0(1,1) = COOEXT(1,1) - ECAMAX
            CADRS0(1,2) = COOEXT(1,2) + ECAMAX
C           ORDONNEES
            CADRS0(2,1) = COOEXT(2,1) - ECAMAX
            CADRS0(2,2) = COOEXT(2,2) + ECAMAX
C           PAS
            NBPASX = 10
            CADRS0(1,3) = ABS( ( CADRS0(1,2) - CADRS0(1,1) ) / NBPASX )
            NBPASY = 10
            CADRS0(2,3) = ABS( ( CADRS0(2,2) - CADRS0(2,1) ) / NBPASY )
         ELSE
C           ABSCISSES
            CADRS0(1,1) =-1.0
            CADRS0(1,2) = 11.0
C           ORDONNEES
            CADRS0(2,1) =-1.0
            CADRS0(2,2) = 11.0
C           PAS
            CADRS0(1,3) = 1.0
            CADRS0(2,3) = 1.0
         ENDIF

C        RECHERCHE DE LA DIMENSION MAXIMALE ECAMAX DE L'OBJET
         IF( CADRS0(1,1) .GT. CADRS0(1,2) ) THEN
            X(1)        = CADRS0(1,1)
            CADRS0(1,1) = CADRS0(1,2)
            CADRS0(1,2) = X(1)
         ENDIF
         ECAMAX = CADRS0(1,2) - CADRS0(1,1)
         IF( ECAMAX .LE. 0. .OR. ECAMAX .GT. 1E30 ) ECAMAX = 10.

         IF( CADRS0(2,1) .GT. CADRS0(2,2) ) THEN
            X(1)        = CADRS0(2,1)
            CADRS0(2,1) = CADRS0(2,2)
            CADRS0(2,2) = X(1)
         ENDIF
         ECAMAY = CADRS0(2,2) - CADRS0(2,1)
         IF( ECAMAY .LE. 0. .OR. ECAMAY .GT. 1E30 ) ECAMAY = 10.

         RAPECR = (LHPXFE / CYMMPX) / (LAPXFE / CXMMPX)
         IF( ECAMAY .LE. ECAMAX*RAPECR )THEN
            ECAMAY = ECAMAX * RAPECR
         ELSE
            ECAMAX = ECAMAY / RAPECR
         ENDIF

C        LE MILIEU
         X(1) = ( CADRS0(1,1) + CADRS0(1,2) ) * 0.5
         CADRSA(1,1) = X(1) - ECAMAX * 0.5
         CADRSA(1,2) = X(1) + ECAMAX * 0.5
         CADRSA(1,3) = CADRS0(1,3)
C
         Y(1) = ( CADRS0(2,1) + CADRS0(2,2) ) * 0.5
         CADRSA(2,1) = Y(1) - ECAMAY * 0.5
         CADRSA(2,2) = Y(1) + ECAMAY * 0.5
         CADRSA(2,3) = CADRS0(2,3)

C        LE CADRSA OBJET EN UNITES UTILISATEUR
         CALL ISOFENETRE( CADRSA(1,1), CADRSA(1,2),
     %                    CADRSA(2,1), CADRSA(2,2) )

C        MISE A JOUR
         CADRSA(1,1) = XOBMIN
         CADRSA(1,2) = XOBMAX
         CADRSA(2,1) = YOBMIN
         CADRSA(2,2) = YOBMAX

         AXOPTV(1) = ( XOBMIN + XOBMAX ) * 0.5
         AXOPTV(2) = ( YOBMIN + YOBMAX ) * 0.5
         AXOLAR = ( XOBMAX - XOBMIN ) * 0.5
         AXOHAU = ( YOBMAX - YOBMIN ) * 0.5
         NOTYVI = 1

C        TRACE AVANT AFFICHAGE DU MENU
         NOMTCL = 1
         NOMTCLTR = NOMTCL

      ELSE

C        TRACE 3D: VISEE PAR DEFAUT
         CALL VISEE0

C        TRACE AVANT AFFICHAGE DU MENU
         NOMTCL = 3
         NOMTCLTR = NOMTCL

      ENDIF


C     TRACE DU REMPLISSAGE DES FACES
      IAVFAC = 1

C     TRACE DES ARETES DES FACES EN NOIR
      IAVARE = 1
      NCOUAF = NCNOIR


C     DEBUT DES MODIFICATIONS DE LA TRIANGULATION
C     -------------------------------------------
 50   IERR = 0

C     L'ECRAN EST EFFACE
      CALL EFFACE
      CALL EFFACEMEMPX

C     REMISE A ZERO DU NOMBRE DES ITEMS
      CALL ITEMS0

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'SURFACE ACTUELLE a MODIFIER:' // KNMSUFI
      ELSE
         KNOM = 'MODIFY the ACTUAL SURFACE:' // KNMSUFI
      ENDIF


C     TRAITEMENT SELON LA DIMENSION DE L'ESPACE DE LA SURFACE
      IF( NDIMLI .GE. 3 ) THEN
C        SURFACE en 3D
         GOTO 300         
      ELSE
C        SURFACE en 2D
         GOTO 80
      ENDIF

C     ==================================================================
C     TRACER en 2D la TRIANGULATION DE LA SURFACE 2D ACTUELLE
C     ==================================================================
C     INITIALISATION DU ZOOM DEPLACEMENT 2D
 65   LORBITE = 1
      CALL ZOOM2D0( NOTYEV )
      GOTO 80

C     ZOOM 2D OU TRANSLATION ACTIFS
 70   CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) THEN
         LORBITE = 0
         GOTO 90
      ENDIF

C     TRACE DES AXES 2D ET DES FACES DE LA TRIANGULATION 2D
 80   CALL TRAXE2
      CALL T21FAC( KNOM, NUSUFI, MNNSEF, MNXYZS )
      LETITRE = 'SURFACE ' // KNOM
      CALL TRFINS( LETITRE )

C     REPRISE DE TRANSLATION ZOOM 2D
      IF( LORBITE .NE. 0 ) GOTO 70


C     LECTURE DU MENU DES MODIFICATIONS 2D POSSIBLES
C     ==============================================
 90   CALL LIMTCL( 'modifm2d' , NOMTCL )
      IF( NOMTCL .LE.  0 ) GOTO 9900
      IF( NOMTCL .EQ. 90 ) GOTO 9900
      IF( NOMTCL .EQ. 80 ) GOTO  280

      GOTO( 65,  65, 130, 140, 150, 160, 170, 180, 190, 200,
     %     210, 220, 230, 240, 250, 260,  90,  90,  90 ), NOMTCL


C     DEPLACER UN SOMMET DE LA TRIANGULATION
C     --------------------------------------
 130  CALL DEPLST2D( MNXYZS, NUMST0, X0, Y0, IERR )
      GOTO 50

C     RETOUR aux XY INITIALES du DERNIER SOMMET DEPLACE
C     -------------------------------------------------
 140  IF( NUMST0 .GT. 0 ) THEN
         MN = MNXYZS + WYZSOM + 3*NUMST0 - 3
         RMCN(MN  ) = X0
         RMCN(MN+1) = Y0
         RMCN(MN+2) = 0
         NUMST0 = 0
      ENDIF
      GOTO 50

C     DETRUIRE UN SOMMET 2D DE LA TRIANGULATION
C     -----------------------------------------
 150  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      CALL SUPPRP2D( NX, NY, MNXYZS, MNNSEF,  NUMST, IERR )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 50

C     AJOUTER UN SOMMET 2D AU TRIANGULATION
C     -------------------------------------
 160  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      CALL AJOUTPT( NX, NY, NDIMLI, MNXYZS, MNNSEF,
     %              MOXYZS, MONSEF, IERR )
      GOTO 50

C     COUPER une ARETE 2D EN AJOUTANT SON MILIEU AUX SOMMETS
C     ------------------------------------------------------
 170  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      CALL AJOUTAR( NX, NY, NDIMLI, MNXYZS, MNNSEF,
     %              MOXYZS, MONSEF, IERR )
      IF( IERR .NE. 0 ) GOTO 50
      GOTO 50

C     ECHANGER une DIAGONALE de 2 TRIANGLES
C     -------------------------------------
C     CLIC D'UN SOMMET SUR L'ECRAN
 180  CALL SAIPTC( NOTYEV,  NX, NY, NOCHAR )
      IF( NOTYEV .LE. 0 ) GOTO 50
      CALL ECHDIAPT( NX, NY, NDIMLI, MNXYZS, MNNSEF, IERR )
      GOTO 50

C     TRIANGULATION DELAUNAY PAR ECHANGE des DIAGONALES de 2 TRIANGLES
C     ----------------------------------------------------------------
 190  CALL TRIDEL( NDIMLI, MNXYZS, MNNSEF, IERR )
      GOTO 50

C     3 SOMMETS CLIQUES AJOUTENT UN TRIANGLE NUTQCR
C     ---------------------------------------------
 200  CALL S3TRIA( NDIMLI, MNXYZS, MONSEF, MNNSEF, NUTQCR, IERR )
      GOTO 50

C     1 POINT EXTERIEUR et 2 SOMMETS CLIQUES AJOUTENT UN TRIANGLE
C     -----------------------------------------------------------
 210  CALL S2P1TRIA( NDIMLI, MOXYZS, MNXYZS, MONSEF, MNNSEF, IERR )
      GOTO 50

C     IDENTIFIER 2 SOMMETS EN UN SEUL POUR FAIRE UNE COUTURE DU
C     TRIANGULATION, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS
C     -------------------------------------------------------------
 220  CALL SAIPTC( NOTYEV, NX,  NY,  NOCHAR )
      CALL SAIPTC( NOTYEV, NX1, NY1, NOCHAR )
      CALL IDENT2ST( 2, MNXYZS, MNNSEF, NUMST )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 50

C     IDENTIFIER 1 SOMMET A SON PLUS PROCHE SOMMET POUR FAIRE UNE COUTURE
C     de la TRIANGULATION, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS
C     -------------------------------------------------------------------
 230  CALL SAIPTC( NOTYEV, NX,  NY,  NOCHAR )
      CALL IDENT1STPP( 2, MNXYZS, MNNSEF, NUMST )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 50

C     SI FONCTION UTILISATEUR TAILLE_IDEALE(x,yz) ou EDGE_LENGTH(x,y,z)
C     ALORS ADAPTER LA TAILLE DES TRIANGLES A CES VALEURS
C     -----------------------------------------------------------------
 240  CALL AJLGARTI( MNXYZS, MNNSEF, MOXYZS, MONSEF, IERR )
      GOTO 50

C     DETRUIRE un TRIANGLE CLIQUE
C     ---------------------------
 250  CALL TUERTQ2D( MNXYZS, MNNSEF, IERR )
      GOTO 50

C     REGENERER le DERNIER TRIANGLE DETRUIT
C     -------------------------------------
 260  CALL REGEN1TQ( MNNSEF, NOSOTQ )
      GOTO 50

C     SAUVER la TRIANGULATION sur le FICHIER xyznsef.s.NomSurface
C     -----------------------------------------------------------
 280  CALL SAVXYZNSEF( 3, KNMSUFI, MNXYZS, MNNSEF, IERR )
      GOTO 50


C     ================================================================
C     CONSTRUCTION DES ARETES et FACES VISIBLES de la TRIANGULATION 3D
C     ================================================================
 300  IERR = 0

C     RECONSTRUCTION DU TABLEAU DES ARETES DES FACES DE LA SURFACE
      MXFAAR = 6
      IF( MNARFA.GT.0 ) CALL TNMCDS('ENTIER',L1ARFA*L2ARFA,MNARFA)
      CALL GEARFA( RMCN(MNXYZS+WYZSOM), MCN(MNNSEF+WBSOEF),
     %              MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF), MXFAAR,
     %             L1ARFA, L2ARFA, MNARFA, NBARXF, IERR )
      IF( IERR .LT. 0 ) THEN
         NOMTCL = 2
         GOTO 310
      ENDIF

C     RECONSTRUCTION DES TABLEAUX DU NO DANS LARETE DES ARETES SELON
C     LE NOMBRE DE FACES AUXQUELLES ELLES SONT ADJACENTES
 302  CALL LAR123F( MCN(MNXYZS+WNBSOM),  RMCN(MNXYZS+WYZSOM),
     %              MCN(MNNSEF+WBSOEF),   MCN(MNNSEF+WBEFOB),
     %              MCN(MNNSEF+WUSOEF),
     %              L1ARFA,    L2ARFA,    MCN(MNARFA),
     %              MXAR1F,    NBAR1F,    MNAR1F,
     %              MXAR2F,    NBAR2F,    MNAR2F,
     %              MXAR3F,    NBAR3F,    MNAR3F,
     %              MXA2SF,    NBA2SF,    MNA2SF,
     %              MXAR2FCOL, NBAR2FCOL, MNAR2FCOL, IERR  )

      IF( NOMTCL .GT. 4 ) THEN

C        TRACE POUR VOIR LES MODIFICATIONS DE LA TRIANGULATION
C        CONSTRUCTION DES ITEMS SOMMETS et FACES VISIBLES dans la FENETRE
C        ----------------------------------------------------------------
         LORBITE = 1
         IF( NOMTCLTR .EQ. 4 ) NOMTCLTR = 3
         CALL TRITEMV0( NOMTCLTR, KNMSUFI, NUSUFI, MNXYZS, MNNSEF,
     %                  L1ARFA,   L2ARFA,  MNARFA, NBAR1F, MNAR1F,
     %                  NBAR2FCOL,MNAR2FCOL, NBA2SF, MNA2SF,
     %                  NBAR2F, MNAR2F, NBAR3F, MNAR3F,
     %                  MXITDI, MNITDI, MXITAN, MNITAN )

      ENDIF


C     LECTURE DU MENU DES MODIFICATIONS 3D POSSIBLES
C     ==============================================
 305  CALL LIMTCL( 'modifm3d' , NOMTCL )
      IF( NOMTCL .LT. 0  ) GOTO 9900
      IF( NOMTCL .EQ. 90 ) GOTO 9900
      IF( NOMTCL .EQ. 80 ) GOTO  800
      IF( NOMTCL .LE. 4  ) NOMTCLTR = NOMTCL

      GOTO( 310, 310, 310, 310, 350, 360, 370, 380, 390, 400,
     %      410, 420, 430, 440, 450, 460, 470, 480, 490, 495,
     %      510, 305, 305, 305, 305, 305, 305, 305, 305, 305,
     %      610, 620, 630, 640, 650, 660, 670, 305, 305, 305 ), NOMTCL


C     1 ou 2 ou 3: DEPLACER(Bouton1) ORBITER(B2) ZOOMER(B3) le TRIANGULATION
C     TRACER LES ITEMS VISIBLES: SOMMETS, ARETES 1F, 3F, 2SF, 2FCOL, FACES
C     ----------------------------------------------------------------------
 310  LORBITE = 1
      CALL TRITEMV0( NOMTCLTR, KNMSUFI, NUSUFI, MNXYZS, MNNSEF,
     %               L1ARFA,   L2ARFA,  MNARFA, NBAR1F, MNAR1F,
     %               NBAR2FCOL,MNAR2FCOL, NBA2SF, MNA2SF,
     %               NBAR2F, MNAR2F, NBAR3F, MNAR3F,
     %               MXITDI, MNITDI, MXITAN, MNITAN )
      GOTO 305

C     DEPLACER UN SOMMET 3D DE LA TRIANGULATION
C     -----------------------------------------
 350  CALL DEPLST3D( MCN(MNXYZS+WNBSOM), RMCN(MNXYZS+WYZSOM),
     %               MCN(MNNSEF+WBEFOB),  MCN(MNNSEF+WUSOEF),
     %               NUMST0, XYZPT0, IERR )
      GOTO 302

C     DEPLACER LE SOMMET NOST AU BARYCENTRE DES BARYCENTRES DES EF
C     DONT IL EST UN SOMMET
C     ------------------------------------------------------------
 360  CALL DEPLSTBA( RMCN(MNXYZS+WYZSOM),
     %                MCN(MNNSEF+WBEFOB),  MCN(MNNSEF+WUSOEF),
     %               NUMST0, XYZPT0, IERR )
      GOTO 302

C     RETOUR aux XYZ INITIALES du DERNIER SOMMET DEPLACE
C     --------------------------------------------------
 370  IF( NUMST0 .GT. 0 ) THEN
         MN = MNXYZS + WYZSOM + 3*NUMST0 - 3

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'suex29  : RESTAURATION du SOMMET',NUMST0,' de',
     %           (RMCN(MN+K),K=0,2),' a',(XYZPT0(K),K=1,3)
         ELSE
            PRINT*,'suex29  : RESTORATION of VERTEX',NUMST0,' from',
     %           (RMCN(MN+K),K=0,2),' to',(XYZPT0(K),K=1,3)
         ENDIF
         RMCN(MN  ) = XYZPT0(1)
         RMCN(MN+1) = XYZPT0(2)
         RMCN(MN+2) = XYZPT0(3)
         NUMST0 = 0
      ENDIF
      GOTO 302

C     CLIC un SOMMET et LE SUPPRIMER et ses EF (NON RECUPERABLE)
C     ----------------------------------------------------------
 380  CALL SUPTQ1ST( MNXYZS, MNNSEF, NUMST, IERR )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 5

C     AJOUTER un SOMMET DANS UN TRIANGLE CLIQUE
C     -----------------------------------------
 390  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      CALL AJOUTPT( NX, NY, NDIMLI, MNXYZS, MNNSEF,
     %              MOXYZS, MONSEF, IERR )
      GOTO 5

C     COUPER une ARETE CLIQUEE COMMUNE a 2 TRIANGLES -> 4 TRIANGLES
C     -------------------------------------------------------------
 400  CALL CL2TR4TR( MNXYZS, MNNSEF, L1ARFA, L2ARFA, MCN(MNARFA),
     %               MOXYZS, MONSEF, IERR )
      GOTO 5

C     ECHANGER une DIAGONALE CLIQUEE de 2 TRIANGLES -> 2 TRIANGLES
C     ------------------------------------------------------------
 410  CALL ECHDIAP3( MNXYZS, MNNSEF, L1ARFA, L2ARFA, MCN(MNARFA),
     %               NUMAR,  IERR )
      GOTO 5

C     MODIFIER vers une TRIANGULATION DELAUNAY
C     ----------------------------------------
 420  CALL TRIDEL( NDIMLI, MNXYZS, MNNSEF, IERR )
      GOTO 5

C     CREER UN TRIANGLE avec 3 SOMMETS CLIQUES
C     ----------------------------------------
 430  CALL S3TRIA( NDIMLI, MNXYZS, MONSEF, MNNSEF, NUTQCR, IERR )
      GOTO 5
C
C     CREER UN TRIANGLE avec 2 SOMMETS CLIQUES + 1 POINT XYZ EXTERNE
C     --------------------------------------------------------------
 440  CALL S2P1TRIA( NDIMLI, MOXYZS, MNXYZS, MONSEF, MNNSEF, IERR )
      GOTO 5

C     IDENTIFIER 2 SOMMETS EN UN SEUL POUR FAIRE UNE COUTURE de la
C     TRIANGULATION, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS
C     -------------------------------------------------------------
 450  CALL IDENT2ST( 3, MNXYZS, MNNSEF, NUMST )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 5

C     IDENTIFIER 1 SOMMET A SON PLUS PROCHE SOMMET POUR FAIRE UNE COUTURE
C     DE LA TRIANGULATION, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS
C     -------------------------------------------------------------------
 460  CALL IDENT1STPP( 3, MNXYZS, MNNSEF, NUMST )
      IF( NUMST .EQ. NUMST0 ) NUMST0 = 0
      GOTO 5

C     SI FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z) ou EDGE_LENGTH(x,y,z)
C     ALORS ADAPTER LA TAILLE DES TRIANGLES A CES VALEURS
C     ------------------------------------------------------------------
 470  CALL AJLGARTI( MNXYZS, MNNSEF,  MOXYZS, MONSEF, IERR )
      GOTO 5

C     DETRUIRE UN QUADRANGLE 3D ou TRIANGLE 3D CLIQUE
C     -----------------------------------------------
 480  CALL TUERTQ3D( MNXYZS, MNNSEF, NOSOTQ, IERR )
      GOTO 5

C     REGENERER le DERNIER QUADRANGLE ou TRIANGLE DETRUIT
C     ---------------------------------------------------
 490  CALL REGEN1TQ( MNNSEF, NOSOTQ )
      GOTO 5

C     CREER TOUT Q-TRIANGLE avec 2 ARETES dans 1 QT et un SOMMET COMMUN
C     -----------------------------------------------------------------
 495  CALL CREETQ1FC( L1ARFA, MNARFA, NBAR1F, MNAR1F,
     %                MNNSEF, MONSEF, NBTRAJ )
      GOTO 5

C     RETIRER et CREER UNE SURFACE de TRIANGLES CLIQUES
C     -------------------------------------------------
 510  CALL CRSFTRCL( KNMSUFI, MNXYZS, MNNSEF )
      GOTO 5

C     DETRUIRE TOUT Q-TRIANGLE AYANT TOUTES SES ARETES DANS 1F LUI-MEME=ISOLE
C     -----------------------------------------------------------------------
 610  CALL TUERTA1F( MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF),
     %               L1ARFA, L2ARFA, MCN(MNARFA), NBFSUP )
      GOTO 5

C     DETRUIRE TOUT Q-TRIANGLE AYANT au MOINS 2ARETES dans 1 QT (NON RECUP)
C     ---------------------------------------------------------------------
 620  CALL TUERTQA1F( 2,  L1ARFA, L2ARFA, MNARFA, NBAR1F, MNAR1F,
     %                MNNSEF, IERR )
      GOTO 5

C     DETRUIRE TOUT Q-TRIANGLE AYANT au MOINS 1 ARETE dans 1 QT (NON RECUP)
C     ---------------------------------------------------------------------
 630  CALL TUERTQA1F( 1, L1ARFA, L2ARFA, MNARFA, NBAR1F, MNAR1F,
     %                MNNSEF, IERR )
      GOTO 5

C     DETRUIRE TOUT QUADRANGLE ou TRIANGLE ayant une ARETE dans +2 QT
C     ---------------------------------------------------------------
 640  CALL TUERTQA3F( L1ARFA, MNARFA, NBAR3F, MNAR3F, MNNSEF, IERR )
      GOTO 5

C     DETRUIRE TOUTES LES FACES des ARETES de RAPPORT des SURFACES
C     ADJACENTES TROP FAIBLE
C     ------------------------------------------------------------
 650  CALL TUER2FCOL( L1ARFA, MNARFA, NBAR2FCOL, MNAR2FCOL, MNNSEF )
      GOTO 5

C     DETRUIRE TOUTES LES FACES des ARETES de RAPPORT des SURFACES
C     ADJACENTES TROP FAIBLE
C     ------------------------------------------------------------
 660  CALL TUERA2SF( L1ARFA, MNARFA, NBA2SF, MNA2SF, MNNSEF )
      GOTO 5

C     DETRUIRE TOUT TRIANGLE de QUALITE < 0.03
C     ----------------------------------------
 670  CALL TUERTRQM( RMCN(MNXYZS+WYZSOM), 0.03,
     %                MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF) )
      GOTO 5

C     SAUVER la TRIANGULATION sur le FICHIER xyznsef.s.NomSurface
C     ===========================================================
C     MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN
C     . IDENTIFIANT LES SOMMETS PROCHES
C     . RENUMEROTANT LES SOMMETS
C     . ELIMINANT LES EF DESACTIVES ou AYANT 2 SOMMETS DE MEME NUMERO
 800  MXNEWS = MCN(MNXYZS+WNBSOM) + 1
      CALL TNMCDC( 'ENTIER', MXNEWS, MNNEWS )
      CALL MAJXYZNSE( MCN(MNXYZS+WNBSOM), RMCN(MNXYZS+WYZSOM),
     %                MCN(MNNEWS), MCN(MNNSEF+WBSOEF),
     %                MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF) )
      IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNEWS, MNNEWS )

      CALL SAVXYZNSEF( 3, KNMSUFI, MNXYZS, MNNSEF, IERR )
      GOTO 300

C     CONFIRMER LA SORTIE DES MODIFICATIONS?
C     ======================================
 9900 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Voulez vous CONTINUER les MODIFICATIONS?'
      ELSE
         KERR(1) = 'Do you want to CONTINUE the MODIFICATIONS?'
      ENDIF
      CALL LERESU
      CALL LIMTCL( 'non_oui', NONOUI )
      IF( NONOUI .EQ. 1 ) THEN
C        OUI: REINITIALISATION
         GOTO 10
      ENDIF

C     SAUVER la TRIANGULATION et QUITTER les MODIFICATIONS
C     ====================================================
      WRITE(NFFRAP,*) '@; '
      LORBITE = 0

C     MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN
C     IDENTIFIANT LES SOMMETS PROCHES
C     RENUMEROTANT LES SOMMETS
C     ELIMINANT LES EF DESACTIVES ou AYANT 2 SOMMETS DE MEME NUMERO
C     =============================================================
      MXNEWS = MCN(MNXYZS+WNBSOM) + 1
      CALL TNMCDC( 'ENTIER', MXNEWS, MNNEWS )
      CALL MAJXYZNSE( MCN(MNXYZS+WNBSOM), RMCN(MNXYZS+WYZSOM),
     %                MCN(MNNEWS), MCN(MNNSEF+WBSOEF),
     %                MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF) )
      IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNEWS, MNNEWS )


 9990 IF( MNITAN .GT. 0 ) CALL TNMCDS( 'ENTIER', MXITAN, MNITAN )
      IF( MNITDI .GT. 0 ) CALL TNMCDS( 'REEL'  , MXITDI, MNITDI )
      IF( MNAR3F .GT. 0 ) CALL TNMCDS( 'ENTIER', MXAR3F, MNAR3F )
      IF( MNA2SF .GT. 0 ) CALL TNMCDS( 'ENTIER', MXA2SF, MNA2SF )
      IF( MNAR2FCOL.GT.0) CALL TNMCDS( 'ENTIER', MXAR2FCOL, MNAR2FCOL )
      IF( MNAR2F .GT. 0 ) CALL TNMCDS( 'ENTIER', MXAR2F, MNAR2F )
      IF( MNAR1F .GT. 0 ) CALL TNMCDS( 'ENTIER', MXAR1F, MNAR1F )
      IF( MNARFA .GT. 0 ) CALL TNMCDS( 'ENTIER', L1ARFA*L2ARFA, MNARFA )

      RETURN
      END
