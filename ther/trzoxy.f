      SUBROUTINE TRZOXY( NDIM,   KNOMOB, MODECO, NBTEMPS, TIME,
     %                   NBTYEL, MNELEM, MNPOGE,
     %                   NDSM,   NTDL,   DEPLAC, AMPLID,
     %                   ZMIN,   ZMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA SURFACE MAILLEE DE L'OBJET 2D
C -----    EN CONSTRUISANT LES TMS XYZSOMMET ET NSEF DES EF SURFACIQUES
C          ET EN PORTANT EN Z LE DEPLACEMENT CALCULE OU L'ERREUR
C          TRACER EN CONTINU LES NDSM CAS POUR VOIR L'ONDE SE PROPAGER
C
C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION DEPLACEMENT_EXACT(t,x,y,z)
C          DANS LE LEXIQUE DES FONCTIONS DE L'UTILISATEUR
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS
C          =1 CE SONT DES DEPLACEMENTS
C          =2 CE SONT DES VECTEURS PROPRES
C NBTEMPS: LE NOMBRE DE TEMPS STOCKES EVENTUELLEMENT ZERO
C TIME   : LES NBTEMPS TEMPS STOCKES
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NDSM   : NOMBRE DE CAS OU SECONDS MEMBRES DU SYSTEME LINEAIRE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C DEPLAC : DEPLACEMENT DES NDSM CAS AVEC NTDL NOEUDS
C AMPLID : COEFFICIENT D'AMPLIFICATION DU DEPLACEMENT VERTICAL
C ZMIN   : DEPLACEMENT MINIMAL DE TOUS LES NOEUDS ET TOUS LES CAS
C ZMAX   : DEPLACEMENT MAXIMAL DE TOUS LES NOEUDS ET TOUS LES CAS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1999
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*90      KNOM
      CHARACTER*4       NOMELE(2)
      INTEGER           NOEUDS(9)
      REAL              TIME(1:NBTEMPS)
      DOUBLE PRECISION  DEPLAC(NTDL,NDSM), PBASE(9), D
C
      IERR   = 0
      MNXYZS = 0
      MNNSEF = 0
      MNZTXY = 0
C
C     SAUVEGARDE DE XYZ MIN MAX ACTUELS
      XMIN0 = COOEXT(1,1)
      XMAX0 = COOEXT(1,2)
      YMIN0 = COOEXT(2,1)
      YMAX0 = COOEXT(2,2)
      ZMIN0 = COOEXT(3,1)
      ZMAX0 = COOEXT(3,2)
      AMPLID = ABS( AMPLID )
C
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
C     ATTENTION: CETTE PROGRAMMATION SOUS ENTEND NOEUDS=POINTS
      NBPOI = MCN(MNPOGE+WNBPOI)
      NBSOM = NBPOI
      CALL DIMCOO( NBPOI, MCN(MNPOGE+WYZPOI), NDIM )
      IF( NDIM .NE. 2 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) ='ERREUR: OBJET NON en 2D'
         ELSE
            KERR(2) ='ERROR: OBJECT NOT in 2D'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     LE NOMBRE DE SOMMETS DES EF DECOUPES EVENTUELLEMENT EN 4
      NBSOM = NBPOI
C
C     CALCUL DU NOMBRE DE TRIANGLES ET QUADRANGLES
C     ============================================
      LDEGRE = 1
      NBTRI  = 0
      NBQUA  = 0
      DO 10 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
         IF( NCOGEL .NE. 3 .AND. NCOGEL .NE. 4 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMOB
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) ='ERREUR: un TYPE EF NON en DIMENSION 2'
            ELSE
               KERR(2) ='ERROR: a FINITE ELEMENT TYPE NOT in 2D'
            ENDIF
            CALL LEREUR
            IERR = 3
            RETURN
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LE NOMBRE D'EF ET DE SOMMETS DU MAILLAGE 2D
         IF( NCOGEL .EQ. 3 ) THEN
            NBTRI = NBTRI + NBELEM
         ELSE
            NBQUA = NBQUA + NBELEM
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LAGRANGE DE DEGRE 1 OU 2?
         IF( NBNOE .GT. 4 ) THEN
            LDEGRE = 2
         ENDIF
C
 10   CONTINUE
C
C     CREATION DU TMC XYZSOMMET DE L'OBJET
C     ====================================
      IF( LDEGRE .EQ. 2 .AND. NBQUA .GT. 0 ) THEN
C        LE BARYCENTRE DU QUAD 2Q2C EST AJOUTE
         NBSOM = NBPOI + NBQUA * (LDEGRE-1)
      ENDIF
      CALL TNMCDC( 'MOTS',  WYZSOM+3*NBSOM, MNXYZS )
      CALL TNMCDC( 'REEL', NBSOM, MNZTXY )
      MNZT = MNZTXY - 1
C
C     LE NOMBRE DE SOMMETS
      MCN( MNXYZS + WNBSOM ) = NBSOM
C
C     LE NOMBRE DE TANGENTES
      MCN( MNXYZS + WNBTGS ) = 0
C
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM
C     ADRESSE MCN DES XYZ DES NOEUDS=POINTS
      MNN = MNPOGE + WYZPOI
C
      XMIN = RINFO('GRAND')
      XMAX =-XMIN
      YMIN = XMIN
      YMAX = XMAX
      DO 70 NS=1,NBPOI
C        L'ABSCISSE
         X   = RMCN( MNN )
         IF( X .LT. XMIN ) XMIN = X
         IF( X .GT. XMAX ) XMAX = X
         RMCN( MNS     ) = X
C        L'ORDONNEE
         X   = RMCN( MNN + 1 )
         IF( X .LT. YMIN ) YMIN = X
         IF( X .GT. YMAX ) YMAX = X
         RMCN( MNS + 1 ) = X
C        LA COTE
         RMCN( MNS + 2 ) = 0.0
C
         MNN = MNN + 3
         MNS = MNS + 3
 70   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNXYZS) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CREATION DU TMS NSEF DE L'OBJET
C     ===============================
C     LE NOMBRE D'EF
      NBEFOB = ( NBTRI + NBQUA ) * LDEGRE**2
      CALL TNMCDC( 'MOTS', WUSOEF+4*NBEFOB, MNNSEF )
C
C     LE TYPE DE L'OBJET  ICI = 3 SURFACE!
      MCN( MNNSEF + WUTYOB ) = 3
C     SURFACE NON FERMEE
      MCN( MNNSEF + WUTFMA ) = 0
C     4 SOMMETS PAR EF
      MCN( MNNSEF + WBSOEF ) = 4
C     0 TANGENTE PAR EF
      MCN( MNNSEF + WBTGEF ) = 0
C     NOMBRE D'EF
      MCN( MNNSEF + WBEFOB ) = NBEFOB
C     NOMBRE D'EF A TGS
      MCN( MNNSEF + WBEFTG ) = 0
C     NOMBRE D'EF AVEC POINTEUR
      MCN( MNNSEF + WBEFAP ) = 0
C     MAILLAGE NON STRUCTURE
      MCN( MNNSEF + WUTYMA ) = 0
C
C     ADRESSE MCN DU NUMERO DES 4 SOMMETS DES EF
      MNN = MNNSEF + WUSOEF - 5
      NBP = NBPOI
C
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM - 4
C
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
         IF( LDEGRE .EQ. 2 .AND. NBNOE .EQ. 8 ) THEN
C
C           RETROUVER LE NUMERO DES INTERPOLATIONS DU DEPLACEMENT
C           ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
            CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &                   NBINVA, NUINVA, NUINTI, NBNDIN )
C           LE DEPLACEMENT EST ICI LA SEULE INCONNUE VARIATIONNELLE
            NOINTE = NUINTI(1)
C           VALEUR DES POLYNOMES DE BASE AU BARYCENTRE DU CARRE UNITE
            CALL INTERP( NOINTE, 0.5D0, 0.5D0, 0D0, K, PBASE )
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
C        ======================================
         DO 95 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, NOEUDS )
C
            IF( LDEGRE .EQ. 1 ) THEN
C
C              1 EF OBJET => 1 EF SURFACE
               MNN = MNN + 4
               DO 80 I=1,NCOGEL
C                 LE NUMERO DE SOMMET DU NOEUD I
                  MCN( MNN + I ) = NOEUDS(I)
  80           CONTINUE
C              LE 4-EME SOMMET D'UN TRIANGLE EST MIS A ZERO
               IF( NCOGEL .EQ. 3 )  MCN( MNN+4 ) = 0
C
            ELSE
C
C              1 EF OBJET => 4 EF SURFACE
               IF( NCOGEL .EQ. 3 ) THEN
C
C                 LE TRIANGLE P2 DONNE 4 SOUS-TRIANGLES P1
C                 ----------------------------------------
C                 TRIANGLE 1 4 6
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(1)
                  MCN( MNN + 2 ) = NOEUDS(4)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 4 2 5
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(4)
                  MCN( MNN + 2 ) = NOEUDS(2)
                  MCN( MNN + 3 ) = NOEUDS(5)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 4 5 6
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(4)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 6 5 3
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(6)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(3)
                  MCN( MNN + 4 ) = 0
C
               ELSE
C
C                 LE QUADRANGLE Q2 DONNE 4 SOUS-QUADRANGLES Q1 AVEC UN BARYCENTR
C                 --------------------------------------------------------------
C                 CONSTRUCTION DES XY DU BARYCENTRE DU QUADRANGLE 2Q2C
                  NBP = NBP + 1
                  NOEUDS(9) = NBP
C
C                 LES 2 COORDONNEES X Y DU BARYCENTRE DU QUADRANGLE 2Q2C
                  MNP = MNS + 3 * NBP
                  DO 84 K=1,2
                     D = 0D0
                     DO 82 L=1,NBPOE
                        D = D + PBASE(L) * RMCN(MNS+3*NOEUDS(L)+K)
 82                  CONTINUE
                     RMCN(MNP+K) = REAL( D )
 84               CONTINUE
                  RMCN(MNP+3) = 0.0
C
C                 QUADRANGLE 1 5 9 8
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(1)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(9)
                  MCN( MNN + 4 ) = NOEUDS(8)
C                 QUADRANGLE 5 2 6 9
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(5)
                  MCN( MNN + 2 ) = NOEUDS(2)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = NOEUDS(9)
C                 QUADRANGLE 8 9 7 4
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(8)
                  MCN( MNN + 2 ) = NOEUDS(9)
                  MCN( MNN + 3 ) = NOEUDS(7)
                  MCN( MNN + 4 ) = NOEUDS(4)
C                 QUADRANGLE 9 6 3 7
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(9)
                  MCN( MNN + 2 ) = NOEUDS(6)
                  MCN( MNN + 3 ) = NOEUDS(3)
                  MCN( MNN + 4 ) = NOEUDS(7)
C
               ENDIF
C
            ENDIF
C
  95     CONTINUE
 100  CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNNSEF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     MISE A JOUR DU TABLEAU COOEXT DU COMMON / XYZEXT /
      COOEXT(1,1) = XMIN
      COOEXT(1,2) = XMAX
      COOEXT(2,1) = XMIN
      COOEXT(2,2) = XMAX
      COOEXT(3,1) = ZMIN * AMPLID
      COOEXT(3,2) = ZMAX * AMPLID
C
C     MISE EN ECRAN POUR LE TRACE DE LA SURFACE=DEPLACEMENT
C     SAUVEGARDE AVEC OU SANS CRITERE DE QUALITE
      LCRIT0 = LCRITR
C     IAVFAC : avec ou non trace des faces
      IAVFA0 = IAVFAC
C     PREDUF : pourcentage de reduction des faces
      PREDU0 = PREDUF
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPAR0 = NEPARF
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUA0 = NCOUAF
C     LE TYPE DE LA VISEE
      NOTYVI = 0
C     PAS DE CRITERE DE QUALITE
      LCRITR = 0
CCC
CCCC     LA PALETTE DES GRIS POUR LE DEPLACEMENT=Z
CCC      CALL PALCDE( 10 )
CCCC     LA PALETTE DES MARRON ROUGE JAUNE POUR LE DEPLACEMENT=Z
CCC      CALL PALCDE( 5 )
CCCC     LES COULEURS PASSENT DU MAGENTA AU ROUGE PUIS AU ROUGE NOIR
CCC      CALL PALCDE( 7 )
C     LES COULEURS PASSENT DU MAGENTA AU BLEU AVEC DE PLUS EN PLUS DE NOIR
      CALL PALCDE( 8 )
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -80.0
      AXOLAT =  40.0
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 0
      CALL VISEE0
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -80.0
      AXOLAT =  40.0
C
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
C
C     LE NUMERO DANS LE LEXIQUE DE L'OBJET KNOMOB
      CALL LXNMNO( NTOBJE, KNOMOB, NUOBJT, I )
C
C     OPTIONS INITIALES DU TRACE
C     IAVFAC : AVEC TRACE des FACES
      IAVFAC = 1
C     PREDUF : pourcentage de reduction des faces
      PREDUF = 0
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPARF = 1
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUAF = NCNOIR
      CALL EFFACE
C
C     ===========================================
C     TRACE DE LA SURFACE AVEC Z=DEPLACEMENT(X,Y)
C     ===========================================
 300  CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9900
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
 310  IF( LORBITE .NE. 0 ) THEN
C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
CCCC     TRACE DU MAILLAGE DANS LE PLAN XOY
CCCC     ----------------------------------
CCCC     ADRESSE MCN DES XYZ DES SOMMETS
CCC      MNS = MNXYZS + WYZSOM
CCC      DO 310 I=1,NBSOM
CCCC        Z EST LE DECALAGE EN Z POUR TRACER LE MAILLAGE DANS UN PLAN // XOY
CCC         RMCN( MNS + 2 ) = DECALZ
CCC         MNS = MNS + 3
CCC 310  CONTINUE
CCC      CALL EFFACE
CCC      CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
CCC      CALL TRAXES
CCC      CALL CLICSO
C
C     LE TRACE DE LA SURFACE AVEC Z=DEPLACEMENT * AMPLID
C     --------------------------------------------------
      DO 900 NCAS=1,NDSM
C
C        LE DEPLACEMENT AMPLIFIE EST REPRESENTE EN Z
         CALL EFFACE
C        ADRESSE MCN DES XYZ DES SOMMETS
         MNS = MNXYZS + WYZSOM
C        LE DEPLACEMENT AMPLIFIE AUX NOEUDS DE L'OBJET
C        SAUF AU BARYCENTRE DES QUADRANGLES DE DEGRE 2
         DO 320 I=1,NBPOI
C           Z EST LE DEPLACEMENT AMPLIFIE
            RMCN( MNS + 2 ) = REAL( DEPLAC(I,NCAS) * AMPLID )
            MNS = MNS + 3
 320     CONTINUE
C
         IF( LDEGRE .EQ. 2 .AND. NBQUA .GT. 0 ) THEN
C
C           CONSTRUCTION DE LE DEPLACEMENT AU BARYCENTRE DES QUADRANGLES 2Q2C
            NBARY = NBPOI
            DO 400 NOTYEL = 1, NBTYEL
C
C              L'ADRESSE DU TABLEAU NPEF"xxxx
               MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C              LE NUMERO DU TYPE DE L'ELEMENT FINI
               NUTYEL = MCN( MNELE + WUTYEL )
C
C              NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
               CALL ELNUCG( NUTYEL, NCOGEL )
               IF( NCOGEL .NE. 4 ) GOTO 400
C
C              LES CARACTERISTIQUES DE L'ELEMENT FINI
               CALL ELNUNM( NUTYEL, NOMELE )
               CALL ELTYCA( NUTYEL )
C
               IF( NBNOE .NE. 8 ) GOTO 400
C
C              RETROUVER LE NUMERO DES INTERPOLATIONS DE LE DEPLACEMENT
C              ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
               CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &                      NBINVA, NUINVA, NUINTI, NBNDIN )
C              LE DEPLACEMENT EST ICI LA SEULE INCONNUE VARIATIONNELLE
               NOINTE = NUINTI(1)
C              VALEUR DES POLYNOMES DE BASE AU BARYCENTRE DU CARRE UNITE
               CALL INTERP( NOINTE, 0.5D0, 0.5D0, 0D0, K, PBASE )
C
C              LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
               NBELEM = MCN( MNELE + WBELEM )
C
               DO 350 NUELEM = 1, NBELEM
C
C                 LES NOEUDS DE L'ELEMENT FINI QUAD 2Q2C
                  CALL EFNOEU( MNELE, NUELEM, NBNDEL, NOEUDS )
C
C                 CONSTRUCTION DE LE DEPLACEMENT AU BARYCENTRE DU QUADRANGLE
                  MNP   = MNXYZS + WYZSOM + 3 * NBARY
                  NBARY = NBARY  + 1
C                 LES 2 COORDONNEES X Y DU BARYCENTRE ONT DEJA ETE CALCULEES
C
C                 LE DEPLACEMENT AU BARYCENTRE
                  D = 0D0
                  DO 330 L=1,NBNOE
                     D = D + PBASE(L) * DEPLAC( NOEUDS(L), NCAS )
 330              CONTINUE
C                 LE DEPLACEMENT AU BARYCENTRE EST PORTEE EN Z ET AMPLIFIE
                  RMCN(MNP+2) = REAL( D * AMPLID )
C
 350           CONTINUE
 400        CONTINUE
         ENDIF
C
CCCC        MISE SUR FICHIER onde**.eps DE CHACUN DES TRACES DE L'ONDE
CCC         CALL XVINITIERPS( 0 )
C
C        TRACE DE L'ONDE POUR CE TEMPS
         CALL TRAXE3
         CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
C
C        LE TRACE DU NOM DE L'OBJET ET DES DEPLACEMENTS EXTREMES
         CALL XVCOULEUR( NCBLEU )
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = 'Z=DEPLACEMENT(X,Y) de l''OBJET: ' // KNOMOB
         ELSE
            KNOM = 'Z=DISPLACEMENT(X,Y) of the OBJECT: ' // KNOMOB
         ENDIF
         I = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
         IF( NBTEMPS .GT. 0 ) THEN
            TEMPS = TIME(NCAS)
         ENDIF
C
         IF( MODECO .EQ. 1 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               KNOM = 'CAS      au TEMPS                : '
            ELSE
               KNOM = 'CASE     at TIME                 : '
            ENDIF
            WRITE( KNOM(5:8),   '(I4)'    ) NCAS
            WRITE( KNOM(19:33), '(G14.6)' ) TEMPS
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
               KNOM = 'VALEUR PROPRE               : '
            ELSE
               KNOM = 'EIGENVALUE                  : '
            ENDIF
            WRITE( KNOM(15:29), '(G15.7)' ) TEMPS
         ENDIF
C
         I = NUDCNB( KNOM )
         IF( LANGAG .EQ. 0 ) THEN
            KNOM(I+1:I+18) = ' DEPLACEMENT MIN= '
         ELSE
            KNOM(I+1:I+19) = ' DISPLACEMENT MIN= '
         ENDIF
         I = NUDCNB( KNOM )
         WRITE( KNOM(I+1:I+14), '(G14.6)' ) ZMIN
C
         I = I + 16
         KNOM(I:I+5) = ' MAX= '
         WRITE( KNOM(I+7:I+20), '(G14.6)' ) ZMAX

C        TRACE DES ANNOTATIONS
         CALL TRFINS( KNOM )
         CALL XVVOIR
         CALL CLICSO
C
CCCC        CONSTRUCTION DU FICHIER onde**.eps
CCCC        ----------------------------------
CCC         KNOM = 'onde    '
CCC         IF( NCAS .LT. 10 ) THEN
CCC            WRITE(KNOM(5:5),'(I1)') NCAS
CCC            I = 5
CCC         ELSE
CCC            WRITE(KNOM(5:6),'(I2)') NCAS
CCC            I = 6
CCC         ENDIF
CCC         CALL XVSAUVERPS( KNOM, I )
C
 900  CONTINUE
C
      IF( NDIM .EQ. 3 .AND. LORBITE .NE. 0 ) GOTO 310
      CALL CLICSO
      GOTO 300
C
C     RESTAURATION DES OPTIONS INITIALES DE TRACE
C     -------------------------------------------
 9900 COOEXT(1,1) = XMIN0
      COOEXT(1,2) = XMAX0
      COOEXT(2,1) = YMIN0
      COOEXT(2,2) = YMAX0
      COOEXT(3,1) = ZMIN0
      COOEXT(3,2) = ZMAX0
      NOTYVI = 0
      ZMIN   = ZMIN0
      ZMAX   = ZMAX0
C     RESTAURATION DES PARAMETRES INITIAUX
      LCRITR = LCRIT0
C     IAVFAC : avec ou non trace des faces
      IAVFAC = IAVFA0
C     PREDUF : pourcentage de reduction des faces
      PREDUF = PREDU0
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPARF = NEPAR0
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUAF = NCOUA0
C
C     DESTRUCTION DES TMC INUTILES
      IF( MNZTXY .GT. 0 ) CALL TNMCDS( 'REEL', NBSOM, MNZTXY )
      IF( MNXYZS .GT. 0 ) CALL TNMCDS( 'MOTS', WYZSOM+3*NBSOM,  MNXYZS )
      IF( MNNSEF .GT. 0 ) CALL TNMCDS( 'MOTS', WUSOEF+4*NBEFOB, MNNSEF )
      RETURN
      END
