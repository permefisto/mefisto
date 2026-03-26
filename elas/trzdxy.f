      SUBROUTINE TRZDXY( NOFOTI, NDIM,   KNOMOB, MODECO,
     %                   NBTYEL, MNELEM, MNPOGE,
     %                   NCAS,   NDSM,   NBNOEU, DEPLAC,
     %                   DEPMIN, NOEMIN, DEPMAX, NOEMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA SURFACE MAILLEE DE L'OBJET 2D
C -----    EN CONSTRUISANT LES TMS XYZSOMMET ET NSEF DES EF SURFACIQUES
C          ET EN PORTANT EN Z L'ERREUR ABSOLUE SUR LA NORME ENTRE
C          LE DEPLACEMENT_EXACT(t,x,y,z,nocomp) ET LE DEPLACEMENT CALCULE
C
C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION DEPLACEMENT_EXACT(t,x,y,z,nocomp)
C          DANS LE LEXIQUE DES FONCTIONS DE L'UTILISATEUR
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS
C          =1 CE SONT DES DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS   : NUMERO DU CAS A TRACER
C NDSM   : NOMBRE DE CAS OU SECONDS MEMBRES DU SYSTEME LINEAIRE
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C DEPLAC : DEPLACEMENT DES NDSM CAS AUX NBNOEU NOEUDS
C DEPMIN : NORME MINIMALE DANS R2 DE L'ERREUR  SUR LE DEPLACEMENT
C NOEMIN : NO DU NOEUD OU LA NORME DE L'ERREUR SUR LE DEPLACEMENT EST MINIMAL
C DEPMAX : NORME MAXIMALE DANS R2 DE L'ERREUR  SUR LE DEPLACEMENT
C NOEMAX : NO DU NOEUD OU LA NORME DE L'ERREUR SUR LE DEPLACEMENT EST MAXIMAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1998
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
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*90      KNOM
      CHARACTER*4       NOMELE(2)
      DOUBLE PRECISION  DEPLAC(NDIM,NBNOEU,NDSM)
C
      DOUBLE PRECISION  PBASE(9), D, DEX(3), DPARAF(5), DBLE
      INTEGER           NOEUDS(9)
      INTRINSIC         REAL
C
      IERR   = 0
      MNXYZS = 0
      MNNSEF = 0
      MNZTXY = 0
      NBEFOB = 0
C
C     SAUVEGARDE DE XYZ MIN MAX ACTUELS
      XMIN0 = COOEXT(1,1)
      XMAX0 = COOEXT(1,2)
      YMIN0 = COOEXT(2,1)
      YMAX0 = COOEXT(2,2)
      ZMIN0 = COOEXT(3,1)
      ZMAX0 = COOEXT(3,2)
      NOEMI0 = NOEMIN
      NOEMA0 = NOEMAX
      DEPMIN0= DEPMIN
      DEPMAX0= DEPMAX
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
            KERR(2) ='ERREUR: OBJET NON EN 2D'
         ELSE
            KERR(2) ='ERROR: NOT 2D OBJECT'
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
               KERR(2) ='ERREUR: OBJET NON EN 2D'
            ELSE
               KERR(2) ='ERROR: NOT 2D OBJECT'
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
      DEPMIN = XMIN
      DEPMAX =-DEPMIN
      DEPEMX = DEPMAX
C
C     BOUCLE SUR LES POINTS=NOEUDS DU MAILLAGE
C     CALCUL DE LA NORME DANS R2 DE L'ERREUR SUR LE DEPLACEMENT
C     =========================================================
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
C        CALCUL DU DEPLACEMENT_EXACT EN CE NOEUD
         DPARAF(1) = TEMPS
         DPARAF(2) = RMCN(MNS)
         DPARAF(3) = RMCN(MNS+1)
         DPARAF(4) = 0D0
         DO 65 K=1,NDIM
C           DEPLACEMENT_EXACT(TEMPS,X,Y,Z,NC)
            DPARAF(5) = K
            CALL FONVAL( NOFOTI, 5, DPARAF, NCODEV, DEX(K) )
            IF( NCODEV .LE. 0 ) GOTO 9999
 65      CONTINUE
C
C        LA NORME R2 DU DEPLACEMENT MAXIMAL
         DEPEMX = REAL( MAX( DBLE(DEPEMX), SQRT(DEX(1)**2+DEX(2)**2) ) )
C
C        LA NORME R2 DE L'ERREUR SUR LE DEPLACEMENT
         D = SQRT( (DEX(1)-DEPLAC(1,NS,NCAS))**2
     %           + (DEX(2)-DEPLAC(2,NS,NCAS))**2 )
         RMCN( MNZT + NS ) = REAL( D )
C
C        LE MINIMUM ET MAXIMUM
         IF( D .LT. DEPMIN ) THEN
            NOEMIN = NS
            DEPMIN = REAL( D )
         ENDIF
         IF( D .GT. DEPMAX ) THEN
            NOEMAX = NS
            DEPMAX = REAL( D )
         ENDIF
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
C     MISE A L'ECHELLE DE Z POUR AVOIR LE MAX DE L'ECART EN X ET Y
      X = DEPMAX - DEPMIN
      IF( X .EQ. 0 ) THEN
C        DEPLACEMENT CONSTANTE
         X = 1.0
      ENDIF
      DECALZ  = MAX( XMAX-XMIN, YMAX-YMIN ) / 1.8
      ECHELLE = DECALZ / X
C
C     COTE DU MAILLAGE PAR RAPPORT A LA SURFACE Z=DEPLACEMENT(X,Y)
      DECALZ  = DECALZ * (-0.17)
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
C
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
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
C                 CONSTRUCTION DU BARYCENTRE DU QUADRANGLE
                  NBP = NBP + 1
                  NOEUDS(9) = NBP
C                 LES 2 COORDONNEES X Y DU BARYCENTRE
                  MNP = MNS + 3 * NBP
                  DO 84 K=1,2
                     D = 0D0
                     DO 82 L=1,NBPOE
                        D = D + PBASE(L) * RMCN(MNS+3*NOEUDS(L)+K)
 82                  CONTINUE
                     RMCN(MNP+K) = REAL( D )
 84               CONTINUE
C
C                 LE DEPLACEMENT AU BARYCENTRE NOEUD=POINT AJOUTE
                  DPARAF(1) = TEMPS
                  DPARAF(2) = RMCN(MNP+1)
                  DPARAF(3) = RMCN(MNP+2)
                  DPARAF(4) = 0D0
                  DO 87 K=1,NDIM
C                    DEPLACEMENT EXACT
                     DPARAF(5) = K
                     CALL FONVAL( NOFOTI, 5, DPARAF, NCODEV, DEX(K) )
                     IF( NCODEV .LE. 0 ) GOTO 9999
C                    DEPLACEMENT APPROCHE
                     D = 0D0
                     DO 86 L=1,NBPOE
                        D = D + PBASE(L) * DEPLAC(K,NOEUDS(L),NCAS)
 86                  CONTINUE
C                    ERREUR SUR CETTE COMPOSANTE DU DEPLACEMENT
                     DEX(K) = DEX(K) - D
 87               CONTINUE
                  D = SQRT( DEX(1)**2 + DEX(2)**2 )
                  RMCN(MNZT+NBP) = REAL( (D-DEPMIN) * ECHELLE )
C
C                 CALCUL DES 2 COORDONNEES X Y DU BARYCENTRE DU QUADRANGLE 2Q2C
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
C     MISE A L'ECHELLE DE L'ERREUR AUX POINTS=NOEUDS
C     ==============================================
      ERRMX = 0
      DO 111 I=1,NBPOI
         ERR = RMCN(MNZT+I)-DEPMIN
         RMCN(MNZT+I) = (RMCN(MNZT+I)-DEPMIN) * ECHELLE
         IF( ERRMX .LT. ERR ) ERRMX = ERR
 111  CONTINUE
C
C     MISE A JOUR DU TABLEAU COOEXT DU COMMON / XYZEXT /
      COOEXT(1,1) = XMIN
      COOEXT(1,2) = XMAX
      COOEXT(2,1) = XMIN
      COOEXT(2,2) = XMAX
      COOEXT(3,1) = (DEPMIN-DEPMIN) * ECHELLE + DECALZ
      COOEXT(3,2) = (DEPMAX-DEPMIN) * ECHELLE
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
C     IAVELO sans eloignement
      IAVEL0 = IAVELO
C     LE TYPE DE LA VISEE
      NOTYVI = 0
C     PAS DE CRITERE DE QUALITE
      LCRITR = -1
      IAVELO = 0
C     TRACE DES AXES 3D
      NTRAXZ = 1
      ZMIAXZ = 0
      ZMXAXZ = ERRMX
C     LA PALETTE ORANGE ROUGE POUR L'ERREUR ABSOLUE
      CALL PALCDE( 5 )
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -70.0
      AXOLAT =  20.0
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 0
      CALL VISEE0
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -70.0
      AXOLAT =  20.0
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
C
C     ===============================================
C     TRACE DE LA SURFACE AVEC Z=||DEPLACEMENT(X,Y)||
C     ===============================================
C
C     OPTIONS DE LA VISEE
 300  CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9900
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
 500  IF( LORBITE .NE. 0 ) THEN
C           ORBITE OU ZOOM OU TRANSLATION ACTIFS
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
C     LE NUMERO DANS LE LEXIQUE DE L'OBJET KNOMOB
      CALL LXNMNO( NTOBJE, KNOMOB, NUOBJT, I )
C
C     TRACE DES AXES 3D
      CALL TRAXE3
C
C     TRACE FIL DE FER DES ARETES DU MAILLAGE DANS LE PLAN XOY
C     --------------------------------------------------------
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM
      DO 310 I=1,NBSOM
C        Z EST LE DECALAGE EN Z POUR TRACER LE MAILLAGE DANS UN PLAN // XOY
         RMCN( MNS + 2 ) = DECALZ
         MNS = MNS + 3
 310  CONTINUE
cccC     IAVFAC : avec trace des faces
ccc      IAVFAC = 1
cccC     PREDUF : pourcentage de reduction des faces
ccc      PREDUF = 0
cccC     NEPARF : nombre d'epaisseurs des ARETES des FACES
ccc      NEPARF = 1
cccC     NCOUAF : COULEUR des ARETES des FACES
ccc      NCOUAF = NCBLAN
      CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
C
C     LE TRACE DE LA SURFACE AVEC Z=DEPLACEMENT * ECHELLE POUR NCAS
C     -------------------------------------------------------------
C     L'ERREUR SUR LE DEPLACEMENT EST IMPOSE EN Z AVEC UNE
C     MISE A L'ECHELLE DE LA PLUS GRANDE DIFFERENCE DE COORDONNEE
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM
      DO 320 I=1,NBSOM
C        Z EST LE DEPLACEMENT MIS A L'ECHELLE
         RMCN( MNS + 2 ) = RMCN( MNZT + I )
         MNS = MNS + 3
 320  CONTINUE
C     NCOUAF : COULEUR des ARETES des FACES SI TRACE DEMANDE
      NC     = NCOUAF
      NCOUAF = NCOAPL
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPARF = 0
      CALL XVEPAISSEUR( NEPARF )
C     IAVFAC : avec trace des faces
      IAVFAC = 1
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUAF = NCNOIR
      CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
      NCOUAF = NC
C
C     LE TRACE DE LA LEGENDE : COULEURS => VALEURS
C     --------------------------------------------
      IF( IAVFAC .GT. 0 ) THEN
C        TRACE SI LES FACES SONT TRACEES
         NBCOUL = NDCOUL - N1COUL
         NCPAS  = NBCOUL / 10
         TPAS   = (DEPMAX-DEPMIN) / 10
         T      =  DEPMIN
C
C        TRACE DE 11 VALEURS
         NCOUL = N1COUL
         NX    = LAPXFE - 150
         NY    = LHPXFE - 30
         DO 600 I=0,10
            CALL XVCOULEUR( NCOUL )
            CALL XVRECTANGLE( NX, NY, 30, 10 )
            WRITE( KNOM(1:10), '(G10.3)' ) T
            CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
            NCOUL = NCOUL + NCPAS
            T     = T  + TPAS
            NY    = NY - 15
 600     CONTINUE
      ENDIF
C
C     LE TRACE DU NOM DE L'OBJET ET DES DEPLACEMENTS EXTREMES
C     -------------------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'NORME DU DEPLACEMENT MAXIMAL = '
      ELSE
         KNOM = 'NORM of MAXIMUM DISPLACEMENT = '
      ENDIF
      WRITE( KNOM(31:44), '(G14.6)' ) DEPEMX
      CALL XVCOULEUR( NCNOIR )
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 50, 80 )
C
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'Z=NORME DE L''ERREUR DU DEPLACEMENT DE L''OBJET:'
     %       // KNOMOB
      ELSE
         KNOM = 'Z=NORM of DISPLACEMENT ERROR of the OBJECT:'
     %       // KNOMOB
      ENDIF
      CALL XVCOULEUR( NCROUG )
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
      IF( MODECO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = 'CAS      au TEMPS                : '
            WRITE( KNOM(5:8),   '(I4)'    ) NCAS
            WRITE( KNOM(19:33), '(G14.6)' ) TEMPS
         ELSE
            KNOM = 'CASE      at TIME                : '
            WRITE( KNOM(6:9),   '(I4)'    ) NCAS
            WRITE( KNOM(20:34), '(G14.6)' ) TEMPS
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = 'FREQUENCE PROPRE                : '
         ELSE
            KNOM = 'EIGENVALUE                      : '
         ENDIF
         WRITE( KNOM(18:32), '(G15.7)' ) SQRT(TEMPS)/(ATAN(1D0)*8D0)
      ENDIF
      I = NUDCNB( KNOM )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM(I+1:I+31) = ' ERREUR NORME DEPLACEMENT MIN= '
      ELSE
         KNOM(I+1:I+31) = ' NORM of DISPLACEMENT ERROR MIN= '
      ENDIF
      I = NUDCNB( KNOM )
      WRITE( KNOM(I+1:I+14), '(G14.6)' ) DEPMIN
      I = I + 14
      KNOM(I:I+5) = ' MAX= '
      WRITE( KNOM(I+7:I+20), '(G14.6)' ) DEPMAX
C     TRACE FORCE DU TITRE
      IAVTIT = 1
      CALL TRFINS( KNOM )
C
C     RETOUR POUR UNE NOUVELLE VISEE
      IF( LORBITE .NE. 0 ) GOTO 500
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
      NOEMIN = NOEMI0
      NOEMAX = NOEMA0
      DEPMIN = DEPMIN0
      DEPMAX = DEPMAX0
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
C     IAVELO : avec ou sans eloignement
      IAVELO = IAVEL0
C     TRACE SPECIAL EN Z
      NTRAXZ = 0
      NDIM   = 2
C
C     DESTRUCTION DES TMC INUTILES
 9999 IF( MNZTXY .GT. 0 ) CALL TNMCDS( 'REEL', NBSOM, MNZTXY )
      IF( MNXYZS .GT. 0 ) CALL TNMCDS( 'MOTS', WYZSOM+3*NBSOM,  MNXYZS )
      IF( MNNSEF .GT. 0 ) CALL TNMCDS( 'MOTS', WUSOEF+4*NBEFOB, MNNSEF )
      RETURN
      END
