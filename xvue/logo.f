      SUBROUTINE LOGO( MEFISTO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER LE LOGO DE MEFISTO
C -----
C
C ENTREE :
C --------
C MEFISTO : LE TEXTE A AFFICHER ENCADR'E
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    JUILLET 1998
C AUTEUR : PERRONNET ALAIN Lab Jacques-Louis LIONS UPMC     JUILLET 2006
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/mecoit.inc"
      INTRINSIC      INT
      CHARACTER*(*)  MEFISTO
      CHARACTER*48   KNOM
      INTEGER        HAUTEU, HAUTEC
      REAL           TETRA(3,10)
C
      A = 0.8
      TETRA(1,1)=0
      TETRA(2,1)=-A/2
      TETRA(3,1)=0
C
      TETRA(1,5)=A* SQRT(3.0) / 2
      TETRA(2,5)=0
      TETRA(3,5)=0
C
      TETRA(1,7)=0
      TETRA(2,7)=A/2
      TETRA(3,7)=0
C
      TETRA(1,10)=A* SQRT(3.0) / 6
      TETRA(2,10)=0
      TETRA(3,10)=A * SQRT(3.0) / 2 * SIN( 70.5  / 180 * 6.28312 )
C
      DO 1 K=1,3
         TETRA(K,2) = (TETRA(K,1)+TETRA(K, 5)) * 0.5
         TETRA(K,6) = (TETRA(K,5)+TETRA(K, 7)) * 0.5
         TETRA(K,3) = (TETRA(K,1)+TETRA(K, 7)) * 0.5
         TETRA(K,4) = (TETRA(K,1)+TETRA(K,10)) * 0.5
         TETRA(K,8) = (TETRA(K,5)+TETRA(K,10)) * 0.5
         TETRA(K,9) = (TETRA(K,7)+TETRA(K,10)) * 0.5
 1    CONTINUE
C
C     LA MEMOIRE PIXELS EST EFFACEE
C     =============================
      CALL EFFACEMEMPX
CCCC
CCCC     DEMANDE DE TRACE POSTSCRIPT
CCCC     ===========================
CCC      CALL XVINITIERPS( 1 )
C
C     LA COULEUR DU FOND
C     ==================
      NCOFON = NCFOND()
      CALL XVFOND( NCOFON )
C
C     REMPLISSAGE AVEC LA COULEUR DE FOND DE TOUTE LA FENETRE
C     =======================================================
C     LE COIN SUPERIEUR GAUCHE DU RECTANGLE
      CALL XVCOULEUR( NCOFON )
      CALL XVRECTANGLE( 0, 0, LAPXFE, LHPXFE )
C
cccC     COPIE DE LA FENETRE AVEC SON FOND DANS LE PIXMAP MEMPX
ccc      CALL FENETREMEMPX
C
C     LE CHOIX DE LA PALETTE 11 => ARC EN CIEL
      CALL PALCDE( 11 )
C
C     CHARGEMENT DE LA FONTE DE PLUS GRANDE HAUTEUR EN PIXELS
C     =======================================================
C     CETTE FONTE EST CHARGEE POUR LE LOGO
      CALL CHOIXFONTE( 20 )
C
C     LE CADRE MIN MAX EN 3D
C     ======================
      A = A *1.3
      COOEXT(1,1) = -A/10
      COOEXT(2,1) = -A/1.5
      COOEXT(3,1) = -A/10
      COOEXT(1,2) =  A/3
      COOEXT(2,2) =  A/3
      COOEXT(3,2) =  A/3
      CALL VISEE0
C
C     LONGITUDE ET LATITUDE EN RADIANS DE LA POSITION DE L'OEIL
      R      = ATAN(1.0)  / 45.0
      AXOLON = -25.0
      AXOLAT =  20.0
      RADLON = AXOLON * R
      RADLAT = AXOLAT * R
      COSLON = COS( RADLON )
      SINLON = SIN( RADLON )
      COSLAT = COS( RADLAT )
      SINLAT = SIN( RADLAT )
C     RAYON DANS LE PLAN XY
      AXODIS    = 0
      RXY       = AXODIS * COSLAT
      AXOEIL(1) = AXOPTV(1) + RXY    * COSLON
      AXOEIL(2) = AXOPTV(2) + RXY    * SINLON
      AXOEIL(3) = AXOPTV(3) + AXODIS * SINLAT
C
C     LA MATRICE DE L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
C
cccC     TRACE DU MOBILE  TAPIS DE TRIANGLES COLORIES
cccC     ============================================
ccc      CALL MOBILE
C
C     TRACE DES TRIANGLES
C     ===================
C     AU DELA LES TRAITS SONT CONTINUS ET DE NEPARF EPAISSEURS
      CALL XVTYPETRAIT( 0 )
      CALL XVEPAISSEUR( 2 )
C     NOMBRE DES EPAISSEURS DES ARETES DES FACES
      NEPARF = 2
C     PAS DE TRACE DES TGS
      IAVTGF = 0
C     TRACE DES TGS EN BLANC
      NCOTGF = NCBLAN
C     TRACE CONTINU DES TGS
      NTRTGF = 0
C
C     TETRAEDRE 1234
      CALL LOGOTE( NCJAUN, NCJAUN, NCBLAN, NCJAUN, NCNOIR,
     %             1, 2, 3, 4, TETRA )
C
C     TETRAEDRE 2568
      CALL LOGOTE( NCROUG, NCROUG, NCSAUM, NCROUG, NCNOIR,
     %             2, 5, 6, 8, TETRA )
C
C     TETRAEDRE 3679
      CALL LOGOTE( NCBLEU, NCBLEU, NCCYAN, NCBLEU, NCNOIR,
     %             3, 6, 7, 9, TETRA )
C
C     TETRAEDRE 489 10
      CALL LOGOTE( NCVERT, NCVERT, NCTURQ, NCVERT, NCNOIR,
     %             4, 8, 9, 10, TETRA )
C
C     TRACE DU REMPLISSAGE DU RECTANGLE ENCADRANT LE NOM
C     ==================================================
      NBEPAI  = 5
      NBBORD  = 20
      NBEPL   = 5
C     LE NOMBRE DE PIXELS EN LARGEUR ET HAUTEUR DE 'MEFISTO'
      NBCAR   = NUDCNB( MEFISTO )
      CALL XVNBPIXELTEXTE( MEFISTO, NBCAR, LARGEC, HAUTEC )
      LARGEU = NBEPAI + NBBORD + LARGEC + NBBORD + NBEPAI
      HAUTEU = NBEPAI + NBBORD + HAUTEC + NBBORD + NBEPAI + NBEPL
C
C     LE COIN SUPERIEUR GAUCHE DU RECTANGLE
      NXR = LAPXFE - 2*LARGEU - NBBORD
      NYR = 150
      CALL XVCOULEUR( NCCYAN )
      CALL XVRECTANGLE( NXR, NYR, LARGEU, HAUTEU )
C
C     TRACE DE L'OMBRAGE DU RECTANGLE
C     ===============================
      NBEPAIO = 10
      NBEPA2O = NBEPAIO / 2
      CALL XVEPAISSEUR( NBEPAIO )
      CALL XVCOULEUR( NCMAGE )
C     LES COORDONNEES DES 2 SOMMETS EXTREMAUX DU RECTANGLE
      NPXMI = NXR + NBEPAIO
      NPYMI = NYR + HAUTEU + NBEPA2O
      NPXMA = NXR + LARGEU + NBEPA2O
      NPYMA = NYR + NBEPAIO
      CALL XVTRAIT( NPXMI, NPYMI, NPXMA+NBEPA2O, NPYMI )
      CALL XVTRAIT( NPXMA, NPYMI, NPXMA,         NPYMA )
C
C     TRACE DE LA LIGNE DU BORD DU RECTANGLE
C     ======================================
      CALL XVCOULEUR( NCBLEU )
      CALL XVEPAISSEUR( NBEPL )
      CALL XVBORDRECTANGLE( NXR, NYR, LARGEU, HAUTEU )
C
C     TRACE DU TEXTE
C     ==============
      NXT = NXR + NBBORD + NBEPAI
      NYT = NYR + NBBORD + NBEPAI + HAUTEC
      CALL XVTEXTE( MEFISTO, NBCAR, NXT, NYT )
C
C     TRACE DU TRAIT ROUGE SOUS LE TEXTE
C     ==================================
      CALL XVCOULEUR( NCROUG )
      NYL = NYT + NBEPL * 2
      CALL XVTRAIT( NXT, NYL, NXT+LARGEC, NYL )
C
C     TRACE DU TEXTE Laboratoire...
C     ==============
      CALL XVCOULEUR( NCNOIR )
      KNOM = 'Laboratoire Jacques-Louis LIONS'
      CALL XVNBPIXELTEXTE( KNOM, 31, LARGEC, HAUTEC )
      NXT = INT( ( LAPXFE - LARGEC ) * 2. / 3. )
      NYT = INT( LHPXFE - 3.5*HAUTEC )
      CALL XVTEXTE( KNOM, 31, NXT, NYT )
C
      KNOM = 'UNIVERSITE Pierre et Marie CURIE Paris FRANCE '
      NN = NUDCNB( KNOM )
      CALL XVNBPIXELTEXTE( KNOM, NN, LARGEC, HAUTEC )
      NYT = INT( NYT + 1.5 * HAUTEC )
      CALL XVTEXTE( KNOM, NN, NXT, NYT )
C
C     NOM DE LA VERSION
      KNOM = ' '
      CALL VRSION( KNOM )
      NN = NUDCNB( KNOM )
      CALL CHOIXFONTE( 20 )
      CALL XVNBPIXELTEXTE( KNOM, NN, LARGEC, HAUTEC )
      NYT = INT( NYT + 1.5 * HAUTEC )
      CALL XVCOULEUR( NCROUG )
      CALL XVTEXTE( KNOM, NN, NXT, NYT )
C
C     COPIE DU PIXMAP MEMPX DANS FENETRE
      CALL MEMPXFENETRE
CCCC
CCCC     PARCOURS DES DIFFERENTES PALETTES DE COULEURS
CCCC     =============================================
CCCC     LE CHOIX DE LA PALETTE 11 => ARC EN CIEL
CCC      DO 50 K=1,11
CCC         CALL PALCDE( K )
CCC 50   CONTINUE
C
C     REMISE A ZERO DU TRACE
      NOTYVI = 0
CCCC
CCCC     SAUVEGARDE DU TRACE DU LOGO DANS UN FICHIER POSTSCRIPT
CCCC     ======================================================
CCC      CALL XVSAUVERPS( 'logo', 4 )
CCCC
CCCC     ENVOI SUR L'IMPRIMANTE DU FICHIER POSTSCRIPT
CCCC     ============================================
CCC      CALL XVIMPRIMERPS( 'logo', 4 )
C
      RETURN
      END
