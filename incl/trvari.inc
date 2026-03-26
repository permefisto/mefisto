C ===========================
C FICHIER ~/incl/trvari.inc
C ===========================

C     IDENTIFICATEURS NECESSAIRES POUR LE TRACE GRAPHIQUE AVEC LA LIBRAIRIE xvue

C     -------------------------------------------------------------------------
      DOUBLE PRECISION  TEMP2TRAC, TEMPUSER
      COMMON / TRVAPS / TEMP2TRAC, TEMPUSER, LASOPS, LVIDEO
C     -------------------------------------------------------------------------
C     TEMP2TRAC : TEMPS UTILISATEUR EN SECONDES ENTRE 2 TRACES
C     LASOPS    : DEMANDE OU NON DE SORTIE DU TRACE DU MAILLAGE EN POSTSCRIPT
C              0     : PAS DE SORTIE POSTSCRIPT
C            +-1     : TRACE PS SANS LES MENUS (- pour desactiver localement)
C            +-2     : TRACE PS AVEC LES MENUS (- pour desactiver localement)
C            >|2|    : ETAT MOMENTANE DE GESTION DES RECTANGLES
C     LVIDEO    : DEMANDE OU NON DE SORTIE DU TRACE SUR FICHIER VIDEO.gif
C             =0     : PAS DE SORTIE DEMANDEE de VIDEO.gif
C            <>0     : SORTIE DEMANDEE de VIDEO.gif

C     -------------------------------------------------------------------------
      COMMON / TRVAR0 / NCVISU
C     -------------------------------------------------------------------------
C     NCVISU : NUMERO DE LA visual->class SOUS X11
C              0 PseudoColor
C              1 GrayScale
C              2 DirectColor
C              3 StaticColor
C              4 StaticGray
C              5 TrueColor

C     -------------------------------------------------------------------------
      COMMON / TRVAR1 / NOPACL, N1CORE, NDCORE, N1COEL, NDCOEL,
     %                     N1COUL, NDCOUL,
     %                     NCNOIR, NCROUG, NCVERT, NCBLEU,
     %                     NCCYAN, NCJAUN, NCMAGE, NCBLAN,
     %                     NCGRIS, NCGRIM, NCGRIC, NCBEIG,
     %                     NCORAN, NCSAUM, NCROSE, NCTURQ
C     -------------------------------------------------------------------------
C     NOPACL : NUMERO DE LA PALETTE DE COULEURS ACTIVE 
C     N1CORE : NUMERO DE LA PREMIERE COULEUR RESERVEE DE L'ECRAN
C     NDCORE : NUMERO DE LA DERNIERE COULEUR RESERVEE DE L'ECRAN
C     N1COEL : NUMERO DE LA PREMIERE COULEUR ELEMENTAIRE (NOIRE)
C     NDCOEL : NUMERO DE LA DERNIERE COULEUR ELEMENTAIRE (BLANCHE)
C     N1COUL : NUMERO DE LA PREMIERE COULEUR NON RESERVEE DE L'ECRAN
C     NDCOUL : NUMERO DE LA DERNIERE COULEUR NON RESERVEE DE L'ECRAN
C     SI N1COUL+1 = NDCOUL L'ECRAN EST CONSIDERE NOIR ET BLANC
C     NCNOIR,...,NCMAGE : NUMERO DE CES COULEURS DANS LA PALETTE DES COULEURS
C     NCNOIR=0 LES AUTRES =1 SUR UN ECRAN NOIR ET BLANC

C     -------------------------------------------------------------------------
      COMMON / TRVAR2 / NOPOCA, NPHACA, NPLACA, RALALH  
C     -------------------------------------------------------------------------
C     NOPOCA : NUMERO DE LA POLICE DES CARACTERES
C     NPHACA : NOMBRE DE PIXELS DANS LA HAUTEUR DES CARACTERES
C     NPLACA : NOMBRE DE PIXELS DANS LA LARGEUR DES CARACTERES
C     RALALH : RAPPORT EN MM ECRAN ENTRE LARGEUR / HAUTEUR D'UN CARACTERE 

C     -------------------------------------------------------------------------
      COMMON / TRVAR3 / CXMMPX,CYMMPX,
     %                   LAMMEC,LHMMEC, LAPXEC,LHPXEC,
     %                   LAPXFE,LHPXFE, 
     %                   XOBMIN,XOBMAX, YOBMIN,YOBMAX,
     %                   AXOBPX,BXOBPX, AYOBPX,BYOBPX
C     -------------------------------------------------------------------------
C     CXMMPX       : COEFFICIENT DE TRANSFORMATION DES MM EN PIXELS EN X
C                    MM * CXMMPX => NB_PIXELS EN X
C     1.0 / CXMMPX : COEFFICIENT DE CONVERSION DES PIXELS EN MM EN X
C                    PX en X / CXMMPX => MM EN X          
C     CYMMPX       : COEFFICIENT DE TRANSFORMATION DES MM EN PIXELS EN Y
C                    MM * CYMMPX => NB_PIXELS EN Y
C     1.0 / CYMMPX : COEFFICIENT DE CONVERSION DES PIXELS EN MM EN Y
C                    PX en Y / CYMMPX => MM EN Y
      REAL     CXMMPX,CYMMPX

C     LAMMEC LHMMEC : LARGEUR ET HAUTEUR DE L'ECRAN EN MM
      INTEGER  LAMMEC,LHMMEC
C     LAPXEC LHPXEC : LARGEUR ET HAUTEUR DE L'ECRAN EN PIXELS
      INTEGER  LAPXEC,LHPXEC 
C     LAPXFE LHPXFE : LARGEUR ET HAUTEUR DE LA FENETRE EN PIXELS   
      INTEGER  LAPXFE,LHPXFE  
C
C     XOBMIN : ABSCISSE OBJET DE LA GAUCHE DE LA FENETRE
C     XOBMAX : ABSCISSE OBJET DE LA DROITE DE LA FENETRE
C     YOBMIN : ORDONNEE OBJET DU    BAS    DE LA FENETRE
C     YOBMAX : ORDONNEE OBJET DU    HAUT   DE LA FENETRE
      REAL     XOBMIN,XOBMAX, YOBMIN,YOBMAX

C     XOB -> X PIXELS     NUPXEX(x) = AXOBPX * x + BXOBPX  
C                         XOB2PX(nx) = ( nx - BXOBPX ) / AXOBPX
C     AXOBPX = LAPXFE / ( XOBMAX - XOBMIN )
C     BXOBPX = - 1.0  / ( XOBMAX - XOBMIN )

C     YOB -> Y PIXELS     NUPXEY(y)  = AYOBPX * y + BYOBPX 
C                         YOB2PX(ny) = ( ny - BYOBPX ) / AYOBPX
C     AYOBPX = - LHPXFE / ( YOBMAX - YOBMIN )
C     BYOBPX = 1.0 / ( YOBMAX - YOBMIN )
      REAL    AXOBPX,BXOBPX, AYOBPX,BYOBPX

C     -------------------------------------------------------------------------
      COMMON / TRVAR4 / NDIMLI, NOTYVI,
     %                  AXOLAR, AXOHAU, AXOARR, AXOAVA,
     %                  AXLMIN, AXLMAX, AXHMIN, AXHMAX, 
     %                  AXOPTV(6), AXOEIL(6), DIREVI(6),
     %                  AXOMAT(6,6), AXOOEIL(3), DPTVOEIL,
     %                  NBLAMP, AXOLAM(6,8)
C     -------------------------------------------------------------------------
C     NDIMLI : 1 VARIABLE  ENTIERE AUXILIAIRE 1=>1D 2=>2D 3=>3D        
C     NOTYVI :   NUMERO DU TYPE DE LA VISEE
C               0 SI NON INITIALISE 
C               1 SI INITIALISE EN 2D
C              11 SI AXONOMETRIE PAR PTV OEIL AXOARR ET AXOAVA
C     AXOLAR : LARGEUR/2 VUE SUR L'ECRAN  (EN VALEUR ABSOLUE )
C     AXLMIN, AXLMAX: EVENTUEL MIN ET MAX EN LARGEUR
C     AXOHAU : HAUTEUR/2 VUE SUR L'ECRAN
C     AXHMIN, AXHMAX: EVENTUEL MIN ET MAX EN HAUTEUR    
C     AXOARR : DISTANCE A PTV DU PLAN ARRIERE DE SECTION 
C     AXOAVA : DISTANCE A PTV DU PLAN AVANT   DE SECTION  
C     AXOPTV : XYZ DU POINT VISE
C     AXOEIL : XYZ DE L'OEIL
C     AXOOEIL: XYZ DE L'AXONOMETRIE DE L'OEIL
C     DPTVOEIL:DISTANCE PTV -> OEIL   
C     DIREVI : TABLEAU AUXILIAIRE = AXOEIL - AXOPTV
C              ECONOMISE DES CALCULS LORS DES TRACES AVEC CALCUL DE LA COULEUR
C              DES FACES
C     AXOMAT : MATRICE DE LA TRANSFORMATION (XYZAXO = AXOMAT * (xyz - PTV))

C     NBLAMP : NOMBRE DE LAMPES
C              0 => L'OEIL SERT D'UNIQUE LAMPE
C     AXOLAM : XYZ DE CHACUNE DES NBLAMP (8 AU PLUS)

C     -------------------------------------------------------------------------
      COMMON / TRVAR5 /  RLONGITUDE0, RLATITUDE0, 
     %                   LORBITE, NOBOUTON, NOETATBOUTON,
     %                   NOPX0, NOPY0, NOPX1, NOPY1, DIAMOB,
     %                   NORBITE
C     -------------------------------------------------------------------------
C     LORBITE:  =1 TEMOIN DE TRACE EN MODE ORBITE ZOOM TRANSLATION
C               =0 SINON
C     NOBOUTON: =1 a 3 NUMERO DU BOUTON ACTIF DE LA SOURIS
C               =0 SI AUCUN ACTIF
C     NOETATBOUTON: 1 SI LE BOUTON NON NUL ACTUEL EST ENFONCE NON RELACHE
C                   0 SI LE BOUTON NON NUL ACTUEL A ETE ENFONCE PUIS RELACHE
C     NOPX0, NOPY0: NUMERO EN X ET Y DU PIXELS OU LE BOUTON ACTUEL A ETE ENFONCE
C     RLONGITUDE0: LONGITUDE INITIALE EN RADIANS OU LE BOUTON ACTUEL A ETE ENFONCE
C     RLATITUDE0 : LATITUDE  INITIALE EN RADIANS OU LE BOUTON ACTUEL A ETE ENFONCE
C     DIAMOB     : DIAMETRE DE L'HEXAEDRE ENGLOBANT L'OBJET
C     NORBITE    : 0 AU DEPART
C                  1 APRES EXECUTION DE ZOOM[123]D0 ou ORBITE0
C                >=2 APRES EXECUTION D'AU MOINS UNE FOIS  ZOOM[123]D1 ou ORBITE[123]
C                  0 A L'ABANDON DE  ZOOM[123]D1 ou ORBITE[123]
