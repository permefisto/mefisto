      SUBROUTINE THEPU3( KNOMOB, NBJEUX, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES SOLUTIONS NON NULLES DU PROBLEME NON LINEAIRE
C -----   -DELTA u + a1 u + a2 u**2 + a3 u**3 = 0  a3<0
C          u sur Gamma = 0
C          ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2 EN 2D OU 3D
C          ATTENTION: JEU 1: CONDUCTIVITE = 1  pour -DELTA u
C                     JEU 1: u INITIALE (1 SUR UN SOUS-DOMAINE, 0 AILLEURS)
C                     JEU 2: a1
C                     JEU 3: a2
C                     JEU 4: a3
C                     JEU 1: CONDITIONS AUX LIMITES
C                     (=> CONTACT=0 IMPOSE    DANS LE JEU 1
C                     LES FLUX SONT CALCULES AVEC  LE JEU 1
C                     (=> CONDUCTIVITE DONNEE DANS LE JEU 1)
C                     JEU 2 3 4: RIEN
C
C          COMPUTE THE NOT ZERO SOLUTIONS of THE NON LINEAR PROBLEM
C          -DELTA u + a1 u + a2 u**2 + a3 u**3 = 0  a3<0
C          u on Gamma = 0
C          FOR LAGRANGE FINITE ELEMENTS of DEGREE 1 or 2 in 1D or 2D or 3D
C          ATTENTION: THE BOUNDARY CONDITIONS ARE DONE WITH THE INPUT DATA
C                     OF THE FIRST GAME OF THERMAL INPUT DATA
C                     THE FLUXES ARE COMPUTED WITH THE CONDUCTIVITY
C                     OF THE FIRST GAME OF THERMAL INPUT DATA
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET DE VECTEURS A CALCULER
C          NAME OF THE OBJECT TO COMPUTE THE VECTORS
C NBJEUX : NOMBRE DE JEUX DE DONNEES = DEGRE+1 DU POLYNOME P(Lambda)
C          >1 ICI POUR AVOIR UN POLYNOME DE DEGRE>=2
C          NUMBER OF INPUT DATA GAMES = DEGREE + 1 OF THE POLYNOMIAL P
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C          O WITHOUT ERROR, >0 ELSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN          DECEMBRE 2009
C                SAINT PIERRE DU PERRAY & LJLL UPMC PARIS   JANVIER 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION   PENALI, EPSITE
      PARAMETER         (PENALI=0D0)
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      PARAMETER         (ITERMX=250)
      PARAMETER         (EPSITE=1D-4)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      EXTERNAL          ETTAEL
      DOUBLE PRECISION  RELMIN, D2PI, DELTA
      DOUBLE PRECISION  DINFO,  DCPU, DPREP, DFACT, DSOLU
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERKGALLOC
      INTRINSIC         ALLOCATED
C
      INTEGER           MNBG(4)
      DOUBLE PRECISION  CCBA(4), LAMBDA0, LAMBDAK, LAMBDAK2
      DOUBLE PRECISION  SQRT, PROSCD, D0, D
      DOUBLE PRECISION  ENERGYk,  NORMUK,  MAXUK, MINUK, MINUK0
      DOUBLE PRECISION  ENERGYk0, NORMUK0, MAXUK0
      DOUBLE PRECISION  ENERGYNQ, ENERGYQ
      DOUBLE PRECISION  NORRES,   MAXRES, MINRES, MAXREL
      CHARACTER*4       NOMELE(2)
      CHARACTER*(*)     KNOMOB
      LOGICAL           COMPGC
      DATA              RELMIN/-1D28/
C     TABLEAU NON UTILISE (VITESSE D'UN FLUIDE NON ICI CALCULE)
      DOUBLE PRECISION  VITEGt(5,3)
C
      IERR   = 0
      NBCOOR = 0
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB, NBJEUX
      ELSE
         WRITE(IMPRIM,20000) KNOMOB, NBJEUX
      ENDIF
10000 FORMAT(/80('=')/'OBJET: ',A/
     %'-DELTA u + a1 u + a2 u**2 + a3 u**3 = 0 a3<0 avec ',
     % i2,' JEUX de DONNEES'/80('='))
20000 FORMAT(/80('=')/'OBJECT: ',A/
     %'-DELTA u + a1 u + a2 u**2 + a3 u**3 = 0 a3<0 with ',
     % i2,' INPUT DATA GAMES'/80('='))
C
C     LE DEGRE DU POLYNOME P(u) DOIT ETRE 3 => NBJEUX=4
      IF( NBJEUX .NE. 4 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = '4 JEUX DE DONNEES A FOURNIR ' // KNOMOB
         ELSE
            KERR(1) = '4 GAMES of INPUT DATA are REQUIRED ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     QUELQUES INITIALISATIONS
      TESTNL = 0
      DCPU   = DINFO( 'CPU' )
      DPREP  = 0D0
      DFACT  = 0D0
      DMODR  = 0D0
      DSOLU  = 0D0
C     TEMPS or TIME
      TEMPS  = 0.0
C     2 PI
      D2PI   = ATAN( 1D0 ) * 8D0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
C     NUMBER OF WORDS FOR A REAL DOUBLE PRECISION
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
C     PROTECTION OF ADDRESSES TO AVOID PROBLEMS DURING DELETION OF ARRAYS
      NTVECT = 0
      MNVECT = 0
      MNVEC01= 0
      MNNPEF = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNNODL = 0
      MNTHER = 0
      MNX    = 0
      MNLPLI = 0
      MNLPCO = 0
      NIVMAX = 5
      MNLPLC = 0
      MNLPCC = 0
      MNLPLU = 0
      MNAGC  = 0
      MNMG   = 0
      DO I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
      ENDDO
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NBDLMX = 0
      NTDL   = 0
      MOFLTO = 0
      MOFLPT = 0
      MONDLX = 0
      MNNDLX = 0
      MNVDLX = 0
      MNAUXGC= 0
      MNBDIR = 0
      MNVEC1 = 0
      MNVEC01= 0
      MNPOL  = 0
      MNPOID = 0
      MNPDEL = 0
      MNF2   = 0
      MNDPOL = 0
      IERKGALLOC = 1
C
C     AFFICHAGE ET VERIFICATION DU NOM_DE_L'OBJET
C     DISPLAY and VERIFICATION ABOUT THE OBJECT
C     ===========================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
C     RESEARCH OF THE OBJECT NAME IN THE DIRECTORY OF OBJECTS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
C     IF IT DOES NOT EXIST, RETURN
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1)='ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 2
         RETURN
      ENDIF
C
C     RECHERCHE DU TMS TABLEAU DEFINITION
C     SEARCH THE TMS DEFINITION
      CALL LXTSOU( NTLXOB,'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: DEFINITION INCONNUE de l''OBJET ' //KNOMOB
         ELSE
            KERR(1)='ERROR: UNKNOWN DEFINITION for the OBJECT '//KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     RECHERCHE TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" de l'OBJET
C     SEARCH THE TMS ASSOCIATED TO THIS OBJECT
      CALL MIMAOB( NBJEUX, NTLXOB, MXDOTH, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 4
         GOTO 9950
      ENDIF
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES PLSV de l'OBJET
C         NUMAOB          LES 4 NUMEROS MAXIMA DES PLSV de l'OBJET
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT les DONNEES THERMIQUES DE L'OBJET COMPLET
C     HERE NUMIOB CONTAINS THE 4 NUMBERS MINIMA of PLSV of the OBJECT
C          NUMAOB CONTAINS THE 4 NUMBERS MAXIMA of PLSV of the OBJECT
C     MNDOEL the 4 MCN ADRESSES of ARRAYS of ADDRESSES of THERMAL INPUT DATA
C
C     OUVERTURE DES TABLEAUX DES DONNEES THERMIQUES DES PLSV DE L'OBJET
C     OPENING of TMS of THERMAL INPUT DATA OF PLSV of THE OBJECT
C     =================================================================
      CALL THEDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             NBJEUX, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             IETEIN, IEVIIN, IEVIANT,IECOBO,
     %             IERR )
      IF( IERR .GT. 0 ) THEN
         IERR = 5
         GOTO 9950
      ENDIF
C
C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     INITIALIZATIONS of VARIABLES and DISPLAY
C     ==========================================
C     NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS=POINTS (3 ou 6)
C              NUMBER OF COORDINATES of  NODES =POINTS
      NBCOOR = MCN( MNXYZP + WBCOOP )
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
C     POINT NUMBER of the OBJECT MESH
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM DIMENSION 1 OU 2 OU 3 OU 6 DE L'ESPACE DES COORDONNEES
C          DIMENSION 1 or 2 or 3 or 6 of the SPACE of POINTS or NODES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
C     PARTICULARITE DES 6-CUBES (3Q1C NON A TRAITER EN THERMIQUE)
      IF( NBCOOR .EQ. 6 ) THEN
         NBTYEL = 1
         NDIM   = NBCOOR
C        ABANDON EN 6D
         GOTO 9950
      ENDIF
C     NOMBRE NOEUDS DU MAILLAGE DE L'OBJET=NOMBRE TOTAL DEGRES DE LIBERTE
C     TOTAL NUMBER of DEGREES OF FREEDOM = NODE NUMBER of THE OBJECT MESH
      NTDL = MCN( MNXYZN + WNBNOE )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM, NTDL, NBJEUX
      ELSE
         WRITE(IMPRIM,20210) NDIM, NTDL, NBJEUX
      ENDIF
10210 FORMAT(/' DIMENSION 1 ou 2 ou 3 de l''ESPACE',T42,'=',I9/
     %' NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR',   T42,'=',I9/
     %' NOMBRE DE JEUX DE DONNEES THERMIQUES',      T42,'=',I9)
20210 FORMAT(/' SPACE DIMENSION (1 or 2 or 3)',     T42,'=',I9/
     %' NUMBER of COMPONENTS of each EIGENVECTOR',  T42,'=',I9/
     %' NUMBER of THERMAL INPUT DATA GAMES',        T42,'=',I9)
C
C     RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA CONSTRUCTION DES
C     TABLEAUX ELEMENTAIRES
C     RETRIEVE THE ARRAYS OF POLYNOMIALS VALUES AT QUADRATURE FORMULA
C     POINTS WHICH ARE NECESSARY TO COMPUTE THE ELEMENT ARRAYS
C     ================================================================
      CALL TAPOBA( NBTYEL, MNNPEF, ETTAEL,
     %             MNTPOB, NBDLMX, MOAUX, NBTTEF, NOAXIS, NCODSM, IERR )
      NCODSK = 1
      IF( IERR .NE. 0 ) THEN
         IERR = 6
         GOTO 9950
      ENDIF
C
C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
C     ADDRESSING of AUXILIARY and ELEMENT TMC ARRAYS
C     ==================================================
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C
C     LES 2 MATRICES ELEMENTAIRES ET LE VECTEUR ELEMENTAIRE
C     The 2 ELEMENT MATRICES and VECTOR
      MOTAEL = NBDLMX * (NBDLMX+1) + NBDLMX
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
C
C     LE NUMERO DES DEGRES DE LIBERTE GLOBAUX DES DL D'UN EF
C     DOF GLOBAL NUMBER of the ELEMENT DEGREES of FREEDOM
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
C
C     LE TENSEUR DE CONDUCTIVITE ou COEFFICIENT TEMPERATURE(Pts INTEGRATION)
C     CONDUCTIVITY TENSOR at QUADRATURE FORMULA POINTS
      CALL TNMCDC( 'REEL2', 128, MNTHER )
C
C     LE TABLEAU DES NBCOOR COORDONNEES DES NBDLMX NOEUDS D'UN ELEMENT FINI
C     NBCOOR COORDINATES of  NBDLMX NODES of AN FINITE ELEMENT
      CALL TNMCDC( 'REEL', NBDLMX*NBCOOR, MNX )
      IF( MNX .LE. 0 ) GOTO 9920
C
C     DECLARATION DES NBJEUX MATRICES GLOBALES MORSES CONDENSEES
C     DECLARATION of the NBJEUX GLOBAL CONDENSED MATRICES
C     ==========================================================
C     CALCUL DU SQUELETTE DE LA MATRICE (LA MATRICE EST ICI SYMETRIQUE)
C     CONSTRUCTION of POINTERS on THE DIAGONAL COEFFICIENT of EVERY LINE
      CALL PRGCMC( MNTOPO, MCN(MNNPEF), MNXYZN, 1, NCODSK,
     %             MNLPLI, MNLPCO, IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 7
         GOTO 9905
      ENDIF
C
C     ADRESSAGE du TMC d'UNE AG SYMETRIQUE NON DIAGONALE SYMETRIQUE
C     ADDRESSING of TMC of ONE SYMMETRIC CONDENSED MATRIX
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE PROFIL KG
      NBRDKG = MCN( MNLPLI + NBDLIB )
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9905
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
      WRITE(IMPRIM,*)
C
C     ICI PAS DE MG.  HERE NO MASS MATRIX
      NCODSM = 0
      MNMG   = 0
      NBRDMG = 0
C
C     ICI PAS DE SECOND MEMBRE GLOBAL POUR THEMKB
C     HERE NO SECOND MEMBER VECTOR FOR THE SUBROUTINE THEMKB
      MNVBG = 0
C
C     CALCUL DE LA MATRICE DE CONDUCTIVITE SUPPOSEE
C     INDEPENDANTE DU TEMPS AVEC LE JEU 1 DES DONNEES
C     COMPUTATION OF the CONDUCTIVITY MATRIX
C     SUPPOSED INDEPENDENT OF TIME WITH THE FIRST DATA GAME
C     ------------------------------------------------------
      IEMG = 0
      IEKG = 1
      IEBG = 0
      CALL THEMKB( NBJEUX, IEMG,   IEKG,    IEBG,   PENALI,
     &             D2PI,   NDIM,   NTDL,    VITEGt,
     &             NBTYEL, MNNPEF, NDPGST,
     &             MNTPOB, MXPOBA, MNTAUX,
     &             MNXYZP, NUMIOB, NUMAOB,  MNDOEL,
     &             MNTHER, MNTAEL, MNX,     MNNODL,
     &             2,      MNLPLI, MNLPCO,
     &             0,      MCN(1), NBRDKG,  KG,     MNVBG,
     &             NCODSM, NCODSK, NBPTAF,  IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 9
         GOTO 9905
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'UNE MATRICE CONDENSEE AG a ete CALCULEE'
      ELSE
         WRITE(IMPRIM,*) 'ONE AG CONDENSED MATRIX has been COMPUTED'
      ENDIF
C
C     CONSTRUCTION DES TABLEAUX DU NUMERO ET VALEUR DES DL FIXES
C     A PARTIR DU PREMIER JEU DE DONNEES (COEFFICIENT 0 du POLYNOME)
C     FROM THE FIRST GAME OF THERMAL INPUT DATA CONSTRUCTION OF THE
C     LIST OF THE DEGREE OF FREEDOM IMPOSED TO ZERO
C     --------------------------------------------------------------
      CALL THDLFX( 0,      NTDL,   NDIM,
     &             NBTYEL, MNNPEF, NDPGST,
     &             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     &             NBDLFX, MONDLX, MNNDLX, MNVDLX, IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 10
         GOTO 9905
      ENDIF
C
      IF( NBDLFX .GT. 0 ) THEN
C
C        RENUMEROTATION DES DEGRES DE LIBERTE LIBRES NODLIB=MCN(MNDLIB)
C        RENUMBERING of DOF REALLY FREE
C        NODLIB(I) = NUMERO DU DEGRE DE LIBERTE S'IL EST LIBRE
C                   -INDICE DANS LA LISTE DES DL BLOQUES S'IL EST BLOQUE
C        NBDLIB    = NOMBRE TOTAL DE DEGRES DE LIBERTE LIBRES
C                    TOTAL NUMBER OF FREE DOF
         CALL REDLIB( NBDLFX, MCN(MNNDLX), NTDL,  NBDLIB, MNDLIB )
C
C        PRISE EN COMPTE DES TEMPERATURES IMPOSEES (DIRICHLET NON PENALISEE)
C        SUPPRESSION LIGNE ET COLONNE DES DEGRES DE LIBERTE IMPOSE A ZERO
C        TAKE IN ACCOUNT THE IMPOSED TEMPERATURES  (DIRICHLET NOT PENALISED)
C        SUPPRESSION OF LINE AND COLUMN OF DOF IMPOSED TO ZERO IN THE MATRICE
C        --------------------------------------------------------------------
         CALL BLDL0GC( NTDL, MCN(MNLPLI+NTDL), 1,
     %                 MCN(MNLPLI), MCN(MNLPCO), KG, MCN(MNDLIB))
C
C        REDUCTION DU STOCKAGE DE LA MATRICE
C        REDUCTION of THE MATRIX STORAGE
         NBRDAG = MCN( MNLPLI + NBDLIB )
         CALL TNMCRA( 'ENTIER', NBRDKG, NBRDAG, MNLPCO )
C
      ELSE
C
C        PAS DE DL FIXE
         NBDLIB = NTDL
C        ADRESSAGE MCN DU TABLEAU NODLIB(K)=K
         CALL TNMCDC( 'ENTIER', NTDL, MNDLIB )
         IF( MNDLIB .LE. 0 ) THEN
            IERR = 11
            GOTO 9920
         ENDIF
         DO K=1,NTDL
            MCN(MNDLIB-1+K) = K
         ENDDO
C
      ENDIF
C
C     NOMBRE DE COEFFICIENTS D'UNE MATRICE MORSE APRES SUPPRESSION DES DL FIXES
C     COEEFICIENT NUMBER of ONE CONDENSED MATRIX AFTER FIXED DOF DELETION
      NBRDAG = MCN( MNLPLI + NBDLIB )
C
C     AFFICHAGE DE LA MATRICE MORSE ou CONDENSEE APRES COMPRESSION
C     DISPLAY OF THE CONDENSED MATRIX AFTER LINES COLUMNS after COMPRESSION
      CALL AFMORSE( 25, NBDLIB, MCN(MNLPLI), MCN(MNLPCO), KG )
C
C     CALCUL DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT DU GC
C     -------------------------------------------------------------
      NIVEAU = 2
      ISTAB  = 0
C     CALCUL DES POINTEURS DE LA MATRICE DE PRECONDITIONNEMENT
C     (FACTORISATION INCOMPLETE DE A SUIVANT LE NIVEAU)
 80   MNLPLC = 0
      MNLPCC = 0
      MNLPLU = 0
      CALL CALPNT( NBDLIB, NIVEAU, NCODSK, MNLPLI, MNLPCO,
     %             MNLPLC, MNLPCC, MNLPLU, COMPGC,   IERR  )
C     VERIFICATION DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
      IF( IERR .NE. 0 ) THEN
         WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='PROBLEME lors de la FACTORISATION de NIVEAU '
     %           // KERR(MXLGER)(1:8)
            KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
         ELSE
            KERR(1)='PROBLEM of INCOMPLETE FACTORIZATION of LEVEL '
     %           // KERR(MXLGER)(1:8)
            KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
         ENDIF
         CALL LEREUR
         IERR = 12
         GOTO 9905
      ENDIF
C
C     DECLARATION INITIALISATION DE LA MATRICE MORSE DE PRECONDITIONNEMENT
C     --------------------------------------------------------------------
      LOLPCC = MCN( MNLPLC + NBDLIB )
      MNAGC  = 0
      CALL TNMCDC( 'REEL2', LOLPCC, MNAGC )
      IF( MNAGC .LE. 0 ) GOTO 9920
      CALL AZEROD( LOLPCC, MCN(MNAGC) )
C
C     CONSTRUCTION EFFECTIVE DE LA MATRICE DE PRECONDITIONNEMENT PAR
C     PAR FACTORISATION INCOMPLETE DE CHOLESKY AVEC NIVEAU
C     --------------------------------------------------------------
      NBLGRC(NRERR) = 1
      WRITE(KERR(2)(1:2),'(I2)') NIVEAU
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'FACTORISATION INCOMPLETE DE NIVEAU' //
     &              KERR(2)(1:2)
      ELSE
         KERR(1) = 'INCOMPLETE CHOLESKY FACTORIZATION of LEVEL' //
     &              KERR(2)(1:2)
      ENDIF
      CALL LERESU
      CALL INCHGC( NBDLIB,
     &             MCN(MNLPLI), MCN(MNLPCO), KG,
     &             MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC), IERR )
C
C     VERIFICATION DE LA STABILITE DE LA FACTORISATION INCOMPLETE
C     QUI DEFINIT LA MATRICE DE PRECONDITIONNEMENT DU GRADIENT CONJUGUE
      IF( IERR .EQ. 1 ) THEN
C        AU MOINS UN PIVOT<=0 LE PIVOT QUI PRECEDE EST PRIS A SA PLACE
         ISTAB = 1
         IERR  = 0
      ELSE IF( IERR .EQ. 2 ) THEN
C        LE NOMBRE DE PIVOTS INCORRECTS<=0 EST DEPASSE
C        LA FACTORISATION EST DECLAREE INSTABLE
         IF( COMPGC ) GOTO 9905
         IF( NIVEAU .LT. NIVMAX ) THEN
C           TENTATIVE D'AUGMENTER LE NIVEAU DE FACTORISATION INCOMPLETE
            NIVEAU = NIVEAU + 1
            CALL TNMCDS( 'ENTIER', NBDLIB+1, MNLPLC )
            CALL TNMCDS( 'ENTIER', LOLPCC,   MNLPCC )
            CALL TNMCDS( 'REEL2',  LOLPCC,   MNAGC  )
            IERR = 0
            GOTO 80
         ELSE
C           LA FACTORISATION EST VRAIMENT TRES INSTABLE
C           PLUS DE NIVMAX NIVEAUX DEPASSE => ABANDON DU GC
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     %              // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
            ELSE
               KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %              // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
            ENDIF
            CALL LEREUR
C           DESTRUCTION DES TABLEAUX DES MATRICES DU GC
            IERR = 13
            GOTO 9905
         ENDIF
      ENDIF
C
C     LE TEMPS CALCUL DE LA FACTORISATION INCOMPLETE
      DFACT = DINFO( 'DELTA CPU' )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10021) NBDLIB+LOLPCC*(1+MOREE2)
      ELSE
         WRITE(IMPRIM,20021) NBDLIB+LOLPCC*(1+MOREE2)
      ENDIF
10021 FORMAT(' PRECONDITIONNEMENT MATRICE MORSE =',I15,' MOTS MEMOIRE' )
20021 FORMAT(' PRECONDITIONED MORSE MATRIX =',I15,' MEMORY WORDS' )
C
C     DECLARATION DES VECTEURS POUR LA RESOLUTION DU SYSTEME PAR
C     GRADIENT CONJUGUE PRECONDITIONNE
C
C     LE PREMIER TABLEAU AUXILIAIRE DE GCPRCH
      LOAUXGC = NBDLIB * 3
      CALL TNMCDC( 'REEL2', LOAUXGC, MNAUXGC )
      IF( MNAUXGC .LE. 0 ) GOTO 9920
C     REPARTITION INTERNE EN SOUS-TABLEAUX
      MNAUXGC1 = MNAUXGC  + NBDLIB * MOREE2
      MNAUXGC2 = MNAUXGC1 + NBDLIB * MOREE2
C
C     STABILISATION DU GC PAR RE-ORTHOGONALISATION DES DIRECTIONS
      IF( ISTAB .EQ. 0 ) THEN
         NBDIR = 1
      ELSE
         NBDIR = 10
      ENDIF
C     LE SECOND TABLEAU AUXILIAIRE DE GCPRCH
      LODIR = (NBDLIB+1) * NBDIR * 2
      MNBDIR = 0
      CALL TNMCDC( 'REEL2', LODIR, MNBDIR )
      IF( MNBDIR .LE. 0 ) GOTO 9920
C     MISE A ZERO DE CE TABLEAU
      CALL AZEROD( LODIR, MCN(MNBDIR) )
C     REPARTITION INTERNE EN SOUS-TABLEAUX
      MNADIR = MNBDIR + NBDLIB * NBDIR * MOREE2
      MNDAD  = MNADIR + NBDLIB * NBDIR * MOREE2
      MNBET  = MNDAD  + NBDIR * MOREE2
C
      MOTSGC = MOREE2*(LOAUXGC+LODIR)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'TABLEAUX AUXILIAIRES DU GC =',
     %        MOTSGC,' MOTS MEMOIRE'
      ELSE
         WRITE(IMPRIM,*)'AUXILIARY ARRAYS of CG =',
     %        MOTSGC,' MEMORY WORDS'
      ENDIF
C
C     RESERVATION DES VECTEURS POUR L'APRES DECOMPRESSION DES DL
C     RESERVATION FOR THE AFTER DECOMPRESSION OF THE DOF
C     ----------------------------------------------------------
      MOVECA = 3 * ( NTDL + NBDLIB )
      CALL TNMCDC( 'REEL2', MOVECA, MNVEC01 )
      IF( MNVEC01 .LE. 0 ) GOTO 9920
C     ADRESSES DES VECTEURS AUXILIAIRES
C     ADDRESSES OF AUXILIARY VECTORS
      MNVEC0  = MNVEC01
      MNVEC1  = MNVEC01 + NTDL   * MOREE2
      MNBG(1) = MNVEC1  + NTDL   * MOREE2
      MNBG(2) = MNBG(1) + NBDLIB * MOREE2
      MNBG(3) = MNBG(2) + NBDLIB * MOREE2
      MNBG(4) = MNBG(3) + NBDLIB * MOREE2
      CALL AZEROD( 2*NTDL, MCN(MNVEC0) )
C
C     TEMPS CALCUL pour CONSTRUIRE la MATRICE GLOBALE MORSE et VECTEURS
C     CPU TIME TO CONSTRUCT THE GLOBAL CONDENSED MATRIX and VECTORS
      DPREP = DINFO( 'DELTA CPU' )
C
C     =========================================================
C     SOLUTION of -DELTA u + a1 u + a2 u**2 + a3 u**3 = 0  a3<0
C     =========================================================
C
C     RECUPERATION DU VECTEUR INITIAL POUR LE PB NON LINEAIRE
C     VEC1 SERA LE VECTEUR GLOBAL DES TEMPERATURES INITIALES
      CALL TEMINI( KNOMOB, NTLXOB, MOREE2, NTDL, TEMPS, IETEIN,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL,
     %             NTVET0, MNVET0, NBTEM0, MNVEC0, IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 14
         GOTO 9905
      ENDIF
C
C     SUPPRESSION DES DL FIXES A ZERO DANS VEC0
      IF( NBDLFX .GT. 0 ) THEN
         CALL COMPVg( NTDL, MCN(MNDLIB), MCN(MNVEC0) )
      ENDIF
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ITERATION   0: NORME**2 U0=',
     %         PROSCD( MCN(MNVEC0), MCN(MNVEC0), NBDLIB )
C
      LAMBDA0 = 1D0
      NORMUK0 = 0D0
      MINUK0  = 0D0
      MAXUK0  = 0D0
      ENERGYk0= 0D0
      ITERK = 0
C
C     =======================================================
C     ITERATION K
C     =======================================================
 5    ITERK = ITERK + 1
C
C     LA SOLUTION A L'ITERATION PRECEDENTE
      MNTHET = MNVEC0
C
C     LES COEFFICIENTS DE CALCUL DE LAMBDAK
      DO I=1,4
         CCBA(I) = 0D0
      ENDDO
C
C     MISE A ZERO DES (NBJEUX-1) SECONDS MEMBRES
      CALL AZEROD( (NBJEUX-1)*NBDLIB, MCN(MNBG(1)) )
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES
C     ---------------------------------------
C     LE NOMBRE DE POINTS D'INTEGRATION PAR ARETE EN 2D
      NBPTAF = 0
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO 500 NOTYEL = 1 , NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
         IF( NUTYEL .LE. 4 ) THEN
C           EF AXISYMETRIQUE
            NOAXIS = 1
         ELSE
            NOAXIS = 0
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
C           POINTS DIFFERENTS DES NOEUDS
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        DU TABLEAU DES POLY(POINTS INTEGRATION), ...
C        ----------------------------------------------
C        SELON LE TYPE DE L'ELEMENT FINI
         GOTO( 10,10,10,10, 1, 1,1, 1, 1, 1,
     %          1, 1,13, 1,10,10,1,10,13,10,
     %         10,10,10,10, 1, 1,1,10,10,10,
     %         10,10,10, 1, 1), NUTYEL
C
C        ERREUR
 1       NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR THEPU3: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR THEPU3: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 15
         GOTO 9905
C
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
 10      L = MNTPOB + (NOTYEL-1) * MXPOBA
C
         IF( NUTYEL .EQ. 30 ) THEN
C           6CUBE 6Q1C N'A POUR L'INSTANT PAS DE FACE...
            NBNSOM = 0
            NARET  = 0
            NFACE  = 0
            GOTO 12
         ENDIF
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION - 1
C        EN 2D SEGMENT DE REFERENCE, EN 3D SURFACE DE REFERENCE
C        ......................................................
         IF( NDIM .EQ. 1 ) GOTO 12
C
         IA     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIMA  = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLA = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPIA   = MCN( IA + 2 )
         NBPTAF = MAX( NBPTAF, NPIA )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOIA = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOLA = MNPOIA + MCN( IA + 3 ) * MOREE2 * NPIA
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
         MNDPOA = MNPOLA + MCN( IA + 4 ) * MOREE2 * NBPOLA * NPIA
C
         IF( NBNSOM .EQ. 5 .OR. NBNSOM .EQ. 6 ) THEN
C           EF AVEC 2 TYPES DE FACES : PAR EXEMPLE LA PYRAMIDE et LE PENTAEDRE
C                                      POLA => TRIANGLE
C                                      POLQ => QUADRANGLE
C
C           LES POLYNOMES DE L'EF DE DIMENSION - 1 POUR LE SECOND TYPE DE
C           FACE=QUADRANGLE POSITIONNE EN 4-EME POSITION (CF SP ETTAEL)
C           ...........................................................
            IA     = MCN( L + 3  )
C           DIMENSION DE L ESPACE
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = MCN( IA + 1 )
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = MCN( IA + 2 )
            NBPTAF = MAX( NBPTAF, NPIQ )
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = IA + 8
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOIQ + MCN( IA + 3 ) * MOREE2 * NPIQ
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNPOLQ + MCN( IA + 4) * MOREE2 * NBPOLQ * NPIQ
C
         ELSE
C
C           POLYNOMES DES FACES DE TYPE 2 = CEUX DE LA FACE DE TYPE 1
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = NBPOLA
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = NPIA
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = MNPOIA
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOLA
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNDPOA
C
         ENDIF
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION TOTALE
C        EN 1D SEGMENT UNITE, EN 2D SURFACE DE REFERENCE,
C        EN 3D ou 6D VOLUME DE REFERENCE
C        .....................................................
         L      = L + 1
 12      IA     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIM   = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLY = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPI    = MCN( IA + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOID = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOL  = MNPOID + MCN( IA + 3 ) * MOREE2 * NPI
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
         MNDPOL = MNPOL + MCN( IA + 4 ) * MOREE2 * NBPOLY * NPI
C
C        LES TABLEAUX AUXILIAIRES
         MNF1   = MNTAUX
         MNF2   = MNF1   + MOREE2 * NPI
ccc        MNPDEL = MNF1 + MOREE2 * NPI * NDIM + MOREE2 * NPI * (NDIM-2) modif l
         MNPDEL = MNF1   + MOREE2 * NPI * NDIM
         MNDP   = MNPDEL + MOREE2 * NPI
         MNDFM1 = MNDP   + MOREE2 * NPI * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NPI * NDIM * NDIM
C        MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY) )
C
C        ADRESSE DU SECOND MEMBRE ELEMENTAIRE
C        DERRIERE LA MATRICE SYMETRIQUE ELEMENTAIRE DE CONDUCTIVITE
         NBDL = NBPOLY
         GOTO 30
C
C        TRIANGLE 2P1D  ET  TETRAEDRE 3P1D
C        =================================
 13      NBPOLY = NDIM + 1
         NBDL   = NBPOLY
         NPI    = 1
         NPIQ   = 1
         NPIA   = 1
         NBPTAF = MAX( NBPTAF, NPIA )
         MNF1   = MNTAUX
         MNDP   = MNF1 + MOREE2 * NDIM
         MNDFM1 = MNDP + MOREE2 * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM
C
C        LES TABLEAUX ELEMENTAIRES: KE, ME ET BE
C        =======================================
 30      IF( IEMG .NE. 0 ) THEN
            MNCAEL = MNTAEL + NBDL * (NBDL+1) / 2  * MOREE2
         ELSE
            MNCAEL = MNTAEL
         ENDIF
         MNSE = MNCAEL + NBDL * (NBDL+1) / 2  * MOREE2
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 200 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           -----------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNNODL) )
C
CCCC           LES POINTS GEOMETRIQUES DE L'ELEMENT
CCCC           ------------------------------------
CCC            CALL EFPOGE( MNELE, NUELEM, NBPGEF, MCN(MNPOEF) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LES COORDONNEES DES NBPOE POINTS DE L'ELEMENT FINI NUELEM
C           ---------------------------------------------------------
            CALL EFXYZP( NDIM, MNXYZP, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   RMCN(MNX) )
C
C           ===============================================================
C           LE CALCUL DES TABLEAUX AUXILIAIRES ET DES TABLEAUX ELEMENTAIRES
C           ===============================================================
            GOTO( 41, 41, 41, 41,  1,  1, 1,  1,  1,  1,
     %             1,  1, 42,  1, 41, 41, 1, 41, 49, 51,
     %            51, 51, 51, 51,  1,  1, 1, 40, 41,  1,
     %            51, 51, 40,  1  ), NUTYEL
C
C           *********************************************
C           1D LAGRANGE F:SEGMENT UNITE -> SEGMENT EST P1
C           *********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 1D
C           ----------------------------------------
 40         IF( RMCN(MNX) .GT. RMCN(MNX+1) ) THEN
C              PERMUTATION DES 2 SOMMETS POUR AVOIR UN EF DE MESURE>0
               XS = RMCN(MNX)
               RMCN(MNX) = RMCN(MNX+1)
               RMCN(MNX+1 ) = XS
C
               NST               = MCN( MNNODL )
               MCN( MNNODL     ) = MCN( MNNODL + 1 )
               MCN( MNNODL + 1 ) = NST
C
               NST       = NOOBPS(1)
               NOOBPS(1) = NOOBPS(2)
               NOOBPS(2) = NST
            ENDIF
            CALL E11LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), DELTA )
C
C           BOUCLE SUR LES JEUX DE DONNEES 2 A 4 POUR LES COEF_TEMP Ajeu-1
            DO JEU = 2, NBJEUX
C
C              {Ve}=INTEGRALE t[P(X)] COEF_TEMP ([P(X)]{Uke})**(JEU-1) dX
C              ----------------------------------------------------------
               CALL TV1LAG( NBJEUX, JEU, MCN(MNNODL), MCN(MNDLIB),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBLA, NUMIOB(2), NUMAOB(2),MCN(MNDOEL(2)),
     %                      MCN(MNF1),   MCN(MNPDEL),
     %                      MCN(MNVEC0), JEU-1, MCN(MNSE) )
C
C              PRODUIT SCALAIRE t{Ue} {Ve} POUR CALCUL A, B, C de LAMBDAk
C              ----------------------------------------------------------
               CALL tVgVe( NBDL, MCN(MNNODL), MCN(MNDLIB), MCN(MNVEC0),
     %                     MCN(MNSE), D )
               CCBA(JEU) = CCBA(JEU) + D
C
C              ASSEMBLAGE DU SECOND MEMBRE JEU POUR CALCUL Uk+1
C              ------------------------------------------------
               CALL ASVeVg( NTDL, NBDLIB, NBDL, MCN(MNNODL),MCN(MNDLIB),
     %                      MCN(MNSE), MCN(MNBG(JEU-1)) )
C
            ENDDO
            GOTO 200
C
C           ********************************************
C           2D OU AXISYMETRIQUE LAGRANGE ISOPARAMETRIQUE
C           ********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 2D OU AXISYMETRIE
C           -------------------------------------------------------
 41         CALL E12LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1), MCN(MNF2),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           BOUCLE SUR LES JEUX DE DONNEES 2 A 4 POUR LES COEF_TEMP Ajeu-1
            DO JEU = 2, NBJEUX
C
C              {Ve}=INTEGRALE t[P(X)] COEF_TEMP ([P(X)]{Uke})**(JEU-1) dX
C              ----------------------------------------------------------
               CALL TV2LAG( D2PI, NOAXIS, NBJEUX, JEU,
     %                      MCN(MNNODL), MCN(MNDLIB),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBSF, NUMIOB(3), NUMAOB(3),MCN(MNDOEL(3)),
     %                      MCN(MNF1), MCN(MNF2), MCN(MNPDEL),
     %                      MCN(MNVEC0), JEU-1, MCN(MNSE) )
C
C              PRODUIT SCALAIRE t{Ue} {Ve} POUR CALCUL A, B, C de LAMBDAk
C              ----------------------------------------------------------
               CALL tVgVe( NBDL, MCN(MNNODL), MCN(MNDLIB), MCN(MNVEC0),
     %                     MCN(MNSE), D )
               CCBA(JEU) = CCBA(JEU) + D
C
C              ASSEMBLAGE DU SECOND MEMBRE JEU POUR CALCUL Uk+1
C              ------------------------------------------------
               CALL ASVeVg( NTDL, NBDLIB, NBDL, MCN(MNNODL),MCN(MNDLIB),
     %                      MCN(MNSE), MCN(MNBG(JEU-1)) )
C
            ENDDO
            GOTO 200
C
C           ***************************
C           3D LAGRANGE ISOPARAMETRIQUE
C           ***************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 3D
C           ----------------------------------------
 51         CALL E13LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           BOUCLE SUR LES JEUX DE DONNEES 2 A 4 POUR LES COEF_TEMP Ajeu-1
            DO JEU = 2, NBJEUX
C
C              {Ve}=INTEGRALE t[P(X)] COEF_TEMP ([P(X)]{Uke})**(JEU-1) dX
C              ----------------------------------------------------------
               CALL TV3LAG( NBJEUX, JEU, MCN(MNNODL), MCN(MNDLIB),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBVC, NUMIOB(4), NUMAOB(4),MCN(MNDOEL(4)),
     %                      MCN(MNF1), MCN(MNPDEL),
     %                      MCN(MNVEC0), JEU-1, MCN(MNSE) )
C
C              PRODUIT SCALAIRE t{Ue} {Ve} POUR CALCUL A, B, C de LAMBDAk
C              ----------------------------------------------------------
               CALL tVgVe( NBDL, MCN(MNNODL), MCN(MNDLIB), MCN(MNVEC0),
     %                     MCN(MNSE), D )
               CCBA(JEU) = CCBA(JEU) + D
C
C              ASSEMBLAGE DU SECOND MEMBRE JEU POUR CALCUL Uk+1
C              ------------------------------------------------
               CALL ASVeVg( NTDL, NBDLIB, NBDL, MCN(MNNODL),MCN(MNDLIB),
     %                      MCN(MNSE), MCN(MNBG(JEU-1)) )
C
            ENDDO
            GOTO 200
C
C           ************************************
C           2D TRIANGLE TRIA 2P1D LAGRANGE DROIT
C           ************************************
C           BOUCLE SUR LES JEUX DE DONNEES 2 A 4 POUR LES COEF_TEMP Ajeu-1
 42         DO JEU = 2, NBJEUX
C
C              {Ve}=INTEGRALE t[P(X)] COEF_TEMP ([P(X)]{Uke})**(JEU-1) dX
C              ----------------------------------------------------------
               CALL TV2P1D( NBJEUX, JEU, RMCN(MNX),
     %                      MCN(MNNODL), MCN(MNDLIB),
     %                      NOOBSF, NUMIOB(3), NUMAOB(3),MCN(MNDOEL(3)),
     %                      MCN(MNVEC0), JEU-1, MCN(MNSE) )
C
C              PRODUIT SCALAIRE t{Ue} {Ve} POUR CALCUL A, B, C de LAMBDAk
C              ----------------------------------------------------------
               CALL tVgVe( NBDL, MCN(MNNODL), MCN(MNDLIB), MCN(MNVEC0),
     %                     MCN(MNSE), D )
               CCBA(JEU) = CCBA(JEU) + D
C
C              ASSEMBLAGE DU SECOND MEMBRE JEU POUR CALCUL Uk+1
C              ------------------------------------------------
               CALL ASVeVg( NTDL, NBDLIB, NBDL, MCN(MNNODL),MCN(MNDLIB),
     %                      MCN(MNSE), MCN(MNBG(JEU-1)) )
C
            ENDDO
            GOTO 200
C
C           *************************************
C           3D TETRAEDRE TETR 3P1D LAGRANGE DROIT
C           *************************************
 49         CALL E13P1D( RMCN(MNX), MCN(MNF1),
     %                   DELTA, MCN(MNDFM1), MCN(MNDP) )
C
C           BOUCLE SUR LES JEUX DE DONNEES 2 A 4 POUR LES COEF_TEMP Ajeu-1
            DO JEU = 2, NBJEUX
C
C              {Ve}=INTEGRALE t[P(X)] COEF_TEMP ([P(X)]{Uke})**(JEU-1) dX
C              ----------------------------------------------------------
               CALL TV3P1D( NBJEUX, JEU, RMCN(MNX), DELTA,
     %                      MCN(MNNODL), MCN(MNDLIB),
     %                      NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDOEL(4)),
     %                      MCN(MNVEC0), JEU-1, MCN(MNSE) )
C
C              PRODUIT SCALAIRE t{Ue} {Ve} POUR CALCUL A, B, C de LAMBDAk
C              ----------------------------------------------------------
               CALL tVgVe( NBDL, MCN(MNNODL), MCN(MNDLIB), MCN(MNVEC0),
     %                     MCN(MNSE), D )
               CCBA(JEU) = CCBA(JEU) + D
C
C              ASSEMBLAGE DU SECOND MEMBRE JEU POUR CALCUL Uk+1
C              ------------------------------------------------
               CALL ASVeVg( NTDL, NBDLIB, NBDL, MCN(MNNODL),MCN(MNDLIB),
     %                      MCN(MNSE), MCN(MNBG(JEU-1)) )
C
            ENDDO
C
 200     CONTINUE
 500  CONTINUE
C
C     CALCUL DE INTEGRALE tGrad Uk Grad Uk dx = CCBA(1)
C     -------------------------------------------------
      CALL MAGCVE( 0, 1D0, NBDLIB, NCODSK,
     %             MCN(MNLPLI), MCN(MNLPCO), KG, MCN(MNVEC0),
     %             MCN(MNVEC1) )
      CCBA(1)= PROSCD( MCN(MNVEC0), MCN(MNVEC1), NBDLIB )
      IF( CCBA(4) .EQ. 0D0 ) THEN
         NBLGRC(NRERR) = 2
         WRITE( KERR(4)(1:5),'(I5)') K
         KERR(1) = 'THEPU3: ITERATION' // KERR(4)(1:5)
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'THEPU3: PROBLEME CCBA(4)=0D0'
         ELSE
            KERR(2) = 'THEPU3: PROBLEM CCBA(4)=0D0'
         ENDIF
         CALL LEREUR
         IERR = 17
         GOTO 9905
      ENDIF
C
C     CALCUL DE LAMBDAK = (-b+-SQRT(b**2-4ac))/(2a) >0
C     ------------------------------------------------
C     LAMBDA K+1 EST LA RACINE POSITIVE DE L'EQUATION DU 2-EME DEGRE
      D = SQRT( CCBA(3)**2 - 4D0 * CCBA(4) * (CCBA(1)+CCBA(2)) )
      LAMBDAK  = ( -CCBA(3) - D ) / ( 2D0 * CCBA(4) )
      LAMBDAK2 = ( -CCBA(3) + D ) / ( 2D0 * CCBA(4) )
C
C     VALEURS DE CCBA(1:4)=> i=2,3,4 Integrale( ai-1 Uk**i )dx POUR Uk
C     ET LES 2 RACINES LAMBDAK et LAMBDAK2
      WRITE(IMPRIM,10500) ITERK, CCBA, LAMBDAK, LAMBDAK2
10500 FORMAT(/'ITERATION',I4,': Int(gradUk**2)dx=',1Pg14.6,
     %' Int(a1Uk2)=',1Pg14.6,' Int(a2Uk3)=',1Pg14.6,
     %' Int(a3Uk4)=',1Pg14.6,
     %'  Lk1=',1Pg14.6,'  Lk2=',1Pg14.6)
C
      IF( LAMBDAK .LE. 0 ) THEN
         LAMBDAK = LAMBDAK2
      ENDIF
C
      IF( LAMBDAK .LE. 0D0 ) THEN
         NBLGRC(NRERR) = 2
         WRITE( KERR(4)(1:5),'(I5)') K
         KERR(1) = 'THEPU3: ITERATION' // KERR(4)(1:5)
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'THEPU3: PROBLEME LAMBDAK<=0'
         ELSE
            KERR(2) = 'THEPU3: PROBLEM LAMBDAK<=0'
         ENDIF
         CALL LEREUR
         IERR = 18
         GOTO 9905
      ENDIF
C
C     CALCUL DES ENERGIES de Lambdak Uk
C     ---------------------------------
C     ENERGIE CINETIQUE de la SOLUTION = INTEGRALE (Lambdak grad Uk)**2/2 dx
      ENERGYQ = CCBA(1) * LAMBDAK**2 /2
C
C     ENERGIE "NON QUADRATIQUE" ou POTENTIELLE? DE LA SOLUTION (LAMBDAk*Uk)
C     INTEGRALE a1/2 (LK UK)**2 + a2/3 (LK UK)**3 + a3/4 (LK UK)**4
      ENERGYNQ = CCBA(2)/2 * LAMBDAK**2
     %         + CCBA(3)/3 * LAMBDAK**3
     %         + CCBA(4)/4 * LAMBDAK**4
C
C     ENERGIE = ENERGIE CINETIQUE + ENERGIE POTENTIELLE NON QUADRATIQUE
      ENERGYk = ENERGYQ + ENERGYNQ
C
      WRITE(IMPRIM,10502) ITERK, LAMBDAK, ENERGYQ, ENERGYNQ, ENERGYk
10502 FORMAT('ITERATION',I4,': Lk=',1Pg14.6,
     %'  EQ=Int(GradLkUk)**2/2dx=',1Pg14.6,
     %'  EN=Int(a1LkUk**2/2+a2LkUk**3/3+a3LkUk**4/4)dx=',1Pg14.6,
     %'  Energyk=EQ+EN=',1Pg14.6)
C
C     CALCUL DES NORMES de Lambdak Uk Lk VERITABLE SOLUTION DU PB
C     -----------------------------------------------------------
      MN = ( MNVEC0 - 1 ) / 2
      NORMUK = 0D0
      MINUK  = DMCN(MN+1) * LAMBDAK
      MAXUK  = DMCN(MN+1) * LAMBDAK
      DO I=1,NBDLIB
         D = DMCN(MN+I) * LAMBDAK
         MINUK  = MIN( MINUK, D )
         MAXUK  = MAX( MAXUK, D )
         NORMUK = NORMUK + D ** 2
      ENDDO
      NORMUK = SQRT( NORMUK / NBDLIB )
      WRITE(IMPRIM,11503) ITERK, MINUK, MAXUK, NORMUK
11503 FORMAT('ITERATION',I4,
     %': MIN LkUk(N)=',1PG14.6,
     %'  MAX LkUk(N)=',1PG14.6,
     %'  ||LkUk||=',1PG14.6)
C
C
C     CONSTRUCTION DU SECOND MEMBRE POUR LE CALCUL DE Uk+1
C     SM = - ( a1 Uk + Lambdak * a2 * Uk**2 + Lambdak**2 * a3 * Uk**3 )
C     COMBINAISON LINEAIRE DES 3 VECTEURS GLOBAUX ASSEMBLES
C     DANS LA BOUCLE SUR LES EF DU MAILLAGE
C     -----------------------------------------------------------------
      CALL CL3VED( NBDLIB, -1D0, MCN(MNBG(1)),-LAMBDAK, MCN(MNBG(2)),
     %            -(LAMBDAK**2), MCN(MNBG(3)),  MCN(MNBG(4)) )
C
C     CALCUL DE Uk+1 PAR GC PRECONDITIONNE SUR LE SYSTEME LINEAIRE
C     -DELTA Uk+1 = - ( a1 Uk + Lambdak * a2 Uk**2 + Lambdak**2 * a3 Uk**3 )
C     ----------------------------------------------------------------------
      CALL GCPRCH( NBDLIB, 1, NBDIR,
     %             MCN(MNLPLI), MCN(MNLPCO), KG, MCN(MNBG(4)),
     %             MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC), MCN(MNVEC0),
     %             MCN(MNAUXGC), MCN(MNAUXGC1), MCN(MNAUXGC2),
     %             MCN(MNBDIR),  MCN(MNADIR), MCN(MNDAD), MCN(MNBET),
     %             MCN(MNVEC1), IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 19
         GOTO 9905
      ENDIF
C
C     ARRET DES ITERATIONS K DE POINT FIXE ?
C     --------------------------------------
      WRITE(IMPRIM,11504) ITERK, ABS(NORMUK-NORMUK0)/NORMUK,
     %                    ABS(MAXUK-MAXUK0)/MAXUK,
     %                    ABS(ENERGYk- ENERGYk0)/ABS(ENERGYk)
11504 FORMAT('ITERATION',I4,
     %': ||LkUk||-||Lk-1Uk-1||/||LkUk||=',1PG14.6,
     %'  |MAX LkUk-MAX Lk-1Uk-1|/MAX LkUk=',1PG14.6,
     %'  |ENERGYk-ENERGYk0|/|ENERGYk|=',1PG14.6)
C
      MAXUK = MAX( ABS(MAXUK), ABS(MINUK) )
      IF( ABS( NORMUK - NORMUK0  ) .GT. EPSITE * NORMUK .OR.
     %    ABS( MAXUK  - MAXUK0   ) .GT. EPSITE * MAXUK  .OR.
     %    ABS( ENERGYk- ENERGYk0 ) .GT. EPSITE * ABS(ENERGYk) ) THEN
C
C        CONTROLE DE LA DIVERGENCE SUR L'ENERGIE CINETIQUE DE Uk
         IF( CCBA(1) .GT. 1D50 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:6),'(I6)') ITERK
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: DIVERGENCE de Int(gradUk**2)dx apres'
     %                 // KERR(5)(1:6) // ' ITERATIONS'
               KERR(2) ='SAUVEGARDE de la SOLUTION Uk ACTUELLE DIVERGEE'
            ELSE
               KERR(1) = 'ERROR: DIVERGENCE of Int(gradUk**2)dx after'
     %                 // KERR(5)(1:6) // ' ITERATIONS'
               KERR(2) = 'The ACTUAL DIVERGED SOLUTION Uk is SAVED'
            ENDIF
            CALL LEREUR
            IERR = 20
            GOTO 800
         ENDIF
C
C        CONTROLE DE LA CONVERGENCE LENTE
         IF( ITERK .GE. ITERMX ) THEN
            NBLGRC(NRERR) = 5
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: CONVERGENCE NON ATTEINTE APRES'
               WRITE(KERR(5)(1:6),'(I6)') ITERK
               KERR(2) = KERR(5)(1:6) // ' ITERATIONS'
               KERR(3) = 'COEFFICIENTS TROP GRANDS A REDUIRE?'
               KERR(4) = 'PROBLEME TROP FORTEMENT NON LINEAIRE?'
               KERR(5) = 'SAUVEGARDE DE LA SOLUTION Uk NON CONVERGEE'
            ELSE
               KERR(1) = 'ERROR: NO CONVERGENCE after'
               WRITE(KERR(5)(1:6),'(I6)') ITERK
               KERR(2) = KERR(5)(1:6) // ' ITERATIONS'
               KERR(3) = 'COEFFICIENTS TOO GREAT HAVE TO BE REDUCED?'
               KERR(4) = 'or TOO GREAT NO LINEARITY?'
               KERR(5) = 'ACTUAL NOT CONVERGED SOLUTION Uk is SAVED'
            ENDIF
            CALL LEREUR
            IERR = 21
C           ATTENTION: LA VERITABLE SOLUTION DU PB EST  LAMBDAk * Uk !
            MN = ( MNVEC0 - 1 ) / 2
            DO I=1,NBDLIB
               DMCN(MN+I) = DMCN(MN+I) * LAMBDAK
            ENDDO
            MN      = MNVEC1
            MNVEC1  = MNVEC0
            MNVEC0  = MN
            GOTO 800
         ENDIF
C
C        UNE ITERATION DE PLUS A FAIRE
C        ATTENTION: LA VERITABLE SOLUTION DU PB EST  LAMBDAk * Uk !
C        CETTE MULTIPLICATION STABILISE LES VALEURS DE LAMBDAk AUTOUR DE 0.5
         MN = ( MNVEC1 - 1 ) / 2
         DO I=1,NBDLIB
cccC           ESSAI D'OBTENIR UNE SOLUTION AVEC DES VALEURS POSITIVES
ccc            IF( DMCN(MN+I) .LT. 0D0 ) DMCN(MN+I) = 0D0
            DMCN(MN+I) = DMCN(MN+I) * LAMBDAK
         ENDDO
         MN      = MNVEC1
         MNVEC1  = MNVEC0
         MNVEC0  = MN
         LAMBDA0 = LAMBDAK
         NORMUK0 = NORMUK
         MAXUK0  = MAXUK
         ENERGYk0= ENERGYk
         GOTO 5
C
      ELSE
C
C        CONVERGENCE DES ITERATIONS DE POINT FIXE
         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10700) ITERK
         ELSE
            WRITE(IMPRIM,20700) ITERK
         ENDIF
10700    FORMAT('ITERATION',I4,': CONVERGENCE des ITERATIONS')
20700    FORMAT('ITERATION',I4,': CONVERGENCE of ITERATIONS')
C
C        ATTENTION: LA VERITABLE SOLUTION DU PB EST  LAMBDAk * Uk !
         MN = ( MNVEC0 - 1 ) / 2
         DO I=1,NBDLIB
            DMCN(MN+I) = DMCN(MN+I) * LAMBDAK
         ENDDO
         MN      = MNVEC1
         MNVEC1  = MNVEC0
         MNVEC0  = MN
C        L'ADRESSE DU VECTEUR SOLUTION Uk est MNVEC1
         GOTO 800
      ENDIF
C
C     VERIFICATION QUE Lambdak Uk EST SOLUTION DE -Delta u + a1u+a2u2+a3u3=0
C     ----------------------------------------------------------------------
C     CALCUL DU VECTEUR a1 Uk + a2 Uk**2 + a3 Uk**3
C     COMBINAISON LINEAIRE DES 3 VECTEURS GLOBAUX ASSEMBLES
C     DANS LA BOUCLE SUR LES EF DU MAILLAGE CONSTRUITS AVEC Uk et non LAMBDAk Uk
 800  CALL CL3VED( NBDLIB, LAMBDAK,    MCN(MNBG(1)),
     %                     LAMBDAK**2, MCN(MNBG(2)),
     %                     LAMBDAK**3, MCN(MNBG(3)),
     %             MCN(MNBG(4)) )
C
C     CALCUL DU VECTEUR -DELTA UK + a1 Uk + a2 Uk**2 + a3 Uk**3
      CALL MAGCVE( 1, 1D0, NBDLIB, NCODSK,
     %             MCN(MNLPLI), MCN(MNLPCO), KG, MCN(MNVEC1),
     %             MCN(MNBG(4)) )
C
C     NORMES DU VECTEUR (-DELTA Uk + a1 Uk + a2 Uk**2 + a3 Uk**3)(N)
      MNV = ( MNVEC1  - 1 ) / 2
      MNR = ( MNBG(4) - 1 ) / 2
      NORRES = 0D0
      MINRES = DMCN(MNR+1)
      MAXRES = DMCN(MNR+1)
      MAXREL = 0D0
      DO I=1,NBDLIB
C        RESIDU AU NOEUD I
         D = DMCN(MNR+I)
         MINRES = MIN( MINRES, D )
         MAXRES = MAX( MAXRES, D )
         NORRES = NORRES + D ** 2
C        SOLUTION AU NOEUD I
         D0 = ABS( DMCN(MNV+I) )
C        RESIDU RELATIF
         IF( D0 .GT. 1D-6*MAXUK ) THEN
            MAXREL = MAX( MAXREL, ABS(D) / D0 )
         ENDIF
      ENDDO
      NORRES = SQRT( NORRES / NBDLIB )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,11800) ITERK, MINRES, MAXRES, MAXREL, NORRES
      ELSE
         WRITE(IMPRIM,21800) ITERK, MINRES, MAXRES, MAXREL, NORRES
      ENDIF
      GOTO 9905
11800 FORMAT('ITERATION',I4,': VERIFICATION RESIDU R(N)=(-LAPLACIEN LkUk
     % + a1 LkUk + a2 (LkUk)**2 + a3 (LkUk)**3)(Noeud)=0 ? :'/
     %'MIN R(N)=',1PG14.6,
     %'  MAX R(N)=',1PG14.6,
     %'  MAX |R(N)/LkUk(N)|=',1PG14.6,
     %'  ||R||=',1PG14.6)
21800 FORMAT('ITERATION',I4,': VERIFICATION RESIDUE R(N)=(-LAPLACIAN LkU
     %k + a1 LkUk + a2 (LkUk)**2 + a3 (LkUk)**3)(Node)=0 ? :'/
     %'MIN R(N)=',1PG14.6,
     %'  MAX R(N)=',1PG14.6,
     %'  MAX |R(N)/LkUk(N)|=',1PG14.6,
     %'  ||R||=',1PG14.6)
C
cccC     ERREUR: EF DEGENERE RENCONTRE
ccc 9900 NBLGRC(NRERR) = 1
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         KERR(1) = 'ERREUR THEPU3: un EF '// NOMELE(1)
ccc     %        // NOMELE(2) // ' est DEGENERE'
ccc      ELSE
ccc         KERR(1) = 'ERROR THEPU3: 1 FE '// NOMELE(1)
ccc     %        // NOMELE(2) // ' is DEGENERATED'
ccc      ENDIF
ccc      CALL LEREUR
ccc      IERR = 2
C
C     DESTRUCTION DE LA MATRICE POUR REDONNER DE LA PLACE EN MC
C     DELETION OF the MATRIX TO RECUPERATE THE MEMORY FOR MCN ARRAYS
 9905 IF( IERKGALLOC .EQ. 0 ) THEN
         DEALLOCATE( KG )
         IERKGALLOC = 1
      ENDIF
      IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBRDAG,  MNLPCO )
      IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1,  MNLPLI )
      IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLIB+1,MNLPLC )
      IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC,  MNLPCC )
      IF( MNAGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC,  MNAGC  )
      IF( MNAUXGC.GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUXGC, MNAUXGC)
      IF( MNBDIR .GT. 0 ) CALL TNMCDS( 'REEL2',  LODIR,   MNBDIR )
      IF( IERR   .GT. 0 .AND. IERR .LT. 20 ) GOTO 9950
C
C     RETOUR A LA NUMEROTATION DE 1 A NTDL DES COMPOSANTES
C     DU VECTEUR SOLUTION GRACE AU TABLEAU NUDLIB(1:NTDL)
C     RETURN TO THE NUMBERING OF INITIAL DOF FROM 1 TO NTDL
      IF( NBDLFX .GT. 0 ) THEN
         CALL RENTDL( NBDLIB, NTDL, MCN(MNDLIB), 1,
     %                MCN(MNVEC1),  MCN(MNVEC1) )
         CALL RENTDL( NBDLIB, NTDL, MCN(MNDLIB), 1,
     %                MCN(MNBG(4)), MCN(MNBG(4)) )
      ENDIF
C
C     CONSTRUCTION DU TMS 'VECTEUR"TEMPERATURE' DE LA SOLUTION CALCULEE
      CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
      IF( NTVECT .GT. 0 ) THEN
C        LE VECTEUR EST DETRUIT POUR ETRE REDECLARE
         CALL LXTSDS( NTLXOB, 'VECTEUR"TEMPERATURE' )
      ENDIF
      L = WECTEU + 2 * NTDL * MOREE2
      CALL LXTNDC( NTLXOB, 'VECTEUR"TEMPERATURE', 'MOTS', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
C     NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR
      MCN( MNVECT + WBCOVE ) = NTDL
C     NOMBRE DE VECTEURS 2: Uk et le Residu V(N)=
C     (-DELTA Uk + a1 Uk + a2 Uk**2 + a3 Uk**3)(N)
      NBVECT = 2
      MCN( MNVECT + WBVECT ) = NBVECT
C     NOMBRE DE VALEURS REELLES SIMPLE PRECISION
      MCN( MNVECT + WBCPIN ) = 0
C     COPIE DES VECTEURS
      MNVECP = MNVECT + WECTEU
      CALL TRTATA( MCN(MNVEC1), MCN(MNVECP), NTDL * MOREE2 )
      MN = MNVECP + MOREE2 * NTDL
      CALL TRTATA( MCN(MNBG(4)), MCN(MN), NTDL * MOREE2 )
C
C     COUT CALCUL DES VECTEURS
C     VECTOR COMPUTATION COST
      DSOLU = DINFO( 'DELTA CPU' )
C
C **************************************************************************
C --------------------------------------------------------------------------
C **************************************************************************
C
C     ATTENTION: SEULEMENT AVEC LE JEU 1 DES DONNEES
C     CALCUL DES FLUX EN CHAQUE POINT D INTEGRATION DES FACES DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     LA CONDUCTIVITE DES MATERIAUX EST SUPPOSEE DANS LE JEU 1 DES DONNEES
C     THE FLUXES ARE COMPUTED WITH THE CONDUCTIVITY OF THE FIRST GAME
C     OF THE NBJEUX THERMAL INPUT DATA
C     ====================================================================
      IF( NBCOOR .LE. 3 .AND. IERR .EQ. 0 ) THEN
         CALL THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     %                NDIM,   MOREE2, NBVECT, NTDL,
     %                NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %                MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                MOTAEL, MNTAEL, MNX,    MNVECT )
      ENDIF
      GOTO 9950
C
C     DIAGNOSTIC POUR PAS ASSEZ DE MEMOIRE MCN
C     ========================================
 9920 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'PAS ASSEZ DE MEMOIRE DANS MCN'
      ELSE
         KERR(1) = 'NOT ENOUGH MEMORY IN MCN'
      ENDIF
      CALL LEREUR
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     DESTRUCTION OF USELESS TEMPORARY ARRAYS
C     =======================================
 9950 IF( IERKGALLOC .EQ. 0 ) THEN
         DEALLOCATE( KG )
         IERKGALLOC = 1
      ENDIF
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO 11000 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
11000 CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA,MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2' , MOAUX ,  MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2' , MOTAEL,  MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX,  MNNODL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2' , 128   ,  MNTHER )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDLMX*NBCOOR, MNX )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONDLX,  MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MONDLX,  MNVDLX )
      IF( MNVEC01.GT. 0 ) CALL TNMCDS( 'REEL2',  MOVECA,  MNVEC01)
      IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBRDAG,  MNLPCO )
      IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1,  MNLPLI )
      IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLIB+1,MNLPLC )
      IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC,  MNLPCC )
      IF( MNAGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC,  MNAGC  )
      IF( MNAUXGC.GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUXGC, MNAUXGC)
      IF( MNBDIR .GT. 0 ) CALL TNMCDS( 'REEL2',  LODIR,   MNBDIR )
      IF( IERR .NE. 0 ) RETURN
C
C     COUT CALCUL DES FLUX
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,12001) DPREP,DFACT,DSOLU,DCPU,DPREP+DFACT+DSOLU+DCPU
      ELSE
       WRITE(IMPRIM,22001) DPREP,DFACT,DSOLU,DCPU,DPREP+DFACT+DSOLU+DCPU
      ENDIF
12001 FORMAT(/
     %' TEMPS CALCUL DE LA PREPARATION   =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DE LA FACTORISATION =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DE LA SOLUTION      =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES FLUX            =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL TOTAL               =',F12.2,' SECONDES CPU')
22001 FORMAT(/
     %' PREPARATION   TIME =',F12.2,' CPU SECONDS'/
     %' FACTORIZATION TIME =',F12.2,' CPU SECONDS'/
     %' SOLUTION      TIME =',F12.2,' CPU SECONDS'/
     %' NORMAL FLUX   TIME =',F12.2,' CPU SECONDS'/
     %' TOTAL         TIME =',F12.2,' CPU SECONDS')
C
      RETURN
      END
