      SUBROUTINE THEINL( KNOMOB, MOREE2, D2PI,
     %                   NDIM,   NDPGST, MNXYZP, MNXYZN,
     %                   NBTYEL, MNNPEF,
     %                   MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                   PENALI, RELMIN,
     %                   MNTHER, MNTAEL, MNX,
     %                   NORESO, MNLPLI, MNLPCO, NIVEAU, NBRDKG,
     %                   NCODSM, MNUG0,
     %                   BETA,   DT,     DTSTOC, TPSINI, TPSFIN,
     %                   NBTEMP, MXTEMP, NTDL,   NTVECT, MNVECT,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES TEMPERATURES ET FLUX DANS UN DOMAINE
C ----- 2D OU 3D OU AXISYMETRIQUE EN THERMIQUE INSTATIONNAIRE LINEAIRE
C       SELON LA THETA METHODE AVEC PAS DE TEMPS CONSTANT
C       CAPACITE, CONDUCTIVITE, ECHANGE, SOURCE, CONTACT
C       PEUVENT DEPENDRE DU TEMPS ET DE LA TEMPERATURE
C       DES ITERATIONS DE POINT FIXE SONT NECESSAIRES
C       LA TEMPERATURE IMPOSEE (CONDITION DE DIRICHLET) EST TRAITEE
C       PAR PENALISATION AVEC ECHANGE=1/EPSILON et SOURCE=TEMPERATURE/EPSILON
C       ET LE CHOIX FAIT EST 1/EPSILON=PENALI=1D20
C
C******************************************************************************
C A L'ETAPE n+1, LE PROBLEME CONSISTE pour
C n l'iteration en temps et m l'iteration de point fixe avec Un+1,0 = Un
C [Mn,m] est la matrice de CAPACITE CALORIFIQUE   a l'instant tn et iteration m
C [Kn,m] est la matrice de CONDUCTIVITE THERMIQUE a l'instant tn et iteration m
C  Fn,m  est le vecteur second membre SOURCE      a l'instant tn et iteration m
C
C A TROUVER Un+1,m+1 SOLUTION DE (1)
C [M n+1,m + DT BETA(1) K n+1,m] U n+1,m+1 = DT BETA(1) Fn+1,m +
C [M n+1,m] ( Un + DT BETA(0) [M n]**-1 ( Fn - [K n] Un) )
C
C JUSQU'A CONVERGENCE des iterations m de POINT FIXE ce qui donne Un+1 puis n=>n
C
C IL A ETE DECIDE DE CALCULER LES MATRICES GLOBALES
C [M n+1] et [M n+1 + DT BETA(1) K n+1]
C en utilisant la TEMPERATURE Un+1 de l'iteration m du point fixe
C
C Apres convergence du point fixe, on a l'egalite sans l'indice m dans (1)
C et cela donne en posant Vn = Un + DT BETA(0) [M n]**-1 ( Fn - [K n] Un )
C [M n+1] ( Un+1 - Vn ) = DT * BETA(1) ( Fn+1 - [K n+1] Un+1 ) ou
C BETA(0)/BETA(1) (Un+1 - Vn) = DT * BETA(0) [M n+1]**-1 ( Fn+1 - [K n+1] Un+1 )
C si bien que
C Vn+1 = Un+1 + DT BETA(0) [M n+1]**-1 ( Fn+1 - [K n+1] Un+1 ) =
C      = Un+1 + BETA(0)/BETA(1) (Un+1 - Vn)
C Vn+1 = (1 + BETA(0)/BETA(1) ) Un+1 - BETA(0)/BETA(1) Vn
C Vn+1 = (BETA(1) + BETA(0))/BETA(1) Un+1 - BETA(0)/BETA(1) Vn
C (2) Vn+1 = 1/BETA(1) Un+1 - BETA(0)/BETA(1) Vn car BETA(1)+BETA(0)=1
C
C Donc a partir du calcul initial de
C V0 = U0 + DT BETA(0) [M 0]**-1 ( F0 - [K 0] U0 )
C les iterations en temps peuvent demarrer et
C le passage de n a n+1 se fera avec (2)
C
C d'ou la necessite de stocker Un puis Un+1,m dans U0, Un+1,m+1 dans U1
C Fn dans F0, Fn+1 dans F1 et Vn dans V0
C
C*******************************************************************************
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET DE THERMIQUE INSTATIONNAIRE A TRAITER
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C D2PI   : NOMBRE 2 x PI DANS UNE VARIABLE REELLE DOUBLE PRECISION
C NDIM   : DIMENSION DES COORDONNEES DES POINTS ( 2 OU 3 )
C NDPGST : CODE D'IDENTIFICATION DES SOMMETS POINTS NOEUDS
C MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DE L'OBJET KNOMOB
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET KNOMOB
C
C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE DE CET OBJET
C MNNPEF : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TMS NPEF"TYPE EF
C NBTTEF : NOMBRE TOTAL D'EF DU MAILLAGE
C MNTPOB : ADRESSE MCN DU TABLEAU POINTEUR SUR LES TABLEAUX POBA DES EF
C NBDLMX : NOMBRE MAXIMAL DE DEGRES DE LIBERTE (TEMPERATURES) D'UN EF
C MOAUX  : NOMBRE DE OTS AUXILIAIRES NECESSAIRES AU CALCUL DES EF
C MNTAUX : ADRESSE MCN DU TABLEAU AUXILIAIRE POUR LES EF
C
C NUMIOB : NUMERO MINIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NUMAOB : NUMERO MAXIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NBOBIN : NOMBRE DE VOLUMES en 3D, SURFACES en 2D DE L'OBJET
C MNOBIN : ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN (PARTIE DE TOPOLOGIE)
C NBOBCL : NOMBRE DE PLS en 3D, PL en 2D DE L'OBJET
C MNOBCL : ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL (PARTIE DE TOPOLOGIE)
C MXTYEL : NOMBRE MAXIMAL DE TYPES D'EF (7)
C MXDOEL : NOMBRE DE MOTS DECLARES DES TABLEAUX D'ADRESSE DANS MNDOEL
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES THERMIQUES DE L'OBJET COMPLET
C
C IEMAST : NOMBRE DE TMS MASSE            DES SV DE L'OBJET RETROUVES
C IECHMA : NOMBRE DE TMS CHALEURMASSIQUE  DES SV DE L'OBJET RETROUVES
C IECOND : NOMBRE DE TMS CONDUCTIVITE     DES SV DE L'OBJET RETROUVES
C IEDILA : NOMBRE DE TMS DILATATION       DES SV DE L'OBJET RETROUVES
C IEVIFL : NOMBRE DE TMS VITESSEFLUIDE    DES SV DE L'OBJET RETROUVES
C IECOET : NOMBRE DE TMS COEFTEMPERATURE  DES SV DE L'OBJET RETROUVES
C IESOIN : NOMBRE DE TMS SOURCE "INTERNE" DES SV DE L'OBJET RETROUVES
C
C IECONT : NOMBRE DE TMS CONTACT              DES PLS DE L'OBJET RETROUVES
C IEECHA : NOMBRE DE TMS ECHANGE              DES PLS DE L'OBJET RETROUVES
C IESOCL : NOMBRE DE TMS SOURCE "AUX LIMITES" DES PLS DE L'OBJET RETROUVES
C IESOPO : NOMBRE DE TMS SOURCE               DES P   DE L'OBJET RETROUVES
C PENALI : COEFFICIENT DE PENALISATION DES TEMPERATURES FIXEES
C          ICI PENALI VAUT 1D20 POUR LE PRENDRE EN COMPTE
C RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
C
C MNTHER : 128 REELS DOUBLE PRECISION POUR LA MATRICE DE CONDUCTIVITE
C MOTAEL : NOMBRE DE MOTS DECLARES DU TABLEAU DES TABLEAUX ELEMENTAIRES
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DU TABLEAU DES 3 COORDONNEES DES NOEUDS=POINTS DE L'EF
C
C NORESO : CODE RESOLUTION DES SYSTEMES LINEAIRES
C          1 FACTORISATION DE CHOLESKY ET MATRICES PROFILS
C          2 GRADIENT CONJUGUE PRECONDITIONNE PAR NIVEAUX ET MATRICES MORSES
C MNLPLI : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPCO : ADRESSE MCN DU NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE
C NIVEAU : NIVEAU DE FACTORISATION INCOMPLETE DE LA MATRICE DE PRECONDITIONNEMEN
C          A CHOISIR PARMI 0 1 2
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DES MATRICES M ET K
C NCODSM : 1 SI MATRICE DE CAPACITE SYMETRIQUE
C          0 SI MATRICE DE CAPACITE DIAGONALE
C         -1 SI MATRICE DE CAPACITE NON SYMETRIQUE
C MNUG0  : ADRESSE MCN DE U0 TEMPERATURE A L'INSTANT INITIAL
C
C BETA   : COEFFICIENT DE LA THETA METHODE
C DT     : PAS CONSTANT DU TEMPS
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"TEMPERATURE
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C NBTEMP : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE
C MXTEMP : NOMBRE DE VECTEUR"TEMPERATURE A STOCKER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1999
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray          Decembre 2022
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      PARAMETER         (ITERMX=16)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__erreurth.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___conductivite.inc"
      include"./incl/a___dilatation.inc"
      include"./incl/a___source.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/trvari.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      CHARACTER*(*)     KNOMOB
C     TABLEAU NON UTILISE (VITESSE D'UN FLUIDE INEXISTANT ICI)
      DOUBLE PRECISION  VITEGt(5,3)

      DOUBLE PRECISION  RELMIN, D2PI
      DOUBLE PRECISION  TECMOY, TECMIN, TECMAX,TEXMIN,TEXMAX
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4)
      LOGICAL           COMP
      DOUBLE PRECISION  PENALI
      DOUBLE PRECISION  BETA(0:1)
      DOUBLE PRECISION  NORMDF, NORMUM, D
C
      DOUBLE PRECISION, allocatable, dimension(:) :: MG
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERMGALLOC, IERKGALLOC
      INTRINSIC         ALLOCATED

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: temp

C     AFFICHAGE
C     =========
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010)
      ELSE
         WRITE(IMPRIM,20010)
      ENDIF
10010 FORMAT(/' CAPACITE CONDUCTIVITE ECHANGE SOURCE CONTACT PEUVENT DEP
     %ENDRE DU TEMPS ET DE LA TEMPERATURE')
20010 FORMAT(/' CAPACITY and CONDUCTIVITY and COEFFICIENTS are DEPENDENT
     % of TIME and/or TEMPERATURE')

C     DECLARATION de temp (BIDON) POUR LE TRACE DES TEMPERATURES
      IALtemp = 1
      allocate( temp( 1:1 ), STAT=IALtemp )
      IF(IALtemp .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'ERREUR en ALLOCATION de temp(1:1)'
         ELSE
            PRINT*,'ALLOCATION ERROR of temp(1:1)'
         ENDIF
      ELSE
         NULLIFY( temp(1)%dptab )
      ENDIF

      NIVMAX = 8

      MDUG0 = (MNUG0+1) / MOREE2
      MNUG1 = 0
      MNFG0 = 0
      MDFG0 = 0
      MNFG1 = 0
      MDFG1 = 0
      MNVG0 = 0
      MDVG0 = 0

      NBTEFX = 0
      MNNTEFX= 0
      MNVTEFX= 0
      MNNFNX = 0
      MNVFNX = 0
      MNDIR  = 0
      MNDAD  = 0
      MNBET  = 0

      MNAUX2 = 0
      MNAUX3 = 0
      MNAUX4 = 0
      MNADIR = 0
      ISTAB  = 0

      IERMGALLOC = 1
      IERKGALLOC = 1
C
C     LES 4 VECTEURS GLOBAUX UG1 FG0 FG1 VG0 QUI S'AJOUTENT A UG0
C     ===========================================================
      CALL TNMCDC( 'REEL2', NTDL, MNUG1 )
      CALL TNMCDC( 'REEL2', NTDL, MNFG0 )
      CALL TNMCDC( 'REEL2', NTDL, MNFG1 )
      CALL TNMCDC( 'REEL2', NTDL, MNVG0 )
      MDUG1 = (MNUG1+1)/MOREE2
      MDFGO = (MNFG0+1)/MOREE2
      MDFG1 = (MNFG1+1)/MOREE2
      MDVGO = (MNVG0+1)/MOREE2

C     DECLARATION DE LA MATRICE GLOBALE DE CAPACITE [MG] ET CONDUCTIVITE [KG]
C     =======================================================================
      IF( NBRDKG .LE. 0 ) THEN
C         PLACE MEMOIRE INSUFFISANTE
          NBLGRC(NRERR) = 3
          WRITE(KERR(MXLGER-1)(1:25),'(G25.0)') DMAT
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) ='ERREUR: PLACE MEMOIRE INSUFFISANTE'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' MOTS NECESSAIRES pour [M] [K] [R]'
             KERR(3) = 'REDUIRE le MAILLAGE'
          ELSE
             KERR(1) ='ERROR: NOT ENOUGH MEMORY'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' NECESSARY WORDS to store [M] [K] [R]'
             KERR(3) ='REDUCE the MESH'
          ENDIF
          CALL LEREUR
          IF( INTERA .LE. 1 ) CALL ARRET( 100 )
          GOTO 9999
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10290) NBRDKG, NBRDKG/NTDL
      ELSE
         WRITE(IMPRIM,20290) NBRDKG, NBRDKG/NTDL
      ENDIF
10290 FORMAT(' 2 MATRICES PROFIL CHACUNE DE',I15,
     %' REELS DOUBLE PRECISION'/
     %' 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT(' 2 SKYLINE MATRICES EACH of',I15,' DOUBLE REALS'/
     %' HALF WIDTH AVERAGE=',I9)
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DES MATRICES PROFIL MG et KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      ALLOCATE ( MG(1:NBRDKG), STAT=IERMGALLOC )
      IF( IERMGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] MATRIX'
         IERR = IERMGALLOC
         GOTO 9999
      ENDIF
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9999
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      WRITE(IMPRIM,*)
C
C     KG & MG SYMETRIQUES NON DIAGONALES PROFIL
      NCODSK = 1
      NCODSM = 1
C
C     CONSTRUCTION DES MATRICE [M0] [K0] ET DU SECOND MEMBRE {FG0}
C     ============================================================
      IEBG = 1
      IEMG = 1
      IEKG = 1
      CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %             D2PI,   NDIM,   NTDL,   VITEGt,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNTHER, MNTAEL, MNX,    MNNODL,
     %             NORESO, MNLPLI, MNLPCO,
     %             NBRDKG, MG,     NBRDKG, KG,   MNFG0,
     %             NCODSM, NCODSK, NBPTAF, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     SI    LA MATRICE DE CAPACITE EST DIAGONALE
C     ALORS ELLE EST RENDUE SYMETRIQUE NON DIAGONALE
C     ==============================================
      IF( NCODSM .EQ. 0 ) THEN
         NCODSM = 1
         DO 20 I=NTDL,1,-1
C           SAUVEGARDE DU COEFFICIENT DIAGONAL AVANT TRANSFERT
            D = MG( I )
C           LE COEFFICIENT EST REMIS A ZERO
            MG( I ) = 0D0
C           LE COEFFICIENT DIAGONAL PROFIL TROUVE SA VALEUR
            MG( MCN(MNLPLI+I) ) = D
 20      CONTINUE
      ENDIF
C
C     CONSTRUCTION NO ET VALEUR DES SOURCES PONCTUELLES POUR {FG0}
C     ============================================================
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C        DES SOURCES PONCTUELLES EXISTENT
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES SOURCES PONCTUELLES DANS FG0
         CALL ASFONO( NTDL, 1, NBFNFX, MCN(MNNFNX), MCN(MNVFNX), RELMIN,
     %                DMCN(MDFG0) )
      ENDIF
C
C     CONSTRUCTION DE V0 = U0 + DT BETA(0) [M0]**-1 ( F0 - [K0] U0 )
C     ==============================================================
C     1. FG1 <= [K0] U0
      IF( NORESO .EQ. 1 ) THEN
         CALL MAPRVE( 0, 1D0, NTDL,
     %                NCODSK, MCN(MNLPLI), KG,  DMCN(MDUG0),
     %                DMCN(MDFG1) )
      ELSE IF ( NORESO .EQ. 2 ) THEN
         CALL MAGCVE( 0, 1D0, NTDL,
     %                NCODSK, MCN(MNLPLI), MCN(MNLPCO), KG,
     %                DMCN(MDUG0),  DMCN(MDFG1) )
      ENDIF
C
C     2. FG0 <= F0 - [K0] U0
      CALL CL2VED ( NTDL,  1D0,  DMCN(MDFG0),
     %                    -1D0,  DMCN(MDFG1),  DMCN(MDFG0) )
C
C     3. [M0]**-1 ( F0 - [K0] U0 )
      IF( NORESO .EQ. 1 .OR. NCODSM .EQ. 0 ) THEN
C
C        FACTORISATION DE CHOLESKY DE [M0]
C        ---------------------------------
         CALL CHOLPR( NTDL, NCODSM, MCN(MNLPLI), MG,
     %                MG, IERR )
         IF( IERR .NE. 0 ) THEN
C            MATRICE NON INVERSIBLE
             NBLGRC(NRERR) = 1
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR: MATRICE CAPACITE NON INVERSIBLE'
             ELSE
                KERR(1) = 'ERROR: CAPACITY MATRIX NOT INVERSIBLE'
             ENDIF
             CALL LEREUR
             IERR = 7
             GOTO 9999
         ENDIF
C
C        RESOLUTION DU SYSTEME FACTORISE PAR CHOLESKY
C        --------------------------------------------
         CALL DRCHPR( NTDL, NCODSM, MCN(MNLPLI), MG, DMCN(MDFG0), 2,
     %                DMCN(MDFG0) )
C
      ELSE IF ( NORESO .EQ. 2 ) THEN
C
C        GC => MATRICE MORSE CONSTRUCTION DE LA MATRICE DE PRECONDITIONNEMENT
C        --------------------------------------------------------------------
C        CALCUL DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
         ISTAB  = 0
C        CALCUL DES POINTEURS DE LA MATRICE DE PRECONDITIONNEMENT
C        (FACTORISATION INCOMPLETE DE A SUIVANT LE NIVEAU)
 40      IERR   = 0
         MNLPLC = 0
         MNLPCC = 0
         MNLPLU = 0
         CALL CALPNT( NTDL,   NIVEAU, NCODAG, MNLPLI, MNLPCO,
     %                MNLPLC, MNLPCC, MNLPLU, COMP,   IERR  )
C        VERIFICATION DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
         IF( IERR .NE. 0 ) THEN
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='PROBLEME DE FACTORISATION INCOMPLETE DE NIVEAU '
     %               // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
               KERR(3) = '=> AUGMENTER NIVEAU ou CHOISIR CHOLESKY'
            ELSE
               KERR(1)='PROBLEM of INCOMPLETE FACTORIZATION of LEVEL '
     %               // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               KERR(3) = '=> AUGMENT LEVEL or CHOOSE CHOLESKY'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
C
C        DECLARATION INITIALISATION DE LA MATRICE MORSE DE CONDITIONNEMENT
C        -----------------------------------------------------------------
         LOLPCC = MCN(MNLPLC+NTDL)
         MNAGC  = 0
         CALL TNMCDC( 'REEL2', LOLPCC, MNAGC )
         CALL AZEROD( LOLPCC, MCN(MNAGC) )
C
C        CONSTRUCTION EFFECTIVE DE LA MATRICE DE PRECONDITIONNEMENT
C        PAR FACTORISATION INCOMPLETE DE CHOLESKY AVEC NIVEAU
C        ----------------------------------------------------------
         NBLGRC(NRERR) = 1
         WRITE(KERR(2)(1:2),'(I2)') NIVEAU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FACTORISATION INCOMPLETE de [M] de NIVEAU'
     %              // KERR(2)(1:2)
         ELSE
            KERR(1) = 'INCOMPLETE FACTORIZATION of [M] of LEVEL'
     &              // KERR(2)(1:2)
         ENDIF
         CALL LERESU
         CALL INCHGC( NTDL,
     %                MCN(MNLPLI), MCN(MNLPCO), MG,
     %                MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC), IERR )
C
C        VERIFICATION DE LA STABILITE DE LA FACTORISATION INCOMPLETE
C        QUI DEFINIT LA MATRICE DE PRECONDITIONNEMENT DU GRADIENT CONJUGUE
C        -----------------------------------------------------------------
         IF( IERR .EQ. 1 ) THEN
C           AU MOINS UN PIVOT<=0 LE PIVOT QUI PRECEDE EST PRIS A SA PLACE
            ISTAB = 1
            IERR  = 0
         ELSE IF( IERR .EQ. 2 ) THEN
C           LE NOMBRE DE PIVOTS INCORRECTS<=0 EST DEPASSE
C           FACTORISATION EST DECLAREE INSTABLE
            IF( COMP ) GO TO 9999
            IF( NIVEAU .LT. NIVMAX ) THEN
C              TENTATIVE D'AUGMENTER LE NIVEAU DE FACTORISATION INCOMPLETE
               NIVEAU = NIVEAU + 1
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
               CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
               CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
               GO TO 40
            ELSE
C              LA FACTORISATION EST VRAIMENT TRES INSTABLE
C              PLUS DE NIVMAX NIVEAUX DEPASSE => ABANDON DU GC
               WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
               NBLGRC(NRERR) = 3
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     %                    // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
                  KERR(3) = '=> AUGMENTER NIVEAU ou CHOISIR CHOLESKY'
               ELSE
                  KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %               // KERR(MXLGER)(1:8)
                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
                  KERR(3) = '=> AUGMENT LEVEL or CHOOSE CHOLESKY'
               ENDIF
C              DESTRUCTION DES TABLEAUX DES MATRICES DU GC
               CALL LEREUR
               IERR = 8
               GOTO 9999
            ENDIF
         ENDIF
C
C        LE PREMIER TABLEAU AUXILIAIRE DE GCPRCH
         MNAUGC = 0
         LOAUGC = NTDL * 3
         CALL TNMCDC( 'REEL2', LOAUGC, MNAUGC )
C        REPARTITION INTERNE EN SOUS-TABLEAUX
         MNAUX2 = MNAUGC
         MNAUX3 = MNAUX2 + NTDL * MOREE2
         MNAUX4 = MNAUX3 + NTDL * MOREE2
C
C        LE SECOND TABLEAU AUXILIAIRE DE GCPRCH SELON LA
C        STABILISATION DU GC PAR RE-ORTHOGONALISATION DES DIRECTIONS
         IF( ISTAB .EQ. 0 ) THEN
             NBDIR = 1
         ELSE
             NBDIR = 10
         ENDIF
         LODIR = (NTDL+1) * NBDIR * 2
         MNBDIR = 0
         CALL TNMCDC( 'REEL2', LODIR, MNBDIR )
C        REPARTITION INTERNE EN SOUS-TABLEAUX
         MNDIR  = MNBDIR
         MNADIR = MNBDIR + NTDL * NBDIR * MOREE2
         MNDAD  = MNADIR + NTDL * NBDIR * MOREE2
         MNBET  = MNDAD  + NBDIR * MOREE2
         MOTSGC = MOREE2 * (LOAUGC+LODIR)
         WRITE(IMPRIM,10150) NTDL+LOLPCC*(1+MOREE2), MOTSGC
C
C        RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE
C        -------------------------------------------
C        LE VECTEUR INITIAL DU GC EST LE VECTEUR UG0
C        MISE A ZERO DU TABLEAU DES DIRECTIONS
         CALL AZEROD( LODIR, MCN(MNBDIR) )
         CALL GCPRCH( NTDL,        1,           NBDIR,
     %                MCN(MNLPLI), MCN(MNLPCO), MG,         DMCN(MDFG0),
     %                MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC),
     %               DMCN(MDUG0),  MCN(MNAUX2), MCN(MNAUX3),
     %                MCN(MNAUX4), MCN(MNDIR),  MCN(MNADIR),
     %                MCN(MNDAD),  MCN(MNBET),  DMCN(MDFG0),  IERR )
C
C        VERIFICATION DE LA CONVERGENCE DU GC
         IF( IERR .NE. 0 ) THEN
C           IL N'Y A PAS EU CONVERGENCE DU GC
C           LA METHODE DU GC NE CONVERGE PAS => ABANDON DES CALCULS
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='NON CONVERGENCE du GC au NIVEAU '
     %               // KERR(MXLGER)(1:8)
               KERR(2)='ABANDON de la METHODE du GRADIENT CONJUGUE'
               KERR(3)='=> ESSAYER METHODE de CHOLESKY MATRICE PROFIL'
            ELSE
               KERR(1) = 'NO CONVERGENCE of CG at LEVEL '
     &                 // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               KERR(3) = '=> TRY direct CHOLESKY SKYLINE METHOD'
            ENDIF
            IERR = 9
            GOTO 9999
         ENDIF
      ENDIF
C
C     4. V0 <= U0 + DT BETA(0) [M0]**-1 ( F0 - [K0] U0 )
      CALL CL2VED ( NTDL, 1D0,        DMCN(MDUG0),
     %                    DT*BETA(0), DMCN(MDFG0),  DMCN(MDVG0) )
C
C     BILAN DE L'INSTANT INITIAL
C     ==========================
C     MNFG0  CONTIENT  F(t0) SECOND MEMBRE A L'INSTANT INITIAL
C     MNUG0  CONTIENT  U0    TEMPERATURE   A L'INSTANT INITIAL
C     MNVG0  CONTIENT  U0 + DT BETA(0) [M0]**-1 ( F0 - [K0] U0 )U0
C
C     LE PROCHAIN TEMPS POUR STOCKER LE VECTEUR TEMPERATURE
      TSTOC  = TPSINI + DTSTOC
      MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBTEMP-1)
      MNTIME = MNVECT + WECTEU + MOREE2 * NTDL * MXTEMP - 1
C     LE DEBUT DU VECTEUR TEMPERATURE INITIALE
      MNTHET = MNTEMP

C     AFFICHAGE DE LA TEMPERATURE AUX NOEUDS AU TEMPS INITIAL
      NUMCAS = 1
      CALL AFTEMP( 3,      NUMCAS, MNXYZN,
     %             NTDL,   NUMCAS, MCN(MNUG0),
     %             TECMOY, TECMIN, NOTMIN, TECMAX, NOTMAX,
     %             NOFOTI, TEXMIN, TEXMAX )


C     ##############################################################
C     ##                                                          ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DT ##
C     ##                                                          ##
C     ##############################################################
C
C     LE NOUVEAU TEMPS OU SE FAIT LE CALCUL
C     =====================================
 100  TEMPS = TEMPS + DT
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) TEMPS
      ELSE
         WRITE(IMPRIM,20100) TEMPS
      ENDIF
10100 FORMAT(/' theinl: AU TEMPS',G14.6,' CALCUL DES TEMPERATURES')
20100 FORMAT(/' theinl: At TIME',G14.6,' COMPUTATION of TEMPERATURES')
C
C     LE NOMBRE D'ITERATIONS DU POINT FIXE
      ITERPF = 0
C
C     =======================================================================
C     A TROUVER Un+1,m+1 SOLUTION DE (1)
C     [Mn+1,m + DT BETA(1) Kn+1,m] Un+1,m+1 = DT BETA(1) Fn+1,m + [Mn+1,m] Vm
C     =======================================================================
C     L'ADRESSE DE LA TEMPERATURE ACTUELLE (n+1,m)
 110  MNTHET = MNUG0
C
C     1. CONSTRUCTION DE [Mn+1,m] [Kn+1,m] Fn+1,m
C     -------------------------------------------
      IEBG = 1
      IEMG = 1
      IEKG = 1
      CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %             D2PI,   NDIM,   NTDL,   VITEGt,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNTHER, MNTAEL, MNX,    MNNODL,
     %             NORESO, MNLPLI, MNLPCO,
     %             NBRDKG, MG,     NBRDKG, KG,     MNFG1,
     %             NCODSM, NCODSK, NBPTAF, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     SI    LA MATRICE DE CAPACITE EST DIAGONALE
C     ALORS ELLE EST RENDUE SYMETRIQUE NON DIAGONALE
      IF( NCODSM .EQ. 0 ) THEN
         NCODSM = 1
         DO 115 I=NTDL,1,-1
C           SAUVEGARDE DU COEFFICIENT DIAGONAL AVANT TRANSFERT
            D = MG( I )
C           LE COEFFICIENT EST REMIS A ZERO
            MG( I ) = 0D0
C           LE COEFFICIENT DIAGONAL PROFIL TROUVE SA VALEUR
            MG( MCN(MNLPLI+I) ) = D
 115     CONTINUE
      ENDIF
C
C     CONSTRUCTION NO ET VALEUR DES SOURCES PONCTUELLES POUR {FG0}
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C        DES SOURCES PONCTUELLES EXISTENT
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES SOURCES PONCTUELLES DANS FG0
         CALL ASFONO( NTDL, 1, NBFNFX, MCN(MNNFNX), MCN(MNVFNX), RELMIN,
     %                DMCN(MDFG0) )
      ENDIF
C
C     2. FG0 <= [Mn+1,m] Vm
C     ---------------------
      IF( NORESO .EQ. 1 ) THEN
         CALL MAPRVE( 0, 1D0, NTDL,
     %                NCODSM, MCN(MNLPLI), MG, DMCN(MDVG0),
     %                DMCN(MDFG0) )
      ELSE IF ( NORESO .EQ. 2 ) THEN
         CALL MAGCVE( 0, 1D0, NTDL,
     %                NCODSM, MCN(MNLPLI), MCN(MNLPCO), MG,
     %                DMCN(MDVG0),  DMCN(MDFG0) )
      ENDIF
C
C     3. FG0 <= DT BETA(1) Fn+1,m + [Mn+1,m] Vm
C     -----------------------------------------
      CALL CL2VED( NTDL, DT*BETA(1), DMCN(MDFG1),
     %                   1D0,        DMCN(MDFG0), DMCN(MDFG0) )
C
C     4. KG <= [Mn+1,m + DT BETA(1) Kn+1,m]
C     -------------------------------------
      CALL CL2VED ( NBRDKG, 1D0,        MG,
     %                      DT*BETA(1), KG,  KG )
C
C     5. UG1 <= [Mn+1,m + DT BETA(1) Kn+1,m]**-1  DT BETA(1) Fn+1,m + [Mn+1,m] V
C     --------------------------------------------------------------------------
      IF( NORESO .EQ. 1 .OR. NCODSM .EQ. 0 ) THEN
C
C        FACTORISATION DE CHOLESKY DE [Mn+1,m + DT BETA(1) Kn+1,m]
         CALL CHOLPR( NTDL, NCODSK, MCN(MNLPLI), KG,
     %                KG, IERR )
         IF( IERR .NE. 0 ) THEN
C            MATRICE NON INVERSIBLE
             NBLGRC(NRERR) = 1
             IF( LANGAG .EQ. 0 ) THEN
             KERR(1)='ERREUR: MATRICE [M]+DT BETA(1) [K] NON INVERSIBLE'
             ELSE
             KERR(1)='ERROR: MATRIX [M]+DT BETA(1) [K] NOT INVERSIBLE'
             ENDIF
             CALL LEREUR
             IERR = 7
             GOTO 9999
         ENDIF
C
C        RESOLUTION DU SYSTEME FACTORISE PAR CHOLESKY
         CALL DRCHPR( NTDL, NCODSK, MCN(MNLPLI), KG, DMCN(MDFG0), 2,
     %                DMCN(MDUG1) )
C
      ELSE IF ( NORESO .EQ. 2 ) THEN
C
C        CONSTRUCTION EFFECTIVE DE LA MATRICE DE PRECONDITIONNEMENT
C        PAR FACTORISATION INCOMPLETE DE CHOLESKY DE [M] +DT BETA(1)[K] AVEC NIV
C        SEULEMENT UNE FOIS PAR PAS DE TEMPS LORS DE LA PREMIERE ITERATION DE PO
C        -----------------------------------------------------------------------
 140     IF( ITERPF .EQ. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(2)(1:2),'(I2)') NIVEAU
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='FACTORISATION INCOMPLETE de [M] + DT BETA(1)[K]
     %de NIVEAU' // KERR(2)(1:2)
            ELSE
               KERR(1)='INCOMPLETE FACTORIZATION of [M] + DT BETA(1)[K]
     %of LEVEL' // KERR(2)(1:2)
            ENDIF
            CALL LERESU
            CALL AZEROD( LOLPCC, MCN(MNAGC) )
            CALL INCHGC( NTDL,
     %                   MCN(MNLPLI),MCN(MNLPCO),KG,
     %                   MCN(MNLPLC),MCN(MNLPCC),MCN(MNAGC), IERR )
C
C           VERIFICATION DE LA STABILITE DE LA FACTORISATION INCOMPLETE
C           QUI DEFINIT LA MATRICE DE PRECONDITIONNEMENT DU GRADIENT CONJUGUE
            IF( IERR .EQ. 1 ) THEN
C              AU MOINS UN PIVOT<=0 LE PIVOT QUI PRECEDE EST PRIS A SA PLACE
               ISTAB = 1
               IERR  = 0
            ELSE IF( IERR .EQ. 2 ) THEN
C              LE NOMBRE DE PIVOTS INCORRECTS<=0 EST DEPASSE
C              FACTORISATION EST DECLAREE INSTABLE
               IF( COMP ) GO TO 9999
               IF( NIVEAU .LT. NIVMAX ) THEN
C                 TENTATIVE D'AUGMENTER LE NIVEAU DE FACTORISATION INCOMPLETE
                  NIVEAU = NIVEAU + 1
                  CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
                  CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
                  CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
                  IERR   = 0
                  MNLPLU = 0
                  CALL CALPNT( NTDL,   NIVEAU, NCODAG, MNLPLI, MNLPCO,
     %                         MNLPLC, MNLPCC, MNLPLU, COMP,   IERR  )
C                 VERIFICATION DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
                  IF( IERR .NE. 0 ) THEN
                     WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
                     NBLGRC(NRERR) = 3
               KERR(1)='PROBLEME DE FACTORISATION INCOMPLETE DE NIVEAU '
     %               // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
                     KERR(3) = 'AUGMENTER NIVEAU OU CHOISIR CHOLESKY'
                     CALL LEREUR
                     GOTO 9999
                  ENDIF
C                 DECLARATION INITIALISATION DE LA MATRICE MORSE DE CONDITIONNEM
                  LOLPCC = MCN(MNLPLC+NTDL)
                  CALL TNMCDC( 'REEL2', LOLPCC, MNAGC )
                  GO TO 140
               ELSE
C                 LA FACTORISATION EST VRAIMENT TRES INSTABLE
C                 PLUS DE NIVMAX NIVEAUX DEPASSE => ABANDON DU GC
                  WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
                  NBLGRC(NRERR) = 3
                  IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     %                    // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
                  KERR(3) = '=> AUGMENTER NIVEAU ou CHOISIR CHOLESKY'
                  ELSE
                  KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %               // KERR(MXLGER)(1:8)
                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
                  KERR(3) = '=> AUGMENT LEVEL or CHOOSE CHOLESKY'
                  ENDIF
                  CALL LEREUR
C                 DESTRUCTION DES TABLEAUX DES MATRICES DU GC
                  IERR = 8
                  GOTO 9999
               ENDIF
            ENDIF
C
C           LE SECOND TABLEAU AUXILIAIRE DE GCPRCH SELON LA
C           STABILISATION DU GC PAR RE-ORTHOGONALISATION DES DIRECTIONS
            IF( ISTAB .EQ. 0 ) THEN
                NBDIR = 1
            ELSE
                NBDIR = 10
            ENDIF
            I = (NTDL+1) * NBDIR * 2
            IF( LODIR .LT. I ) THEN
C              LE TABLEAU DEMANDE A ETRE AUGMENTE
               CALL TNMCAU( 'REEL2', LODIR, I, 0, MNBDIR )
               LODIR = I
C              REPARTITION INTERNE EN SOUS-TABLEAUX
               MNDIR  = MNBDIR
               MNADIR = MNBDIR + NTDL * NBDIR * MOREE2
               MNDAD  = MNADIR + NTDL * NBDIR * MOREE2
               MNBET  = MNDAD  + NBDIR * MOREE2
               MOTSGC = MOREE2*(LOAUGC+LODIR)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10150) NTDL+LOLPCC*(1+MOREE2), MOTSGC
               ELSE
                  WRITE(IMPRIM,20150) NTDL+LOLPCC*(1+MOREE2), MOTSGC
               ENDIF
            ENDIF
         ENDIF
10150 FORMAT(' PRECONDITIONNEMENT MATRICE MORSE =',I15,' MOTS'/
     %       ' TABLEAUX SUPPLEMENTAIRES GC      =',I15,' MOTS')
20150 FORMAT(' PRECONDITIONED MORSE MATRIX =',I15,' MEMORY WORDS'/
     %       ' AUXILIARY ARRAYS CG =',I15,' MEMORY WORDS')
C
C        RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE
C        LE VECTEUR INITIAL DU GC EST LE VECTEUR UG0
C        MISE A ZERO DU TABLEAU DES DIRECTIONS
         CALL AZEROD( LODIR, MCN(MNBDIR) )
         CALL GCPRCH( NTDL,        1,           NBDIR,
     %                MCN(MNLPLI), MCN(MNLPCO), KG,         DMCN(MDFG0),
     %                MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC),
     %               DMCN(MDUG0),  MCN(MNAUX2), MCN(MNAUX3),
     %                MCN(MNAUX4), MCN(MNDIR),  MCN(MNADIR),
     %                MCN(MNDAD),  MCN(MNBET),  DMCN(MDUG1),  IERR )
C
C        VERIFICATION DE LA CONVERGENCE DU GC
         IF( IERR .NE. 0 ) THEN
C           IL N'Y A PAS EU CONVERGENCE DU GC
C           LA METHODE DU GC NE CONVERGE PAS => ABANDON DES CALCULS
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='NON CONVERGENCE du GC au NIVEAU '
     %               // KERR(MXLGER)(1:8)
               KERR(2)='ABANDON de la METHODE du GRADIENT CONJUGUE'
               KERR(3)='=> ESSAYER METHODE de CHOLESKY MATRICE PROFIL'
            ELSE
               KERR(1) = 'NO CONVERGENCE of CG at LEVEL '
     &                 // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               KERR(3) = '=> TRY direct CHOLESKY SKYLINE METHOD'
            ENDIF
            IERR = 9
            GOTO 9999
         ENDIF
      ENDIF

C     LISTE DES NUMEROS ET VALEURS DES TEMPERATURES FIXEES AU TEMPS tn+1
C     ==================================================================
      CALL THDLFX( 1,      NTDL,   NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBTEFX, MONTEFX, MNNTEFX, MNVTEFX, IERR )

C     PRISE EN COMPTE DES CONDITIONS AUX LIMITES de TEMPERATURE SUR MNUG1
C     NBTEFX TEMPERATURES SONT FIXEES AUX VALEURS du TABLEAU CONTACT des PLSV
C     =======================================================================
      IF( NBTEFX .GT. 0 ) THEN
         CALL BLDLFX( NTDL, 1, NBTEFX, MCN(MNNTEFX), MCN(MNVTEFX),
     %                1D0, DMCN(MDUG1) )
      ENDIF

C     Y A T IL CONVERGENCE DU POINT FIXE?
C     ===================================
C     CALCUL DE || Um+1 - Um || et ||Um+1||
      MN0 = ( MNUG0 - 1 ) / 2
      MN  = ( MNUG1 - 1 ) / 2
      NORMDF = 0D0
      NORMUM = 0D0
      DO 160 I=1,NTDL
         NORMDF = NORMDF + ( DMCN(MN+I) - DMCN(MN0+I) ) ** 2
         NORMUM = NORMUM +   DMCN(MN+I) ** 2
 160  CONTINUE
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10160) TEMPS, ITERPF+1,
     %                       NORMUM, NORMDF, NORMDF/NORMUM
      ELSE
         WRITE(IMPRIM,20160) TEMPS, ITERPF+1,
     %                       NORMUM, NORMDF, NORMDF/NORMUM
      ENDIF
10160 FORMAT(' Au TEMPS',G14.6,' ITERATION',I3,
     %       ' ||Un||=',G12.4,' ||Un-Un-1||=',G12.4,
     %       ' ||Un-Un-1||/||Un||=',G10.2)
20160 FORMAT(' At TIME',G14.6,' ITERATION',I3,
     %       ' ||Un||=',G12.4,' ||Un-Un-1||=',G12.4,
     %       ' ||Un-Un-1||/||Un||=',G10.2)
C
      IF( NORMDF .GT. 1D-4 * NORMUM ) THEN
C
C        NON CONVERGENCE => UNE ITERATION DE POINT FIXE DE PLUS
         IF( ITERPF .GE. ITERMX ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: POINT FIXE NON ATTEINT apres'
            ELSE
               KERR(1) = 'ERROR: FIX POINT NOT CONVERGED after'
            ENDIF
            KERR(2) = KERR(5)(1:6) // ' ITERATIONS'
            CALL LEREUR
            IERR = 29
            GOTO 9999
         ENDIF
         ITERPF = ITERPF + 1

C        UG1 CONTIENT UGn+1,m+1 => PERMUTATION AVEC UG0
         MN    = MNUG0
         MNUG0 = MNUG1
         MNUG1 = MN
         MDUG0 = (MNUG0+1)/MOREE2
         MDUG1 = (MNUG1+1)/MOREE2
         GOTO 110

      ENDIF
C
C     CONVERGENCE DU POINT FIXE
C     =========================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) TEMPS, ITERPF+1
      ELSE
         WRITE(IMPRIM,20170) TEMPS, ITERPF+1
      ENDIF
10170 FORMAT(' AU TEMPS',G14.6,': CONVERGENCE APRES',I5,
     %       ' ITERATION(S) DE POINT FIXE')
20170 FORMAT(' At TIME',G14.6,' CONVERGENCE after',I5,
     %       ' FIX POINT ITERATIONS')
C
      IF( TEMPS .GE. TSTOC*0.9999 ) THEN
C        STOCKAGE DE LA TEMPERATURE A CET INSTANT TEMPS
         MNTEMP = MNTEMP + NTDL * MOREE2
         CALL TRTATD( MCN(MNUG1), MCN(MNTEMP), NTDL )
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBTEMP = NBTEMP + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBTEMP) = TEMPS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10180) TEMPS,NTDL,NBTEMP
         ELSE
            WRITE(IMPRIM,20180) TEMPS,NTDL,NBTEMP
         ENDIF
10180 FORMAT(' Au TEMPS',G14.6,' STOCKAGE DE',I7,' TEMPERATURES',
     %   ' en COLONNE',I5,' du TMS VECTEUR"TEMPERATURE')
20180 FORMAT(' At TIME',G14.6,' STORAGE of',I9,' TEMPERATURES',
     %   ' in the COLUMN',I5,' of TMS VECTEUR"TEMPERATURE')
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC
      ENDIF

C     AFFICHAGE DE LA TEMPERATURE AUX NOEUDS AU TEMPS
      CALL AFTEMP( 3,      NUMCAS, MNXYZN,
     %             NTDL,   NUMCAS, MCN(MNUG1),
     %             TECMOY, TECMIN, NOTMIN, TECMAX, NOTMAX,
     %             NOFOTI, TEXMIN, TEXMAX )

C     TRACE DES ZONES DE COULEURS des TEMPERATURES en 2D et ISOTHERMES en 3D
C     ======================================================================
      IF( INTERA .GE. 1 ) THEN
C        MODE GRAPHIQUE AVEC X11: LA MEMOIRE PIXELS EST EFFACEE
         CALL EFFACEMEMPX
         CALL VISEE0

C        TRACE DES ARETES DES FACES AVEC LA COULEUR
         NCOUAF = NCGRIM
C        POURCENTAGE DE REDUCTION DES FACES
         PREDUF = 0.
C        COULEUR PAR DEFAUT DES ARETES DES FACES FRONTIERE
         NCOAFR = NCGRIC
         LORBITE = 0
         NCAS = 1
         TMIN = REAL( TECMIN )
         TMAX = REAL( TECMAX )

         IF( NDIM .EQ. 2 ) THEN
ccc            PRINT *,'trzont Temps=',TEMPS,'  ---------->'
            CALL TRZONT( 0,      NDIM,    KNOMOB, -11,
     %                   NBTYEL, MNNPEF,  MNXYZN, MNXYZN, NDPGST,
     %                   NCAS,   NCAS,    NTDLTE,
     %                   0,      DMCN(MDUG1),  temp,
     %                   TMIN,   NOTMIN, NCAS, TMAX, NOTMAX, NCAS,
     %                   TEMPS )
         ELSE
C           LA VISEE EN 3D
C           LONGITUDE et LATITUDE
ccc            CALL LONLAT( -80.0, 10. )
ccc            CALL LONLAT( -82.0, 8. )
ccc            CALL LONLAT( -110., 20. )
            CALL LONLAT( -85., 3. )
C           LOUPE GROSSISSANTE
            GROSSI = 0.75
            AXOLAR = AXOLAR / GROSSI
            AXOHAU = AXOHAU / GROSSI

ccc            PRINT *,'trisot Temps=',TEMPS,'  ---------->'
            CALL TRISOT( NDIM,   KNOMOB, -11,
     %                   NBTYEL, MNNPEF, MNXYZN, NDPGST,
     %                   NCAS,   NCAS,   NTDLTE,
     %                   0,      DMCN(MDUG1),  temp,
     %                   TMIN,   NOTMIN, NCAS, TMAX, NOTMAX, NCAS,
     %                   TEMPS )
         ENDIF
      ENDIF


C     CALCUL DE VG0 <= (2) Vn+1 = 1/BETA(1) Un+1 - BETA(0)/BETA(1) Vn
C     ===============================================================
      CALL CL2VED ( NTDL,      1D0/BETA(1),  DMCN(MDUG1),
     %                    -BETA(0)/BETA(1),  DMCN(MDVG0), DMCN(MDVG0) )

C     MISE A JOUR DE MNUG0
C     ====================
      MN    = MNUG0
      MNUG0 = MNUG1
      MNUG1 = MN
      MDUG0 = (MNUG0+1)/MOREE2
      MDUG1 = (MNUG1+1)/MOREE2

      IF( TEMPS + DT .LT. TPSFIN*1.00001 ) GOTO 100
C
C    ##############################################################
C    ##                                                          ##
C    ##                FIN DE LA BOUCLE EN TEMPS                 ##
C    ##                                                          ##
C    ##############################################################

      IF( RMCN(MNTIME+NBTEMP) .NE. TEMPS ) THEN
C        STOCKAGE DES TEMPERATURES A CET INSTANT
         MNTEMP = MNTEMP + NTDL * MOREE2
         CALL TRTATD( MCN(MNUG0), MCN(MNTEMP), NTDL )
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBTEMP = NBTEMP + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBTEMP) = TEMPS
         WRITE(IMPRIM,10180) TEMPS,NTDL,NBTEMP
      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"TEMPERATURE'
C     ========================================
      MCN( MNVECT + WBCOVE ) = NTDL
      MCN( MNVECT + WBVECT ) = NBTEMP
      MCN( MNVECT + WBCPIN ) = NBTEMP
      IF( NBTEMP .LT. MXTEMP ) THEN
C        LE TMS EST RACOURCI
         L  = MNVECT + WECTEU + NTDL * MXTEMP * MOREE2 - 1
         L1 = MNVECT + WECTEU + NTDL * NBTEMP * MOREE2 - 1
         DO 200 I=1,NBTEMP
            RMCN(L1+I) = RMCN(L+I)
200      CONTINUE
         CALL TAMSRA( NTVECT, WECTEU+NTDL*NBTEMP*MOREE2+NBTEMP )
      ENDIF
C     LA DATE
      CALL ECDATE( MCN(MNVECT) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVECT + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
C
C     LES TEMPERATURES A L'INSTANT FINAL SONT A L'ADRESSE MNTEMP
      MNTEMP = MNVECT + WECTEU
C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS TEMPERATURE
      MNTIME = MNTEMP + NTDL * NBTEMP * MOREE2 -1
C
C     AFFICHAGE DES TEMPERATURES
C     ==========================
      CALL AFTEMP( 3,      NBTEMP, MNXYZN,
     %             NTDL,   NBTEMP, MCN(MNTEMP),
     %             TECMOY, TECMIN, NOTMIN, TECMAX, NOTMAX,
     %             NOFOTI, TEXMIN, TEXMAX )
C
C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
 9999 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( IERMGALLOC .EQ. 0 ) DEALLOCATE( MG )
      IF( MNVG0  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNVG0  )
      IF( MNUG1  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNUG1  )
      IF( MNFG1  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNFG1  )
      IF( MNUG0  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNUG0  )
      IF( MNFG0  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNFG0  )
      IF( MNNFNX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONFNX, MNNFNX )
      IF( MNVFNX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOVFNX, MNVFNX )
      IF( NBTEFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTEFX, MNNTEFX )
      IF( NBTEFX .GT. 0 ) CALL TNMCDS( 'REEL2',  NBTEFX, MNVTEFX )

      IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
      IF( NORESO .EQ. 2 ) THEN
          IF( MNAUGC .GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUGC, MNAUGC )
          IF( MNBDIR .GT. 0 ) CALL TNMCDS( 'REEL2',  LODIR,  MNBDIR )
          IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCO )
          IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
          IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
          IF( MNAGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
      ENDIF
C
C     BILAN SUR LA PLACE MEMOIRE CENTRALE OCCUPEE PAR LES MATRICES ...
C     ================================================================
      IF( NORESO .EQ. 1 ) THEN
C        CHOLESKY PROFIL
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19002) NTDL + MOREE2 * NBRDKG * 2
         ELSE
            WRITE(IMPRIM,29002) NTDL + MOREE2 * NBRDKG * 2
         ENDIF
19002    FORMAT(/' STOCKAGE MATRICES PROFIL =',I15,' MOTS'/ )
29002    FORMAT(/' MATRIX SKYLINE STORAGE =',I15,' MEMORY WORDS' )
      ELSE IF( NORESO .EQ. 2 ) THEN
C        GRADIENT CONJUGUE MORSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19004) NIVEAU, NTDL+NBRDKG*(1+MOREE2)*2,
     %                          NTDL+LOLPCC*(1+MOREE2),MOTSGC,
     %      NTDL+NBRDKG*(1+MOREE2)*2+NTDL+LOLPCC*(1+MOREE2)+MOTSGC
         ELSE
            WRITE(IMPRIM,29004) NIVEAU, NTDL+NBRDKG*(1+MOREE2)*2,
     %                          NTDL+LOLPCC*(1+MOREE2),MOTSGC,
     %      NTDL+NBRDKG*(1+MOREE2)*2+NTDL+LOLPCC*(1+MOREE2)+MOTSGC
         ENDIF
      ENDIF
19004    FORMAT(/' FACTORISATION INCOMPLETE de NIVEAU ',I2/
     %           ' MATRICES MORSE GC ',         T30,I15,' MOTS'/
     %           ' MATRICE PRECONDITIONNEMENT', T30,I15,' MOTS'/
     %           ' TABLEAUX SUPPLEMENTAIRES GC',T30,I15,' MOTS'/
     %           ' AU TOTAL LE GC DEMANDE',     T30,I15,' MOTS'/)
29004    FORMAT(/' INCOMPLETE FACTORIZATION of LEVEL ',I2/
     %           ' CG MORSE MATRIX  ',      T30,I15,' MEMORY WORDS'/
     %           ' PRECONDITIONED MATRIX',  T30,I15,' MEMORY WORDS'/
     %           ' CG AUXILIARY ARRAYS',    T30,I15,' MEMORY WORDS'/
     %           ' TOTAL: The CG REQUIRES', T30,I15,' MEMORY WORDS'/)
      RETURN
      END
