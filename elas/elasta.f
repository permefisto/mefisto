      SUBROUTINE ELASTA( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES DEPLACEMENTS ET CONTRAINTES DANS UN DOMAINE
C -----    2D OU 3D OU AXISYMETRIQUE EN ELASTICITE LINEAIRE STATIONNAIRE
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET D'ELASTICITE STATIONNAIRE A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION
C          NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1989
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___young.inc"
      include"./incl/a___dilatation.inc"
      include"./incl/a___coefdeplacement.inc"
      include"./incl/a___force.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/a___contrinit.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERKGALLOC
      INTRINSIC         ALLOCATED
C
      EXTERNAL          EETAEL
      DOUBLE PRECISION  RELMIN,D2PI,PENALI
      DOUBLE PRECISION  DEXMAX,DECMAX
      DOUBLE PRECISION  DINFO,DCPU,DFABG,DFACT,DDEPL
      INTEGER           NUMIOB(4),NUMAOB(4),MNDOEF(4),MXDOEF(4)
      CHARACTER*(*)     KNOMOB
      LOGICAL           COMP
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('*')/
     %' RESOLUTION DE L''ELASTICITE STATIONNAIRE DE L''OBJET: ',
     %A/1X,80('*')/)
20000 FORMAT(/1X,80('*')/
     %' COMPUTATION of the STEADY ELASTICITY of the OBJECT: ',
     %A/1X,80('*')/)
C
C     QUELQUES INITIALISATIONS
      DCPU   = DINFO( 'DELTA CPU' )
      TEMPS  = 0.0
C     PENALI : EPSILON DE LA PENALISATION DE LA FIXATION
C              0 => DIRICHLET PUR SANS PENALISATION
      PENALI = 0D0
      DDEPL  = 0
      IERR   = 0
      D2PI   = ATAN( 1D0 ) * 8D0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNNPEF = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNNODL = 0
      MNIP   = 0
      MNELAS = 0
      MNFORC = 0
      MNX    = 0
      MNB    = 0
      MNNFNX = 0
      MNVFNX = 0
      MNNDLX = 0
      MNVDLX = 0
      MNAUX  = 0
      MNAUXE = 0
      DO 5 I=1,4
         MNDOEF(I) = 0
         MXDOEF(I) = 0
 5    CONTINUE
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NDSM   = 1
      NBDLMX = 0
      NTDL   = 0
      NBRFNX = 0
      NBRDLX = 0
      MOFLTO = 0
      MOFLPT = 0
      MONFNX = 0
      MONDLX = 0
      IERKGALLOC = 1
      NTDLT  = 0
      MNLPLI = 0
      MNLPCO = 0
      MNLPLC = 0
      MNLPCC = 0
      MNKGC  = 0
      DCPU   = 0D0
      DFABG  = 0D0
      DFACT  = 0D0
      DDEPL  = 0D0
C
C     VERIFICATION DU NOM_DE_L'OBJET
C     ==============================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: DEFINITION INCONNUE de l''OBJET ' // KNOMOB
         ELSE
            KERR(1)='ERROR: UNKNOWN DEFINITION of the OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     L'ANCIEN HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
C     PROBLEME A COEFFICIENTS INDEPENDANTS OU NON DU DEPLACEMENT?
C     ===========================================================
      CALL LIMTCL( 'resoelst', TESTNL )
      IF( TESTNL .LE.  0 ) RETURN
      TESTNL = TESTNL - 1
C
C     L'OBJET EST-IL STRUCTURE CLASSIQUE, SOUS-DOMAINES OU JOINTS?
C     ------------------------------------------------------------
      NDOUNO = MCN(MNDFOB+WDOUNO)
      NDOUNO = 0
C
      IF (NDOUNO .EQ. 1) THEN
C
C        RESOLUTION PAR LA METHODE DES SOUS-DOMAINES
C        ===========================================
         WRITE(IMPRIM,10013) KNOMOB
10013    FORMAT(/,1X,80('=')/,
     % ' ELASTICITE STATIONNAIRE PAR SOUS-DOMAINES DE L''OBJET ',A,/,
     %   1X,80('='))
C         CALL SDELAS( KNOMOB, IERR )
         DDEPL = DINFO( 'DELTA CPU' )
         RETURN
C
      ELSE IF (NDOUNO .EQ. 2) THEN
C
C        RESOLUTION PAR LA METHOOE DES JOINTS
C        ====================================
         WRITE(IMPRIM,10113) KNOMOB
10113    FORMAT(/,1X,80('=')/,
     % ' ELASTICITE STATIONNAIRE PAR JOINTS DE L''OBJET ',A,/,
     %   1X,80('='))
         WRITE(IMPRIM,10114)
10114    FORMAT(1X,' OPTION EN COURS DE PROGRAMMATION')
         NBLGRC(NRERR) = 2
         KERR(1) = 'RESOLUTION PAR LA METHODE DES JOINTS '
         KERR(2) = 'OPTION EN COURS DE PROGRAMMATION '
         CALL LEREUR
C         CALL JOELAS( KNOMOB, IERR )
         IERR = 1
         RETURN
C
      ENDIF
C
C     RESOLUTION CLASSIQUE
C     ====================
 10   CONTINUE
CCC      WRITE(IMPRIM,10009) KNOMOB
CCC10009 FORMAT(/,1X,80('*')/,
CCC     %' RESOLUTION ELASTICITE STATIONNAIRE CLASSIQUE DE L''OBJET ',A/
CCC     %  1X,80('*'))
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF"
C     ASSOCIES A L'OBJET
C     ========================================================
      CALL MIMAOB( 1,      NTLXOB, MXDOEL, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEF, MNDOEF, IERR )
      IF( IERR .NE. 0 ) RETURN
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES OBJETS
C         NUMAOB          LES 4 NUMEROS MAXIMA DES OBJETS
C     MNDOEF LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT L'ELASTICITE DE L'OBJET COMPLET
C
C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     ==========================================
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM LA DIMENSION 1 OU 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )
C     LE NOMBRE TOTAL DE DEGRES DE LIBERTE
      NTDL   = NDIM * NBNOEU
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) NDIM,NBNOEU,NTDL,NDSM
      ELSE
         WRITE(IMPRIM,20010) NDIM,NBNOEU,NTDL,NDSM
      ENDIF
10010 FORMAT(/' DIMENSION 1 ou 2 ou 3 de l''ESPACE',T37,'=',I9/
     %' NOMBRE de NOEUDS'                     ,T37,'=',I9/
     %' NOMBRE de DL des DEPLACEMENTS'        ,T37,'=',I9/
     %' NOMBRE de JEUX de DONNEES'            ,T37,'=',I9)
20010 FORMAT(/' SPACE DIMENSION (1 or 2 or 3)',T37,'=',I9/
     %' NUMBER of NODES'                      ,T37,'=',I9/
     %' NUMBER of DoF of DISPLACEMENTS'       ,T37,'=',I9/
     %' NUMBER of DATA CASES'                 ,T37,'=',I9)
C
C     CHOISIR THERMO-ELASTICITE OU ELASTICITE SEULE
C     =============================================
      CALL LIMTCL( 'therelas', NOTHEL )
C     NOTHEL = 1  POUR LA RESOLUTION DE L'ELASTICITE SEULE
C              2  POUR LA RESOLUTION DE LA THERMO-ELASTICITE
      IF( NOTHEL .LE. 0 ) RETURN
      IF( NOTHEL .EQ. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'THERMO-ELASTICITE DEMANDEE'
         ELSE
            WRITE(IMPRIM,*) 'THERMO-ELASTICITY REQUIRED'
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ELASTICITE SANS THERMIQUE'
         ELSE
            WRITE(IMPRIM,*) 'ONLY ELASTICITY WITHOUT HEAT TRANSFER'
         ENDIF
      ENDIF
C
C     LA TEMPERATURE A T ELLE ETE CALCULEE ?
      NTVETE = 0
      MNVETE = 0
      MNTEMP = 0
      IF( NOTHEL .EQ. 2 ) THEN
C
C        THERMO-ELASTICITE: RECHERCHE DES TEMPERATURES
         CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVETE, MNVETE )
         IF( NTVETE .LE. 0 ) THEN
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: THERMO-ELASTICITE DEMANDEE'
               KERR(2) = 'SANS VECTEUR DES TEMPERATURES'
               KERR(2) = 'EXECUTER THERMICER AUPARAVANT'
            ELSE
               KERR(1) = 'ERROR: THERMO-ELASTICITT REQUIRED'
               KERR(2) = 'WITHOUT the TEMPERATURE VECTOR'
               KERR(3) = 'EXECUTE THERMICER PREVIOUSLY'
            ENDIF
            CALL LEREUR
            IERR = 1
            RETURN
         ENDIF
C
C        ADRESSE DU DEBUT DES VECTEURS TEMPERATURE
         MNTEMP = MNVETE + WECTEU
         NTDLT  = MCN( MNVETE + WBCOVE )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10012) MCN( MNVETE + WBVECT ), NTDLT
            WRITE(IMPRIM,10014)
         ELSE
            WRITE(IMPRIM,20012) MCN( MNVETE + WBVECT ), NTDLT
            WRITE(IMPRIM,20014)
         ENDIF
10012    FORMAT(' NOMBRE DE CARTES DE TEMPERATURE',I8/
     %          ' NOMBRE DE NOEUDS DE TEMPERATURE',I8)
10014 FORMAT(' LA TEMPERATURE DEJA CALCULEE EST PRISE EN COMPTE')
20012    FORMAT(' NUMBER of CASES of TEMPERATURE',I8/
     %          ' NUMBER of NODES of TEMPERATURE',I8)
20014 FORMAT(' The PREVIOUS COMPUTED TEMPERATURES are USED')
      ENDIF
C
C     CHOIX DE LA METHODE DE RESOLUTION ET DU STOCKAGE DES MATRICES
C     =============================================================
      CALL LIMTCL( 'methreso', NORESO )
      IF( NORESO .LE. 0 ) GO TO 9999
      IF( NORESO .LE. 0 .OR. NORESO .GE. 3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='ERREUR ELASTA: METHODE de RESOLUTION NON PROGRAMMEE'
         ELSE
           KERR(1)='ERROR ELASTA: METHOD of COMPUTATION NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF
C
C     RETROUVER LES ADRESSES MCN DES DONNEES ELASTICITE
C     DES   SV "OBJETS INTERNES"    DE L'OBJET
C     DES PLS  "OBJETS AUX LIMITES" DE L'OBJET
C     =================================================
      CALL ELADON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL, MNDOEF,
     %             IEMASS, IEYOUN, IEDILA, IECOED, IECOIN, IEFOIN,
     %             IEFIXA, IEFOCL, IEFOPO,
     %             IEDEIN, IEVIIN, IERR )
      IF( IEYOUN .NE. NBOBIN ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            IF( NDIM .EQ. 3 ) THEN
               KERR(1) = 'Au MOINS 1 VOLUME SANS YOUNG'
            ELSE IF( NDIM .EQ. 2 ) THEN
               KERR(1) = 'Au MOINS 1 SURFACE SANS YOUNG'
            ELSE IF( NDIM .EQ. 1 ) THEN
               KERR(1) = 'Au MOINS 1 LIGNE SANS YOUNG'
            ENDIF
            KERR(2) = '=> PROBLEME SANS SOLUTION'
         ELSE
            IF( NDIM .EQ. 3 ) THEN
               KERR(1) = 'At LEAST 1 VOLUME WITHOUT YOUNG'
            ELSE IF( NDIM .EQ. 2 ) THEN
               KERR(1) = 'At LEAST 1 SURFACE WITHOUT  YOUNG'
            ELSE IF( NDIM .EQ. 1 ) THEN
               KERR(1) = 'At LEAST 1 LINE WITHOUT  YOUNG'
            ENDIF
            KERR(2) = '=> PROBLEM WITHOUT SOLUTION'
         ENDIF
         CALL LEREUR
         IERR = 2
      ENDIF
      IF( IEFIXA .LE. 0 ) THEN
          NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET SANS FIXATION'
            KERR(2) = '=> PROBLEME SANS SOLUTION'
         ELSE
            KERR(1) = 'OBJECT WITHOUT FIXATION'
            KERR(2) = '=> PROBLEM WITHOUT SOLUTION'
         ENDIF
         CALL LEREUR
         IERR = 2
      ENDIF
      IF( IERR .NE. 0 ) RETURN
C
C     RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA
C     CONSTRUCTION DES TABLEAUX ELEMENTAIRES
C     ===============================================
      CALL TAPOBA( NBTYEL, MNNPEF, EETAEL,
     %             MNTPOB, NBDLMX, MOAUX, NBTTEF, NOAXIS, NCODMG, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
C     ===================================================
      MNTAUX = 0
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C     LA MATRICE ELEMENTAIRE ET LES NDSM SECONDS MEMBRES
      MOTAEL = NBDLMX * (NBDLMX+1) / 2 + NBDLMX * NDSM
      MNTAEL = 0
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
      MNIP   = 0
      CALL TNMCDC( 'ENTIER', NBDLMX, MNIP   )
      MNELAS = 0
      CALL TNMCDC( 'REEL2',  21,      MNELAS )
      MNFORC = 0
C     6*NDSM A CAUSE DES CONTRAINTES INITIALES
      MOFORC = 6 * NDSM
      CALL TNMCDC( 'REEL2', MOFORC, MNFORC )
C     LE TABLEAU DES COORDONNEES DES NOEUDS D'UN ELEMENT FINI
      MNX    = 0
      CALL TNMCDC( 'REEL', NBDLMX*NDIM, MNX )
C
C     CONSTRUCTION DES TABLEAUX DES NUMEROS ET VALEURS DES FORCES NODALES
C     ===================================================================
      CALL ELFNFX( NTDL,   NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEF, RELMIN,
     %             NBRFNX, MONFNX, MNNFNX, MNVFNX,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     CONSTRUCTION DES TABLEAUX DES NO ET VALEURS DES DEPLACEMENTS NODAUX FIXES
C     =========================================================================
      CALL ELDLFX( NTDL,   NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEF, RELMIN,
     %             NBRDLX, MONDLX, MNNDLX, MNVDLX,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     PREPARATION DE LA RESOLUTION
C     ============================
C     LA MATRICE PROFIL EST ICI SYMETRIQUE NON DIAGONALE
      NCODKG = 1
      IF ( NORESO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10015)
         ELSE
            WRITE(IMPRIM,20015)
         ENDIF
10015    FORMAT(/' RESOLUTION par FACTORISATION de CHOLESKY avec STOCKAG
     %E PROFIL de [K]')
20015    FORMAT(/' RESOLUTION by CHOLESKY FACTORIZATION with SKYLINE STO
     %RAGE of [K]')
C
C        METHODE DE CHOLESKY :
C        =====================
C        CALCUL DU PROFIL DE LA MATRICE
C        (LA MATRICE PROFIL EST ICI SYMETRIQUE)
         MNLPLI = 0
         CALL TNMCDC( 'ENTIER', NTDL+1, MNLPLI )
         CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZN ,
     %                NDIM  , NCODKG ,
     %                MCN(MNLPLI), IERR )
         IF( IERR .GT. 0 ) GOTO 9999
C
C        NOMBRE DE COEFFICIENTS DE LA MATRICE KG
         NBRDKG = MCN( MNLPLI + NTDL )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10290) NBRDKG
         ELSE
            WRITE(IMPRIM,20290) NBRDKG
         ENDIF
10290    FORMAT(' MATRICE PROFIL DE',I15,' REELS DOUBLE PRECISION')
20290    FORMAT(' SKYLINE MATRIX of',I15,' DOUBLE PRECISION REAL')
C
      ELSE IF( NORESO .EQ. 2 ) THEN
C
C        METHODE DU GRADIENT CONJUGUE
C        ============================
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10016)
         ELSE
            WRITE(IMPRIM,20016)
         ENDIF
10016    FORMAT(/' RESOLUTION PAR GRADIENT CONJUGUE MATRICE CONDENSEE')
20016    FORMAT(/' SOLUTION by SPARSE MATRIX CONJUGATE GRADIENT')
C        CALCUL DU SQUELETTE DE LA MATRICE
C        (LA MATRICE EST ICI SYMETRIQUE)
         MNLPLI = 0
         MNLPCO = 0
         CALL PRGCMC( MNTOPO, MCN(MNNPEF), MNXYZN ,
     %                NDIM  , NCODKG ,
     %                MNLPLI, MNLPCO, IERR )
         IF( IERR .NE. 0 ) THEN
            IERR = 20
            GOTO 9999
         ENDIF
C
C        DECLARATION INITIALISATION DE LA MATRICE MORSE SYMETRIQUE
         NBRDKG = MCN(MNLPLI+NTDL)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10291) NBRDKG
         ELSE
            WRITE(IMPRIM,20291) NBRDKG
         ENDIF
10291    FORMAT(' MATRICE MORSE ',I15,' REELS DOUBLE PRECISION')
20291    FORMAT(' SPARSE MATRIX ',I15,' DOUBLE PRECISION REAL')
C
      ENDIF
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE KG
C     ---------------------------------------------------
      IF( NBRDKG .LE. 0 ) THEN
C        PLACE MEMOIRE INSUFFISANTE
         NBLGRC(NRERR) = 3
         WRITE(KERR(MXLGER-1)(1:25),'(G25.0)') NBRDKG
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: PLACE MEMOIRE INSUFFISANTE'
            KERR(2) = KERR(MXLGER-1)(1:25) //
     %                 ' MOTS NECESSAIRES pour KG'
            KERR(3) = 'REDUIRE le MAILLAGE'
         ELSE
            KERR(1) ='ERROR: NOT ENOUGH MEMORY'
            KERR(2) = KERR(MXLGER-1)(1:25) //
     %           ' NECESSARY WORDS to store KG'
            KERR(3) ='REDUCE the MESH'
         ENDIF
         CALL LEREUR
         GOTO 9999
      ENDIF
C
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
C
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9999
      ENDIF
C
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
      WRITE(IMPRIM,*)
C
C     ADRESSAGE INITIALISATION DU VECTEUR DEPLACEMENT INITIAL
C     =======================================================
      CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTVEDE, MNVEDE )
      IF( NTVEDE .GT. 0 ) THEN
C        LE VECTEUR EST DETRUIT POUR ETRE REDECLARE
         CALL LXTSDS( NTLXOB, 'VECTEUR"DEPLACT' )
      ENDIF
      L = ( WECTEU + NDSM * NTDL * MOREE2 ) / MOREE2
      CALL LXTNDC( NTLXOB, 'VECTEUR"DEPLACT', 'REEL2', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTVEDE, MNVEDE )
      MNU0 = MNVEDE + WECTEU
      MNB  = 0
      CALL TNMCDC( 'REEL2', NDSM*NTDL, MNB )
      CALL AZEROD( NDSM*NTDL, MCN(MNB) )
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES ET DU SYSTEME LINEAIRE
C     ==============================================================
      CALL ELAMKB( 0,      1,      1,
     %             PENALI, D2PI,   NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEF,
     %             MNELAS, MNFORC, MNTAEL, MNX,    MNIP,   MNNODL,
     %             MNTEMP, NTDLT,
     %             NORESO, MNLPLI, MNLPCO,
     %             0,   0, NBRDKG, KG,     NTDL, MNU0,
     %             NCODMG, NCODKG, NPIMAX )
C
C     ASSEMBLAGE DES FORCES NODALES
C     =============================
      CALL ASFONO(NTDL,NDSM,NBRFNX,MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %            MCN(MNU0))
CCC
CCCC     SAUVEGARDE SUR FICHIER NFCONT DES SECONDS MEMBRES ET DE LA MATRICE
CCCC     AVANT PRISE EN COMPTE DES CONDITIONS AUX LIMITES POUR LES REACTIONS
CCC      REWIND NFCONT
CCC      WRITE(NFCONT) NDSM,NTDL
CCC      KK = 2 * NDSM * NTDL
CCC      WRITE(NFCONT) (MCN(MNU0-1+I),I=1,KK)
CCC      KK = NTDL + 1
CCC      WRITE(NFCONT) (MCN(MNLPLI-1+I),I=1,KK)
CCC      KK = 2 * MCN(MNLPLI+NTDL)
CCC      WRITE(NFCONT) (MCN(MNKG-1+I),I=1,KK)
C
C     CALCUL DE LA RESULTANTE TOTALE DES FORCES AVANT C.L.
C     ====================================================
      MNFALL = 0
      CALL TNMCDC( 'REEL2', NDSM*NDIM, MNFALL )
      CALL RESFOR( NBNOEU,NDIM,NTDL,NDSM,MCN(MNU0),'FORCE',MCN(MNFALL) )
      CALL TNMCDS( 'REEL2', NDSM*NDIM, MNFALL )
C
C     TEMPS CALCUL DE FORMATION DU SYSTEME LINEAIRE
      DFABG = DINFO( 'DELTA CPU' )
C
      IF( NORESO .EQ. 1 ) THEN
C
C        PRISE EN COMPTE DES CONDITIONS AUX LIMITES. DL FIXES
C        ====================================================
         CALL BLDLPC(NTDL,NDSM,NBRDLX,MCN(MNNDLX),MCN(MNVDLX),NCODKG,
     %               MCN(MNLPLI),KG,MCN(MNU0))
C
C        FACTORISATION DE CHOLESKY
C        =========================
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FACTORISATION DE CHOLESKY. PATIENCE...'
         ELSE
            KERR(1) = 'CHOLESKY FACTORIZATION. WAIT PATIENTLY...'
         ENDIF
         CALL LERESU
         CALL CHOLPR( NTDL, NCODKG, MCN(MNLPLI), KG,
     %                KG, IERR )
C        TEMPS CALCUL DE FORMATION DU SYSTEME LINEAIRE
         DFACT = DINFO( 'DELTA CPU' )
         IF( IERR .NE. 0 ) THEN
C            MATRICE NON INVERSIBLE
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR MATRICE NON INVERSIBLE'
                KERR(2) = 'REVOYEZ LES CONDITIONS AUX LIMITES'
             ELSE
                KERR(1) = 'ERROR NOT INVERSIBLE MATRIX'
                KERR(2) = 'SEE AGAIN the BOUNDARY CONDITIONS'
             ENDIF
             CALL LEREUR
             IERR = 7
             GOTO 9999
         ENDIF
C
C        RESOLUTION DU SYSTEME FACTORISE
C        ================================
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'DESCENTE REMONTEE DU SYSTEME FACTORISE'
         ELSE
            KERR(1) = 'SOLUTION by L (Lt x) = b'
         ENDIF
         CALL LERESU
         CALL DRCH1D(NTDL,NCODKG,MCN(MNLPLI),KG,NDSM,MCN(MNU0),2,
     %               MCN(MNU0) )
C
CCCC        MISE A JOUR DU TMS 'PROFIL"RAIDEUR'
CCCC        ===================================
CCC         MCN( MNPROF + WBLIPR ) = NTDL
CCC         MCN( MNPROF + WPPROF ) = LPPROF
CCCC        LA DATE
CCC         CALL ECDATE( MCN(MNPROF) )
CCCC        LE NUMERO DU TABLEAU DESCRIPTEUR
CCC         MCN( MNPROF + MOREE2 ) = NONMTD( '~>>>PROFIL' )
C
      ELSE IF ( NORESO .EQ. 2 ) THEN
C
C        PRISE EN COMPTE DES CONDITIONS AUX LIMITES. DL FIXES
C        ====================================================
         CALL BLDLGC(NTDL,NDSM,NBRDLX,MCN(MNNDLX),MCN(MNVDLX),NCODKG,
     %               MCN(MNLPLI),MCN(MNLPCO),KG,MCN(MNU0))
C
C        COMPRESSION DE LA MATRICE  KG : SUPPRESSION DES ZEROS
C        =====================================================
         CALL COMORS( NTDL, MCN(MNLPLI), MCN(MNLPCO), KG )
C        IL NE FAUT PAS MODIFIER LPMORS !
C
CCCC        MISE A JOUR DU TMS 'MORSE"RAIDEUR'
CCCC        =================================
CCC         MCN( MNMORS + WBLIMO ) = NTDL
CCC         MCN( MNMORS + WPMORS ) = LPMORS
CCCC        LA DATE
CCC         CALL ECDATE( MCN(MNMORS) )
CCCC        LE NUMERO DU TABLEAU DESCRIPTEUR
CCC         MCN( MNMORS + MOREE2 ) = NONMTD( '~>>>MORSE' )
C
C        CALCUL DE LA MATRICE DE PRECONDITIONNEMENT
C        ==========================================
C        LECTURE DU NIVEAU DE FACTORISATION INCOMPLETE
         NACTUA=0
11       CALL INVITE( 36 )
         NCVALS = 0
         NIVMAX = 10
         NIVEAU = 0
         CALL LIRENT( NCVALS, NIVEAU )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'FACTORISATION INCOMPLETE de NIVEAU',NIVEAU
         ELSE
            WRITE(IMPRIM,*) 'INCOMPLETE FACTORIZATION of LEVEL',NIVEAU
         ENDIF
C
C        CALCUL DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
         ISTAB = 0
 123     IERR  = 0
         CALL CALPNT( NTDL, NIVEAU, NCODKG, MNLPLI, MNLPCO,
     %                MNLPLC, MNLPCC, MNLPLU, COMP, IERR )
C        VERIFICATION DE LA FACTORISATION
         LOLPCC = MCN(MNLPLC+NTDL)
         IF( IERR .NE. 0 ) THEN
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PROBLEME lors de la FACTORISATION de NIVEAU '
     %              // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
             ELSE
               KERR(1)='PROBLEM of INCOMPLETE FACTORIZATION of LEVEL '
     %               // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
            ENDIF
            CALL LEREUR
            IF( MNLPLI .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
            IF( MNLPCO .GT. 0) CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCO )
            IF( MNLPLC .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
            IF( MNLPCC .GT. 0) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
            IF( MNKGC  .GT. 0) CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
            GO TO 9999
         ENDIF
C
C        DECLARATION INITIALISATION DE LA MATRICE MORSE DE CONDITIONNEMENT
         LOLPCC = MCN(MNLPLC+NTDL)
CCC
CCC         CALL LXTSOU( NTLXOB, 'MORSE"RAIDEUR_C', NTMORC, MNMORC )
CCC         IF( NTMORC .GT. 0 ) THEN
CCCC           LA MATRICE MORSE EST DETRUITE POUR ETRE REDECLAREE
CCC            CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR_C' )
CCC         ENDIF
CCC         LPMORC = WPLIGN + 1 + NTDL + LOLPCC
CCC         IF( MOD(LPMORC,MOREE2) .EQ. 1 ) LPMORC = LPMORC + 1
CCC         LO = LPMORC / MOREE2 + LOLPCC
CCC
C        TEST DE LA PLACE MEMOIRE DISPONIBLE
         CALL TNMCMX( 'REEL2', MAXVAR )
         LO = LOLPCC
         IF( MAXVAR .LT. LO ) THEN
C            PLACE MEMOIRE INSUFFISANTE
             NBLGRC(NRERR) = 2
             WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR PLACE MEMOIRE INSUFFISANTE'
                KERR(2) = 'DIMINUEZ la VALEUR du NIVEAU '
     %                 // KERR(MXLGER)(1:8)
             ELSE
                KERR(1) = 'ERROR NOT ENOUGH MEMORY'
                KERR(2) = 'REDUCE the LEVEL VALUE'
     %                 // KERR(MXLGER)(1:8)
             ENDIF
             CALL LEREUR
             IF( MNLPLC .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
             IF( MNLPCC .GT. 0) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
             IF( MNKGC  .GT. 0) CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
             IF( NACTUA .EQ. 0 ) THEN
                GOTO 11
             ELSE
                GOTO 9999
             ENDIF
         ENDIF
C
         MNKGC = 0
         CALL TNMCDC( 'REEL2', LOLPCC, MNKGC )
C
CCC         CALL LXTNDC( NTLXOB, 'MORSE"RAIDEUR_C', 'REEL2', LO )
CCC         CALL LXTSOU( NTLXOB, 'MORSE"RAIDEUR_C', NTMORC, MNMORC )
CCCC        COPIE DU TABLEAU LPLIGN
CCC         CALL TRTATA( MCN(MNLPLC), MCN(MNMORC+WPLIGN), NTDL+1 )
CCCC        COPIE DU TABLEAU LPCOLO
CCC         CALL TRTATA( MCN(MNLPCC), MCN(MNMORC+WPLIGN+NTDL+1), LOLPCC )
CCCC        DESTRUCTION DES TABLEAUX TEMPORAIRES
CCC         CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
CCC         CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
CCCC        MISE A JOUR DES ADRESSES
CCC         MNLPLC = MNMORC + WPLIGN
CCC         MNLPCC = MNLPLC + NTDL + 1
CCCC        INITIALISATION DE LA MATRICE
CCC         MNKGC = MNMORC + LPMORC
CCC         CALL AZEROD( LOLPCC, MCN(MNKGC) )
C
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FACTORISATION INCOMPLETE PAR NIVEAUX'
         ELSE
            KERR(1) = 'INCOMPLETE FACTORIZATION by LEVELS'
         ENDIF
         CALL LERESU
         CALL INCHGC( NTDL ,
     %                MCN(MNLPLI),MCN(MNLPCO),KG,
     %                MCN(MNLPLC),MCN(MNLPCC),MCN(MNKGC), IERR )
C        TEMPS CALCUL DE FORMATION DU SYSTEME LINEAIRE
         DFACT = DINFO( 'DELTA CPU' )
C
CCCC        MISE A JOUR DU TMS 'MORSE"RAIDEUR_C'
CCCC        ===================================
CCC         MCN( MNMORC + WBLIMO ) = NTDL
CCC         MCN( MNMORC + WPMORS ) = LPMORC
CCCC        LA DATE
CCC         CALL ECDATE( MCN(MNMORC) )
CCCC        LE NUMERO DU TABLEAU DESCRIPTEUR
CCC         MCN( MNMORC + MOREE2 ) = NONMTD( '~>>>MORSE' )
C
C        VERIFICATION DE LA STABILITE
         IF( IERR .EQ. 1 ) THEN
            ISTAB = 1
            IERR  = 0
         ELSE IF( IERR .EQ. 2 ) THEN
            IF ( COMP ) GO TO 9999
            IF ( NIVEAU .LT. NIVMAX ) THEN
               NIVEAU = NIVEAU + 1
               NACTUA=1
CCCC              DESTRUCTION DE LA STRUCTURE MORSE DE CONDITIONNEMENT
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR_C',NTMORC,MNMORC)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR_C' )
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
               CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
               CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
               IERR  = 0
               GO TO 123
            ELSE
C           LA FACTORISATION EST INSTABLE : ABANDON DES CALCULS
               WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     %                    // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
               ELSE
                  KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %               // KERR(MXLGER)(1:8)
                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               ENDIF
               CALL LEREUR
CCCC              DESTRUCTION DE LA STRUCTURE MORSE
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR',NTMORS,MNMORS)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR' )
CCCC              DESTRUCTION DE LA STRUCTURE MORSE DE CONDITIONNEMENT
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR_C',NTMORC,MNMORC)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR_C' )
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
               CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCO )
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
               CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
               CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
               IERR = 8
               GO TO 9999
            END IF
         END IF
C
C        RESOLUTION DU SYSTEME FACTORISE
C        ===============================
C        LES TABLEAUX AUXILIAIRES DE GCPRCH
         LOAUX = NTDL * NDSM * 4
         CALL TNMCDC( 'REEL2', LOAUX, MNAUX )
         MNAUX1 = MNAUX
         MNAUX2 = MNAUX1 + NTDL * NDSM * MOREE2
         MNAUX3 = MNAUX2 + NTDL * NDSM * MOREE2
         MNAUX4 = MNAUX3 + NTDL * NDSM * MOREE2
C        MISE A ZERO DU TABLEAU AUXILIAIRE
         CALL AZEROD( LOAUX, MCN(MNAUX) )
C        PRISE EN COMPTE DES CONDITIONS AUX LIMITES
         IF (NBRDLX.NE.0) CALL BLDLX0(NTDL,NDSM,NBRDLX,
     %       MCN(MNNDLX),MCN(MNVDLX),MCN(MNAUX1))
C
C        STABILISATION PAR RE-ORTHOGONALISATION
         IF (ISTAB.EQ.0) THEN
            NBDIR = 1
         ELSE
            NBDIR = 10
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,12011)
            ELSE
               WRITE(IMPRIM,22011)
            ENDIF
12011 FORMAT(/' RE-ORTHOGONALISATION DES DIRECTIONS')
22011 FORMAT(/' RE-ORTHOGONALISATION of DIRECTIONS')
         ENDIF
         LODIR = (NTDL+1) * NBDIR * 2
         MNBDIR = 0
         CALL TNMCDC( 'REEL2', LODIR, MNBDIR )
C        MISE A ZERO
         CALL AZEROD( LODIR, MCN(MNBDIR) )
         MNDIR  = MNBDIR
         MNADIR = MNBDIR + NTDL * NBDIR * MOREE2
         MNDAD  = MNADIR + NTDL * NBDIR * MOREE2
         MNBET  = MNDAD  + NBDIR * MOREE2
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'DEBUT DES ITERATIONS DE GRADIENT CONJUGUE'
         ELSE
            KERR(1) = 'START of CONJUGATE GRADIENT ITERATIONS'
         ENDIF
         CALL LERESU
         CALL GCPRCH( NTDL, NDSM, NBDIR ,
     %                MCN(MNLPLI),MCN(MNLPCO),KG,MCN(MNU0),
     %                MCN(MNLPLC),MCN(MNLPCC),MCN(MNKGC),
     %                MCN(MNAUX1),MCN(MNAUX2),MCN(MNAUX3),
     %                MCN(MNAUX4),MCN(MNDIR),MCN(MNADIR),
     %                MCN(MNDAD),MCN(MNBET),MCN(MNU0),IERR)
         CALL TNMCDS( 'REEL2', LOAUX, MNAUX )
         CALL TNMCDS( 'REEL2', LODIR, MNBDIR )
C        VERIFICATION DE LA CONVERGENCE
         IF( IERR .NE. 0 ) THEN
            IF ( NIVEAU .LT. NIVMAX ) THEN
               NIVEAU = NIVEAU + 1
CCCC              DESTRUCTION DE LA STRUCTURE MORSE DE CONDITIONNEMENT
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR_C',NTMORC,MNMORC)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR_C' )
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
               CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
               CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
               GO TO 123
            ELSE
C           LA METHODE NE CONVERGE PAS : ABANDON DES CALCULS
               WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'NON CONVERGENCE AU NIVEAU '
     %                   // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
               ELSE
                  KERR(1) = 'NO CONVERGENCE at LEVEL '
     %                   // KERR(MXLGER)(1:8)
                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               ENDIF
C
CCCC              DESTRUCTION DE LA STRUCTURE MORSE
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR',NTMORS,MNMORS)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR' )
CCCC              DESTRUCTION DE LA STRUCTURE MORSE DE CONDITIONNEMENT
CCC               CALL LXTSOU( NTLXOB,'MORSE"RAIDEUR_C',NTMORC,MNMORC)
CCC               CALL LXTSDS( NTLXOB, 'MORSE"RAIDEUR_C' )
CCC
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
               CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCO )
               CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
               CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
               CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
               IERR = 9
               GOTO 9999
            END IF
         END IF
C
      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"DEPLACT'
C     ====================================
      MCN( MNVEDE + WBCOVE ) = NTDL
      MCN( MNVEDE + WBVECT ) = NDSM
      MCN( MNVEDE + WBCPIN ) = 0
C     LA DATE
      CALL ECDATE( MCN(MNVEDE) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVEDE + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
C
C     COUT CALCUL DES DEPLACEMENTS
C     ============================
      DDEPL = DINFO( 'DELTA CPU' )
C
C     AFFICHAGE DU DERNIER VECTEUR"DEPLACT CALCULE
C     ============================================
      CALL AFDEPL( 10,     NDSM,   MNXYZN,
     %             NDIM,   NBNOEU, NDSM,   MCN(MNU0),
     %             DECMAX, NOFOTI, DEXMAX )
C
C
CCCC     RESTAURATION   FICHIER NFCONT DES SECONDS MEMBRES ET DE LA MATRICE
CCCC     AVANT PRISE EN COMPTE DES CONDITIONS AUX LIMITES POUR LE
CCCC     CALCUL DES REACTIONS
CCCC     ==================================================================
CCC      REWIND NFCONT
CCC      READ(NFCONT) NDSM,NTDL
CCC      KK = 2 * NDSM * NTDL
CCC      READ(NFCONT) (MCN(MNB-1+I),I=1,KK)
CCC      KK = NTDL + 1
CCC      READ(NFCONT) (MCN(MNLPLI-1+I),I=1,KK)
CCC      KK = 2 * MCN(MNLPLI+NTDL)
CCC      READ(NFCONT) (MCN(MNKG-1+I),I=1,KK)
C
CCCC     B - A * U0 => B
CCCC     ===============
CCC      CALL PRBMAX(1,NDSM,NTDL,MCN(MNLPLI),KG,MCN(MNU0),MCN(MNB))
CCC      IF(NFDEPL.GT.0) WRITE(NFDEPL) LB4,(MCN(MNB-1+I),I=1,LB4)
C
CCCC     IMPRESSION DES REACTIONS
CCCC     ========================
CCC      WRITE(IMPRIM,10540)
CCC10540 FORMAT(/'LES REACTIONS'/1X,79(1H=))
CCC      MN = ( MNB - 1 ) / 2
CCC           DO 545 I=1,NBNOEU
CCC                DO 542 J=1,NDIM
CC            ATTENTION B(NTDL,NDSM) MAINTENANT !
CCC                WRITE(IMPRIM,10542) I,J,(K,DMCN(MN+K),K=1,NDSM)
CCC                MN = MN + NDSM
CCC 542            CONTINUE
CCC 545       CONTINUE
CCC10542 FORMAT(' NOEUD',I5,' D.L.',I2,
CCC     %(T25,4(' CAS',I2,' B-AX=',G12.4)))
C
C .....................................................................
C
C     FERMETURE DE LA MATRICE POUR REDONNER DE LA PLACE EN MC
C     =======================================================
      IF( MNLPLI .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
      IF( NORESO .EQ. 2 ) THEN
          IF( MNLPCO .GT. 0) CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCO )
          IF( MNLPLC .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
          IF( MNLPCC .GT. 0) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
          IF( MNKGC  .GT. 0) CALL TNMCDS( 'REEL2',  LOLPCC, MNKGC  )
      ENDIF
C
C **************************************************************************
C //////////////////////////////////////////////////////////////////////////
C **************************************************************************
C
C **************************************************************************
C --------------------------------------------------------------------------
C **************************************************************************
C
C     CALCUL DES CONTRAINTES EN CHAQUE POINT D INTEGRATION DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     =======================================================
      IF( IECOIN .GT. 0 ) THEN
         IECOIN = 1
      ELSE
         IECOIN = 0
      ENDIF
      CALL ELASTR( NTLXOB, MNTOPO, NOAXIS, NDIM,   MOREE2,
     %             NPIMAX, NDSM,   IECOIN, NTDL,
     %             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %             MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEF,
     %             MNELAS, MOTAEL, MNTAEL, MNX,
     %             MNVEDE, NOTHEL, MNVETE, DEXMAX, DECMAX )
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
 9999 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER' ,2*MXTYEL, MNNPEF )
      DO 19000 I=1,4
         IF( MNDOEF(I) .GT. 0 .AND. MXDOEF(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEF(I), MNDOEF(I) )
         ENDIF
19000 CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA ,MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX,  MNTAUX )
      IF( MNAUXE .GT. 0 ) CALL TNMCDS( 'REEL2',  9*NDSM, MNAUXE )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL, MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNNODL )
      IF( MNIP   .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNIP )
      IF( MNELAS .GT. 0 ) CALL TNMCDS( 'REEL2',  21   ,  MNELAS )
      IF( MNFORC .GT. 0 ) CALL TNMCDS( 'REEL2',  MOFORC, MNFORC )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDLMX*NDIM, MNX )
      IF( MNB    .GT. 0 ) CALL TNMCDS( 'REEL2',  NDSM*NTDL,   MNB )
      IF( MNNFNX .GT. 0 .AND. MONFNX .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', MONFNX, MNNFNX )
      ENDIF
      IF( MNVFNX .GT. 0.AND. MONFNX .GT. 0 ) THEN
         CALL TNMCDS( 'REEL2', MONFNX, MNVFNX )
      ENDIF
      IF( MNNDLX .GT. 0 .AND. MONDLX .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
      ENDIF
      IF( MNVDLX .GT. 0 .AND. MONDLX .GT. 0 ) THEN
         CALL TNMCDS( 'REEL2', MONDLX, MNVDLX )
      ENDIF
      IF( MNAUX .GT. 0 )  CALL TNMCDS( 'REEL2', LOAUX, MNAUX )
C
C     GESTION DES ERREURS
C     ===================
      IF( IERR .EQ. 7 ) THEN
C     RETOUR SI MATRICE NON INVERSIBLE
         IERR = 0
         RETURN
      ELSE IF( IERR .EQ. 8 ) THEN
C     RETOUR SI FACTORISAION INCOMPLETE INSTABLE ET TRAVAIL INTERACTIF
         IERR = 0
         GOTO 10
      ELSE IF( IERR .EQ. 9 ) THEN
C     RETOUR SI NON CONVERGENCE DU G.C. ET TRAVAIL INTERACTIF
         IERR = 0
         GOTO 10
      ELSE IF( IERR .NE. 0 ) THEN
         RETURN
      ENDIF
C
C     AFFICHAGE DES TEMPS CALCUL
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) DFABG,DFACT,DDEPL,DCPU,
     %                       DFABG+DFACT+DDEPL+DCPU
      ELSE
         WRITE(IMPRIM,22001) DFABG,DFACT,DDEPL,DCPU,
     %                       DFABG+DFACT+DDEPL+DCPU
      ENDIF
12001 FORMAT(/
     % 'TEMPS FORMATION     DE LA MATRICE =',F12.2,' SECONDES CPU'/,
     % 'TEMPS FACTORISATION DE LA MATRICE =',F12.2,' SECONDES CPU'/,
     % 'TEMPS CALCUL DES DEPLACEMENTS     =',F12.2,' SECONDES CPU'/
     % 'TEMPS CALCUL DES CONTRAINTES      =',F12.2,' SECONDES CPU'/
     % 'TEMPS CALCUL TOTAL                =',F12.2,' SECONDES CPU'/)
22001 FORMAT(/
     % 'MATRIX FORMATION         TIME =',F12.2,' CPU SECONDS'/,
     % 'MATRIX FACTORIZATION     TIME =',F12.2,' CPU SECONDS'/,
     % 'DISPLACEMENT COMPUTATION TIME =',F12.2,' CPU SECONDS'/,
     % 'STRESS       COMPUTATION TIME =',F12.2,' CPU SECONDS'/,
     % 'SOLUTION     TOTAL       TIME =',F12.2,' CPU SECONDS'/)
      RETURN
      END
