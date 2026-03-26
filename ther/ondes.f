      SUBROUTINE ONDES( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE DEPLACEMENT VERTICAL DANS UN DOMAINE 2D ou 3D
C -----    SOLUTION DE L'EQUATION DES ONDES
C          INTEGREE EN TEMPS SELON
C          Soit, le SCHEMA de NEWMARK
C          Soit, un SCHEMA DECENTRE IMPLICITE
C          (Soit, un SCHEMA EXPLICITE avec condition DT/h<Cte
C           aout 2006: PENALI => DIVERGENCE meme pour petit DT
C           => suppression de td/m/scheonde )
C          A partir de U0 et U1=U0 + DT * V0 donnes
C          et avec des COEFFICIENTS Rot C, [A], g INDEPENDANTS DU TEMPS
C          ET DU DEPLACEMENT. Les SOURCES PEUVENT DEPENDRE DU TEMPS
C
C ATTENTION:
C          Une condition de DIRICHLET U=UD est penalisee sous forme FOURIER
C          avec un coefficient d'ECHANGE g=1/Epsilon et F=UD/Epsilon
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1999
C23456---------------------------------------------------------------012
      DOUBLE PRECISION   PENALI
      PARAMETER         (PENALI=1D20)
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
      include"./incl/a___deplactinit.inc"
      include"./incl/a___vitesseinit.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE       (MCN(1), RMCN(1), DMCN(1))
C
      EXTERNAL          ETTAEL
      DOUBLE PRECISION  RELMIN, D2PI,  DT0,   DT, DTDT
      DOUBLE PRECISION  XPOIN,  YPOIN, ZPOIN
      DOUBLE PRECISION  DINFO,  DCPU,  DPREP, DDEPL
      DOUBLE PRECISION  DECMAX, DEXMAX
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      CHARACTER*(*)     KNOMOB
      LOGICAL           COMP
C     TABLEAU NON UTILISE (VITESSE D'UN FLUIDE INEXISTANT ICI)
      DOUBLE PRECISION  VITEGt(5,3)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: MG
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERMGALLOC, IERKGALLOC
      INTRINSIC         ALLOCATED
C
      DOUBLE PRECISION  DEPLA0, VITES0(1), D, DMAT
ccc       DOUBLE PRECISION comxbg
      CHARACTER*80      NOMTS, KNOMTD
C
      DATA              RELMIN/-1D28/
C
C     QUELQUES INITIALISATIONS DE VARIABLES ET ADRESSES MCN
C     =====================================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('*')/
     %' CALCUL du DEPLACEMENT de l''ONDE sur l''OBJET ',A/1X,80('*'))
20000 FORMAT(/1X,80('*')/
     %' COMPUTATION of the WAVE of the OBJECT ',A/1X,80('*'))
C
C     PROBLEME LINEAIRE
      TESTNL = 0
C     A PRIORI LE TEMPS INITIAL DE CALCUL
      TEMPS  = 0.0
      MNTIME = 0
C     LES TEMPS CALCUL
      DCPU   = DINFO( 'DELTA CPU' )
      DPREP  = 0D0
      DDEPL  = 0D0
      IERR   = 0
C     LA CONSTANTE 2 Pi
      D2PI   = ATAN( 1D0 ) * 8D0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      NBDLMXEF = 0
      NBCOOR = 0
      MNDIR  = 0
      MNDAD  = 0
      MNBET  = 0
      MNAUX2 = 0
      MNAUX3 = 0
      MNAUX4 = 0
      MNADIR = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNTHER = 0
      MNNODL = 0
      MNX    = 0
C
      NBFNFX = 0
      MONFNX = 0
      MNNFNX = 0
      MOVFNX = 0
      MNVFNX = 0
      MNNDLX = 0
      MNVDLX = 0
C
      MNAUX  = 0
      MNFLPT = 0
      MNFLTO = 0
      MNERTH = 0
      MNTAUX = 0
      MNNPEF = 0
C
      MNLPLK = 0
      MNLPCK = 0
      MNAUGC = 0
      MNBDIR = 0
      MNLPLC = 0
      MNLPCC = 0
      MNAGC  = 0
      MOTMAT = 0
C
      MNSG   = 0
      MNUG   = 0
      MNUG0  = 0
      MNUG1  = 0
      MNUG2  = 0
      MNFG   = 0
      MNFG0  = 0
      MNFG1  = 0
      MNFG2  = 0
C
      MNTHET = 0
      MNTHDL = 0
      DO 5 I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
 5    CONTINUE
      TEMPEL = 0D0
      IERMGALLOC = 1
      IERKGALLOC = 1
C
C     AFFICHAGE ET VERIFICATION DE L'OBJET
C     ====================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'OBJET INCONNU'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'UNKNOWN OBJECT'
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'TMS DEFINITION INCONNU'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'UNKNOWN DEFINITION TMS'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     ENTREE DU TEMPS FINAL: TPSINI
C     =============================
      CALL INVITE( 95 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, TPSINI )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        IERR = 2
        RETURN
      ENDIF
      TEMPS = TPSINI
C
C     ENTREE DU DEPLACEMENT INITIAL U0 SUR L'OBJET:
C     LECTURE DU TMS DEPLACTINIT
C     =============================================
      L1     = NUDCNB(KNOMOB)
      NOMTS  = '~>OBJET>' // KNOMOB(1:L1) // '>DEPLACTINIT'
      KNOMTD = '~>>>DEPLACTINIT'
      L1 = NUDCNB(NOMTS)
      CALL MOTSTD( KNOMTD, NOMTS(1:L1), IERR )
      IF( IERR .NE. 0 ) THEN
        NBLGRC(NRERR) = 1
        IF( LANGAG .EQ. 0 ) THEN
           KERR(1) = 'ERREUR: DEPLACEMENT INITIAL INCONNU'
        ELSE
           KERR(1) = 'ERROR: INITIAL DISPLACEMENT UNKNOWN'
        ENDIF
        CALL LEREUR
        RETURN
      ENDIF
C     AFFICHAGE DU TMS
      CALL AFTSTD( NOMTS )
C
C     CODE DE CALCUL DU DEPLACEMENT INITIAL U0
      CALL LXTSOU( NTLXOB,'DEPLACTINIT', NTDEIN, MNDEIN )
      LTDEP0 = MCN( MNDEIN + WTDEP0 )
C
      IF( LTDEP0 .NE. 3 ) THEN
C
C        ENTREE DE LA VITESSE INITIALE V0: LECTURE DU TMS VITESSEINIT
C        ============================================================
         L1     = NUDCNB(KNOMOB)
         NOMTS  = '~>OBJET>' // KNOMOB(1:L1) // '>VITESSEINIT'
         KNOMTD = '~>>>VITESSEINIT'
         L1 = NUDCNB(NOMTS)
         CALL MOTSTD( KNOMTD, NOMTS(1:L1), IERR )
         IF( IERR .NE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: VITESSE INITIALE INCONNUE'
            ELSE
               KERR(1) = 'ERROR: UNKNOWN INITIAL SPEED'
            ENDIF
            CALL LEREUR
            RETURN
         ENDIF
C        AFFICHAGE DU TMS
         CALL AFTSTD( NOMTS )
C
      ENDIF
C
C     ENTREE DU TEMPS FINAL: TPSFIN
C     =============================
      CALL INVITE( 94 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, TPSFIN )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        IERR = 2
        RETURN
      ENDIF
      IF( TPSFIN .LE. 0 ) THEN
        NBLGRC(NRERR) = 1
        IF( LANGAG .EQ. 0 ) THEN
           KERR(1) = 'ERREUR: TEMPS FINAL INCORRECT (<=0)'
        ELSE
           KERR(1) = 'ERROR: INCORRECT FINAL TIME (<=0)'
        ENDIF
        CALL LEREUR
        IERR = 1
        RETURN
      ENDIF
C
C     ENTREE DU PAS DE TEMPS CONSTANT: DT
C     ===================================
      CALL INVITE( 86 )
      NCVALS = 0
      CALL LIRRDP( NCVALS, DT )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        IERR = 2
        RETURN
      ENDIF
      IF( DT .LE. 0 ) THEN
        NBLGRC(NRERR) = 1
        IF( LANGAG .EQ. 0 ) THEN
           KERR(1) = 'ERREUR: PAS DE TEMPS <=0 INCORRECT'
        ELSE
           KERR(1) = 'ERROR: INCORRECT STEP of TIME (<=0)'
        ENDIF
        CALL LEREUR
        IERR = 1
        RETURN
      ENDIF
      DTDT = DT * DT
C
C     ENTREE DU NOMBRE MAXIMAL DE VECTEURS SOLUTIONS A STOCKER
C     ========================================================
      CALL INVITE( 75 )
      MXTEMP = 10
      NCVALS = 4
      CALL LIRENT( NCVALS, MXTEMP )
      IF( NCVALS .LT. 0 ) THEN
C        ABANDON DE LA LECTURE DES DONNEES
         IERR = 2
         RETURN
      ENDIF
      IF( MXTEMP .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'STOCKAGE DES DEPLACEMENTS A T=0 ET MAX'
         ELSE
            KERR(1) = 'STORAGE of DISPLACEMENTS at T=0 and MAX'
         ENDIF
         CALL LEREUR
         MXTEMP = 2
         RETURN
      ENDIF
C
C     CHOISIR LA METHODE DE RESOLUTION DU SYSTEME LINEAIRE
C     ====================================================
      CALL LIMTCL( 'methreso', NORESO )
      IF( NORESO .LE. 0 ) GO TO 9000
      IF( NORESO .GE. 3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'METHODE DE RESOLUTION NON PROGRAMMEE'
         ELSE
            KERR(1) = 'METHOD NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
      IF( NORESO .EQ. 2 ) THEN
C        LECTURE DU NIVEAU DE FACTORISATION INCOMPLETE
         CALL INVITE( 36 )
         NCVALS = 0
         NIVMAX = 10
         NIVEAU = 0
         CALL LIRENT( NCVALS, NIVEAU )
      ENDIF
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF"
C     ASSOCIES A L'OBJET
C     ========================================================
      CALL MIMAOB( 1, NTLXOB, MXDOTH, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES OBJETS
C         NUMAOB          LES 4 NUMEROS MAXIMA DES OBJETS
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     RETROUVER LES ADRESSES MCN DES DONNEES THERMIQUES
C     DES   SV "OBJETS INTERNES"    DE L'OBJET
C     DES PLS  "OBJETS AUX LIMITES" DE L'OBJET
C     =================================================
      CALL THEDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             1, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             IETEIN, IEVIIN, IEVIANT,IECOBO,
     %             IERR )
      IF( IEMAST .NE. NBOBIN ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1)='SURFACE ou VOLUME SANS MASSE'
             KERR(2)='=> PROBLEME SANS SOLUTION'
          ELSE
             KERR(1)='SURFACE or VOLUME WITHOUT MASS'
             KERR(2)='=> PROBLEM WITHOUT SOLUTION'
          ENDIF
          CALL LEREUR
          IERR = 11
      ENDIF
      IF( IECOND .NE. NBOBIN ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'SURFACE OU VOLUME SANS CONDUCTIVITE'
             KERR(2) = '=> PROBLEME SANS SOLUTION'
          ELSE
             KERR(1)='SURFACE or VOLUME WITHOUT CONDUCTIVITY'
             KERR(2)='=> PROBLEM WITHOUT SOLUTION'
          ENDIF
          CALL LEREUR
          IERR = 12
      ENDIF
      IF( IERR .NE. 0 ) GOTO 9000
C
C     INITIALISATIONS DE TABLEAUX ET AFFICHAGES
C     =========================================
C     LE NOMBRE TOTAL DE DEGRES DE LIBERTE EST EGAL ICI
C     AU NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET CAR ON A
C     UN DEGRE DE LIBERTE (DEPLACEMENT VERTICAL) PAR NOEUD
      NBNOEU = MCN( MNXYZN + WNBNOE )
      NTDL   = NBNOEU
C
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
      IF( NBCOOR .GT. 3 ) NDIM = NBCOOR
C
C     NOMBRE DE TEMPS OU LE CALCUL SE FAIT
      NBINST = NINT( (TPSFIN-TPSINI) / DT )
      MXTEMP = MIN( NBINST, MXTEMP )
C     LE PAS DE TEMPS POUR LES STOCKAGES DES VECTEURS DEPLACEMENTS
      DTSTOC = ( TPSFIN - TPSINI ) / MXTEMP
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) NDIM, NBNOEU, NTDL, TPSINI, TPSFIN, DT,
     %                       NBINST, MXTEMP
      ELSE
         WRITE(IMPRIM,20010) NDIM, NBNOEU, NTDL, TPSINI, TPSFIN, DT,
     %                       NBINST, MXTEMP
      ENDIF
10010 FORMAT(/' DIMENSION 2 ou 3 de l''ESPACE'  ,T41,'=',I6,/,
     %' NOMBRE de NOEUDS'                       ,T41,'=',I6,/,
     %' NOMBRE de DEGRES de LIBERTE'            ,T41,'=',I6,/,
     %' TEMPS INITIAL de PROPAGATION de L''ONDE',T41,'=',G14.6,/,
     %' TEMPS FINAL   de PROPAGATION de L''ONDE',T41,'=',G14.6,/,
     %' VALEUR du PAS de TEMPS CONSTANT'        ,T41,'=',G14.6/,
     %' NOMBRE de PAS de TEMPS DU CALCUL'       ,T41,'=',I6,/
     %' NOMBRE STOCKAGES du DEPLACEMENT EN PLUS',T41,'=',I6,/)
20010 FORMAT(/' SPACE DIMENSION 2 or 3'      ,T41,'=',I6,/,
     %' NUMBER of NODES'                     ,T41,'=',I6,/,
     %' NUMBER of DEGREES of FREEDOM (DoF)'  ,T41,'=',I6,/,
     %' INITIAL TIME of WAVE PROPAGATION'    ,T41,'=',G14.6,/,
     %' FINAL   TIME of WAVE PROPAGATION'    ,T41,'=',G14.6,/,
     %' CONSTANT VALUE of STEP of TIME'      ,T41,'=',G14.6/,
     %' NUMBER of STEPS of TIME'             ,T41,'=',I6,/
     %' NUMBER of STORAGES of DISPLACEMENTS' ,T41,'=',I6,/)
C
C     CHOISIR LE SCHEMA D'INTEGRATION EN TEMPS
C     ========================================
      CALL LIMTCL( 'scheonde', NOSCHM )
      IF( NOSCHM .LE. 0 ) GO TO 9000
      IF( LANGAG .EQ. 0 ) THEN
         IF( NOSCHM .EQ. 1 ) THEN
            WRITE(IMPRIM,*)'SCHEMA DE NEWMARK IMPLICITE EN TEMPS'
         ELSE IF( NOSCHM .EQ. 2 ) THEN
            WRITE(IMPRIM,*)'SCHEMA DECENTRE IMPLICITE EN TEMPS'
         ELSE
            WRITE(IMPRIM,*)'SCHEMA EXPLICITE CONDITIONNELLEMENT STABLE'
         ENDIF
      ELSE
         IF( NOSCHM .EQ. 1 ) THEN
            WRITE(IMPRIM,*)'NEWMARK SCHEME IMPLICIT in TIME'
         ELSE IF( NOSCHM .EQ. 2 ) THEN
            WRITE(IMPRIM,*)'DECENTRED SCHEME IMPLICIT in TIME'
         ELSE
            WRITE(IMPRIM,*)'EXPLICIT SCHEME CONDITIONALY STABLE'
         ENDIF
      ENDIF
C
C     RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA
C     CONSTRUCTION DES TABLEAUX ELEMENTAIRES
C     ===============================================
      CALL TAPOBA( NBTYEL, MNNPEF, ETTAEL,
     %             MNTPOB, NBDLMX, MOAUX, NBTTEF, NOAXIS, NCODSM, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C
C     DECLARATION DES 3 VECTEURS GLOBAUX DES DEPLACEMENTS EN 3 TEMPS SUCCESSIFS
C     =========================================================================
      CALL TNMCDC( 'REEL2', NTDL*3, MNUG )
      MNUG0 = MNUG
      MNUG1 = MNUG0 + MOREE2 * NTDL
      MNUG2 = MNUG1 + MOREE2 * NTDL
C
C     DECLARATION DES 3 VECTEURS GLOBAUX DES SECONDS MEMBRES EN 3 TEMPS SUCCESSI
C     ==========================================================================
      CALL TNMCDC( 'REEL2', NTDL*3, MNFG )
      MNFG0 = MNFG
      MNFG1 = MNFG0 + MOREE2 * NTDL
      MNFG2 = MNFG1 + MOREE2 * NTDL
C
C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
C     ===================================================
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C
C     MOTAEL = 2 MATRICES ELEMENTAIRES (DE CAPACITE + DE CONDUCTIVITE)
C              ET MXTEMP SECOND MEMBRES
      MOTAEL = NBDLMX * (NBDLMX+1) + NBDLMX * (MXTEMP+1)
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
C
C     LE TABLEAU DU NUMERO DES DEGRES DE LIBERTE D'UN ELEMENT FINI
      MNNODL = 0
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
C
C     CAPACITE OU CONDUCTIVITE AUX POINTS D'INTEGRATION NUMERIQUE DE L'EF
      CALL TNMCDC( 'REEL2', 128, MNTHER )
C
C     LE TABLEAU DES NBCOOR COORDONNEES DES NBDLMX NOEUDS DE L'EF MAXIMAL
      CALL TNMCDC( 'REEL', NBDLMX*NBCOOR, MNX )
C
C     DECLARATION POUR LE TERME DE TRANSPORT (VITESSEFLUIDE*GRADIENT DEPLACEMT)
C     DECLARATION DU VECTEUR DEPLACEMENT LOCAL SUR UN ELEMENT FINI
      CALL NLDATADC( NBDLMX )
C
C     ===========================================================
C     CALCUL DES MATRICES DE CAPACITE [M] ET CONDUCTIVITE [C]
C     SUPPOSEES INDEPENDANTES DU TEMPS ET DU DEPLACEMENT VERTICAL
C     ===========================================================
C
C     PREPARATION DE LA MATRICE PROFIL OU MORSE
C     =========================================
      NCODSK = 1
      IF ( NORESO .EQ. 1 ) THEN
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10015)
         ELSE
            WRITE(IMPRIM,20015)
         ENDIF
10015    FORMAT(' RESOLUTION DIRECTE PAR FACTORISATION DE CHOLESKY')
20015    FORMAT(' DIRECT SOLUTION by CHOLESKY FACTORIZATION')
C
C        STOCKAGE PROFIL ET FACTORISATION DE CHOLESKY
C        --------------------------------------------
C        CALCUL DU PROFIL DE LA MATRICE
C        (LA MATRICE PROFIL EST ICI SYMETRIQUE)
         CALL TNMCDC( 'ENTIER', NTDL+1     , MNLPLK )
         CALL PRPRMC( MNTOPO  , MCN(MNNPEF), MNXYZN,
     %                1       , NCODSK,
     %                MCN(MNLPLK), IERR )
         IF( IERR .NE. 0 ) GOTO 9000
C
C        DECLARATION INITIALISATION DE LA MATRICE PROFIL SYMETRIQUE
         NBRDKG = MCN( MNLPLK + NTDL )
         DMAT   = NBRDKG * 2D0
         IF( NBRDKG .LE. 0 ) THEN
C           PLACE MEMOIRE INSUFFISANTE
            NBLGRC(NRERR) = 3
            WRITE(KERR(MXLGER-1)(1:25),'(G25.0)') DMAT
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='ERREUR: PLACE MEMOIRE INSUFFISANTE'
               KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' MOTS NECESSAIRES pour [M] [K]'
               KERR(3) = 'REDUIRE le MAILLAGE'
            ELSE
               KERR(1) ='ERROR: NOT ENOUGH MEMORY'
               KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' NECESSARY WORDS to store [M] [K]'
               KERR(3) ='REDUCE the MESH'
            ENDIF
            CALL LEREUR
            IF( INTERA .LE. 1 ) CALL ARRET( 100 )
            GOTO 9000
         ENDIF
C
         MOTMAT = NTDL + MOREE2 * NBRDKG
         IF( NCODSM .NE. 0 ) THEN
            NBRDMG = NBRDKG
         ELSE
            NBRDMG = NTDL
         ENDIF
         MOTMAT = MOTMAT + MOREE2 * NBRDMG
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10016) MOTMAT
         ELSE
            WRITE(IMPRIM,20016) MOTMAT
         ENDIF
10016    FORMAT(' STOCKAGE DES MATRICES [M] [K] PROFIL =',I15,' MOTS MEM
     %OIRE' )
20016    FORMAT(' STORAGE of SKYLINE MATRICES [M] [K] =',I15,' MEMORY WO
     %RDS' )
C
      ELSE IF( NORESO .EQ. 2 ) THEN
C
C        POINTEURS DU STOCKAGE MORSE POUR METHODE DU GRADIENT CONJUGUE
C        -------------------------------------------------------------
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10018)
         ELSE
            WRITE(IMPRIM,20018)
         ENDIF
10018    FORMAT(' RESOLUTION ITERATIVE PAR GRADIENT CONJUGUE')
20018    FORMAT(' ITERATIVE SOLUTION by CONJUGATE GRADIENT')
C
C        CALCUL DU SQUELETTE DE LA MATRICE
C        (LA MATRICE EST ICI SYMETRIQUE)
         CALL PRGCMC( MNTOPO, MCN(MNNPEF), MNXYZN,
     %                1     , NCODSK,
     %                MNLPLK, MNLPCK, IERR )
         IF( IERR .NE. 0 ) GOTO 9000
C
C        DECLARATION INITIALISATION DE LA MATRICE MORSE SYMETRIQUE
         NBRDKG = MCN( MNLPLK + NTDL )
         MOTMAT = NTDL + (1+MOREE2) * NBRDKG
         IF( NCODSM .NE. 0 ) THEN
            NBRDMG = NBRDKG
         ELSE
            NBRDMG = NTDL
         ENDIF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10020) MOTMAT
         ELSE
            WRITE(IMPRIM,20020) MOTMAT
         ENDIF
10020    FORMAT(' STOCKAGE DES MATRICES [M] [K] MORSE=',I15,' MOTS MEMOI
     %RE' )
20020    FORMAT(' STORAGE of SPARSE MATRICES [M] [K]=',I15,' MEMORY WORD
     %S' )
      ENDIF
C
C     DECLARATION DE LA MATRICE DE MASSE
      IF( NCODSM .NE. 0 .OR. (NCODSM.EQ.0 .AND. NOSCHM.EQ.1) ) THEN
         NBRDMG = NBRDKG
         NCODSM = NCODSK
      ELSE
         NBRDMG = NTDL
      ENDIF
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DES MATRICES PROFIL MG et KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      ALLOCATE ( MG(1:NBRDMG), STAT=IERMGALLOC )
      IF( IERMGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDMG,
     %                ' DOUBLE PRECISION of the [MG] MATRIX'
         IERR = IERMGALLOC
         GOTO 9000
      ENDIF
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9000
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      WRITE(IMPRIM,*)
C
C     CONSTRUCTION DES MATRICES PROFIL [M] ET [K]
C     ===========================================
      IF( LTDEP0 .EQ. 3 ) TEMPS = REAL( TEMPS - DT )
      IEMG = 1
      IEKG = 1
      IEBG = 0
      CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %             D2PI,   NDIM,   NTDL,   VITEGt,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNTHER, MNTAEL, MNX,    MNNODL,
     %             NORESO, MNLPLK, MNLPCK,
     %             NBRDMG, MG,     NBRDKG, KG,     MNFG0,
     %             NCODSM, NCODSK, NBPTAF, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C
C     ICI 2 CAS POSSIBLES
C     NCODSM<>0 ALORS MEME PROFIL OU MORSE POUR MG ET KG
C            =0 ALORS DIAGONALE POUR MG ET PROFIL OU MORSE POUR KG
C
C     SI    LA MATRICE DE MASSE EST DIAGONALE ET METHODE DE NEWMARK
C     ALORS ELLE EST RENDUE SYMETRIQUE PROFIL
      IF( NCODSM .EQ. 0 .AND. NOSCHM .EQ. 1 ) THEN
         DO 26 I=NTDL,1,-1
C           SAUVEGARDE DU COEFFICIENT DIAGONAL AVANT TRANSFERT
            D = MG( I )
C           LE COEFFICIENT EST REMIS A ZERO
            MG( I ) = 0D0
C           LE COEFFICIENT DIAGONAL PROFIL TROUVE SA VALEUR
            MG( MCN(MNLPLK+I) ) = D
 26      CONTINUE
C        LA MATRICE MG N'EST PLUS DIAGONALE
         NCODSM = 1
         NBRDMG = NBRDKG
      ENDIF
C
C     La matrice de CAPACITE CALORIFIQUE [M] et de CONDUCTIVITE [K]
C     SONT ICI CENSEES ETRE INDEPENDANTES DU TEMPS
C     => ELLES SONT CALCULEES UNE SEULE FOIS
C
C     *******************************************************************
C     NOSCHM = 1 : NEWMARK  IMPLICIT inconditionnellement STABLE
C     -----------  ----------------------------------------------------
C     M ( Un+1 - 2 Un + Un-1 ) / DT**2 + K ( Un+1 + 2 Un + Un-1 ) / 4 =
C                                          ( Fn+1 + 2 Fn + Fn-1 ) / 4
C
C     U0 donn'e
C     U1 = U0 + DT * V0 (donn'e)
C
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [ M + DT**2/4 * K ] ( Un+1 - 2Un + Un-1 ) =
C     - DT**2 [K] Un + DT**2/4 ( Fn+1 + 2Fn + Fn-1 )
C     soit encore dans ce sous-programme
C     [AG](Un+1- 2Un + Un-1) = -DT**2 [BG]Un + DT**2/4 (Fn+1 + 2Fn + Fn-1)
C     avec {U0} = Un-1    {U1} = Un   {U2} = Un+1
C     et   {F0} = Fn-1    {F1} = Fn   {F2} = Fn+1
C     ******************************************************************
C     NOSCHM = 2 : SCHEMA DECENTRE IMPLICIT inconditionnellement STABLE
C     -----------  ----------------------------------------------------
C     M ( Un+1 - 2 Un + Un-1 ) / DT**2 + K Un+1 = Fn+1
C
C     U0 donn'e
C     U1 = U0 + DT * V0 (donn'e)
C
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [ M + DT**2 * K ] Un+1 =  [M] ( 2 Un - Un-1 ) + DT**2 * Fn+1
C     soit encore dans ce sous-programme
C     [AG] Un+1 = [BG] ( 2 Un - Un-1 ) + DT**2 * Fn+1
C     avec {U0} = Un-1    {U1} = Un   {U2} = Un+1  {F2}=Fn+1
C     ******************************************************************
C     NOSCHM = 3 : EXPLICITE=>CONDITION CFL impose un petit pas de temps
C     -----------  -----------------------------------------------------
C     M ( Un+1 - 2 Un + Un-1 ) / DT**2 + K Un = Fn
C
C     U0 donn'e
C     U1 = U0 + DT * V0 (donn'e)
C
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [ M ] (Un+1 + Un-1) =  [ 2M - DT**2 K ] Un + DT**2 * Fn
C     puis Un+1 = (Un+1 + Un-1) - Un-1
C     avec {U0} = Un-1     {U1} = Un    {U2} = Un+1
C     et   {F0} = Fn-1     {F1} = Fn    {F2} = Fn+1
C     ******************************************************************
C
      IF( NOSCHM .EQ. 1 ) THEN
C
C        NOSCHM = 1 : NEWMARK  IMPLICIT inconditionnellement STABLE
C        MG <= CAPACITE + DT**2/4 * CONDUCTIVITE
C        ----------------------------------------------------------
         NCODSA = NCODSK
         MNLPLA = MNLPLK
         MNLPCA = MNLPCK
C
         NCODSB = NCODSK
         MNLPLB = MNLPLK
         MNLPCB = MNLPCK
         CALL MUA2PD( NTDL,      1D0, NCODSM, MCN(MNLPLK), MG,
     %                      DTDT/4D0, NCODSK, MCN(MNLPLK), KG,
     %                                        MCN(MNLPLA), MG )
C
      ELSE IF( NOSCHM .EQ. 2 ) THEN
C
C        NOSCHM = 2 : DECENTRE IMPLICIT inconditionnellement STABLE
C        MG <= CAPACITE + DT**2 * CONDUCTIVITE
C        ----------------------------------------------------------
         NCODSA = NCODSK
         MNLPLA = MNLPLK
         MNLPCA = MNLPCK
C
         NCODSB = NCODSM
         MNLPLB = MNLPLK
         MNLPCB = MNLPCK
         CALL MUA2PD( NTDL,  1D0, NCODSM, MCN(MNLPLK), MG,
     %                      DTDT, NCODSK, MCN(MNLPLK), KG,
     %                                    MCN(MNLPLA), MG )
C
      ELSE
C
C        NOSCHM = 3 : EXPLICITE => CONDITION CFL pas de temps petit
C        KG <= 2 CAPACITE - DT**2 * CONDUCTIVITE
C        ----------------------------------------------------------
         NCODSA = NCODSM
         MNLPLA = MNLPLK
         MNLPCA = MNLPCK
C
         NCODSB = NCODSK
         MNLPLB = MNLPLK
         MNLPCB = MNLPCK
         CALL MUA2PD( NTDL,  2D0, NCODSM, MCN(MNLPLK), MG,
     %                     -DTDT, NCODSK, MCN(MNLPLK), KG,
     %                                    MCN(MNLPLB), KG )
C
      ENDIF
C
      IF( NCODSA .NE. 0 ) THEN
C
C        MATRICE AG NON DIAGONALE
         IF( NORESO .EQ. 1 ) THEN
C
C           FACTORISATION DE CHOLESKY DE LA MATRICE PROFIL MG
C           =================================================
            CALL CHOLPR( NTDL, NCODSA, MCN(MNLPLA), MG,   MG, IERR )
            IF( IERR .NE. 0 ) THEN
C              MATRICE NON INVERSIBLE
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'MATRICE MG NON INVERSIBLE'
                  KERR(2) = 'REVOYEZ LES CONDITIONS AUX LIMITES'
               ELSE
                  KERR(1) = 'NOT INVERSIBLE MATRIX MG'
                  KERR(2) = 'MODIFY BOUNDARY CONDITIONS'
               ENDIF
               CALL LEREUR
               IERR = 7
               GOTO 9000
            ENDIF
C
         ELSE
C
C           GC => MATRICE MORSE CONSTRUCTION DE LA MATRICE DE PRECONDITIONNEMENT
C           ====================================================================
C           CALCUL DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
            ISTAB  = 0
C           CALCUL DES POINTEURS DE LA MATRICE DE PRECONDITIONNEMENT
C           (FACTORISATION INCOMPLETE DE A SUIVANT LE NIVEAU)
 30         IERR   = 0
            MNLPLC = 0
            MNLPCC = 0
            MNLPLU = 0
            CALL CALPNT( NTDL,   NIVEAU, NCODSA, MNLPLA, MNLPCA,
     %                   MNLPLC, MNLPCC, MNLPLU, COMP,   IERR  )
C           VERIFICATION DU SQUELETTE DE LA MATRICE DE PRECONDITIONNEMENT
            IF( IERR .NE. 0 ) THEN
               WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
               NBLGRC(NRERR) = 3
               IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='PROBLEME DE FACTORISATION INCOMPLETE DE NIVEAU '
     %               // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
                  KERR(3) = 'AUGMENTER NIVEAU OU CHOISIR CHOLESKY'
               ELSE
                 KERR(1)='PROBLEM of INCOMPLETE FACTORIZATION of LEVEL '
     %               // KERR(MXLGER)(1:8)
                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
                  KERR(3) = 'AUGMENT LEVEL or CHOOSE CHOLESKY'
               ENDIF
               CALL LEREUR
               GOTO 9000
            ENDIF
C
C           DECLARATION INITIALISATION DE LA MATRICE MORSE DE CONDITIONNEMENT
            LOLPCC = MCN(MNLPLC+NTDL)
            MNAGC  = 0
            CALL TNMCDC( 'REEL2', LOLPCC, MNAGC )
            CALL AZEROD( LOLPCC, MCN(MNAGC) )
C
C           CONSTRUCTION EFFECTIVE DE LA MATRICE DE PRECONDITIONNEMENT
C           PAR FACTORISATION INCOMPLETE DE CHOLESKY AVEC NIVEAU
            NBLGRC(NRERR) = 1
            WRITE(KERR(2)(1:2),'(I2)') NIVEAU
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'FACTORISATION INCOMPLETE DE NIVEAU' //
     &                    KERR(2)(1:2)
            ELSE
               KERR(1) = 'INCOMPLETE CHOLESKY FACTORIZATION of LEVEL' //
     &                    KERR(2)(1:2)
            ENDIF
            CALL LERESU
            CALL INCHGC( NTDL,
     %                   MCN(MNLPLA),MCN(MNLPCA),MG,
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
               IF( COMP ) GO TO 9000
               IF( NIVEAU .LT. NIVMAX ) THEN
C                 TENTATIVE D'AUGMENTER LE NIVEAU DE FACTORISATION INCOMPLETE
                  NIVEAU = NIVEAU + 1
                  CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
                  CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
                  CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
                  GO TO 30
               ELSE
C                 LA FACTORISATION EST VRAIMENT TRES INSTABLE
C                 PLUS DE NIVMAX NIVEAUX DEPASSE => ABANDON DU GC
                  WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
                  NBLGRC(NRERR) = 3
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     %                    // KERR(MXLGER)(1:8)
                  KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
                     KERR(3) = 'AUGMENTER NIVEAU OU CHOISIR CHOLESKY'
                  ELSE
                   KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %                    // KERR(MXLGER)(1:8)
                     KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
                     KERR(3) = 'AUGMENT LEVEL or CHOOSE CHOLESKY'
                  ENDIF
C                 DESTRUCTION DES TABLEAUX DES MATRICES DU GC
                  CALL LEREUR
                  IERR = 8
                  GOTO 9000
               ENDIF
            ENDIF
C
C           LE PREMIER TABLEAU AUXILIAIRE DE GCPRCH
            MNAUGC = 0
            LOAUGC = NTDL * 3
            CALL TNMCDC( 'REEL2', LOAUGC, MNAUGC )
C           REPARTITION INTERNE EN SOUS-TABLEAUX
            MNAUX2 = MNAUGC
            MNAUX3 = MNAUX2 + NTDL * MOREE2
            MNAUX4 = MNAUX3 + NTDL * MOREE2
C
C           LE SECOND TABLEAU AUXILIAIRE DE GCPRCH SELON LA
C           STABILISATION DU GC PAR RE-ORTHOGONALISATION DES DIRECTIONS
            IF( ISTAB .EQ. 0 ) THEN
                NBDIR = 1
            ELSE
                NBDIR = 10
            ENDIF
            LODIR = (NTDL+1) * NBDIR * 2
            MNBDIR = 0
            CALL TNMCDC( 'REEL2', LODIR, MNBDIR )
C           REPARTITION INTERNE EN SOUS-TABLEAUX
            MNDIR  = MNBDIR
            MNADIR = MNBDIR + NTDL * NBDIR * MOREE2
            MNDAD  = MNADIR + NTDL * NBDIR * MOREE2
            MNBET  = MNDAD  + NBDIR * MOREE2
            MOTSGC = MOREE2*(LOAUGC+LODIR)
            MOTMAT = MOTMAT + NTDL+LOLPCC*(1+MOREE2) + MOTSGC
C
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10035) NTDL+LOLPCC*(1+MOREE2), MOTSGC,MOTMAT
            ELSE
               WRITE(IMPRIM,20035) NTDL+LOLPCC*(1+MOREE2), MOTSGC,MOTMAT
            ENDIF
10035       FORMAT(' PRECONDITIONNEMENT MATRICE MORSE   =',I15,' MOTS'/
     %             ' TABLEAUX SUPPLEMENTAIRES GC        =',I15,' MOTS'/
     %             ' STOCKAGE TOTAL DE LA METHODE DU GC =',I15,' MOTS'/)
20035     FORMAT(' SPARSE PRECONDITIONED MATRIX =',I15,' MEMORY WORDS'/
     %           ' AUXILIARY CG ARRAYS          =',I15,' MEMORY WORDS'/
     %           ' TOTAL STORAGE of CG METHOD   =',I15,' MEMORY WORDS'/)
         ENDIF
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'MATRICES [M] [K] CONSTANTES CALCULEES'
      ELSE
         WRITE(IMPRIM,*) 'CONSTANT MATRICES [M] [K] COMPUTED'
      ENDIF
C     TEMPS CALCUL DE LA FORMATION DE [AG]=[L] t[L]  et  [M]
      DPREP = DINFO( 'DELTA CPU' )
C
C     INITIALISATION DU VECTEUR UG0 DU DEPLACEMENT INITIAL
C     AUX NOEUDS DU MAILLAGE ET A L'INSTANT  TEMPS INITIAL
C     ====================================================
      IF( LTDEP0 .EQ. 3 ) TEMPS = REAL( TEMPS + DT )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10040) TEMPS
      ELSE
         WRITE(IMPRIM,20040) TEMPS
      ENDIF
10040 FORMAT(/' Au TEMPS',G14.6,' CALCUL des DEPLACEMENTS VERTICAUX')
      CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
20040 FORMAT(/' At TIME',G14.6,' COMPUTATION of VERTICAL DISPLACEMENTS')
      NBTEMP = 0
      DT0    = 0D0
      IF( LTDEP0 .GE. 2 ) THEN
C
C        LE DEPLACEMENT INITIAL EST DEJA CALCULE DANS VECTEUR"DEPLACT
C        ------------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTVECT, MNVECT )
         IF( MNVECT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PAS DE VECTEUR"DEPLACT INITIAL'
            ELSE
               KERR(1) = 'UNKNOWN INITIAL VECTEUR"DEPLACT'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9000
         ENDIF
         IF( NTDL .NE. MCN(MNVECT+WBCOVE) ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NTDL
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOMBRE DE COMPOSANTES DU VECTEUR"TEMPERATURE'
               KERR(2) = 'DIFFERENT NOMBRE '//KERR(5)(1:8) // ' NOEUDS'
            ELSE
               KERR(1) = 'NUMBER of COMPONENTS of VECTEUR"TEMPERATURE'
               KERR(2) = 'DIFFERENT to'//KERR(5)(1:8)// ' NODE NUMBER'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9000
         ENDIF
C
C        RECHERCHE DU TEMPS INITIAL TPSINI COMME PLUS PROCHE TEMPS
C        PARMI LES TEMPS STOCKES DU CALCUL PRECEDENT
C        NOMBRE DE TEMPS STOCKES
         N      = MCN( MNVECT + WBCPIN )
         MNTIME = MNVECT + WECTEU + MOREE2*NTDL*MCN(MNVECT+WBVECT) - 1
C        NBTEMP LE NUMERO DU VECTEUR STOCKE AU TEMPS LE PLUS PROCHE DE TPSINI
         NBTEMP = 1
         R      = RINFO( 'GRAND' )
         DO 45 I=1,N
            RR = ABS( RMCN(MNTIME+I) - TPSINI )
            IF( RR .LT. R ) THEN
               NBTEMP = I
               R      = RR
            ENDIF
 45      CONTINUE
C
C        TRANSFERT DU DEPLACEMENT INITIAL NBTEMP DANS LE TABLEAU MNUG0
         MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBTEMP-1)
         TEMPS  = RMCN(MNTIME+NBTEMP)
         TPSINI = TEMPS
         WRITE(IMPRIM,*) 'TEMPS INITIAL EFFECTIF =',TEMPS
         CALL TRTATD( MCN(MNTEMP), MCN(MNUG0), NTDL )
C
         IF( N .GT. 0 .AND. NBTEMP .GT. 0 ) THEN
C           LE PAS DE TEMPS DE L'ETUDE PRECEDENTE
            DT0 = RMCN(MNTIME+NBTEMP) - RMCN(MNTIME+NBTEMP-1)
         ENDIF
C
         IF( LTDEP0 .EQ. 3 .AND. DT0 .GT. 0D0 ) THEN
C           SAUVEGARDE DE UG0 ET UG1
            CALL TRTATD( MCN(MNTEMP), MCN(MNUG1), NTDL )
            CALL TRTATD( MCN(MNTEMP-MOREE2*NTDL), MCN(MNUG0), NTDL )
         ELSE
            CALL TRTATD( MCN(MNTEMP), MCN(MNUG1), NTDL )
         ENDIF
C
      ELSE IF( LTDEP0 .EQ. 1 ) THEN
C
C        LE DEPLACEMENT UG0 EST OBTENU PAR CONSTANTE
C        -------------------------------------------
         DEPLA0 = RMCN( MNDEIN + WADEP0 )
         MN = (MNUG0-1) / 2
         DO 47 I=1,NBNOEU
C           LE VECTEUR UG0 EST INITIALISE
            DMCN( MN+I ) = DEPLA0
 47      CONTINUE
C
      ELSE IF( LTDEP0 .EQ. -1 ) THEN
C
C        LE DEPLACEMENT UG0 EST OBTENU PAR FONCTION
C        ------------------------------------------
         MN = (MNUG0-1) / 2
         MM = MNXYZN + WYZNOE - 3
         DO 50 I=1,NBNOEU
            MM    = MM + 3
            XPOIN = RMCN( MM )
            YPOIN = RMCN( MM + 1 )
            ZPOIN = RMCN( MM + 2 )
            CALL REDEIN( 5,NOOB,1,XPOIN,YPOIN,ZPOIN,MNDEIN, DEPLA0 )
C           LE VECTEUR UG0 EST INITIALISE
            DMCN( MN+I ) = DEPLA0
 50      CONTINUE
      ENDIF
c
ccc      mnn = ( MNUG0 - 1 ) / 2
ccc      print *,('  onde',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
C     ICI LE VECTEUR UG0 EST INITIALISE
C
      IF( NOSCHM .EQ. 1 ) THEN
C
C        SCHEMA DE NEWMARK INCONDITIONNELEMENT STABLE
C        CALCUL DU VECTEUR GLOBAL SECOND MEMBRE FG0=FG(T0)
         IEMG = 0
         IEKG = 0
         IEBG = 1
         CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %                D2PI,   NDIM,   NTDL,   VITEGt,
     %                NBTYEL, MNNPEF, NDPGST,
     %                MNTPOB, MXPOBA, MNTAUX,
     %                MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                MNTHER, MNTAEL, MNX,    MNNODL,
     %                NORESO, MNLPLK, MNLPCK,
     %                NBRDKG, MG,     NBRDKG, KG,     MNFG0,
     %                NCODSM, NCODSK, NBPTAF, IERR )
         IF( IERR .NE. 0 ) GOTO 9000
C
         IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C           CONSTRUCTION NO ET VALEUR DES SOURCES PONCTUELLES
            CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                   NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C           ASSEMBLAGE DES SOURCES PONCTUELLES
            CALL ASFONO( NTDL, 1, NBFNFX,MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %                   MCN(MNFG0) )
C
         ENDIF
C
      ENDIF
C
C     CONSTRUCTION DU VECTEUR UG1 A l'INSTANT T + DT
C     ==============================================
      IF( LTDEP0 .EQ. 3 .AND. DT0 .GT. 0D0 ) THEN
C
C        UG0 EST DEJA EN FAIT LE VECTEUR"DEPLACT QUI PRECEDE UG1=U(TPSINI)
C        L'ANCIEN PAS DE TEMPS EST DT0
C        UG0 EST LA COMBINAISON LINEAIRE DE UG0 et UG1 POUR LE NOUVEAU PAS DT
C        DE TEMPS QUI A PU CHANGER
         IF( ABS( DT - DT0 ) .GT. 1D-5 * DT ) THEN
            CALL CL2VED( NTDL, 1D0-DT/DT0, MCN(MNUG1),
     %                             DT/DT0, MCN(MNUG0),  MCN(MNUG0) )
         ENDIF
C
      ELSE
C
C        OUVERTURE DU TMS VITESSEINIT
C        ============================
         CALL LXTSOU(NTLXOB,'VITESSEINIT',NTVIIN,MNVIIN)
         LTVIT0=MCN(MNVIIN+WTVIT0)
C        L'ADRESSE DMCN DE UG0
         MN   = (MNUG0-1) / 2
C        L'ADRESSE DMCN DE UG1
         MNU1 = (MNUG1-1) / 2
         MM   = MNXYZN + WYZNOE - 3
         DO 60 I=1,NBNOEU
            MM    = MM + 3
            XPOIN = RMCN( MM     )
            YPOIN = RMCN( MM + 1 )
            ZPOIN = RMCN( MM + 2 )
            CALL REVIIN( 5,NOOB,1,XPOIN,YPOIN,ZPOIN,MNVIIN, VITES0 )
C           LE VECTEUR UG1 EST INITIALISE
            DMCN( MNU1+I ) = DMCN( MN+I ) + DT * VITES0(1)
 60      CONTINUE
C
C        TRAITEMENT DE LA VITESSE INITIALE PARTICULARISEE SUR DES PLSV
C        -------------------------------------------------------------
         IF( IEVIIN .GT. 0 ) THEN
C           AU MOINS UN PLSV SUPPORTE UN TMS 'VITESSEINIT'
            CALL DEPLT1( MOREE2, MNXYZN, NDPGST, MXDOTH, LPVIIN, MNVIIN,
     %                   1,      NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %                   DT,     MNUG0,  MNUG1 )
         ENDIF
C
      ENDIF
C     ICI LE VECTEUR UG1 EST INITIALISE
C
C     ADRESSAGE INITIALISATION DU VECTEUR"DEPLACT SUR LE MAILLAGE TOTAL
C     =================================================================
      CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTVECT, MNVECT )
      IF( NTVECT .GT. 0 ) THEN
C
         IF( NBTEMP .LE. 0 ) THEN
C
C           LE VECTEUR EST DETRUIT POUR ETRE REDECLARE
            CALL LXTSDS( NTLXOB, 'VECTEUR"DEPLACT' )
C
         ELSE
C
C           LE VECTEUR EXISTE DEJA. MXTEMP NOUVEAUX VECTEURS LUI
C           SONT AJOUTES A PARTIR DU VECTEUR NBTEMP QUI EST VECTEUR INITIAL
            MXTEMP = NBTEMP + MXTEMP
            IF( MXTEMP .GT. MCN(MNVECT+WBVECT) ) THEN
C              LE TMS EST AUGMENTE
               L = WECTEU + NTDL * MXTEMP * MOREE2 + MXTEMP
               CALL TAMSAU( NTVECT, L )
               CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT',
     %                      NTVECT, MNVECT )
            ENDIF
            MNTEMP = MNVECT + WECTEU
            MN     = MNTEMP + MOREE2 * NTDL * MCN(MNVECT+WBVECT)
            MNTIME = MNTEMP + MOREE2 * NTDL * MXTEMP
C           COPIE DES NBTEMP PREMIERS TEMPS EXISTANTS
            IF( NBTEMP .EQ. 0 ) GOTO 70
            DO 62 I=0,NBTEMP-1
               RMCN( MNTIME+I ) = RMCN( MN+I )
 62         CONTINUE
            GOTO 70
         ENDIF
C
      ENDIF
C
C     LE VECTEUR N'EXISTE PAS. IL EST CREE AVEC 1+MXTEMP VECTEURS A STOCKER
      NBTEMP = 1
      MXTEMP = 1 + MXTEMP
      L      = WECTEU + NTDL * MXTEMP * MOREE2 + MXTEMP
      CALL LXTNDC( NTLXOB, 'VECTEUR"DEPLACT', 'MOTS', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTVECT, MNVECT )
C
C     LE STOCKAGE DU TEMPS POUR CHAQUE VECTEUR STOCKE
 70   MNTIME = MNVECT + WECTEU + MOREE2 * NTDL * MXTEMP - 1
C
C     INITIALISATION DES TEMPS DE STOCKAGE DES VECTEUR"
      TSTOC = TEMPS + DTSTOC
C
C     STOCKAGE DU VECTEUR UG0 INITIAL DANS MNTEMP
      MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBTEMP-1)
      MNTHET = MNTEMP
      IF( LTDEP0 .EQ. 3 ) GOTO 80
C
C     STOCKAGE DU DEPLACEMENT INITIAL
      CALL TRTATD( MCN(MNUG0), MCN(MNTEMP), NBNOEU )
C
C     LE TEMPS DU VECTEUR TEMPERATURE INITIAL STOCKE
      RMCN(MNTIME+NBTEMP) = TEMPS
C
C     AFFICHAGE DU STOCKAGE DU VECTEUR INITIAL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10070) TEMPS,NTDL,NBTEMP
      ELSE
         WRITE(IMPRIM,20070) TEMPS,NTDL,NBTEMP
      ENDIF
C
10070 FORMAT(' Au TEMPS',G14.6,' STOCKAGE de',I7,' DEPLACEMENTS',
     %   ' COLONNE',I5,' du TMS VECTEUR"DEPLACT')
20070 FORMAT(' At TIME',G14.6,' STORAGE of',I7,' DISPLACEMENTS',
     %   ' COLUMN',I5,' of TMS VECTEUR"DEPLACT')
C
C     ==========================================
C     ACTUALISER LE TEMPS 1 = TEMPS INITIAL + DT
C     ==========================================
      TEMPS = REAL( TEMPS + DT )
      IF( TEMPS .GE. TSTOC*0.9999 ) THEN
C        STOCKAGE DU DEPLACEMENT A CET INSTANT TEMPS 1
         MNTEMP = MNTEMP + MOREE2 * NTDL
         CALL TRTATD( MCN(MNUG1), MCN(MNTEMP), NBNOEU )
C        LE NOMBRE DE VECTEURS DEPLACEMENT STOCKES
         NBTEMP = NBTEMP + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBTEMP) = TEMPS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10070) TEMPS,NTDL,NBTEMP
         ELSE
            WRITE(IMPRIM,20070) TEMPS,NTDL,NBTEMP
         ENDIF
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10040) TEMPS
      ELSE
         WRITE(IMPRIM,20040) TEMPS
      ENDIF
c
c     deplacements verticaux au temps DT
ccc      mnn = ( MNUG1 - 1 ) / 2
ccc      print *,('  onde',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
C
 80   IEMG = 0
      IEKG = 0
      IEBG = 1
      IF( NOSCHM .EQ. 1 .OR. NOSCHM .EQ. 3 ) THEN
C
C        CALCUL DU VECTEUR GLOBAL SECOND MEMBRE FG1=FG(T1)
C        =================================================
         CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %                D2PI,   NDIM,   NTDL,   VITEGt,
     %                NBTYEL, MNNPEF, NDPGST,
     %                MNTPOB, MXPOBA, MNTAUX,
     %                MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                MNTHER, MNTAEL, MNX,    MNNODL,
     %                NORESO, MNLPLK, MNLPCK,
     %                NBRDKG, MG,     NBRDKG, KG,     MNFG1,
     %                NCOM,   NCOK,   NBPTAF, IERR )
         IF( IERR .NE. 0 ) GOTO 9000
C
         IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C           CONSTRUCTION NO ET VALEUR DES SOURCES PONCTUELLES
            CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                   NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C           ASSEMBLAGE DES SOURCES PONCTUELLES
            CALL ASFONO( NTDL, 1, NBFNFX,MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %                   MCN(MNFG1) )
C
         ENDIF
C
      ENDIF
C
C     BILAN DE L'INSTANT TEMPS1 = TEMPS INITIAL + DT
C     ==============================================
C     ##############################################################
C     ##                                                          ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DT ##
C     ##                                                          ##
C     ##############################################################
C
C     MNUG0  CONTIENT  U(0) DEPLACEMENT   A L'INSTANT t n-1
C     MNUG1  CONTIENT  U(1) DEPLACEMENT   A L'INSTANT t n
C     MNFG0  CONTIENT  F(0) SECOND MEMBRE A L'INSTANT t n-1
C     MNFG1  CONTIENT  F(1) SECOND MEMBRE A L'INSTANT t n
C
C     LE NOUVEAU TEMPS OU SE FAIT LE CALCUL
 100  TEMPS = REAL( TEMPS + DT )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10040) TEMPS
      ELSE
         WRITE(IMPRIM,20040) TEMPS
      ENDIF
C
C     *******************************************************************
C     NOSCHM = 1 : NEWMARK  IMPLICIT inconditionnellement STABLE
C     -----------  ----------------------------------------------------
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [MG](Un+1- 2Un + Un-1) = -DT**2 [KG]Un + DT**2/4 (Fn+1 + 2Fn + Fn-1)
C     avec [MG] = [ M + DT**2/4 * K ] et [KG] = [ K ]
C     et se RESOUD en 5 PARTIES DISTINCTES :
C     1.  calcul FG2 = SECOND MEMBRE A L'INSTANT T(N+1)
C     2.  FG0 <=  DT**2/4 ( FG2 + 2*FG1 + FG0 )
C     3.  FG0 <= -DT**2 [KG] UG1 + FG0
C     4.  UG2 <= [MG]**-1 FG0
C     5.  UG2 <= UG2 + 2 UG1 - UG0
C     ******************************************************************
C     NOSCHM = 2 : SCHEMA DECENTRE IMPLICIT inconditionnellement STABLE
C     -----------  ----------------------------------------------------
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [MG] Un+1 =  [KG] ( 2 Un - Un-1 ) + DT**2 * Fn+1
C     avec [MG] = [ M + DT**2 * K ] et [KG] = [ M ]
C     et se RESOUD en 5 PARTIES DISTINCTES :
C     1.  calcul FG2 = SECOND MEMBRE A L'INSTANT T(N+1)
C     2.  UG0 <= 2*UG1 - UG0
C     3.  UG2 <= [KG] UG0
C     4.  UG0 <= UG2 + DT**2 * FG2
C     5.  UG2 <= [MG]**-1 UG0
C     ******************************************************************
C     NOSCHM = 3 : EXPLICITE=>CONDITION CFL impose un petit pas de temps
C     -----------  -----------------------------------------------------
C     A L'ETAPE n+1, LE PROBLEME CONSISTE a TROUVER Un+1 SOLUTION DE
C     [ M ] (Un+1 + Un-1) = [ 2M - DT**2 K ] Un + DT**2 Fn
C     puis Un+1 = (Un+1+Un-1) - Un-1
C     et se RESOUD en 5 PARTIES DISTINCTES :
C     1.  calcul FG2 = SECOND MEMBRE A L'INSTANT T(N+1)
C     2.  FG0 <= [KG]  UG1
C     3.  FG0 <=  FG0 + DT**2 FG1 <--! FG1=SECOND MEMBRE A L'INSTANT T(N)
C     4.  UG2 <= [MG]**-1 FG0
C     5.  UG2 <=  UG2 - UG0
C     ******************************************************************
C
C     ETAPE 1: CALCUL DU VECTEUR GLOBAL SECOND MEMBRE FG2=FG(T2)
C     ======== =================================================
      CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     %             D2PI,   NDIM,   NTDL,   VITEGt,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNTHER, MNTAEL, MNX,    MNNODL,
     %             NORESO, MNLPLK, MNLPCK,
     %             NBRDKG, MG,     NBRDKG, KG,     MNFG2,
     %             NCOM,   NCOK,   NBPTAF, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C        CONSTRUCTION NO ET VALEUR DES SOURCES PONCTUELLES
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES SOURCES PONCTUELLES
         CALL ASFONO( NTDL, 1, NBFNFX, MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %                MCN(MNFG2) )
C
      ENDIF
C
C     ETAPES SUIVANTES
C     ================
      IF( NOSCHM .EQ. 1 ) THEN
C
C        NOSCHM = 1 : NEWMARK  IMPLICIT inconditionnellement STABLE
C        -----------  ---------------------------------------------
C        2. FG0 <=  DT**2/4 ( FG2 + 2*FG1 + FG0 )
         CALL CL3VED( NTDL, DTDT/4D0, MCN(MNFG2),
     %                      DTDT/2D0, MCN(MNFG1),
     %                      DTDT/4D0, MCN(MNFG0),  MCN(MNFG0) )
C
C        3. FG0 <= -DT**2 [KG] UG1 + FG0
         IF( NORESO .EQ. 1 ) THEN
            CALL MAPRVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), KG, MCN(MNUG1),
     %                   MCN(MNUG2) )
         ELSE
            CALL MAGCVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), MCN(MNLPCB), KG,
     %                   MCN(MNUG1),  MCN(MNUG2) )
         ENDIF
         CALL CL2VED( NTDL, -DTDT,  MCN(MNUG2),
     %                        1D0,  MCN(MNFG0),  MCN(MNFG0) )
C
C        4. UG2 <= [MG]**-1 FG0
         IF( NORESO .EQ. 1 ) THEN
C           RESOLUTION DU SYSTEME FACTORISE
            CALL DRCHPR( NTDL,NCODSA,MCN(MNLPLA),MG,
     %                   MCN(MNFG0),2, MCN(MNUG2) )
         ELSE
C
C           RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE PRECONDITIONNE
C           LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR CALCULE UG1
C           MISE A ZERO DU TABLEAU DES DIRECTIONS
            CALL AZEROD( LODIR, MCN(MNBDIR) )
            CALL GCPRCH( NTDL,        1,           NBDIR,
     %                   MCN(MNLPLA), MCN(MNLPCA), MG,       MCN(MNFG0),
     %                   MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC),
     %                   MCN(MNUG1),  MCN(MNAUX2), MCN(MNAUX3),
     %                   MCN(MNAUX4), MCN(MNDIR),  MCN(MNADIR),
     %                   MCN(MNDAD),  MCN(MNBET),  MCN(MNUG2),  IERR )
C
C           VERIFICATION DE LA CONVERGENCE DU GC
            IF( IERR .NE. 0 ) GOTO 9990
         ENDIF
C
C        5. UG2 <= UG2 + 2 UG1 - UG0
         CALL CL3VED( NTDL, 1D0, MCN(MNUG2),
     %                      2D0, MCN(MNUG1),
     %                     -1D0, MCN(MNUG0),  MCN(MNUG2) )
C
      ELSE IF( NOSCHM .EQ. 2 ) THEN
C
C        NOSCHM = 2 : SCHEMA DECENTRE IMPLICIT inconditionnellement STABLE
C        -----------------------------------------------------------------
C        2. UG0 <= 2*UG1 - UG0
         CALL CL2VED( NTDL,  2D0, MCN(MNUG1),
     %                      -1D0, MCN(MNUG0),  MCN(MNUG0) )
C
C        3. UG2 <= [KG] UG0
         IF( NORESO .EQ. 1 ) THEN
            CALL MAPRVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), KG, MCN(MNUG0),
     %                   MCN(MNUG2) )
         ELSE
            CALL MAGCVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), MCN(MNLPCB), KG,
     %                   MCN(MNUG0),  MCN(MNUG2) )
         ENDIF
C
C        4. UG0 <= UG2 + DT**2 * FG2
         CALL CL2VED( NTDL,  1D0, MCN(MNUG2),
     %                      DTDT, MCN(MNFG2),  MCN(MNUG0) )
C
C        5. UG2 <= [MG]**-1 UG0
         IF( NORESO .EQ. 1 ) THEN
C
C           RESOLUTION DU SYSTEME FACTORISE
            CALL DRCHPR( NTDL,NCODSA,MCN(MNLPLA),MG,
     %                   MCN(MNUG0),2, MCN(MNUG2) )
         ELSE
C
C           RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE PRECONDITIONNE
C           LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR CALCULE UG1
C           MISE A ZERO DU TABLEAU DES DIRECTIONS
            CALL AZEROD( LODIR, MCN(MNBDIR) )
            CALL GCPRCH( NTDL,        1,           NBDIR,
     %                   MCN(MNLPLA), MCN(MNLPCA), MG,       MCN(MNUG0),
     %                   MCN(MNLPLC), MCN(MNLPCC), MCN(MNAGC),
     %                   MCN(MNUG1),  MCN(MNAUX2), MCN(MNAUX3),
     %                   MCN(MNAUX4), MCN(MNDIR),  MCN(MNADIR),
     %                   MCN(MNDAD),  MCN(MNBET),  MCN(MNUG2),  IERR )
C
C           VERIFICATION DE LA CONVERGENCE DU GC
            IF( IERR .NE. 0 ) GOTO 9990
         ENDIF
C
      ELSE
C
C        NOSCHM = 3 : EXPLICITE=>CONDITION CFL impose un petit pas de temps
C        -----------  -----------------------------------------------------
C     1.  calcul FG2 = SECOND MEMBRE A L'INSTANT T(N+1)
C     2.  FG0 <= [KG]  UG1
C     3.  FG0 <=  FG0 + DT**2 FG1 <--! FG1=SECOND MEMBRE A L'INSTANT T(N)
C     4.  UG2 <= [MG]**-1 FG0
C     5.  UG2 <=  UG2 - UG0
C
C        2. FG0 <= [KG]  UG1
         IF( NORESO .EQ. 1 ) THEN
            CALL MAPRVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), KG, MCN(MNUG1),
     %                   MCN(MNFG0) )
ccc      call prcomx( ntdl, ncodsb,  mcn(mnlplb), kg, comxbg )
ccc      print *,'comxbg=',comxbg
         ELSE
            CALL MAGCVE( 0, 1D0, NTDL,
     %                   NCODSB, MCN(MNLPLB), MCN(MNLPCB), KG,
     %                   MCN(MNUG1),  MCN(MNFG0) )
         ENDIF
ccc      print *,('  KG ',kkl,'=',kg(MCN(MNLPLB+kkl)),kkl=1,6)
C
ccc      mnn = ( MNFG0 - 1 ) / 2
ccc      print *,('  FG0',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
C        3. FG0 <=  FG0 + DT**2 FG1 <--! FG1=SECOND MEMBRE A L'INSTANT T(N)
         CALL CL2VED( NTDL,  1D0, MCN(MNFG0),
     %                      DTDT, MCN(MNFG1),  MCN(MNFG0) )
ccc      mnn = ( MNFG0 - 1 ) / 2
ccc      print *,('  FG0',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
C
C        4. UG2 <= [MG]**-1 FG0
         IF( NCODSA .NE. 0 ) THEN
C
C           MATRICE MG NON DIAGONALE
            IF( NORESO .EQ. 1 ) THEN
C              MATRICE MG PROFIL SYMETRIQUE
ccc      call prcomx( ntdl, ncodsa,  mcn(mnlpla), mg, comxbg )
ccc      print *,'comxag=',comxbg
               CALL DRCHPR( NTDL, NCODSA,  MCN(MNLPLA), MG,
     %                      MCN(MNFG0), 2, MCN(MNUG2) )
            ELSE
C
C              RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE PRECONDITIONNE
C              LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR CALCULE UG1
C              MISE A ZERO DU TABLEAU DES DIRECTIONS
               CALL AZEROD( LODIR, MCN(MNBDIR) )
               CALL GCPRCH( NTDL,      1,          NBDIR,
     %                     MCN(MNLPLA),MCN(MNLPCA),MG,       MCN(MNFG0),
     %                     MCN(MNLPLC),MCN(MNLPCC),MCN(MNAGC),
     %                     MCN(MNUG1), MCN(MNAUX2),MCN(MNAUX3),
     %                     MCN(MNAUX4),MCN(MNDIR), MCN(MNADIR),
     %                     MCN(MNDAD), MCN(MNBET), MCN(MNUG2),  IERR )
C
C              VERIFICATION DE LA CONVERGENCE DU GC
               IF( IERR .NE. 0 ) GOTO 9990
C
            ENDIF
C
         ELSE
C
C           MATRICE DIAGONALE MG => DIVISION PAR LES COEFFICIENTS DIAGONAUX DE M
            MNB = ( MNFG0 - 1 ) / MOREE2
            MN  = ( MNUG2 - 1 ) / MOREE2
            DO 110 I=1,NTDL
               DMCN(MN+I) = DMCN(MNB+I) / MG(I)
 110        CONTINUE
         ENDIF
ccc      mnn = ( MNUG2 - 1 ) / 2
ccc      print *,('  UG2',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
C
C        5. UG2 = UG2 - UG0
         CALL CL2VED( NTDL,  1D0, MCN(MNUG2),
     %                      -1D0, MCN(MNUG0),  MCN(MNUG2) )
c
ccc         mnn = ( MNUG2 - 1 ) / 2
ccc         print *,('  onde',kkl,'=',dmcn(mnn+kkl),kkl=1,6)
      ENDIF
C
      IF( TEMPS .GE. TSTOC*0.9999 ) THEN
C
C        STOCKAGE DU DEPLACEMENT A CET INSTANT TEMPS
         MNTEMP = MNTEMP + MOREE2 * NTDL
         CALL TRTATD( MCN(MNUG2), MCN(MNTEMP), NBNOEU )
C        LE NOMBRE DE VECTEURS DEPLACEMENT STOCKES
         NBTEMP = NBTEMP + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBTEMP) = TEMPS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10070) TEMPS,NTDL,NBTEMP
         ELSE
            WRITE(IMPRIM,20070) TEMPS,NTDL,NBTEMP
         ENDIF
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC
C
      ENDIF
C
C     PERMUTATION DES VECTEURS UG0 UG1 UG2 et FG0 FG1 FG2
C     ===================================================
      MN    = MNUG0
      MNUG0 = MNUG1
      MNUG1 = MNUG2
      MNUG2 = MN
C
      MN    = MNFG0
      MNFG0 = MNFG1
      MNFG1 = MNFG2
      MNFG2 = MN
C
C     UNE ITERATION EN TEMPS DE PLUS EST ELLE NECESSAIRE?
      IF( TEMPS + DT .LT. TPSFIN*1.00001 ) GOTO 100
C
C    ##############################################################
C    ##                                                          ##
C    ##                FIN DE LA BOUCLE EN TEMPS                 ##
C    ##                                                          ##
C    ##############################################################
C
      IF( RMCN(MNTIME+NBTEMP) .NE. TEMPS ) THEN
C        STOCKAGE DES DEPLACEMENTS A CET INSTANT
         MNTEMP = MNTEMP + MOREE2 * NTDL
         CALL TRTATD( MCN(MNUG1), MCN(MNTEMP), NBNOEU )
C        LE NOMBRE DE VECTEURS DEPLACEMENT STOCKES
         NBTEMP = NBTEMP + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBTEMP) = TEMPS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10070) TEMPS,NTDL,NBTEMP
         ELSE
            WRITE(IMPRIM,20070) TEMPS,NTDL,NBTEMP
         ENDIF
      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"DEPLACT'
C     ====================================
      MCN( MNVECT + WBCOVE ) = NTDL
      MCN( MNVECT + WBVECT ) = NBTEMP
      MCN( MNVECT + WBCPIN ) = NBTEMP
      IF( NBTEMP .LT. MXTEMP ) THEN
C        LE TMS EST RACOURCI
         L  = MNVECT + WECTEU + NTDL * MXTEMP * MOREE2 - 1
         L1 = MNVECT + WECTEU + NTDL * NBTEMP * MOREE2 - 1
         DO 150 I=1,NBTEMP
            RMCN(L1+I) = RMCN(L+I)
 150     CONTINUE
         CALL TAMSRA( NTVECT, WECTEU+NTDL*NBTEMP*MOREE2+NBTEMP )
      ENDIF
C     LA DATE
      CALL ECDATE( MCN(MNVECT) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVECT + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
C
C     LES VECTEURS DEPLACEMENTS COMMENCENT A L'ADRESSE MNTEMP
      MNTEMP = MNVECT + WECTEU
C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS DEPLACEMENT
      MNTIME = MNTEMP + NTDL * NBTEMP * MOREE2 - 1
C
C     AFFICHAGE DES DEPLACEMENTS
C     ==========================
      CALL AFDEPL( 10,    NBTEMP, MNXYZN,
     %              1,    NTDL,   NBTEMP, MCN(MNTEMP),
     %            DECMAX, NOFOTI, DEXMAX )
C
C     DESTRUCTION DES TABLEAUX TMC INUTILES
      IF( MNLPLK .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLK )
C
      IF( NORESO .EQ. 2 ) THEN
          IF( MNAUGC .GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUGC, MNAUGC )
          IF( MNBDIR .GT. 0 ) CALL TNMCDS( 'REEL2',  LODIR,  MNBDIR )
          IF( MNLPCK .GT. 0 ) CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCK )
          IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
          IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
          IF( MNAGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
      ENDIF
C
C     COUT CALCUL DES DEPLACEMENTS ET AFFICHAGE
      DDEPL = DINFO( 'DELTA CPU' )
C
C     ==========================================================
C     CALCUL DES FLUX EN CHAQUE POINT D INTEGRATION DES FACES DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     ==========================================================
      CALL THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     %             NDIM,   MOREE2, NBTEMP, NTDL,
     %             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %             MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MOTAEL, MNTAEL, MNX,    MNVECT )
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
9000  IF( IERKGALLOC .EQ. 0 )  DEALLOCATE( KG )
      IF( IERMGALLOC .EQ. 0 )  DEALLOCATE( MG )
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF )
      DO 9001 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
9001  CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA, MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX,    MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL,   MNTAEL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2',  128,      MNTHER )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX,   MNNODL )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL' ,  NBDLMX*NBCOOR, MNX )
      IF( MNFG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL*3,   MNFG   )
      IF( MNUG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL*3,   MNUG   )
      IF( MNNFNX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONFNX,   MNNFNX )
      IF( MNVFNX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOVFNX,   MNVFNX )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MODLFX,   MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MODLFX,   MNVDLX )
      IF( MNFLTO .GT. 0 ) CALL TNMCDS( 'REEL2',  MOFLTO,   MNFLTO )
      IF( MNFLPT .GT. 0 ) CALL TNMCDS( 'REEL2',  MOFLPT,   MNFLPT )
      IF( MNLPLK .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1,   MNLPLK )
      IF( MNTHDL .GT. 0 ) CALL NLDATADS
C
      IF( NORESO .EQ. 2 ) THEN
          IF( MNAUGC .GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUGC, MNAUGC )
          IF( MNBDIR .GT. 0 ) CALL TNMCDS( 'REEL2',  LODIR,  MNBDIR )
          IF( MNLPCK .GT. 0 ) CALL TNMCDS( 'ENTIER', NBRDKG, MNLPCK )
          IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
          IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
          IF( MNAGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC, MNAGC  )
      ENDIF
      IF( IERR .NE. 0 ) RETURN
C
C     COUT CALCUL DES FLUX ET DESTRUCTIONS
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19900) DPREP,DDEPL,DCPU, DPREP+DDEPL+DCPU
      ELSE
         WRITE(IMPRIM,29900) DPREP,DDEPL,DCPU, DPREP+DDEPL+DCPU
      ENDIF
19900 FORMAT(/' COUT CALCUL PREPARATION MATRICES=',F12.2,' SECONDES CPU'
     % /' COUT CALCUL DES DEPLACEMENTS    =',F12.2,' SECONDES CPU'
     % /' COUT CALCUL DES FLUX            =',F12.2,' SECONDES CPU'
     % /' COUT TOTAL                      =',F12.2,' SECONDES CPU')
29900 FORMAT(/' PREPARATION MATRICES=',F12.2,' CPU SECONDS'
     %       /' DISPLACEMENTS       =',F12.2,' CPU SECONDS'
     %       /' FLUX                =',F12.2,' CPU SECONDS'
     %       /' TOTAL               =',F12.2,' CPU SECONDS')
      RETURN
C
C     LA METHODE DU GC NE CONVERGE PAS => ABANDON DES CALCULS
 9990 WRITE(KERR(MXLGER)(1:8),'(I8)') NIVEAU
      NBLGRC(NRERR) = 3
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'NON CONVERGENCE DU GC AU NIVEAU '
     %       //  KERR(MXLGER)(1:8) // ' de PRECONDITIONNEMENT'
         KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
         KERR(3) = 'ESSAYER la METHODE de CHOLESKY MATRICE PROFIL'
      ELSE
         KERR(1) = 'NO CG CONVERGENCE DU GC at LEVEL '
     %       //  KERR(MXLGER)(1:8)
         KERR(2) = 'EXIT of CONJUGATE GRADIENT METHOD'
         KERR(3) = 'AUGMENT the PRECOND LEVEL or CHOOSE CHOLESKY'
      ENDIF
      IERR    = 9
      GOTO 9000
C
      END
