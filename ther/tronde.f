      SUBROUTINE TRONDE( KNOMOB , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LE DEPLACEMENT DE L'ONDE DANS UN OBJET 2D OU 3D
C -----
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1999
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray              Mars 2021
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7, LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ctemps.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  DTMIN,  DTMAX
      DOUBLE PRECISION  DEXMAX, DECMAX
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(0:1) :: dptemp

      IF( INTERA .LE. 0 ) THEN
C        DEMANDE DE TRACE EN MODE BATCH => ARRET DE MEFISTO
         CALL ARRET( 100 )
      ENDIF
C
C     INITIALISATIONS POUR LA REMANENCE DES VALEURS
      MNNPEF = 0
      MNTIMES= 0
      MNTIME = 0
      NCAS   = 1
      NCAS0  = 1
      NCAS1  = 1
      TEMPS  = 0.0
      NOPT   = 1
      CMFLEC = 2.5
      CMPGRA = 0.0
C     PAS DE PROJECTION DE R6 DANS R3
      NOPROJ = -1
      IERR   = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
      CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         ELSE
            KERR(1)='ERROR: UNKNOWM DEFINITION of the OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9000
      ENDIF
C
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ADRESSAGE DES ADRESSES DES TABLEAUX ELEMENTS DE CET OBJET
      CALL TNMCDC( 'ENTIER' , 2*MXTYEL , MNNPEF )
      MNTELE = MNNPEF + MXTYEL
      CALL NDPGEL( NTLXOB , NTTOPO , MNTOPO ,
     %             NTXYZP , MNXYZP , NTXYZN , MNXYZN ,
     %             NBTYEL , MCN(MNTELE) , MCN(MNNPEF) , IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'    DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'    DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'    DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'    DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS DU MAILLAGE
C     NTNPEF : NUMERO      DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
C     MNNPEF : ADRESSE MCN DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
      IF( IERR .NE. 0 ) GOTO 10
C
C     NDPGST : CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS DU MAILLAGE
C              0 : NOEUDS=POINTS=SOMMETS
C              1 : NOEUDS=POINTS#SOMMETS
C              2 : NOEUDS#POINTS=SOMMETS
C              3 : NOEUDS#POINTS#SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     NDIM DIMENSION EFFECTIVE DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( MCN(MNXYZP+WNBPOI), MCN(MNXYZP+WYZPOI), NDIM )
      NDIMLI = NDIM
C
C     RECHERCHE DES VECTEURS DEPLACEMENTS DE L'OBJET
      CALL  LXLXOU( NTLXOB, 'VECTEUR"DEPLACT', NTTEMP, MNTEMP )
      IF( NTTEMP .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET '// KNOMOB
            KERR(2) = 'VECTEURS DEPLACEMENTS NON CALCULES'
         ELSE
            KERR(1) = 'ERROR: OBJECT '// KNOMOB
            KERR(2) = 'DISPLACEMENTS VECTORS NOT COMPUTED'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9000
      ELSE
C        MODE DE TRACE DES DEPLACEMENTS
         MODECO = 1
      ENDIF
C
C     NOMBRE DE DEGRES DE LIBERTE
      NTDL = MCN( MNTEMP + WBCOVE )
      IF( NTDL .NE. MCN(MNXYZN+WNBNOE) ) THEN
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR : OBJET ' // KNOMOB
            KERR(2) = 'NOMBRE des NOEUDS du MAILLAGE NON EGAL au'
            KERR(3) = 'NOMBRE des NOEUDS de CALCUL des DEPLACEMENTS'
            KERR(4) = 'RESOUDRE a NOUVEAU l''EQUATION des ONDES'
         ELSE
            KERR(1) = 'ERROR : OBJECT ' // KNOMOB
            KERR(2) = 'NUMBER of MESH NODES NOT EQUAL TO'
            KERR(3) = 'NUMBER of NODES WITH DISPLACEMENTS'
            KERR(4) = 'SOLVE AGAIN the WAVE EQUATION'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9000
      ENDIF
C
C     NOMBRE TOTAL DE CAS
      NDSM = MCN( MNTEMP + WBVECT )
C
C     L'ADRESSE DES TEMPS SI PB INSTATIONNAIRE
      MNTIME = MNTEMP + WECTEU + NTDL * NDSM * MOREE2 - 1
C
C     SI PB INSTATIONNAIRE LE DERNIER VECTEUR EST CHOISI
      NBTEMPS = MCN( MNTEMP + WBCPIN )
      IF( NBTEMPS .GT. 0 ) THEN
C        LE TEMPS FINAL
         NCAS0 = 1
         NCAS1 = NDSM
         NCAS  = NDSM
         TEMPS = RMCN(MNTIME+NDSM)
      ENDIF
C
C     CONSTRUCTION DU TABLEAU LESTEMPS A L'ADRESSE MCN MNTIMES
C     NBVECT NOMBRE DE VECTEUR"VITESSEPRESSION STOCKES
      CALL LESTEMPS( KNOMOB, MNTEMP, NDSM, MNTIMES, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C
C     LA FENETRE EST EFFACEE
      CALL EFFACEMEMPX
C
C     LE TITRE EST TRACE
      IAVTIT = 1
C
C     COEFFICIENT D'AMPLIFICATION DU DEPLACEMENT DE L'ONDE
      AMPLID = 1.0
C
C     *********************************************************
C     NO DU CAS OU TRACE DU DEPLACEMENT OU DE SON GRADIENT ?
C     *********************************************************
 10   CALL RECTEF( NRHIST )
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
      CALL LIMTCL( 'traconde' , NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
      IF( NMTCL .EQ. 0 ) GOTO 50
      GOTO ( 100, 200, 300, 400, 500 ) , NMTCL
C
C     COEFFICIENT D'AMPLIFICATION DU DEPLACEMENT DE L'ONDE
C     ====================================================
 50   NCVALS = 5
      AMPLID = 1.0
      CALL INVITE( 3 )
      CALL LIRRSP( NCVALS, AMPLID )
      IF( NCVALS .EQ. -1 ) GOTO 10
      IF( AMPLID .EQ. 0.0 ) AMPLID = 1.0
      GOTO 10
C
C     NUMERO DU CAS A VISUALISER
C     ==========================
 100  CALL LIRNOCAS( NDSM, NCAS0, NCAS1, N )
      IF( N .LT. 0 ) GOTO 10
      IF( NCAS0 .EQ. NCAS1 ) THEN
         NCAS = NCAS1
      ELSE
         NCAS = NDSM
      ENDIF
C     LES TEMPS ONT ILS ETE STOCKES?
      NBTEMPS = MCN( MNTEMP + WBCPIN )
      IF( NBTEMPS .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR"DEPLACEMENT
         TEMPS = RMCN( MNTIME + NCAS )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRACE AU TEMPS ',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'DRAWING at TIME ',TEMPS
         ENDIF
      ELSE
C        TEMPS INITIAL SUPPOSE NUL
         NBTEMPS = 0
         TEMPS   = 0
      ENDIF
      GOTO 10
C
C     TRACE DU DEPLACEMENT
C     ====================
C     CALCUL DES DEPLACEMENTS MIN ET MAX ET DE LEURS NOEUDS
 200  CALL MXVECT( NTDL, NDSM, MCN(MNTEMP+WECTEU),
     %             DTMIN, NOEMIN, NCAMIN, DTMAX, NOEMAX, NCAMAX )
      TMIN = REAL( DTMIN )
      TMAX = REAL( DTMAX )
C
C     OPTIONS DU TRACE DU DEPLACEMENT NCAS
 201  CALL LIMTCL( 'trdepond' , NMTCL0 )
      IF( NMTCL0 .LE. 0 ) GOTO 10
      GOTO ( 210, 220, 230, 240, 250, 260, 201, 201, 201, 290 ) , NMTCL0
C
C     TRACE DES ISOTHERMES (LIGNES EN 2D ET SURFACES EN 3D)
C     =====================================================
 210  CALL TRISOT( NDIM,   KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      MCN(MNTEMP+WECTEU), dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,   NOEMAX, NCAMAX,
     %             RMCN(MNTIMES) )
      GOTO 201
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES
C     ======================================
 220  CALL TRZONT( 0,      NDIM,   KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      MCN(MNTEMP+WECTEU), dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,   NOEMAX, NCAMAX,
     %             RMCN(MNTIMES) )
      GOTO 201
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 230  IF( NDIM .EQ. 2 ) GOTO 220
      CALL TRPLSE( 0,      KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      MCN(MNTEMP+WECTEU), dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,   NOEMAX, NCAMAX,
     %             RMCN(MNTIMES) )
      GOTO 201
C
C     TRACE DES PROFILS DE COULEURS ISOTHERMES PAR SECTIONS X ou Y ou Z=CTE
C     =====================================================================
 240  IF( NDIM .EQ. 2 ) GOTO 220
      CALL TRPLSE( 1,      KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      MCN(MNTEMP+WECTEU), dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,   NOEMAX, NCAMAX,
     %             RMCN(MNTIMES) )
      GOTO 201
C
C     TRACE EN 2D SURFACE(X,Y,DEPLACEMENT(X,Y)) DE TOUS LES PAS DE TEMPS
C     ==================================================================
 250  IF( NDIM .GT. 2 ) GOTO 10
C     ICI PAS DE DEPLACEMENT_EXACT A PRENDRE EN COMPTE
      NOFOTI = 0
      CALL TRZOXY( NDIM,   KNOMOB, MODECO, NBTEMPS, RMCN(MNTIME+1),
     %             NBTYEL, MNNPEF, MNXYZP,
     %             NDSM,   NTDL,   MCN(MNTEMP+WECTEU), AMPLID,
     %             TMIN,   TMAX )
      GOTO 201
C
C     TRACE DE L'ERREUR ABSOLUE EN CHAQUE NOEUD ET EN 2D
C     ==================================================
C     EXISTENCE OU NON DE LA FONCTION 'DEPLACEMENT_EXACT(t,x,y,z,nocomp)'
 260  NOFOTI = NOFODEEX()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NDIM .EQ. 2 .AND. NOFOTI .GT. 0 ) THEN
         CALL TRZTXY( NOFOTI, NDIM,   KNOMOB, MODECO,
     %                NBTYEL, MNNPEF, MNXYZP,
     %                NCAS0,  NCAS1,  NTDL,
     %                0,      MCN(MNTEMP+WECTEU), dptemp,
     %                TMIN,   NOEMIN, NCAMIN, TMAX, NOEMAX, NCAMAX,
     %                RMCN(MNTIMES) )
      ENDIF
      GOTO 201
C
C     AFFICHAGE DU VECTEUR DEPLACEMENT
C     ================================
 290  NOFOTI = NOFODEEX()
      CALL AFDEPL( NTDL,   NCAS,   MNXYZN,
     %             1,      NTDL,   NDSM,   MCN(MNTEMP+WECTEU),
     %             DECMAX, NOFOTI, DEXMAX )
      GOTO 201
C
C     TRACE DU GRADIENT DU DEPLACEMENT
C     ================================
 300  CALL TRGRAD(NOPROJ, NCAS0, NCAS1, NDIM, KNOMOB, NTLXOB, MODECO,
     %            NBTYEL, MNNPEF, MNXYZP, NDPGST, RMCN(MNTIMES))
      GOTO 10
C
C     TRACE DU FLUX NORMAL DU DEPLACEMENT
C     ===================================
 400  CALL TRFLUX(NOPROJ, NCAS0, NCAS1, NDIM, KNOMOB, NTLXOB, MODECO,
     %            NBTYEL, MNNPEF, MNXYZP, NDPGST, RMCN(MNTIMES))
      GOTO 10
C
C     TRACE DES ESTIMATEURS D'ERREUR EN 2D ET EN DEPLACEMENT SEULEMENT
C     ================================================================
 500  IF( NDIM .EQ. 2 .AND. MODECO .EQ. 1 ) THEN
         CALL TRERTH( NCAS1,  NDIM,   NTLXOB,
     %                NBTYEL, MNNPEF, MNXYZP, NDPGST )
      ENDIF
      GOTO 10
C
C     FIN DE L'EXECUTION
C     ==================
C     RETOUR AUX PARAMETRES INITIAUX
 9000 CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER' , 2*MXTYEL , MNNPEF )
      IF( MNTIMES.GT. 0 ) CALL TNMCDS( 'REEL', NDSM, MNTIMES)
      RETURN
      END
