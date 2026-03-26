      SUBROUTINE TRELAS( KNOMOB, NTYSOL, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACER LA DEFORMEE ET/OU LES CONTRAINTES PRINCIPALES
C -----   D'UN OBJET 2D OU 3D DE NOM KNOMOB APRES LE CALCUL
C         SOIT DES VECTEURS DEPLACEMENTS
C         SOIT DES MODES PROPRES
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET
C NTYSOL : NO DU TYPE DE SOLUTION A TRACER
C          1 VECTEURS TEMPERATURES
C          2 VECTEURS PROPRES
C          3 VECTEURS DEPLACEMENTS D'UNE ONDE
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1989
C MODIF  : GUILLAUME GERBER ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2000
C MODIF  : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      PARAMETER (MXTYEL=7, LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
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
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      CHARACTER*(*)     KNOMOB
C
      IF( INTERA .LE. 0 ) THEN
C        DEMANDE DE TRACE EN MODE BATCH => ARRET DE MEFISTO
         CALL ARRET( 1 )
      ELSE IF( INTERA .GE. 3 ) THEN
         LORBITE = 1
         NORBITE = 0
      ENDIF
C
C     INITIALISATIONS POUR LA REMANENCE DES VALEURS
      NCAS   = 1
      CMFLEC = 2.5
      CMPCON = 0.
      NOPT   = 1
      IERR   = 0
      NBISO  = 11
      IAVTIT = 1
      MNTIME = 0
      NTVVPR = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
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
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN DEFINITION for OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     L'OBJET EST-IL STRUCTURE EN SOUS-DOMAINES ?
      NDOUNO = MCN(MNDFOB+WDOUNO)
      NDOUNO = 0
ccc      IF( NDOUNO .EQ. 1 ) THEN
cccC
cccC        TRACE PAR SOUS-DOMAINES
cccC        =======================
ccc         CALL SDTREL( KNOMOB, NTLXOB, MNDFOB, IERR )
ccc        RETURN
cccC
ccc      ENDIF
C
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ===============================================================
C     ADRESSAGE DES ADRESSES DES TABLEAUX ELEMENTS FINIS DE CET OBJET
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
      MNTELE = MNELEM + MXTYEL
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO ,
     %             NTPOGE, MNPOGE, NTNOEU, MNNOEU ,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTPOGE : NUMERO      DU TMS 'XYZPOINT'    DE L'OBJET
C     MNPOGE : ADRESSE MCN DU TMS 'XYZPOINT'    DE L'OBJET
C     NTNOEU : NUMERO      DU TMS 'XYZNOEUD'    DE L'OBJET
C     MNNOEU : ADRESSE MCN DU TMS 'XYZNOEUD'    DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS DU MAILLAGE
C     NTELEM : NUMERO      DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
C     MNELEM : ADRESSE MCN DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
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
      CALL DIMCOO( MCN(MNPOGE+WNBPOI), MCN(MNPOGE+WYZPOI), NDIM )
C
      IF( NTYSOL .EQ. 1 ) THEN
C
C        RECHERCHE DU TABLEAU VECTEUR"DEPLACT
C        ------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"DEPLACT', NTDEPL, MNDEPL )
         IF( NTDEPL .LE. 0 ) THEN
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
         ENDIF
C        MODE DE TRACE DES DEPLACEMENTS
         MODECO = 1
         AMPLID = 1.0
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Trace des DEPLACEMENTS VECTEUR"DEPLACT'
         ELSE
            WRITE(IMPRIM,*) 'Drawings of VECTOR"DISPLACEMENTS'
         ENDIF
C
      ELSE
C
C        RECHERCHE DU TABLEAU VECTEUR"MODEPROPRE
C        ---------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"MODEPROPRE', NTVVPR, MNVVPR )
         IF( NTVVPR .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET '// KNOMOB
               KERR(2) = 'MODES PROPRES NON CALCULES'
            ELSE
               KERR(1) = 'ERROR: OBJECT '// KNOMOB
               KERR(2) = 'EIGENMODES NOT COMPUTED'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9000
         ENDIF
         NTDEPL = NTVVPR
         MNDEPL = MNVVPR
C        TRACE DES MODES PROPRES
         MODECO = 2
         AMPLID = 0.02
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Trace VECTEURS"MODES PROPRES'
         ELSE
            WRITE(IMPRIM,*) 'Drawings of VECTOR"EIGENMODES'
         ENDIF
C
      ENDIF
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES NOEUDS AVEC DEFORMATION AMPLIF
C     NOMBRE DE COMPOSANTES DU VECTEUR DEPLACEMENT = NOMBRE DE DEGRES DE LIBERTE
      NTDL = MCN( MNDEPL + WBCOVE )
      IF( NTDL .NE. MCN(MNNOEU+WNBNOE)*NDIM ) THEN
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'NOMBRE DE DL DU MAILLAGE NON EGAL AU'
            KERR(3) = 'NOMBRE DES DL DES DEPLACEMENTS'
            KERR(4) = 'RESOUDRE A NOUVEAU LE PROBLEME D''ELASTICITE'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'NUMBER of MESH DoF NOT EQUAL to'
            KERR(3) = 'NUMBER of DISPLACEMENT Dof'
            KERR(4) = 'SOLVE AGAIN the ELASTICITY PROBLEM'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9000
      ENDIF
C
C     NOMBRE DE VECTEURS DEPLACEMENT STOCKES
      NBVECT = MCN(MNDEPL+WBVECT)
C     NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOE  = MCN(MNNOEU+WNBNOE)
C     DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
      NDIMLI = NTDL / NBNOE
C     NUMERO DU CAS A VISUALISER : AU DEPART C'EST LE DERNIER CALCULE
      NCAS = NBVECT
      CALL MIMXDE( NDIMLI, NBNOE, MCN(MNNOEU+WYZNOE),
     %             AMPLID, MCN(MNDEPL+WBCOVE), NBVECT,
     %             NCAS,   MCN(MNDEPL+WECTEU),
     %             COOEXT )
C
C     L'INSTANT TEMPS DU CAS NCAS
      NBTIME = MCN(MNDEPL+WBCPIN)
      IF( NBTIME .LE. 0 ) THEN
C
C        LES TEMPS NE SONT PAS STOCKES
         MNTIME = 0
         TEMPS  = 0.0
C
      ELSE
C
C        L'ADRESSE DES TEMPS SI PB INSTATIONNAIRE OU FREQUENCES PROPRES
         MNTIME = MNDEPL + WECTEU + NTDL * NBVECT * MOREE2 - 1
C        LE TEMPS EST CELUI DU DERNIER DEPLACEMENT STOCKE
         TEMPS = RMCN(MNTIME+NCAS)
C
      ENDIF
C
C     LA FENETRE EST EFFACEE
      CALL EFFACEMEMPX
C
C     -----------------------------------------------------------
C     TRACE DES CONTRAINTES, DES DEPLACEMENTS OU DE L'ERREUR 2D ?
C     -----------------------------------------------------------
C     CADRE ENGLOBANT DE l'OBJET
 10   CALL CADEXT( MNPOGE, COOEXT )
C
C     PLUS D'ITEMS VISIBLES
      CALL ITEMS0
C
C     LA PREPARATION DU TRACE: VISEE PAR DEFAUT
      NOTYVI = 0
      CALL VISEE0
C
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
      CALL LIMTCL( 'contdepl', NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
      IF( NMTCL .EQ. 0 ) THEN
C
C        NUMERO DU CAS A VISUALISER
C        ==========================
         NCVALS = 4
         CALL INVITE( 84 )
         CALL LIRENT( NCVALS, NCAS )
         IF( NCVALS .EQ. -1 ) GOTO 9000
C        PROTECTION DU NUMERO DE CAS A TRACER
         NCAS = MAX( 1, ABS(NCAS) )
         NCAS = MIN( NBVECT, NCAS )
         GOTO 10
C
      ENDIF
C
C     PROTECTION DU NUMERO DE CAS A TRACER
      NCAS = MAX( 1, ABS(NCAS) )
      NCAS = MIN( NBVECT, NCAS )
C     LES TEMPS ONT ILS ETE STOCKES?
      IF( MNTIME .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR"DEPLACT
         TEMPS = RMCN( MNTIME + NCAS )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRACE au TEMPS ',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'DRAWING at TIME ',TEMPS
         ENDIF
      ELSE
C        TEMPS INITIAL SUPPOSE NUL
         TEMPS  = 0
         MNTIME = 0
      ENDIF
C
      GOTO ( 100, 200, 300, 400, 10, 600, 10, 800 ), NMTCL
C
C     TRACE DES CONTRAINTES PRICIPALES
C     ================================
 100  CALL TRCONT( KNOMOB, NTLXOB, MODECO, NDIM,
     %             NBTYEL, MNTOPO, MNELEM, MNPOGE,
     %             NBVECT, MNTIME,
     %             NCAS  , NOPT  , CMFLEC, CMPCON )
      GOTO 10
C
C     TRACE DE L'OBJET DEFORME PAR SES DEPLACEMENTS AMPLIFIES
C     =======================================================
 200  CALL TRDEPL( KNOMOB, MODECO, NDIM,
     %             NBTYEL, MNTOPO, MNELEM, MNNOEU,
     %             NBVECT, MNDEPL, MNTIME, NCAS  , AMPLID )
      GOTO 10
C
C     Seuil de PLASTICITE = LIMITE de l'ELASTICITE du MATERIAU
C     CRITERE DE VON MISES avec si = i-eme contrainte principale
C     En 3d: Von MISES = sqrt( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 ) / sqrt(2)
C     En 2d: Von MISES = TRESCA = Abs( s1 - s2 )
C     =========================================================================
 300  MISTRE = 1
      GOTO 405
C
C     Seuil de PLASTICITE = LIMITE de l'ELASTICITE du MATERIAU
C     CRITERE DE TRESCA avec si = i-eme contrainte principale
C     En 3d: TRESCA = MAX( abs(s1-s2), abs(s2-s3), abs(s3-s1) ) / 2
C     En 2d: Von MISES = TRESCA = Abs( s1 - s2 )
C     =============================================================
 400  MISTRE = 2
C
C     TRACE DU CRITERE DE VON MISES ou de TRESCA de ce CAS
C     EN 2D ZONES de COULEURS des TRIANGLES et QUADRANGLES
 405  CALL TRVMTR( MISTRE, NDIM,   KNOMOB, NTLXOB, MODECO,
     %             NBTYEL, MNELEM, MNPOGE, NDPGST, NCAS,   IERR )
      GOTO 10
C
C     TRACE DE L'OBJET DEFORME PAR SES DEPLACEMENTS AMPLIFIES EN FREQUENCE PROPR
C     ==========================================================================
 600  IF( NTVVPR .LE. 0 ) GOTO 10
C      LE TMS VVPR EXISTE
      CALL TRFRPR( KNOMOB, MODECO, NDIM,
     %             NBTYEL, MNTOPO, MNELEM, MNNOEU,
     %             NBVECT, MNDEPL, MNTIME,
     %             NCAS  , AMPLID )
      GOTO 10
C
C     EN 2D ERREUR SI DEPLACEMENT_EXACT()
C     TRACE DE L'ERREUR ABSOLUE EN CHAQUE NOEUD
C     =========================================
C     EXISTENCE OU NON DE LA FONCTION 'DEPLACEMENT_EXACT'
 800  NOFOTI = NOFODEEX()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NDIM .EQ. 2 .AND. NOFOTI .GT. 0 ) THEN
         CALL TRZDXY( NOFOTI, NDIM,   KNOMOB, MODECO,
     %                NBTYEL, MNELEM, MNPOGE,
     %                NCAS,   NBVECT, MCN(MNNOEU+WNBNOE),
     %                MCN(MNDEPL+WECTEU),
     %                DEPMIN, NOEMIN, DEPMAX, NOEMAX )
      ENDIF
      GOTO 10
C
C     FIN DE L'EXECUTION
 9000 RETURN
      END
