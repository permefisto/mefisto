      SUBROUTINE T1OBJE( NOMOBJ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER L'OBJET NOMOBJ AVEC ou SANS INTERPOLATION DEFINIE
C------
C ENTREES:
C --------
C NOMOBJ : NOM DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 1994
C MODIF :PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C.......................................................................
      include"./incl/mecoit.inc"
      PARAMETER        (MXTYEL=7, MXPILE=128)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / MSSFTA / MSSF(28),NTADAM
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*(*)     NOMOBJ
      CHARACTER*24      KNOMOB
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*8       NMSOMM
      CHARACTER*9       KSYMBO
C
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ponoel.inc"
      include"./incl/xyzext.inc"

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'t1obje: TRACE de l''OBJET ',NOMOBJ
      ELSE
         PRINT*,'t1obje: Drawing of the OBJECT ',NOMOBJ
      ENDIF

      IERR    = 0
      MNPILE  = 0
      MNELEM  = 0
      MNFASU  = 0
      MNARLI  = 0
      MNSOPO  = 0
      MNPLS   = 0
      MNNUBA  = 0
      MNBARY  = 0
      MNXYZP  = 0
      MXPLS   = 0
      NBSO    = 0
      NBAR    = 0
      NBFA    = 0
      NBCOOR  = 3
      QUAMIN  = 0

C     SI LA SOURIS EST DISPONIBLE EN INTERACTIF ELLE EST ACTIVEE
C     POUR MODIFIER LE CADRE DES TRACES
ccc      IF( INTERA .GE. 2 ) LORBITE=1
      LORBITE = 0
      CALL XVEPAISSEUR( 0 )
      IAVNSO0 = IAVNSO
      PREDUF0 = PREDUF
      PREDUF  = 50.

C     L'OBJET EXISTE-T-IL ?
C     =====================
      CALL LXLXOU( NTOBJE, NOMOBJ, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // NOMOBJ
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // NOMOBJ
         ENDIF
         CALL LEREUR
         IERR   = 1
         KNOMOB = NOMOBJ
         GOTO 9999
      ENDIF
C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'OBJET: ' //  NOMOBJ
      CALL LHISTO

C     L'OBJET EST IL AVEC OU SANS INTERPOLATION ?
C     ===========================================
C     LE TABLEAU TOPOLOGIE DE L'OBJET
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )
      IF( NTTOPO .GT. 0 ) GOTO 100
C
 1    IF( NDIMLI .LE. 2 ) THEN
C        OBJET 1D ou 2D: TRACE AVEC ZOOM ou NON?
         IF( LORBITE .EQ. 0 ) GOTO 4
C        INITIALISATION DU ZOOM DEPLACEMENT
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 99999
C        POUR TRACER TOUS LES PLS AVANT DE REPRENDRE L'ORBITE
C        SINON ORBITE APRES TRACE DE L'UN DES PLS
         LORBITE = - ABS( LORBITE )
      ENDIF

C     ======================================================================
C     OBJET 1D ou 2D ou 3D SANS INTERPOLATION ou 1D ou 2D AVEC INTERPOLATION
C     ======================================================================
C     DECLARATION DE LA PILE DES PLSVO
 4    IF( MNPILE .EQ. 0 ) CALL TNMCDC( 'ENTIER', MXPILE, MNPILE )
      LHPILE = 1
      CALL NUOBNM( 'OBJET', NOMOBJ, NUO )
      IF( NUO .LE. 0 ) THEN
         KNOMOB = NOMOBJ
         GOTO 9999
      ENDIF
      MCN(MNPILE) = NUO
C
C     ********************************************************
C     TANT QUE LA PILE DES OBJETS EST NON VIDE TRACER LES PLSV
C     ********************************************************
 5    IF( LHPILE .GT. 0 ) THEN
C
C        LE NUMERO DE L'OBJET
         LHPILE = LHPILE - 1
         NUO    = MCN( MNPILE + LHPILE )
C
C        LE NOM KNOMOB DE L'OBJET
         CALL NMOBNU( 'OBJET', NUO, KNOMOB )
C
C        RECHERCHE DE L'OBJET
         CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
         IF( NTLXOB .LE. 0 ) GOTO 9999
C
C        OUVERTURE DU TABLEAU  'DEFINITION' DE L'OBJET
         CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
         IF( NTDFOB .LE. 0 ) THEN
C           TABLEAU 'DEFINITION' INEXISTANT
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'OBJET NON DEFINI:'//KNOMOB
            ELSE
               KERR(1) = 'NOT DEFINED OBJECT: '//KNOMOB
            ENDIF
            CALL LEREUR
            GOTO 20
         ENDIF
C
C        LE TRACE DANS L'ORDRE DES VOLUMES PUIS SURFACES PUIS
C        LIGNES ET ENFIN DES POINTS DES OBJETS EMPILES DE L'OBJET
C        --------------------------------------------------------
         IF( NTTOPO .GT. 0 . AND. IAVNSO0 .GT. 0 ) THEN
C           POUR NE PAS AVOIR LE TRACE DU NUMERO DES SOMMETS
            IAVNSO = 0
         ENDIF
         DO 15 NTYK=5,1,-1
C
C           TRAITEMENT DES OBJETS PUIS VOLUMES PUIS SURFACES ...
            MN = MNDFOB + WTYOBJ
            DO 10 I=1,MCN(MNDFOB+WBDOBJ)
C
C              LE NO DU TYPE DU PLSVO
               NTY = MCN( MN )
               IF( NTY .NE. NTYK ) GOTO 8
C
C              LE NUMERO DU PLSVO
               NUO = MCN( MN + 1 )
C
               IF( NTY .NE. 5 ) THEN
C
C                 LE NOM DU TYPE POINT OU LIGNE OU SURFACE OU VOLUME
                  KNOMTY = NMTYOB( NTY )
C                 LE NOM KNOMOB DE L'OBJET
                  CALL NMOBNU( KNOMTY, NUO, KNOMOB )
C                 TRACE EFFECTIF DU PLSV
                  CALL T1MOBJ( KNOMTY, KNOMOB, NUO )
C
               ELSE
C
C                 OBJET A EMPILER
                  IF( LHPILE .GE. MXPILE ) THEN
                     NBLGRC(NRERR) = 1
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(1) = 'PILE T1OBJE SATUREE'
                     ELSE
                        KERR(1) = 'T1OBJE: STACK SATURATED'
                     ENDIF
                     CALL LEREUR
                     GOTO 20
                  ENDIF
                  MCN( MNPILE + LHPILE ) = NUO
                  LHPILE = LHPILE + 1
               ENDIF
C
  8            MN = MN + 2
 10         CONTINUE
 15      CONTINUE
         GOTO 5
      ENDIF
C
C     TRACE DES POINTS DES EF DU MAILLAGE DE L'OBJET 1D ou 2D AVEC INTERPOLATION
C     --------------------------------------------------------------------------
      IF( NTTOPO .GT. 0 .AND. IAVNSO0 .GT. 0 .AND. MNXYZP .GT. 0 ) THEN
C        DIMENSION DE L'ESPACE DES COORDONNEES POUR LE TMS XYZPOINT
         NBCOOR = MCN(MNXYZP+WBCOOP)
         NBPOMA = MCN(MNXYZP+WNBPOI)
         MNX = MNXYZP + WYZPOI
         DO 18 NP=1,NBPOMA
            WRITE( NMSOMM, '(I7)' ) NP
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            CALL SYMBOLE2D( NCONSO, RMCN(MNX), RMCN(MNX+1), KSYMBO(1:L))
            MNX = MNX + NBCOOR
 18      CONTINUE
      ENDIF

C     TRACE ET AFFICHAGE DE LA QUALITE DU MAILLAGE DE L'OBJET
      CALL IMPQUAOB( NOMOBJ, QUAMIN, NBEFMQ, IERR )

C     TRACE DU TITRE ET FERMETURE
      LORBITE = ABS( LORBITE )
      CALL TRFINS( NOMOBJ )

C     REPRISE DE TRANSLATION ZOOM EN 1D ou 2D
      IF( NDIMLI .LE. 2 .AND. LORBITE .GT. 0 ) THEN
C        ZOOM et TRANSLATION ACTIFS EN 1D ou 2D
         CALL ZOOM2D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
C        POUR TRACER TOUS LES PLS AVANT DE REPRENDRE L'ORBITE
C        SINON ORBITE APRES TRACE DE L'UN DES PLS
         LORBITE = - ABS( LORBITE )
         LHPILE = 1
         GOTO 5
      ENDIF

C     DESTRUCTION DE LA PILE
 20   IF( MNPILE .GT. 0 ) CALL TNMCDS( 'ENTIER', MXPILE, MNPILE )
      GOTO 99999

C     ========================
C     OBJET AVEC INTERPOLATION
C     ========================
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ADRESSAGE DES ADRESSES DES TABLEAUX ELEMENTS DE CET OBJET
 100  MNELEM = 0
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
      MNTELE = MNELEM + MXTYEL
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'    DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'    DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'    DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'    DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS FINIS DU MAILLAGE
C     MNTELE : NUMERO      DU TMS DES NBTYEL TYPES D'ELEMENTS
C     MNELEM : ADRESSE MCN DU TMS DES NBTYEL TYPES D'ELEMENTS
      IF( IERR .NE. 0 ) GOTO 9000
C
C     NDPGST : CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS DU MAILLAGE
C              0 : NOEUDS=POINTS=SOMMETS
C              1 : NOEUDS=POINTS#SOMMETS
C              2 : NOEUDS#POINTS=SOMMETS
C              3 : NOEUDS#POINTS#SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     DIMENSION DE L'ESPACE DES COORDONNEES POUR LE TMS XYZPOINT
      NBCOOR = MCN(MNXYZP+WBCOOP)
C
C     LES COORDONNEES EXTREMES SONT CELLES DE CET OBJET
      NBPOMA = MCN(MNXYZP+WNBPOI)
      CALL MIMXPT( NBCOOR, NBPOMA,
     %             RMCN(MNXYZP+WYZPOI), COOEXT )
C     CADRE AVEC 15% EN PLUS
      DO 30 I=1,NBCOOR
         EC = ( COOEXT(I,2) - COOEXT(I,1) ) * 1.15 / 2
         CM = ( COOEXT(I,1) + COOEXT(I,2) ) / 2
         COOEXT(I,1) = CM - EC
         COOEXT(I,2) = CM + EC
 30   CONTINUE

      IF( NBCOOR .EQ. 6 ) THEN

C        =================================================================
C        TRACE DES 192 ARETES DES 6-CUBES DE L'OBJET 6D AVEC INTERPOLATION
C        =================================================================
C        LA DIRECTION DE VISEE EN 6D
         AXOLAR = 0
         DO 35 N=1,6
            AXOPTV(N) = ( COOEXT(N,1) + COOEXT(N,2) ) * 0.5
            S = COOEXT(N,2) - COOEXT(N,1)
            AXOLAR = MAX( AXOLAR, S )
            AXOEIL(N) = AXOPTV(N) + S * N * 0.4
 35      CONTINUE
         AXOLAR = AXOLAR / 1.2
         AXOHAU = AXOLAR

C        CONSTRUCTION DE LA MATRICE D'AXONOMETRIE 6D
         CALL MATAXO6

C        DEFINITION DE LA FENETRE EN LARGEUR ET HAUTEUR
         CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )

C        TRACE FINAL DANS R3
         NDIM   = 3
         NDIMLI = NDIM

C        TRACE DES ARETES DES 6-CUBES EN PROJECTION XYZ (OUBLI DE UVW)
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
         GOTO 55

C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
 53      CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C        TRACE DES 192 ARETES DES 6-CUBES
 55      CALL T3AR6C( NOMOBJ )
C        TRACE DES AXES DU CUBE ENGLOBANT
         CALL TRAXE3
C        TRACE DU TITRE ET FERMETURE
         CALL TRFINS( NOMOBJ )

C        REPRISE DE L'ORBITE
         IF( LORBITE .NE. 0 ) GOTO 53
C        SORTIE DU TRACE
         NOTYVI = 0
         GOTO 9000

      ENDIF

C     OBJET 1D ou 2D ou 3D AVEC POINTS a 3 COORDONNEES dans XYZPOINT
      CALL DIMCOO( NBPOMA, MCN(MNXYZP+WYZPOI), NDIM )
      NDIMLI = NDIM
      IF( NDIM .LE. 2 ) GOTO 1

C     ======================================
C     TRACE DE L'OBJET 3D AVEC INTERPOLATION
C     ======================================

C     LES ARETES FRONTALIERES DES VOLUMES DE L'OBJET
C     ----------------------------------------------
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE
      CALL HACHOB( NOMOBJ, 4, NTFAOB, MNFAOB, IERR )
      IF( IERR .GT. 0 .OR. NTFAOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR T1OBJE: OBJET ' // NOMOBJ
            KERR(2) = 'CALCUL IMPOSSIBLE DE SES FACES'
         ELSE
            KERR(1) = 'ERROR T1OBJE: OBJECT ' // NOMOBJ
            KERR(2) = 'FACES NOT COMPUTABLE'
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF

C     CREATION DU HACHAGE DES ARETES DES FACES FRONTALIERES DE L'OBJET
      CALL HACHAF( NOMOBJ, 0, NTFAOB, MNFAOB, NTAFOB, MNAFOB, I )

C     LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
      MOARFR = MCN( MNAFOB + WOARFR )
C     LA MAJORATION DU NOMBRE DES ARETES FRONTALIERES
      MXARFR = MCN( MNAFOB + WXARFR )
C     LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
      L1ARFR = MCN( MNAFOB + W1ARFR )
C     LE NOMBRE D'ARETES FRONTALIERES DANS LE CHAINAGE
      NBARFR = MCN( MNAFOB + WBARFR )

C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C     POUR CALCULER LE NOMBRE DE FACES, ARETES, SOMMETS APPARTENANT A
C     DES SURFACES, LIGNES ET POINTS NOMMES PAR L'UTILISATEUR
C     ===============================================================
      NBFA = 0
      NBAR = 0
      NBSO = 0

      DO 110 I = 0, NBTYEL-1

C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )

C        LES SOMMETS IDENTIQUES A DES POINTS UTILISATEUR
         NBSO = NBSO + MCN( MNELE + WBPSEL )

C        LES ARETES APPARTENANT A UNE LIGNE
         NBAR = NBAR + MCN( MNELE + WBLAEL )

C        LES FACES APPARTENANT A UNE SURFACE
         NBFA = NBFA + MCN( MNELE + WBSFEL )

 110  CONTINUE

C     CREATION DU TABLEAU DES NUMEROS DANS XYZPOINT DES SOMMETS
C     DES FACES ARETES SOMMETS APPARTENANT A UNE SURFACE LIGNE POINT
C     --------------------------------------------------------------
      MOTFAS = 5*NBFA + 3*NBAR + 2*NBSO
      PRINT*,'t1obje:',MOTFAS,' MAXIMUM de FACES+ARETES+SOMMETS pour TRA
     %CER ',NOMOBJ
      CALL TNMCMX( 'ENTIER', MAXVAR )
      IF( MAXVAR .LT. MOTFAS ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'t1obje: super-tableau MCN TROP COURT MOTMCN=',MOTMCN
         ELSE
            PRINT*,'t1obje: super-array MCN TOO SHORT MOTMCN=',MOTMCN
         ENDIF
         GOTO 9000
      ENDIF
 
      IF( NBFA .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', 5*NBFA, MNFASU )
      ENDIF
      MNF = MNFASU - 1

      MNARLI = 0
      IF( NBAR .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', 3*NBAR, MNARLI )
      ENDIF
      MNA = MNARLI - 1

      MNSOPO = 0
      IF( NBSO .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', 2*NBSO, MNSOPO )
      ENDIF
      MNS = MNSOPO - 1

C     INITIALISATION DE CES 3 TABLEAUX
      DO 200 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN( MNELE + WBELEM )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         NBNDEL = MCN( MNELE + WBNDEL )
         NBPGEL = MCN( MNELE + WBPGEL )
C        L'ADRESSE DE NUNDEL PUIS NUPGEL
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .LT. 2 .OR. NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 )THEN
C           NOEUDS=POINTS
            NBPGEL = NBNDEL
         ELSE
C           NOEUDS#POINTS
            NBPGEL = MCN( MNELE + WBPGEL )
            MNPGEL = MNPGEL + NBELEM * NBNDEL
         ENDIF
C
C        ATTENTION : RECUL D'UN MOT POUR EVITER LA SOUSTRACTION DANS LES BOUCLES
         MNPGEL = MNPGEL - 1
C
C        LES SOMMETS = POINTS UTILISATEUR
         NBPSEL = MCN( MNELE + WBPSEL )
         MOPSEL = MCN( MNELE + WOPSEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLPSEL
         MNPSEL = MNPGEL + NBELEM * NBPGEL + NBPSEL
         DO 120 J = 1,NBPSEL
C           LE NUMERO DU POINT
            MCN( MNS + 2 ) = MCN( MNPSEL - NBPSEL + J )
C           LE NUMERO LOCAL DU SOMMET DANS L'EF
            N   = MCN( MNPSEL + J )
C           LE NUMERO DE L'EF DE CE SOMMET
            NEF = MCN( MNPSEL + MOPSEL + J )
C           LE NUMERO XYZPOI DE CE SOMMET
            NS  = MCN( MNPGEL + NEF + NBELEM * ( N - 1 ) )
C           CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
            MCN( MNS + 1 ) = NS
            MNS = MNS + 2
 120     CONTINUE
C
C        LES ARETES APPARTENANT A UNE LIGNE
         NBLAEL = MCN( MNELE + WBLAEL )
         MOLAEL = MCN( MNELE + WOLAEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLLAEL
         MNLAEL = MNPSEL + MOPSEL + MOPSEL + NBLAEL
         DO 130 J = 1,NBLAEL
C           LE NUMERO DE LA LIGNE
            MCN( MNA + 3 ) = MCN( MNLAEL - NBLAEL + J )
C           LE NUMERO LOCAL DE L'ARETE DANS L'EF
            N   = MCN( MNLAEL + J )
C           LE NUMERO DE L'EF DE CETTE ARETE
            NEF = MCN( MNLAEL + MOLAEL + J )
C           LE NUMERO XYZPOI DE SES 2 SOMMETS
            DO 125 K=1,2
C              LE NUMERO LOCAL DU SOMMET DANS L'EF
               NS = NOSOAR( K, N )
               NS = MCN( MNPGEL + NEF + NBELEM * ( NS - 1 ) )
C              CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
               MCN( MNA + K ) = NS
 125        CONTINUE
            MNA = MNA + 3
 130     CONTINUE
C
C        LES FACES APPARTENANT A UNE SURFACE
         NBSFEL = MCN( MNELE + WBSFEL )
         MOSFEL = MCN( MNELE + WOSFEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLSFEL
         MNSFEL = MNLAEL + MOLAEL + MOLAEL + NBSFEL
         DO 140 J = 1,NBSFEL
C           LE NUMERO DE LA SURFACE
            MCN( MNF + 5 ) = MCN( MNSFEL - NBSFEL + J )
C           LE NUMERO LOCAL DE LA FACE DANS L'EF
            N   = MCN( MNSFEL + J )
C           LE NUMERO DE L'EF DE CETTE FACE
            NEF = MCN( MNSFEL + MOSFEL + J )
C           LE NOMBRE DE SOMMETS DE CETTE FACE
            NBSTFA = NBSOFA(N)
C           TEMOIN DE TRIANGLE, ECRASE SI CETTE FACE A 4 SOMMETS
            MCN( MNF + 4 ) = 0
C           LE NUMERO XYZPOI DE SES NBSTFA SOMMETS
            DO 135 K=1,NBSTFA
C              LE NUMERO LOCAL DU SOMMET DANS L'EF
               NS = NOSOFA( K, N )
C              LE NUMERO XYZPOI DE CE SOMMET
               NS = MCN( MNPGEL + NEF + NBELEM * ( NS - 1 ) )
C              CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
               MCN( MNF + K ) = NS
 135        CONTINUE
            MNF = MNF + 5
 140     CONTINUE
 200  CONTINUE

C     LE TABLEAU DES Z AXONOMETRIQUES DES BARYCENTRES DES FACES ARETES ET SOMMET
      NBBARY = NBFA + NBAR + NBSO + NBARFR
ccc      PRINT*,'t1obje:',NBBARY,' FACES+ARETES+SOMMETS pour TRACER ',
ccc     %                 NOMOBJ
      CALL TNMCMX( 'ENTIER', MAXVAR )
      IF( MAXVAR .LT. NBBARY ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'t1obje: SUPER-TABLEAU MCN TROP COURT MOTMCN=',MOTMCN
         ELSE
            PRINT*,'t1obje: SUPER-ARRAY MCN TOO SHORT MOTMCN=',MOTMCN
         ENDIF
         GOTO 9000
      ENDIF
      CALL TNMCDC( 'REEL',   NBBARY, MNBARY )
      CALL TNMCDC( 'ENTIER', NBBARY, MNNUBA )

      IF( LORBITE .GT. 0 ) THEN
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
      ENDIF

C     FORMATION TABLEAUX Z-AXONOMETRIQUES ET NUMERO DES BARYCENTRES AVANT TRI
C     -----------------------------------------------------------------------
 219   CALL TBAFAS( NBCOOR, NBPOMA, RMCN(MNXYZP+WYZPOI),
     %              NBSO,   MCN(MNSOPO),
     %              NBAR,   MCN(MNARLI),
     %              NBFA,   MCN(MNFASU),
     %              MOARFR, MXARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %              NBBARY, MCN(MNNUBA), RMCN(MNBARY),
     %              NUMXPO, NUMXLI, NUMXSU )

C     LE TRI PAR TAS DE CETTE COTE AXONOMETRIQUE
      CALL TRITRP( NBBARY, RMCN(MNBARY), MCN(MNNUBA) )

C     TRACE EFFECTIF DES ARETES FRONTALIERES SURFACES LIGNES ET POINTS
C     ----------------------------------------------------------------
      MXPLS = NUMXPO + NUMXLI + NUMXSU
      IF( MXPLS .GT. 0 ) THEN
         PRINT*,'t1obje:',NBBARY,
     %          ' FACES+ARETES+SOMMETS TRACENT ', NOMOBJ
         CALL TNMCMX( 'ENTIER', MAXVAR )
         IF( MAXVAR .LT. NBBARY ) THEN
            IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'t1obje: SUPER-TABLEAU MCN TROP COURT MOTMCN=',MOTMCN
            ELSE
               PRINT*,'t1obje: SUPER-ARRAY MCN TOO SHORT MOTMCN=',MOTMCN
            ENDIF
            GOTO 9000
         ENDIF
         MNPLS = 0
         CALL TNMCDC( 'ENTIER', MXPLS, MNPLS )
         CALL AZEROI( MXPLS, MCN(MNPLS) )
      ELSE
C        ADRESSE BIDON NON UTILISEE ENSUITE
         MNPLS = MNNUBA
      ENDIF
      CALL T3AFAS( NBCOOR, RMCN(MNXYZP+WYZPOI),
     %             NBSO,   MCN(MNSOPO),
     %             NBAR,   MCN(MNARLI),
     %             NBFA,   MCN(MNFASU),
     %             MOARFR, MXARFR, MCN(MNAFOB+WAREFR),
     %             NBBARY, MCN(MNNUBA),
     %             NUMXPO, NUMXLI, NUMXSU, MCN(MNPLS) )

C     TRACE DES POINTS DES EF DU MAILLAGE DE L'OBJET 2D AVEC INTERPOLATION
      IF( NTTOPO .GT. 0 .AND. IAVNSO0 .GT. 0 ) THEN
         MNX = MNXYZP + WYZPOI
         DO 230 NP=1,NBPOMA
            WRITE( NMSOMM, '(I7)' ) NP
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            CALL SYMBOLE3D( NCONSO, RMCN(MNX), KSYMBO(1:L) )
            MNX = MNX + NBCOOR
 230     CONTINUE
      ENDIF

C     TRACE ET AFFICHAGE DE LA QUALITE DU MAILLAGE DE L'OBJET
      CALL IMPQUAOB( NOMOBJ, QUAMIN, NBEFMQ, IERR )

C     TRACE DU TITRE ET FERMETURE
      CALL TRFINS( NOMOBJ )

C     REPRISE DE L'ORBITE SI ACTIVE
C     =============================
      IF( LORBITE .GT. 0 ) THEN
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         GOTO 219
      ENDIF

C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
 9000 IF( MNNUBA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBBARY,   MNNUBA )
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBBARY,   MNBARY )
      IF( MNPLS  .GT. 0 ) CALL TNMCDS( 'ENTIER', MXPLS,    MNPLS  )
      IF( MNSOPO .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*NBSO,   MNSOPO )
      IF( MNARLI .GT. 0 ) CALL TNMCDS( 'ENTIER', 3*NBAR,   MNARLI )
      IF( MNFASU .GT. 0 ) CALL TNMCDS( 'ENTIER', 5*NBFA,   MNFASU )
      IF( MNELEM .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNELEM )
      GOTO 99999

C     OBJET INCONNU => ERREUR AVEC LISTE DES OBJETS ACTUELS
C     =====================================================
 9999 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'OBJET INCONNU:' // KNOMOB
         KERR(2) = 'A choisir parmi :'
         WRITE(IMPRIM,19999) 'OBJET',KNOMOB
      ELSE
         KERR(1) = 'UNKNOWN OBJECT:' // KNOMOB
         KERR(2) = 'CHOOSE AMONG:'
         WRITE(IMPRIM,29999) 'OBJECT',KNOMOB
      ENDIF
      CALL LEREUR
19999 FORMAT(1X,A,1X,A,' NON RETROUVE PARMI')
29999 FORMAT(1X,A,1X,A,' NOT FOUND AMONG')
C     OUVERTURE DES OBJETS
      CALL LXLXOU( NTADAM, 'OBJET', NTOBJE, MNOBJE )
      CALL LXIM0( MNOBJE )
C
99999 IAVNSO = IAVNSO0
      PREDUF = PREDUF0
      RETURN
      END
