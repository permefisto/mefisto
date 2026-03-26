      PROGRAM PPFLUI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LA MEMOIRE SECONDAIRE
C ----- TRAITER UN PROBLEME DE MECANIQUE D'UN FLUIDE INCOMPRESSIBLE
C       FERMER LA MEMOIRE SECONDAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: BENHAMADOUCHE-BOYER-PERRONNET-ZACHARIE           Janvier 2000
C MODIFS: ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris      Juin 2007
C MODIFS: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2008
C MODIFS: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Avril 2010
C MODIFS: ALAIN PERRONNET TEXAS A & M University at QATAR   Fevrier 2012
C MODIFS: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray       Mai 2012
C MODIFS: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray   Fevrier 2013
C MODIFS: ALAIN PERRONNET             St Pierre du Perray      Aout 2020
C MODIFS: ALAIN PERRONNET             St Pierre du Perray     Avril 2022
C2345X7..............................................................012
C$    USE OMP_LIB 
      include"./incl/threads.inc"
      include"./incl/nmproj.inc"
      include"./incl/pp.inc"
      include"./incl/lu.inc"
      include"./incl/ppmck.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/typnoobj.inc"
      include"./incl/mecoit.inc"
      include"./incl/trvari.inc"
      include"./incl/traaxe.inc"

      REAL              RMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      COMMON             MCN(MOTMCN)
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      LOGICAL            EXIST
      CHARACTER*24       KNOMOB
      DOUBLE PRECISION   DATE0, DATE, SECONDES

C     LES ARGUMENTS A L'APPEL DU LOAD MODULE ppflui
      INTEGER            EXISTNF, NBARGS
      CHARACTER*128      NMFILEDO, ARGUMENT


C///////////////////////////////////////////////////////////////////////
      NBTHREADS = 1
C$OMP PARALLEL
C$OMP MASTER
C$    NBTHREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

C     RECONNAISSANCE DE LA LANGUE DES DONNEES DE Mefisto
      CALL LANGUE

C     NOM DU REPERTOIRE DE TRAVAIL
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Nom du REPERTOIRE de TRAVAIL :'
      ELSE
         PRINT *,'WORKING DIRECTORY NAME :'
      ENDIF
      CALL SYSTEM( 'echo `pwd` ' )

C     ESSAI DE RECUPERATION DU NOM DU PROJET
C     ======================================
      CALL NOMJOB( EXIST )
      IF( .NOT. EXIST ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR ppflui: EXECUTEZ INITIER AUPARAVANT'
         ELSE
            KERR(1) ='ppflui ERROR: EXECUTE INITIER BEFORE'
         ENDIF
         CALL LEREUR
         STOP
      ENDIF

C     NOMBRE DES ARGUMENTS DE LA COMMANDE D'EXECUTION DU LOAD-MODULE
C     ==============================================================
      NBARGS = IARGC()
      IF( LANGAG .EQ. 0 ) THEN
         print *,'Mefisto-FLUIDER: NOMBRE DES ARGUMENTS=',NBARGS+1
      ELSE
         print *,'Mefisto-FLUIDER: ARGUMENT NUMBER=',NBARGS+1
      ENDIF
      DO N=0,NBARGS
C        EXPLORATION DES ARGUMENTS
         CALL GETARG( N, ARGUMENT )
         print *,'Mefisto-FLUIDER: ARGUMENT',N,' = ',ARGUMENT
      ENDDO

C     LA COMMANDE Mefisto-FLUIDER est elle suivie d'un nom
C     de fichier de donnees MEFISTO SELON le LANGAGE LU?
C     ====================================================
      CALL NMFIDOCO( 'ppflui', EXISTNF, NMFILEDO )

C     NOMBRE DE CARATERES NON BLANCS DE NMFILEDO
      N = NUDCNB( NMFILEDO )

C     LE MODE DES INTER-ACTIONS ENTRE LE PROGRAMME ET L'UTILISATEUR
C     INTERA=0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS

      IF( EXISTNF .EQ. 0 ) THEN
C        PAS DE NOM DE FICHIER DE DONNEES => INTERACTIF AVEC X11
         INTERA = IINFO( 'INTERACTIVITE INITIALE' )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER EN EXECUTION INTERACTIVE'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER in INTERACTIVE EXECUTION'
         ENDIF
      ELSE
C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
         INTERA = 0
         LAPXFE = 800
         LHPXFE = 600
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER: Nom Fichier Donnees=',
     %                       NMFILEDO(1:N)
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER s''execute en mode BATCH'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER: Data File Name=',
     %                       NMFILEDO(1:N)
            WRITE(IMPRIM,*) 'Mefisto-FLUIDER in BATCH EXECUTION'
         ENDIF
      ENDIF
      LHLECT = 1
      INTERB(LHLECT) = INTERA
      NBPAS  = 0

C     INITIALISATIONS
C     ===============
      CALL INITIA( EXIST )

C     OUVERTURE DE X-WINDOW X11
      IF( INTERA .GE. 1 ) THEN
         CALL XTINIT
C        OUVERTURE DU SERVEUR ET DE LA FENETRE X11 DE NOM Mefisto
         CALL XVINIT
C        TRACE DU LOGO DE Mefisto
         CALL LOGO( 'Mefisto-FLUIDER' )
      ENDIF

C     OUVERTURE DE LA MS
      CALL OUVRMS

C     LE LANGAGE UTILISATEUR RETROUVE SES VARIABLES
      CALL LUOU

C     AFFICHAGE DE LA TAILLE DU SUPER-TABLEAU NUMERIQUE
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) MOTMCN
      ELSE
         WRITE(IMPRIM,20010) MOTMCN
      ENDIF
10010 FORMAT('SUPER-TABLEAU MCN DE',I12,' MOTS-MEMOIRE')
20010 FORMAT('SUPER-ARRAY MCN WITH',I12,' WORDS of MEMORY')

C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
      CALL SECONDES1970( DATE0 )

C     RECUPERATION DU DERNIER OBJET DECLARE
C     -------------------------------------
 20   CALL REDEOB( KNOMOB, NUOB, MN )

      IF( NUOB .LE. 0 ) THEN
C        PAS D'OBJET RETROUVE => ERREUR
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: PAS D''OBJET MAILLE'
            KERR(2) = 'CREER un OBJET avec Mefisto-MAILLER'
         ELSE
            KERR(1) = 'ERROR: NO OBJECT'
            KERR(2) = 'CREATE an OBJECT with Mefisto-MESHER'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      NBPAS = NUOB

C     EXISTE-T-IL UN FICHIER DE DONNEES POUR FLUIDER?
      IF( EXISTNF .NE. 0 ) THEN

C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
         N = NUDCNB( NMFILEDO )
         KLG(1) = 'READF ' // NMFILEDO(1:N) // ' ; '
         LHKLG  = 1
         NLPTV1 = 1
         NCPTV1 = 0
         NOTYPE =-2
         CALL DONNMF( NOTYPE, NOTYPS, NCVALS, DBLVAL, KCHAIN )

      ENDIF

C     LECTURE DES MOTS CLE
C     ====================
C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
 30   CALL SECONDES1970( DATE )
      SECONDES = DATE - DATE0
      DATE0    = DATE
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10030) SECONDES
10030 FORMAT(' ppflui: Temps REEL derniere commande=',G15.6,' secondes')
      ELSE
         WRITE(IMPRIM,20030) SECONDES
20030    FORMAT(' ppflui: Last Command REAL Time=',G15.6,' seconds')
      ENDIF

C     L'ANCIEN HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C     TRACE DES AXES REINITIALISE
      NETAXE = 0

C     ==============================================================
C     SAISIE INTERACTIVE OU SUR FICHIER DES DONNEES DE L'UTILISATEUR
C     ==============================================================
      CALL LIMTCL( 'debuflui', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 30
      IF( NMTCL .EQ. 19 ) GOTO 1900
      IF( NMTCL .EQ. 20 ) GOTO 2000
      IF( NMTCL .EQ. 35 ) GOTO 3500
      IF( NMTCL .EQ. 36 ) GOTO 3600
      IF( NMTCL .EQ. 60 ) GOTO 6000
      IF( NMTCL .EQ. 70 ) GOTO 7000
      IF( NMTCL .EQ. 97 ) GOTO 9700
      IF( NMTCL .EQ. 99 ) GOTO 9900
      GOTO( 100, 200, 300, 30, 500, 600, 700, 800, 900, 1000), NMTCL

C     NOM DE L'OBJET A TRAITER
C     (NECESSAIRE SEULEMENT SI LE PROJET CONTIENT PLUSIEURS OBJETS)
 100  IF( NBPAS .GT. 0 ) THEN
C        FERMETURE DE L'OBJET PRECEDENT
         CALL LXLXFE( NTOBJE, KNOMOB )
         NBPAS = 0
C        HISTORIQUE REMIS A ZERO
         CALL RECTEF( NRHIST )
         NBLGRC( NRHIST ) = 0
      ENDIF

      CALL INVITE( 45 )
      NCVALS = 0
      CALL LIRLEX( NTOBJE, NCVALS, KNOMOB, NUOB )
      IF( NCVALS .EQ. -1 ) GOTO 20
      IF( NUOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: OBJET INCONNU:' // KNOMOB
         ELSE
            KERR(1) ='ERROR: UNKNOWN OBJECT:' // KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 20
      ENDIF

C     REMPLISSAGE DES COMMONS de ./incl/typnoobj.inc
      NUMTYPOBJ = 5
      NUMOBJLX  = NUOB
      IF( LANGAG .EQ. 0 ) THEN
         KNMTYPOBJ = 'Objet'
      ELSE
         KNMTYPOBJ = 'Object'
      ENDIF
      KNMOBJLX  = KNOMOB

C     MISE A JOUR DES XYZ MIN MAX DANS COOEXT de xyzext.inc
      CALL LXNLOU( NTOBJE, NUOB, NTLXOB, MNLXOB )
      CALL LXTSOU( NTLXOB, 'XYZPOINT', NTXYZP, MNXYZP )
      IF( NTXYZP .LE. 0 ) THEN
         CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZP, MNXYZP )
         IF( NTXYZP .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET ' // KNOMOB
               KERR(2) = 'SANS TMS XYZSOMMET et XYZPOINT'
            ELSE
               KERR(1) = 'ERROR: OBJECT ' // KNOMOB
               KERR(2) = 'WITHOUT TMS XYZSOMMET and XYZPOINT'
            ENDIF
            CALL LEREUR
            GOTO 20
         ENDIF
      ENDIF

C     MISE A JOUR DU CADRE MINMAX DES NBCOOR COORDONNEES DES POINTS
      INIEXT = 0
      CALL MAJEXT( MNXYZP )

C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO

      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL VISEE1
         CALL T1OBJE( KNOMOB )
      ENDIF
      NBPAS = 1
      GOTO 30


C     DONNEES des CARACTERISTIQUES PHYSIQUES et CONDITIONS aux LIMITES du FLUIDE
 200  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL DFFLUI( KNOMOB, IERR )
      GOTO 30


C     DONNEES des CARACTERISTIQUES THERMIQUES et CONDITIONS aux LIMITES du FLUIDE
 300  IF( NBPAS .EQ. 0 ) GOTO 100
      NBJEUX = 1
      CALL DFTHER( KNOMOB, NBJEUX, IERR )
      GOTO 30

C     RESOLUTIONS d'un PROBLEME d'un FLUIDE INCOMPRESSIBLE
C     5; STOKES_STATIONNAIRE
C     5:'STEADY STOKES SOLVER',
 500  IF( NBPAS .EQ. 0 ) GOTO 100
      NAVSTO = -1
      CALL STOKESTA( KNOMOB, IERR )
      IF( IERR .NE. 0 ) GOTO 9600
      GOTO 30


C     6; STOKES_INSTATIONNAIRE
C     6:'UNSTEADY STOKES SOLVER',
 600  IF( NBPAS .EQ. 0 ) GOTO 100
      NAVSTO = 0
      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
      IF( IERR .NE. 0 ) GOTO 9600
      GOTO 30


C     7; NAVIER-STOKES_INSTATIONNAIRE EULER IMPLICITE Algorithme1
C     7:'IMPLICIT UNSTEADY NAVIER-STOKES SOLVER Algorithm1'
 700  IF( NBPAS .EQ. 0 ) GOTO 100
      NAVSTO = 1
      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
      IF( IERR .NE. 0 ) GOTO 9600
      GOTO 30


C     8; NAVIER-STOKES_INSTATIONNAIRE PAR PAS FRACTIONNAIRES Algorithme2
C     8:'FRACTIONAL UNSTEADY NAVIER-STOKES SOLVER Algorithm2'
 800  IF( NBPAS .EQ. 0 ) GOTO 100
      NAVSTO = 2
      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
      IF( IERR .NE. 0 ) GOTO 9600
      GOTO 30


C     9; NAVIER-STOKES + THERMIQUE INSTATIONNAIRE PAS FRACTIONNAIRES
C     9; UNSTEADY HEAT NAVIER-STOKES SOLVER
 900  IF( NBPAS .EQ. 0 ) GOTO 100
      NAVSTO = 3
      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
      IF( IERR .NE. 0 ) GOTO 9600
      GOTO 30


cccC     7; NAVIER-STOKES_INSTATIONNAIRE PAR SCHEMA PISO+CARACTERISTIQUE Algorithme3
cccC     7:'PISO+Charac UNSTEADY NAVIER-STOKES SOLVER Algorithm3'
ccc 700  IF( NBPAS .EQ. 0 ) GOTO 100
ccc      NAVSTO = 4
ccc      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
ccc      IF( IERR .NE. 0 ) GOTO 9600
ccc      GOTO 30


cccC     8; NAVIER-STOKES_INSTATIONNAIRE PAR SCHEMA PISO+(V.D)V Algorithme4
cccC     8:'PISO+(V.D)V UNSTEADY NAVIER-STOKES SOLVER Algorithm4',
ccc 800  IF( NBPAS .EQ. 0 ) GOTO 100
ccc      NAVSTO = 5
ccc      CALL FLUIDENS( KNOMOB, NAVSTO, IERR )
ccc      IF( IERR .NE. 0 ) GOTO 9600
ccc      GOTO 30


C     DESSIN DES VITESSES ET/OU PRESSIONS du FLUIDE INCOMPRESSIBLE
 1000 IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TFLUIDE( KNOMOB, IERR )
      GOTO 30

C     TRACE DES MAILLAGES
 1900 CALL TRMAIL
      GOTO 30

C     PRECISION pour INVERSER A x = b
 2000 CALL ZEROGC
      GOTO 30

C     TABLE des RESULTATS 2D du TOURBILLON de TAYLOR-GREEN
 3500 CALL NS2DTGV
      GOTO 30

C     TABLE des RESULTATS 3D du TOURBILLON de TAYLOR-GREEN
 3600 CALL NS3DTGV
      GOTO 30

C     MANAGEMENT de la FENETRE GRAPHIQUE de Mefisto
C     LARGEUR HAUTEUR en PIXELS COULEUR du FOND
 6000 CALL MANAGFEN
      GOTO 30

C     LA GESTION DES RESSOURCES DE Mefisto
C     MANAGEMENT des TMS Files Unites de Mefisto
 7000 CALL MANAGMEF
      GOTO 30

C     ERREUR DANS LE CALCUL:
 9600 IF( INTERA .GE. 3 ) THEN
C         INTERACTIVITE AVEC CLAVIER, TRACES GRAPHIQUES et SOURIS
          LHLECT = 1
          INTERB(LHLECT) = INTERA
          LECTEU = IINFO( 'LECTEUR INITIAL' )
      ELSE
C         SAUVEGARDE DES DONNEES SUR LA MS ET FIN DU TRAITEMENT
          CALL ARRET( 100 )
      ENDIF
      GOTO 30

C     NOM de la VERSION DE Mefisto
 9700 NBLGRC(NRERR) = 1
      KERR(1) = ' '
      CALL VRSION( KERR(1) )
      CALL LERESU
      GOTO 30

C     SAUVEGARDE DES DONNEES SUR LA MS ET FIN DU TRAITEMENT
 9900 CALL ARRET(  0 )

      STOP
      END
