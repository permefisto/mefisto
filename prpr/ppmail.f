      PROGRAM PPMAIL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PAR BATCH ou en INTERACTIF ECRAN CLAVIER SOURIS
C -----    MAILLER DES POINTS, LIGNES, SURFACES, VOLUMES, OBJETS
C          DEFINIR L'INTERPOLATION DE CHACUN DES OBJETS
C
C          PPINIT A DU ETRE EXECUTE AUPARAVANT POUR DEFINIR LES FICHIERS
C          DE LA MEMOIRE SECONDAIRE ET CONSTRUIRE LES TABLES DE GESTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     Janvier 1997
C MODIF  : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray    Avril 2013
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray             Avril 2020
C2345X7..............................................................012
C$    USE OMP_LIB
      include"./incl/threads.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/majmai.inc"
      include"./incl/pp.inc"
      include"./incl/ppmck.inc"
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/mecoit.inc"
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/darete.inc"
      include"./incl/traaxe.inc"
C
C     DECLARATION DU SUPER-TABLEAU NUMERIQUE MCN
      REAL               RMCN(MOTMCN)
      DOUBLE PRECISION   DMCN(MOTMCN/2)
      COMMON              MCN(MOTMCN)
      EQUIVALENCE        (MCN(1),RMCN(1),DMCN(1))
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
C     UNE VARIABLE LOGIQUE LOCALE
      LOGICAL      EXIST
      CHARACTER*24 KNOMOB

C     LES ARGUMENTS A L'APPEL DU LOAD MODULE ppmail ET AUXILIAIRES
      INTEGER          EXISTNF,  NBARGS
      CHARACTER*128    NMFILEDO, ARGUMENT, KCHAIN
      DOUBLE PRECISION DBLVAL

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
C
C     NOM DU REPERTOIRE DE TRAVAIL
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Nom du REPERTOIRE de TRAVAIL :'
      ELSE
         PRINT *,'WORKING DIRECTORY NAME :'
      ENDIF
      CALL SYSTEM( 'echo `pwd` ' )
C
C     ESSAI DE RECUPERATION DU NOM DU PROJET
C     ======================================
      CALL NOMJOB( EXIST )
      IF( .NOT. EXIST ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: EXECUTEZ INITIER AUPARAVANT'
         ELSE
            KERR(1) ='ERROR: EXECUTE INITIER BEFORE'
         ENDIF
         CALL LEREUR
         STOP
      ENDIF
C
C     NOMBRE DES ARGUMENTS DE LA COMMANDE D'EXECUTION DU LOAD-MODULE
C     ==============================================================
      NBARGS = IARGC()
      IF( LANGAG .EQ. 0 ) THEN
         print *,'Mefisto-MAILLER: NOMBRE DES ARGUMENTS=', NBARGS+1
      ELSE
         print *,'Mefisto-MESHER: ARGUMENT NUMBER=',NBARGS+1
      ENDIF
      DO N=0,NBARGS
C        EXPLORATION DES ARGUMENTS de Mefisto-MAILLER
         CALL GETARG( N, ARGUMENT )
         print *,'Mefisto-MAILLER: ARGUMENT',N,' = ',ARGUMENT
      ENDDO
C     RECHERCHE D'UN EVENTUEL NOM DE FICHIERS DE DONNEES POUR Mefisto
      CALL NMFIDOCO( 'ppmail', EXISTNF, NMFILEDO )

C     LE MODE DES INTER-ACTIONS ENTRE LE PROGRAMME ET L'UTILISATEUR
C     =============================================================
C     INTERA=0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS

      IF( EXISTNF .NE. 0 ) THEN

C        UN NOM DE FICHIER EXISTE => BATCH SANS X11 SANS CLAVIER SANS SOURIS
C        INITIALISATIONS POUR EVITER LES TRACES X11
         INTERA = 0
         LAPXFE = 800
         LHPXFE = 600
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'Mefisto-MAILLER: Nom Fichier Donnees=',NMFILEDO
         WRITE(IMPRIM,*)'Mefisto-MAILLER s''execute en mode BATCH'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-MESHER: Data File Name=',NMFILEDO
            WRITE(IMPRIM,*) 'Mefisto-MESHER in BATCH EXECUTION'
         ENDIF

      ELSE

C        PAS DE NOM DE FICHIER DE DONNEES => INTERACTIF AVEC X11
         INTERA = IINFO( 'INTERACTIVITE INITIALE' )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Mefisto-MAILLER: PAS de Fichier Donnees'
            WRITE(IMPRIM,*) 'Mefisto-MAILLER: EXECUTION INTERACTIVE X11'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-MESHER: NO MESHER DATA FILE'
            WRITE(IMPRIM,*) 'Mefisto-MESHER: X11 INTERACTIVE EXECUTION'
         ENDIF

      ENDIF
      LHLECT = 1
      INTERB(LHLECT) = INTERA

C     PAR DEFAUT, PAS DE MISE A JOUR AUTOMATIQUE DES MAILLAGES
C     LORS DE L'EXECUTION DE LA SUBROUTINE MAILEX
      MAJAUT = 0
C     TRACE DES AXES REINITIALISE
      NETAXE = 0
      NOTYVI = 0
      COOEXT(1,1) = RINFO('GRAND')

C     INITIALISATIONS
C     ===============
C     OUVERTURE DU PROJET
      CALL INITIA( EXIST )

C     OUVERTURE DE X-WINDOW X11
      IF( INTERA .GE. 1 ) THEN
         CALL XTINIT
C        OUVERTURE DU SERVEUR ET DE LA FENETRE X11 DE NOM Mefisto
         CALL XVINIT
C        TRACE DU LOGO DE Mefisto-MAILLAGES
         IF( LANGAG .EQ. 0 ) THEN
            CALL LOGO( 'Mefisto-MAILLER' )
         ELSE
            CALL LOGO( 'Mefisto-MESHER' )
         ENDIF
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
10010 FORMAT('SUPER-TABLEAU MCN DE',I12,' MOTS-MEMOIRE'/)
20010 FORMAT('SUPER-ARRAY MCN WITH',I12,' WORDS of MEMORY'/)

C     RECUPERATION DU DERNIER OBJET DECLARE
C     -------------------------------------
      CALL REDEOB( KNOMOB, NUOB, MN )
      IF( NUOB .GT. 0 ) GOTO 22

C     SI PAS d'OBJET TRACE DES POINTS ET DES LIGNES POUR VOIR CE QUI EXISTE
C     ---------------------------------------------------------------------
      IF( INTERA .GE. 1 ) THEN
         CALL VISEE0
         CALL EFFACE
         CALL TRTOBJ( 'LIGNE', 0 )
         CALL TRTOBJ( 'POINT', 0 )
         IF( LANGAG .EQ. 0 ) THEN
            CALL TRFINS( 'POINTS & LIGNES' )
         ELSE
            CALL TRFINS( 'POINTS & LINES' )
         ENDIF
      ENDIF

C     EXISTE-T-IL UN FICHIER DE DONNEES POUR MAILLER?
C     -----------------------------------------------
 22   IF( EXISTNF .NE. 0 ) THEN

C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
C        CONSTRUCTION DE LA DONNEE  "READF NMFILEDO ;" dans KLG
         N      = NUDCNB( NMFILEDO )
         KLG(1) = 'READF ' // NMFILEDO(1:N) // ' ; '
         LHKLG  = 1
         NLPTV1 = 1
         NCPTV1 = 0
         NOTYPE = 0
         NOTYPE =-2
         CALL DONNMF( NOTYPE, NOTYPS, NCVALS, DBLVAL, KCHAIN )

      ENDIF

C     ==============================================================
C     SAISIE INTERACTIVE OU SUR FICHIER DES DONNEES DE L'UTILISATEUR
C     ==============================================================
C     HISTORIQUE REMIS A ZERO
 90   NBLGRC( NRHIST ) = 0

C     TRACE DES AXES REINITIALISE
      NETAXE = 0

C     LECTURE DU MOT-CLE OU OPTION A EXECUTER PARMI TOUTES LES OPTIONS
C     CONTENUES DANS LE MENU DE NOM 'debut'
C     ================================================================
      CALL LIMTCL( 'debut', NMTCL )

C     TRAITEMENT DE L'OPTION NMTCL
      IF( NMTCL .LT. 0 ) GOTO 90

      IF( NMTCL .EQ. 0 ) THEN
C       'La LONGUEUR SOUHAITEE par DEFAUT des ARETES des MAILLAGES'
C        ---------------------------------------------------------
 95      CALL INVITE( 152 )
         NCVALS = 6
         CALL LIRRDP( NCVALS, DARETE )
         IF( NCVALS .LT. 0   ) GOTO 90
         IF( DARETE .LE. 0D0 ) GOTO 95
C        ATTENTION: DARETE DOIT ETRE >0
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10095) DARETE
         ELSE
            WRITE(IMPRIM,20095) DARETE
         ENDIF
         GOTO 90
      ENDIF
10095 FORMAT(G14.6,'=LONGUEUR SOUHAITEE des ARETES des MAILLAGES' )
20095 FORMAT(G14.6,'=WISHED LENGTH of MESH EDGES' )

      IF( NMTCL .EQ. 60 ) GOTO 6000
      IF( NMTCL .EQ. 70 ) GOTO 7000
      IF( NMTCL .EQ. 80 ) GOTO 8000
      IF( NMTCL .EQ. 90 ) GOTO 9000
      IF( NMTCL .EQ. 98 ) GOTO 9800
      IF( NMTCL .EQ. 99 ) GOTO 9900
      GOTO( 100, 200, 300, 400, 500, 600, 700,  90,   90, 1000,
     %     1100,  90,  90,  90,  90,  90,  90,  90, 1900, 2000,
     %     2100,  90,  90,  90,  90,  90,  90,  90,   90,   90), NMTCL

C     DEFINITION DES POINTS
 100  CALL DFPOIN
      GOTO 90

C     DEFINITION DES LIGNES
 200  CALL DFLIGN
      GOTO 90

C     DEFINITION DES SURFACES
 300  CALL DFSURF
      GOTO 90

C     DEFINITION DES VOLUMES
 400  CALL DFVOLU
      GOTO 90

C     DEFINITION DES OBJETS
 500  CALL DFOBJE
      GOTO 90

C     DEFINITION DES TRANSFORMATIONS MATHEMATIQUES R3->R3
 600  CALL DFTRAN
      GOTO 90

C     INTERPOLATION ET RENUMEROTATION DES NOEUDS DES ELEMENTS FINIS
 700  CALL DFTOPO
      GOTO 90

C     TRACE DES MAILLAGES
 1000 CALL TRMAIL
      GOTO 90

C     RENOMMER LA FONCTION TAILLE_IDEALE(x,y,z)
1100  CALL RNTAIL
      GOTO 90

C     CADRE EXTREME DES COORDONNEES DES OBJETS
1900  CALL COEXOB
      CALL EFFACE
      CALL ITEMS0
      CALL VISEE0
      CALL TRAXES
CCC      CALL TRTOBJ( 'POINT', 0 )
CCC      CALL TRTOBJ( 'LIGNE', 0 )
      IF( LANGAG .EQ. 0 ) THEN
         CALL TRFINS( 'AXES' )
      ELSE
         CALL TRFINS( 'AXIS' )
      ENDIF
      GOTO 90

C     PRECISIONS pour IDENTIFIER 2 POINTS ou SOMMETS
 2000 CALL ZEROS
      GOTO 90

C     MISE A JOUR AUTOMATIQUE OU NON DES MAILLAGES LORS DE
C     L'EXECUTION DE MAILEX
 2100 IF( MAJAUT .EQ. 0 ) THEN
         MAJAUT = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='MISE A JOUR AUTOMATIQUE DES MAILLAGES'
         ELSE
            KERR(1) ='AUTOMATIC UPDATE of MESHES'
         ENDIF
         CALL LERESU
      ELSE
         MAJAUT = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='PAS DE MISE A JOUR AUTOMATIQUE DES MAILLAGES'
         ELSE
            KERR(1) ='NO AUTOMATIC UPDATE of MESHES'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 90

C     MANAGEMENT de la FENETRE GRAPHIQUE de Mefisto
C     LARGEUR HAUTEUR en PIXELS COULEUR du FOND
 6000 CALL MANAGFEN
      GOTO 90

C     LA GESTION DES RESSOURCES DE Mefisto
C     MANAGEMENT des TMS Files Unites de Mefisto
 7000 CALL MANAGMEF
      GOTO 90

C     IMPORTER un MAILLAGE de PLSVO
 8000 CALL IMPORTER
      GOTO 90

C     EXPORTER un MAILLAGE de PLSVO
 9000 CALL EXPORTER
      GOTO 90

C     NOM DE LA VERSION de Mefisto
 9800 NBLGRC(NRERR) = 1
      KERR(1) = ' '
      CALL VRSION( KERR(1) )
      CALL LERESU
      GOTO 90

C     DESTRUCTION DES TABLEAUX DES ITEMS A TRACER
 9900 CALL ITEMDS

C     SAUVEGARDE DES DONNEES SUR LA MS ET FIN DU TRAITEMENT
      CALL ARRET(  0 )

      STOP
      END
