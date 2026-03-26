      PROGRAM PPELAS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LA MEMOIRE SECONDAIRE
C ----- TRAITER EVENTUELLEMENT UN PROBLEME D'ELASTICITE
C       FERMER LA MEMOIRE SECONDAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     Janvier 1997
C MODIFS : ALAIN PERRONNET St Pierre du Perray & LJLL UPMC  Mars    2013
C2345X7..............................................................012
C$    USE OMP_LIB
      include"./incl/threads.inc"
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

C     LES ARGUMENTS A L'APPEL DU LOAD MODULE ppelas
      INTEGER        EXISTNF, NBARGS
      CHARACTER*128  NMFILEDO, ARGUMENT

C///////////////////////////////////////////////////////////////////////
      NBTHREADS = 1
C$OMP PARALLEL
C$OMP MASTER
C$    NBTHREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

C     RECONNAISSANCE DE LA LANGUE DES DONNEES DE MEFISTO
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
         print *,'Mefisto-ELASTICER: NOMBRE DES ARGUMENTS',NBARGS+1
      ELSE
         print *,'Mefisto-ELASTICER: ARGUMENT NUMBER=',NBARGS+1
      ENDIF
      DO N=0,NBARGS
C        EXPLORATION DES ARGUMENTS
         CALL GETARG( N, ARGUMENT )
         print *,'Mefisto-ELASTICER: ARGUMENT',N,' = ',ARGUMENT
      ENDDO
C
C     LA COMMANDE Mefisto-ELASTICER est elle suivie d'un nom de fichier de donnees?
C     =============================================================================
      CALL NMFIDOCO( 'ppelas', EXISTNF, NMFILEDO )
C
C     LE MODE DES INTER-ACTIONS ENTRE LE PROGRAMME ET L'UTILISATEUR
C     INTERA=0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS

      IF( EXISTNF .EQ. 0 ) THEN
C        PAS DE NOM DE FICHIER DE DONNEES => INTERACTIF AVEC X11
         INTERA = IINFO( 'INTERACTIVITE INITIALE' )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Mefisto-ELASTICER EN EXECUTION INTERACTIVE'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-ELASTICER in INTERACTIVE EXECUTION'
         ENDIF
      ELSE
C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
         INTERA = 0
         LAPXFE = 800
         LHPXFE = 600
         IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,*)'Mefisto-ELASTICER: Nom Fichier Donnees=',NMFILEDO
       WRITE(IMPRIM,*)'Mefisto-ELASTICER s''execute en mode BATCH'
         ELSE
           WRITE(IMPRIM,*) 'Mefisto-ELASTICER: Data File Name=',NMFILEDO
           WRITE(IMPRIM,*) 'Mefisto-ELASTICER in BATCH EXECUTION'
         ENDIF
      ENDIF
      LHLECT = 1
      INTERB(LHLECT) = INTERA
      NBPAS  = 0
C
C     INITIALISATIONS
C     ===============
      CALL INITIA( EXIST )
C
C     OUVERTURE DE X11-WINDOW
      IF( INTERA .GE. 1 ) THEN
         CALL XTINIT
         CALL XVINIT
C        TRACE DU LOGO DE MEFISTO
         CALL LOGO( 'MEFISTO-ELASTICER' )
      ENDIF
C
C     OUVERTURE DE LA MS
 10   CALL OUVRMS
C
C     LE LANGAGE UTILISATEUR RETROUVE SES VARIABLES
      CALL LUOU
C
C     AFFICHAGE DE LA TAILLE DU SUPER-TABLEAU NUMERIQUE
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) MOTMCN
      ELSE
         WRITE(IMPRIM,20010) MOTMCN
      ENDIF
10010 FORMAT('SUPER-TABLEAU MCN DE',I12,' MOTS-MEMOIRE')
20010 FORMAT('SUPER-ARRAY MCN WITH',I12,' WORDS of MEMORY')
C
C     RECUPERATION DU DERNIER OBJET DECLARE
C     -------------------------------------
 20   CALL REDEOB( KNOMOB, NUOB, MN )

      IF( NUOB .LE. 0 ) THEN
C        PAS D'OBJET RETROUVE => ERREUR
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: PAS D''OBJET MAILLE'
            KERR(2) = 'CREER UN OBJET AVEC MAILLER'
         ELSE
            KERR(1) = 'ERROR: NO OBJECT'
            KERR(2) = 'CREATE AN OBJECT WITH MAILLER'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      NBPAS = NUOB

C     EXISTE-T-IL UN FICHIER DE DONNEES POUR ELASTICER?
      IF( EXISTNF .NE. 0 ) THEN
C
C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
         KLG(1) = 'READF ' // NMFILEDO
         N = NUDCNB( NMFILEDO )
         KLG(1) = 'READF ' // NMFILEDO(1:N) // ' ; '
         LHKLG  = 1
         NLPTV1 = 1
         NCPTV1 = 0
         NOTYPE =-2
         CALL DONNMF( NOTYPE, NOTYPS, NCVALS, DBLVAL, KCHAIN )

      ENDIF
C
C     LECTURE DES MOTS CLE
C     ====================
C     L'ANCIEN HISTORIQUE EST EFFACE
 30   CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C     TRACE DES AXES REINITIALISE
      NETAXE = 0
C
      CALL LIMTCL( 'debuelas' , NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 30
      IF( NMTCL .GT. 70 ) GOTO 7001
      IF( NMTCL .EQ. 38 ) GOTO 3800
      IF( NMTCL .EQ. 39 ) GOTO 3900
      GOTO( 100, 200, 300, 400, 30, 600, 700, 800, 30, 1000,
     %       30,  30,  30,  30, 30,  30,  30,  30, 30, 2000), NMTCL
C
C     NOM DE L'OBJET A TRAITER
 100  IF( NBPAS .GT. 0 ) THEN
C        FERMETURE DE L'OBJET PRECEDENT
         CALL LXLXFE( NTOBJE , KNOMOB )
         NBPAS = 0
C        HISTORIQUE REMIS A ZERO
         CALL RECTEF( NRHIST )
         NBLGRC( NRHIST ) = 0
      ENDIF
C
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
      CALL CADEXT( MNXYZP, COOEXT )
C
C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL VISEE0
         CALL T1OBJE( KNOMOB )
      ENDIF
      NBPAS = 1
      GOTO 30
C
C     DONNEES_ELASTICITE
 200  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL DFELAS( KNOMOB , IERR )
      GOTO 30
C
C     ELASTICITE_STATIONNAIRE
 300  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL ELASTA( KNOMOB , IERR )
      GOTO 30
C
C     ELASTICITE_INSTATIONNAIRE LINEAIRE
 400  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL ELAINS( KNOMOB , IERR )
      GOTO 30
C
C     CALCUL DES MODES PROPRES
 600  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL ELAVVP( KNOMOB , IERR )
      GOTO 30
C
C     DESSIN DES MODES PROPRES
 700  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRELAS( KNOMOB, 2, IERR )
      GOTO 30
C
C     DESSIN VECTEURS DEPLACEMENTS
 800  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRELAS( KNOMOB, 1, IERR )
      GOTO 30
C
C     TRACE DES MAILLAGES
 1000 CALL TRMAIL
      GOTO 30
C
C     PRECISION POUR RESOUDRE Ax=b
 2000 CALL ZEROGC
      GOTO 30
C
C     NOMBRE PIXELS de la LARGEUR HAUTEUR de la FENETRE
 3800 CALL XVPXFE
      GOTO 30
C
C     COULEUR du FOND
 3900 CALL INVITE( 17 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 30
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOFON = 0
      ELSE
         NCOFON = N1COEL + I
      ENDIF
C     LA COULEUR DU FOND EST IMPOSEE
      CALL XVFOND( NCOFON )
      CALL EFFACE
      GOTO 30
C
C     LA GESTION DES RESSOURCES DE MEFISTO
C     ====================================
 7001 NMTCL = NMTCL - 70
      GOTO( 7100, 7200, 7300, 7400, 30, 30,  30,   30,   30, 30,
     %        30,   30,   30,   30, 30, 30,  30,   30,   30, 30,
     %        30,   30,   30,   30, 30, 30,9700, 9800, 9900, 30), NMTCL
C
C     SUIVI_TMS
 7100 CALL SUITMS
      GOTO 30
C
C     SUIVI DES FICHIERS DE LA MS
 7200 CALL SUIFMS
      GOTO 30
C
C     UTILITAIRES DE GESTION DES UNITES DE LECTURE AFFICHAGE
 7300 CALL SUIVES
      GOTO 30
C
C     DESTRUCTION DE TMS POINT, LIGNE, ..., OBJET, ...
 7400 CALL TUER
      GOTO 30
C
C     AFFICHAGE DU NOM DE LA VERSION DE MEFISTO
 9700 NBLGRC(NRERR) = 1
      KERR(1) = ' '
      CALL VRSION( KERR(1) )
      CALL LERESU
      GOTO 30
C
C     SAUVEGARDE DES DONNEES SUR LA MS ET REDEMARRAGE
 9800 CALL ARRET( -2 )
      GOTO 10
C
C     SAUVEGARDE DES DONNEES SUR LA MS ET FIN DU TRAITEMENT
 9900 CALL ARRET(  0 )
C
      STOP
      END
