      PROGRAM PPNLSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  OUVRIR  LA MEMOIRE SECONDAIRE      (OPEN  THE SECONDARY MEMORY)
C -----  TRAITER UN PROBLEME NON LINEAIRE de SHRODINGER (TREAT NLSE )
C        i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(|U|**2) U(t,X)
C                         - i OmegaZ (x dU/dy - y dU/dx) = F
C        FERMER  LA MEMOIRE SECONDAIRE      (CLOSE THE SECONDARY MEMORY)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR  Fevrier 2011
C MODIF  : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Mars    2013
C MODIF  : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Mai     2014
C23456---------------------------------------------------------------012
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
      include"./incl/cnonlin.inc"
      include"./incl/traaxe.inc"
C
      REAL              RMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      COMMON             MCN(MOTMCN)
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      LOGICAL            EXIST
      CHARACTER*24       KNOMOB
      DOUBLE PRECISION   EXPOSANT
C
C     LES ARGUMENTS A L'APPEL DU LOAD MODULE ppnlse
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
C     RECUPERATION DU NOM DU PROJET
C     RECUPERATION OF THE OBJECT NAME
C     ===============================
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
         print *,'Mefisto-NLSER: NOMBRE DES ARGUMENTS=',NBARGS+1
      ELSE
         print *,'Mefisto-NLSER: ARGUMENT NUMBER=',NBARGS+1
      ENDIF
      DO N=0,NBARGS
C        EXPLORATION DES ARGUMENTS
         CALL GETARG( N, ARGUMENT )
         print *,'Mefisto-NLSER: ARGUMENT',N,' = ',ARGUMENT
      ENDDO
C
C     LA COMMANDE Mefisto-NLSER est elle suivie d'un nom de fichier de donnees?
C     ===========================================================================
      CALL NMFIDOCO( 'ppnlse', EXISTNF, NMFILEDO )
C
C     LE MODE DES INTER-ACTIONS ENTRE LE PROGRAMME ET L'UTILISATEUR
C     INTERA=0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C            3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS

      IF( EXISTNF .EQ. 0 ) THEN
C        PAS DE NOM DE FICHIER DE DONNEES => INTERACTIF AVEC X11
         INTERA = IINFO( 'INTERACTIVITE INITIALE' )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Mefisto-NLSER EN EXECUTION INTERACTIVE'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-NLSER in INTERACTIVE EXECUTION'
         ENDIF
      ELSE
C        UN NOM DE FICHIER EXISTE => BATCH SANS X11
C        INITIALISATIONS POUR EVITER LES TRACES X11
         INTERA = 0
         LAPXFE = 800
         LHPXFE = 600
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'Mefisto-NLSER: Nom Fichier Donnees=',NMFILEDO
         WRITE(IMPRIM,*) 'Mefisto-NLSER s''execute en mode BATCH'
         ELSE
            WRITE(IMPRIM,*) 'Mefisto-NLSER: Data File Name=',NMFILEDO
            WRITE(IMPRIM,*) 'Mefisto-NLSER in BATCH EXECUTION'
         ENDIF
      ENDIF
      LHLECT = 1
      INTERB(LHLECT) = INTERA
C
C     INITIALISATIONS
C     ===============
      CALL INITIA( EXIST )
C
C     OUVERTURE DE X11-WINDOW
C     OPEN OF X11-WINDOW
      IF( INTERA .GE. 1 ) THEN
         CALL XTINIT
C        INITIALISATION DE XVUE
         CALL XVINIT
C        TRACE DU LOGO DE MEFISTO
C        DRAWING OF MEFISTO LOGO
         CALL LOGO( 'MEFISTO-NLSE' )
      ENDIF
C
cccC     INITIALISATION of PETSC   11/12/2009 Cf bin.petsc
ccc      call PetscIni()
C
C     OUVERTURE DE LA MS
C     OPEN THE MS
 10   CALL OUVRMS
C
C     LE LANGAGE UTILISATEUR RETROUVE SES VARIABLES
C     THE USER LANGUAGE RECOVERS ITS VARIABLES
      CALL LUOU
C
C     AFFICHAGE DE LA TAILLE DU SUPER-TABLEAU NUMERIQUE
C     PRINTING OF THE MEMORY SIZE FOR COMPUTATION
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) MOTMCN
      ELSE
         WRITE(IMPRIM,20010) MOTMCN
      ENDIF
10010 FORMAT('SUPER-TABLEAU MCN DE',I12,' MOTS-MEMOIRE')
20010 FORMAT('SUPER-ARRAY MCN WITH',I12,' WORDS of MEMORY')
C
C     RECUPERATION DU DERNIER OBJET DECLARE
C     RECUPERATION OF THE LAST DECLARED OBJECT
C     ----------------------------------------
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
C
      NBPAS  = NUOB
      NBJEUX = 1

C     EXISTE-T-IL UN FICHIER DE DONNEES POUR NLSER?
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
C     LECTURE SU MENU DES MOTS CLE
C     READING THE MENU OF KEY-WORDS
C     =============================
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
      CALL LIMTCL( 'debunlse', NMTCL )
      IF( NMTCL .LT.  0 ) GOTO   30
      IF( NMTCL .GT. 70 ) GOTO 7001
      IF( NMTCL .EQ. 38 ) GOTO 3800
      IF( NMTCL .EQ. 39 ) GOTO 3900
      GOTO( 100,  200,  300,  400,  500,  600, 700,  30,   30, 1000,
     %     1100, 1200, 1300, 1400, 1500, 1600,  30,  30, 1900, 2000
     %    ), NMTCL
C
C     NOM de L'OBJET a TRAITER
C     NAME of the OBJECT to TREAT
C     ===========================
 100  IF( NBPAS .GT. 0 ) THEN
C        FERMETURE DE L'OBJET PRECEDENT
C        CLOSING OF THE PREVIOUS OBJECT
         CALL LXLXFE( NTOBJE, KNOMOB )
         NBPAS = 0
C        HISTORIQUE REMIS A ZERO
C        HISTORY IS ANNULED
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
C     UPDATE OF MIN MAX XYZ OF THE OBJECT in COOEXT of xyzext.inc
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
      CALL MAJEXT( MNXYZP )
C
C     TRACE DU NOM DE L'OBJET
C     DRAWING OF THE OBJECT NAME
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
C        DRAWING OF THE OBJECT IF SUFFICIENT INTERACTIVITY
         CALL VISEE1
         CALL T1OBJE( KNOMOB )
      ENDIF
      NBPAS = 1
      GOTO 30
C
C     LECTURE DES DONNEES DE L'OPERATEUR NLSE
C     =======================================
 200  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL DFNLSE( KNOMOB, IERR )
      GOTO 30
C
C     CALCULER les solutions (E,u(X)) de type U(t,X)=u(X) exp(iEt) dans
C     i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(|U|**2) U(t,X) = 0
C     (SANS TERME DE ROTATION OmegaZ=0)
C     c-a-d CALCULER LES VALEURS E et VECTEURS PROPRES u(X) reels tels que
C     -Alfa LAPLACIEN u + N(u**2) u = E u
C     NonLinear SCHRODINGER EQUATIONS for an OPTIC BEAM
C     TRANSFORMED to an Nonlinear EIGENVALUE PROBLEM
C     ====================================================================
 300  IF( NBPAS .EQ. 0 ) GOTO 100
      TESTNL = 5
      CALL NLSEVVP( KNOMOB, IERR )
      GOTO 30
C
C     Schema1: idu(t,X)/dt -Alfa LAPLACIEN u(t,X) +N(u(t,X)**2) u(t,X)=0
C     Propagation d'un rayon optique selon l'EQUATION COMPLEXE
C     NON LINEAIRE de SCHRODINGER avec Matrice GLOBALE nxn SYMETRIQUE
C     NonLinear SCHRODINGER EQUATIONS for an OPTIC BEAM
C     Methode SEMI-IMPLICITE MODIFIEE AVEC MATRICE GLOBALE nxn
C     =======================================================================
 400  IF( NBPAS .EQ. 0 ) GOTO 100
      TESTNL = 6
      CALL NLSEINS( KNOMOB, IERR )
      GOTO 30
C
C     Schema2: idu(t,X)/dt -Alfa LAPLACIEN u(t,X) + N(u(t,X)**2) u(t,X)=0
C     Propagation d'un rayon optique selon l'EQUATION COMPLEXE
C     NON LINEAIRE de SCHRODINGER avec Matrice GLOBALE 2nx2n NON SYMETRIQUE
C     NonLinear SCHRODINGER EQUATIONS for an OPTIC BEAM
C     TRAITEMENT IMPLICITE EN TEMPS AVEC UNE MATRICE GLOBALE 2nx2n
C     =======================================================================
 500  IF( NBPAS .EQ. 0 ) GOTO 100
      TESTNL = 7
      CALL NLSEINS( KNOMOB, IERR )
      GOTO 30

cccC     Schema3: idu(t,X)/dt -Alfa LAPLACIEN u(t,X) + N(u(t,X)**2) u(t,X)=0
cccC     Propagation d'un rayon optique selon l'EQUATION COMPLEXE
cccC     NON LINEAIRE de SCHRODINGER SANS MATRICE GLOBALE
cccC     NonLinear SCHRODINGER EQUATIONS for an OPTIC BEAM
cccC     TRAITEMENT EXPLICITE EN TEMPS SANS MATRICE GLOBALE
cccC     =======================================================================
ccc 600  IF( NBPAS .EQ. 0 ) GOTO 100
ccc      TESTNL = 8
ccc      CALL NLSEINS( KNOMOB, IERR )
ccc      GOTO 30

C     GROSS-PITAEVSKII equation : 
C     i du(t,X)/dt +LAPLACIEN u(t,X)/2 - (V + Beta |u(t,X)|**2) u(t,X)
C                  +i Omega (x du/dy - y du/dx) =0
C     i time scheme => Changement t -> it  
C       du(t,X)/dt -LAPLACIEN u(t,X)/2 + (V + Beta |u(t,X)|**2) u(t,X)
C                  -i Omega (x du/dy - y du/dx) =0
C     Methode i-time SEMI-IMPLICITE MODIFIEE AVEC MATRICE GLOBALE nxn
C     ================================================================
ccc 1400 IF( NBPAS .EQ. 0 ) GOTO 100
 600  IF( NBPAS .EQ. 0 ) GOTO 100
      TESTNL = 9
      CALL NLSEINS( KNOMOB, IERR )
      GOTO 30
C
C     NLSE: i dU(t,X)/dt -Alfa LAPLACIEN U(t,X) -(V(X)+Beta U(t,X)**2) U(t,X)
C     - i OmegaZ (x d/dy - y d/dx) U(t,X) = F(t,X)
C     EQUATION COMPLEXE NON LINEAIRE de GROSS-PITAEVSKII (Alfa=-1/2)
C     INTERPOLATION ELEMENTS FINIS P1 PRODUISANT UNE EQUATION DU
C     3-EME DEGRE EN CHACUN DES SOMMETS et SANS Matrice GLOBALE nxn
C     =======================================================================
 700  IF( NBPAS .EQ. 0 ) GOTO 100
      TESTNL = 10
      CALL NLSEINS( KNOMOB, IERR )
      GOTO 30
C
C     DESSIN DES VALEURS ET VECTEURS PROPRES DE NLSE
C     DRAWING OF EIGENSOLUTIONS of NLSE
C     ==============================================
 1000  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRTHER( KNOMOB, 2, IERR )
      GOTO 30
C
C     DESSIN DU MODULE des SOLUTIONS OndeNLSE, GRADIENT, FLUX
C     DRAWING OF SOLUTIONS, GRADIENTS, FLUXES
C     =======================================================
 1100  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRTHER( KNOMOB, 4, IERR )
      GOTO 30
C
C     DESSIN DE LA PARTIE REELLE des SOLUTIONS OndeNLSE, GRADIENT, FLUX
C     DRAWING OF SOLUTIONS, GRADIENTS, FLUXES
C     =================================================================
 1200  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRTHER( KNOMOB, 5, IERR )
      GOTO 30
C
C     DESSIN DE LA PARTIE IMAGINAIRE des SOLUTIONS OndeNLSE, GRADIENT, FLUX
C     DRAWING OF SOLUTIONS, GRADIENTS, FLUXES
C     =====================================================================
 1300 IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRTHER( KNOMOB, 6, IERR )
      GOTO 30
C
C     DESSIN des VALEURS du TEST sur les ITERATIONS
C     =============================================
 1400 CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) GOTO 30
      CALL TRNLSETST( NTLXOB )
      GOTO 30
C
C     DESSIN des Max|U(Noeud)|(Temps)
C     ==================================================
 1500 CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) GOTO 30
      CALL TRNLSEMXU( NTLXOB )
      GOTO 30
C
C     DESSIN des ERREUR(Temps) = ER(Temps) + i EI(Temps)
C     ==================================================
 1600 CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) GOTO 30
      CALL TRNLSERR( NTLXOB )
      GOTO 30
C
C     TRACE des MAILLAGES
C     DRAWING of MESHES
C     ===================
 1900 CALL TRMAIL
      GOTO 30
C
C     SEUILS DE PRECISION POUR RESOUDRE Ax=b
C     PRECISION CUTOFFS to SOLVE LINEAR SYSTEM Ax=b
 2000 CALL ZEROGC
      GOTO 30
C
C     NOMBRE PIXELS de la LARGEUR HAUTEUR de la FENETRE
C     PIXELS NUMBER of WINDOW WIDTH & HEIGHT
 3800 CALL XVPXFE
      GOTO 30
C
C     COULEUR du FOND
C     BACKGROUND COLOR
 3900 CALL INVITE( 17 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 30
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE   BLACK COLOR
         NCOFON = 0
      ELSE
         NCOFON = N1COEL + I
      ENDIF
C     LA COULEUR DU FOND EST IMPOSEE
C     THE BACKGROUND COLOR IS IMPOSED
      CALL XVFOND( NCOFON )
      CALL EFFACE
      GOTO 30
C
C     LA GESTION DES RESSOURCES DE MEFISTO
C     THE MANAGEMENT OF MEFISTO RESOURCES
C     ====================================
 7001 NMTCL = NMTCL - 70
      GOTO( 7100, 7200, 7300, 7400, 30, 30, 30,   30,   30,  30,
     %        30,   30,   30,   30, 30, 30, 30,   30,   30,  30,
     %        30,   30,   30,   30, 30, 30,9700, 9800, 9900, 30), NMTCL
C
C     SUIVI_TMS
C     UPDATE of TMS
 7100 CALL SUITMS
      GOTO 30
C
C     SUIVI DES FICHIERS DE LA MS
C     UPDATE of FILES of MS
 7200 CALL SUIFMS
      GOTO 30
C
C     UTILITAIRES DE GESTION DES UNITES DE LECTURE AFFICHAGE
C     MANAGEMENT OF READING PRINTING...
 7300 CALL SUIVES
      GOTO 30
C
C     DESTRUCTION DE TMS POINT, LIGNE, ..., OBJET, ...
C     CANCELLATION OF POINT, LINE, ..., OBJECT, ...
 7400 CALL TUER
      GOTO 30
C
C     NOM DE LA VERSION DE MEFISTO
C     MEFISTO VERSION NAME
 9700 NBLGRC(NRERR) = 1
      KERR(1) = ' '
      CALL VRSION( KERR(1) )
      CALL LERESU
      GOTO 30
C
C     SAUVEGARDE DES DONNEES SUR LA MS ET REDEMARRAGE
C     SAVE DATA on MS and RESTART
 9800 CALL ARRET( -2 )
      GOTO 10
C
C     SAUVEGARDE DES DONNEES SUR LA MS ET FIN DU TRAITEMENT
C     SAVE DATA and QUIT
 9900 CALL ARRET(  0 )
C
cccC     END USING PETSC   11/12/2009   Cf bin.petsc
ccc      call PetscFinal()
C
      STOP
      END
