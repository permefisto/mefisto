      PROGRAM PPADAP
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LA MEMOIRE SECONDAIRE
C ----- TRAITER EVENTUELLEMENT L'ADAPTATION D'UN PROBLEME DE THERMIQUE
C       FERMER LA MEMOIRE SECONDAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1995
C2345X7..............................................................012
C$    USE OMP_LIB
      include"./incl/threads.inc"
      include"./incl/pp.inc"
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
      COMMON / EPSSSS /  EPZERO, EPSXYZ
      LOGICAL            EXIST
      CHARACTER*24       KNOMOB

C///////////////////////////////////////////////////////////////////////
      NBTHREADS = 1
C$OMP PARALLEL
C$OMP MASTER
C$    NBTHREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////
C
      NBJEUX = 1
C
C     ESSAI DE RECUPERATION DU NOM DU PROJET
C     ======================================
      CALL NOMJOB( EXIST )
      IF( .NOT. EXIST ) THEN
C        LE PROJET N'EXISTE PAS => ARRET
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: EXECUTEZ INITIER AUPARAVANT'
         ELSE
            KERR(1) ='ERROR: PREVIOUSLY EXECUTE INITIER'
         ENDIF
         CALL LEREUR
         STOP
      ENDIF
C
C     LE MODE DES INTER-ACTIONS ENTRE LE PROGRAMME ET L'UTILISATEUR
C     =============================================================
      INTERA = IINFO( 'INTERACTIVITE INITIALE' )
      LHLECT = 1
      INTERB(LHLECT) = INTERA
      NBPAS  = 0
      NBISO  = 0
C
C     INITIALISATIONS
C     ===============
      CALL INITIA( EXIST )
C
C     OUVERTURE DE XTOOLKIT ET MOTIF
      CALL XTINIT
C
C     INITIALISATION DE XV
      CALL XVINIT
C
C     TRACE DU LOGO DE MEFISTO
      CALL LOGO( 'MEFISTO-ADAPTATION' )
C
C     OUVERTURE DE LA MS
 10   CALL OUVRMS
C
C     LE LANGAGE UTILISATEUR RETROUVE SES VARIABLES
      CALL LUOU
C
C     AFFICHAGE DE LA TAILLE DU SUPER-TABLEAU NUMERIQUE
      WRITE(IMPRIM,10010) MOTMCN
10010 FORMAT(/' NOMBRE DE MOTS DU SUPER-TABLEAU MCN =',I12)
C
C     RECUPERATION DU DERNIER OBJET DECLARE
C     -------------------------------------
 20   CALL REDEOB( KNOMOB, NUOB, MN )
      IF( NUOB .LE. 0 ) THEN
C        PAS D'OBJET RETROUVE
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
C
C     LECTURE DES MOTS CLE
C     ====================
C     L'ANCIEN HISTORIQUE EST EFFACE
 30   CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'OBJET: ' // KNOMOB
      CALL LHISTO
C     TRACE DES AXES REINITIALISE
      NETAXE = 0
C
      CALL LIMTCL( 'debuadap' , NMTCL )
      IF( NMTCL .EQ. -1 ) GOTO 30
      IF( NMTCL .GT. 70 ) GOTO 7001
      IF( NMTCL .EQ. 39 ) GOTO 3900
      GOTO( 100, 200, 300, 400,  30,  30, 700, 800, 30, 1000,
     %       30,  30,  30,  30,  30,  30,  30,  30, 30, 2000), NMTCL
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
         KERR(1) ='ERREUR: OBJET INCONNU:' // KNOMOB
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
      INIEXT = 0
      CALL MAJEXT( MNXYZP )
C
C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'OBJET: ' // KNOMOB
      CALL LHISTO
C
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL VISEE1
         CALL T1OBJE( KNOMOB )
      ENDIF
      NBPAS = 1
      GOTO 30
C
C     INTERPOLATION ET RENUMEROTATION DES NOEUDS DES ELEMENTS FINIS
C     PROTECTION PUIS REGENERATION DE LA PRECISION POUR IDENTIFIER LES POINTS
C     POUR EVITER L'IDENTIFICATION DE 2 SOMMETS D'UNE ARETE APRES ADAPTATION
 200  EPZER0 = EPZERO
      EPSXY0 = EPSXYZ
      EPZERO = 0
      EPSXYZ = 0
      CALL DFTOPO
      EPZERO = EPZER0
      EPSXYZ = EPSXY0
      GOTO 30
C
C     DONNEES_THERMIQUE
 300  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL DFTHER( KNOMOB, NBJEUX, IERR )
      GOTO 30
C
C     THERMIQUE_STATIONNAIRE
 400  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL THESTA( KNOMOB, IERR )
      GOTO 30
C
C     ADAPTATION DU MAILLAGE EN 2D
 700  CALL ADAP2D( KNOMOB, NBISO, IERR )
      GOTO 30
C
C     DESSIN_RESULTATS
 800  IF( NBPAS .EQ. 0 ) GOTO 100
      CALL TRTHER( KNOMOB, 1, IERR )
      GOTO 30
C
C     TRACE DES MAILLAGES
 1000 CALL TRMAIL
      GOTO 30
C
C     PRECISIONS AUTOUR D'UN POINT D'UN OBJET ...
 2000 CALL ZEROGC
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
C     NOM DE LA VERSION DE MEFISTO
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
