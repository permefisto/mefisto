      SUBROUTINE ERR3CUB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES COURBES en fonction des maillages nM
C ----- ERREUR(1,nM) Max |ue-uc| / ( Max ue - Min ue )
C       ERREUR(2,nM) Max ( |ue-uc| / |ue| )
C       ERREUR(3,nM) Sum |ue-uc| / Sum |ue|
C       NORMES RELATIVES DE L'ERREUR AUX NOEUDS ENTRE
C       ue solution exacte   aux noeuds
C       uc solution calculee aux noeuds
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL Paris & St Pierre du Perray Fevrier 2009
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*120 TITRE
C
C     LES DONNEES CALCULEES PAR THERMICER
      PARAMETER (NBM=10)
      INTEGER    NB3CUB(NBM), NBNODE(NBM)
      REAL       ERREUR(3,NBM)
C     ERREUR(1,nM) Max |ue-uc| / ( Max ue - Min ue )
C     ERREUR(2,nM) Max ( |ue-uc| / |ue| )
C     ERREUR(3,nM) Som |ue-uc| / Som |ue|
C
C     INITIALISATIONS DES DONNEES
      DATA NBNODE/1432,2648,4064,4096,5248,7351,7677,9784,9784,10963/
      DATA NB3CUB/1325,2525,3871,3375,5047,7000,7424,9477,9477,10600/
C
      DATA ERREUR/ 0.03,   0.14,   0.0079,
     %             0.016,  0.053,  0.0047,
     %             0.01,   0.097,  0.0036,
     %             0.009,  0.032,  0.012,
     %             0.01,   0.093,  0.0028,
     %             0.0054, 0.079,  0.003,
     %             0.0069, 0.048,  0.0025,
     %             0.0068, 0.096,  0.0019,
     %             0.0057, 0.06,   0.0021,
     %             0.0055, 0.084,  0.0019  /
C
C     LE CADRE DE LA FENETRE DE TRACE
      print *
      HMIN   =  1E20
      HMAX   = -1E20
      ERRMIN =  1E20
      ERRMAX = -1E20
      DO 3 M = 1, NBM
         HH = NBNODE( M )
         HMIN = MIN( HMIN, HH )
         HMAX = MAX( HMAX, HH )
         ER = ERREUR(3,M)
         ERRMIN = MIN( ERRMIN, ER )
         ERRMAX = MAX( ERRMAX, ER )
 3    CONTINUE
      print *,'HMIN=',HMIN,' HMAX=',HMAX,' ERRMIN=',ERRMIN,
     %        '  ERRMAX=',ERRMAX
      HH = ( HMAX - HMIN ) / 10
      RH = ( ERRMAX - ERRMIN ) / 10
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( HMIN-HH, HMAX+HH, ERRMIN-RH, ERRMAX+RH )
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
ccc      CALL TEXTE2D( NCVERT, HMIN, ERRMAX+RH/5,
ccc     %             'Sum |ue-uc| / Sum |ue|')
      CALL TEXTE2D( NCMAGE, HMAX-HH, ERRMIN-RH/2, 'Vertices Nb' )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE DES ( Sum |ue-uc| / Sum |ue| )  en fonction de Nb de Sommets
      DO 20 M = 1, NBM
C
C        CHOIX DE LA COULEUR DE TRACE
         NC = MOD(M-1,4) + 1
C
C        LE POINT A TRACER (  NBNODE( M ), ERREUR(3,M) )
         H  = NBNODE( M )
         ER = ERREUR(3,M)
C        LE NO DU MAILLAGE A GAUCHE
         CALL ENTIER2D(  NC, H, ER, NBNODE( M ) )
         CALL SYMBOLE2D( NC, H, ER, '*' )
C
 20      CONTINUE
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 25 )
      TITRE ='- k Delta u + c u = f on Omega=[-10,10]**3, u=0 on Gamma '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+HH, ERRMAX+RH/2, TITRE(1:L) )
C
      TITRE = 'u exact = 1E-6 (10+x)(10-x)(10+y)(10-y)(10-z)(10+z) '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+HH, ERRMAX, TITRE(1:L) )
C
      TITRE='Sum |u exact - u computed| / Sum |u exact| (Nodes Number) '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+HH, ERRMAX-Rh/2, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
C     MISE SUR FICHIER 3cuberr.eps du TRACE
C     ATTENTION PASSAGE PAR VARIABLE OBLIGATOIRE
      NBC = 7
      CALL xvsauverps( '3cuberr', NBC )
      print *, 'NBC=',NBC
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
      CALL XVEPAISSEUR( 1 )
      RETURN
      END
