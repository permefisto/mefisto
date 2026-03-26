      SUBROUTINE TRNLSEMXU( NTLXOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER Max Module U EN FONCTION DU TEMPS
C -----   'VECTEUR"TESTM_ERREUR' DONNEES CALCULEES PAR nlseimpl.f
C          VECTEUR 1 = VALEUR DU TEST AU TEMPS
C          VECTEUR 2 = VALEUR DU Max|U(Noeud)| AU TEMPS
C ENTREE :
C --------
C NTLXOB : NUMRO DU TMS DE L'OBJET
C
C SORTIE :
C --------
C MISE SUR FICHIER Max Module U_temps.eps du TRACE des TESTSm(Temps)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___vecteur.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DOUBLE PRECISION DMCN(1)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*80      TITRE
      CHARACTER*20      KDT0
      CHARACTER*16      KDT
ccc
ccc      CHARACTER*1       CAR(2)
ccc      DATA              CAR / '+', '.' /
C
C     TABLEAUX AUXILIAIRES
      REAL              PT(2,2)
C
C     OUVERTURE DU TMS 'VECTEUR"TESTM_ERREUR' de l'OBJET
      CALL  LXTSOU( NTLXOB, 'VECTEUR"TESTM_ERREUR', NTVERT, MNVERT )
      IF( MNVERT .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'Le tms VECTEUR"TESTM_ERREUR est INCONNU'
         ELSE
         ENDIF
            WRITE(IMPRIM,*)'The tms VECTEUR"TESTM_ERREUR is UNKNOWN'
         RETURN
      ENDIF
C
C     NOMBRE DE COMPOSANTES D'UN VECTEUR
      NBC = MCN( MNVERT + WBCOVE )
C
C     NOMBRE DE VECTEURS (2: TESTm, Max|U|
C                      ou 4: TESTm, Max|U|, ERREUR PR, ERREUR PI)
      NBV = MCN( MNVERT + WBVECT )
      IF( NBV .LE. 0 ) RETURN
C
C     NOMBRE DE TEMPS FOURNIS = NOMBRE DE COMPOSANTES
CCC   NBT = MCN( MNVERT + WBCPIN ) = NBC
C
C     CALCUL DES DIMENSIONS DE LA FENETRE DE TRACE
C     --------------------------------------------
C     Max|U| EST LE VECTEUR APRES LES TESTm
      MOREE2 = MOTVAR(6)
      MNE  = ( MNVERT + WECTEU + 1 ) / MOREE2 + NBC
      EMAX = REAL( DMCN( MNE ) )
C
C     TEMPS(1) EST LA PREMIERE INFORMATION APRES LES NBV VECTEURS
      MNT  = MNVERT + WECTEU + MOREE2 * NBC * NBV
      TMIN = RMCN( MNT )
      TMAX = RMCN( MNT )
C
      MNT = MNT - 1
      MNE = MNE - 1
C
C     CALCUL DU TEMPS MIN, TEMPS MAX, MaxU MIN, MaxU MAX
      DTMAX = 0
      T0    = RMCN( MNT + 1 )
      DO K = 1, NBC
C
C        TEMPS K
         T = RMCN( MNT + K )
         TMIN = MIN( TMIN, T )
         TMAX = MAX( TMAX, T )
         DTMAX = MAX( DTMAX, T-T0 )
         T0    = T
C
C        MaxU K
         E = REAL( DMCN( MNE + K ) )
         EMIN = MIN( EMIN, E )
         EMAX = MAX( EMAX, E )
C
      ENDDO
      EMIN = 0.0
C
      WRITE(IMPRIM,*) 'trnlsemxu: nbc=',NBC,' nbv=',NBV,
     %                ' Time MIN=',TMIN,' Time MAX=',TMAX,
     %                ' Max|U| MIN=',EMIN,' Max|U| MAX=',EMAX
      WRITE(IMPRIM,*) 'trnlsemxu: Last Time=',T,' Last  Max|U|=',E,
     %                ' k=',K
C
      TH = ( TMAX - TMIN ) / 15
      EH = ( EMAX - EMIN ) / 15
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
C     LA COULEUR DU FOND EST IMPOSEE
C     THE BACKGROUND COLOR IS IMPOSED
      NC0    = NCOFON
      NCOFON = NCTURQ
      CALL XVFOND( NCOFON )
      CALL EFFACE
      CALL FENETRE( TMIN-TH, TMAX+TH, EMIN-EH, EMAX+2*EH )

C     TRACE NOM PROJET UTILISATEUR DATE ...
      CALL TIT1LG
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )

C     TEXTE de l'AXE X
      CALL TEXTE2D( NCBLEU, TMAX, EMIN, 'TIME' )

C     TEXTE de l'AXE Y
ccc      CALL TEXTE2D( NCBLEU, TMIN, EMAX+EH/3, 'Max|U|' )
C
C     LE TRACE DES AXES 2D
C     --------------------
      NETAXE = 0
      CALL TRAXE2
C
C     LA COURBE Max|U|(temps)
C     -----------------------
      CALL XVEPAISSEUR( 0 )
C
C     LE PAS DE TEMPS INITIAL
      T0 = RMCN( MNT+1 )
      K = 1
 3    K = K+ 1
      T  = RMCN( MNT+K )
      DT0 = T - T0
      IF( DT0 .LE. 0 ) GOTO 3
C     AFFICHAGE DU PAS DE TEMPS INITIAL
      WRITE( KDT0(1:14), '(G14.6)' ) DT0
      CALL REELCA( KDT0(1:14), KDT )
      K = NUDCNB( KDT )
      KDT0 = 'dt0=' // KDT(1:K)
      CALL TEXTSB( KDT0, LDT0 )

C     LE 1-ER SOMMET
      DT0 = 0.0
      DT  = 0.0
      HDT0= EMAX
      NCT = 0
      NC  = 0
C     TEMPS(1)
      T0 = RMCN( MNT+1 )
      PT(1,2) = T0
C     Max |U|(1)
      PT(2,2) = REAL( DMCN( MNE+1 ) )
C
ccc   CALL ENTIER2D(  NCBLEU,  PT(1,2), PT(2,2), 1 )
ccc   CALL SYMBOLE2D( NCNOIR, PT(1,2), PT(2,2), CAR(1) )
      NBSTPT = 0
C
C     LE TRACE DES ( TEMPS(K), Max|U|(K) )
      DO 10 K = 2, NBC
C
C        LE POINT INITIAL EST LE POINT FINAL PRECEDENT
         PT(1,1) = PT(1,2)
         PT(2,1) = PT(2,2)
C
C        LE POINT FINAL DE L'ARETE
         T = RMCN( MNT+K )
         IF( T .LT. T0 ) GOTO 10

         PT(1,2) = T
         PT(2,2) = REAL( DMCN( MNE+K ) )
C
C        CHOIX DE LA COULEUR EN UN MEME TEMPS
         IF( T .NE. T0 ) THEN
            NBSTPT = NBSTPT + 1
            NC = 0
            DT = T - T0
            T0 = T
         ELSE
C           LA COULEUR DE TESTm CHANGE A CHAQUE ITERm
            NC = NC + 1
            IF( NC .GT. 8 ) NC = 0
         ENDIF
         CALL TRAIT2D(   NC, PT(1,1), PT(2,1), PT(1,2), PT(2,2) )
ccc      CALL SYMBOLE2D( NCNOIR, PT(1,2), PT(2,2), CAR(1) )

C        TRACE DU SEGMENT Max |U|
         HDT = EMAX + DT / DTMAX * EH
         IF( ABS(DT-DT0) .GT. 1E-3*DT ) THEN
C           TRACE DE BAS EN HAUT POUR MONTRER LE CHANGEMENT DE PAS DE TEMPS
            CALL TRAIT2D( NCGRIS, PT(1,1), EMIN, PT(1,1), EMAX-EH/2 )
            CALL TRAIT2D( NCGRIS, PT(1,1), EMAX, PT(1,1), MAX(HDT0,HDT))
            HDT0 = HDT
            NCT  = NCT + 1
            IF( NCT .GT. 8 ) NCT = 0
            DT0 = DT
         ENDIF
         CALL TRAIT2D( NCT, PT(1,1), HDT, PT(1,2), HDT )

 10   CONTINUE
ccc   CALL ENTIER2D( NCBLEU,  PT(1,2), PT(2,2), 1 )
      CALL TRAIT2D( NCNOIR, TMIN, EMIN, TMAX, EMIN )

C     L'AXE DES PAS DE TEMPS AU COURS DU TEMPS
      TITRE='TIME'
      NB = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLEU, TMAX, EMAX, TITRE(1:NB) )
      CALL TRAIT2D( NCNOIR, TMIN, EMAX, TMAX, EMAX )
      CALL TRAIT2D( NCNOIR, TMIN, EMAX, TMIN, EMAX+EH )
C
C     LE TITRE DU GRAPHIQUE
C     ---------------------
      CALL CHOIXFONTE( 25 )

      KDT = ' '
      WRITE( KDT(1:8), '(I8)' ) NBSTPT
      CALL TEXTSB( KDT, K )
      TITRE = KDT0(1:LDT0)//'  Step Time(time)  STNb='//KDT(1:K)//' '
      NB = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLEU, TMIN+TH/6, EMAX+EH, TITRE(1:NB) )

      TITRE='Max|Wave(Node)|(time)'
      NB = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLEU, TMIN+TH/6, EMAX-EH*2, TITRE(1:NB) )

C     TRACE DU NOMBRE ITERm
      KDT = 'ITERm='
      WRITE( KDT(7:16), '(I9)' ) NBC
      CALL TEXTSB( KDT, NB )
      CALL TEXTE2D( NCBLEU, TMAX-TH*0.5, EMAX-EH*2, KDT(1:NB) )

C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
C     MISE SUR FICHIER maxmoduleu_temps.eps du TRACE
C     ATTENTION NB PASSAGE PAR VARIABLE OBLIGATOIRE
      NB = NUDCNB( 'maxmoduleu_temps' )
      CALL xvsauverps( 'maxmoduleu_temps', NB )
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
C     REMISE EN ETAT
      CALL XVEPAISSEUR( 1 )
      NCOFON = NC0
      CALL XVFOND( NCOFON )
C
      RETURN
      END
