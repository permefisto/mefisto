      SUBROUTINE TRNLSERR( NTLXOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER ERREUR(NB) EN FONCTION DU TEMPS(NB) ECRETE A 30%
C -----   'VECTEUR"TESTM_ERREUR' DONNEES CALCULEES PAR nlseimpl.f
C          VECTEUR 1 = VALEUR DU TESTm AU TEMPS
C          VECTEUR 2 = ERREUR SUR LA PARTIE REELLE
C          VECTEUR 3 = ERREUR SUR LA PARTIE IMAGINAIRE
C ENTREE :
C --------
C NTLXOB : NUMRO DU TMS DE L'OBJET
C
C SORTIE :
C --------
C MISE SUR FICHIER erreur_temps.eps du TRACE des ERREURS(Temps)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Juillet 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/traaxe.inc"
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
C     NOMBRE DE COMPOSANTES DES 2 OU 4 VECTEURS=NOMBRE DES ITERATIONS m
      NBC = MCN( MNVERT + WBCOVE )
C
C     NOMBRE DE VECTEURS: Testm, Max |U|, eventuellement ERREUR PR, ERREUR PI
      NBV = MCN( MNVERT + WBVECT )
C
      IF( NBV .LE. 2 ) RETURN
C
C     LES ERREURS ONT ETE STOCKEES
C     ----------------------------
C     NOMBRE DE TEMPS FOURNIS = NOMBRE DE COMPOSANTES
CCC   NBT = MCN( MNVERT + WBCPIN ) = NBC
C
C     CALCUL DES DIMENSIONS DE LA FENETRE DE TRACE
C     --------------------------------------------
C     ERREUR(1) EST LE VECTEUR APRES LES TESTm et Max|U(Noeud)|
      MOREE2 = MOTVAR(6)
      MNE  = ( MNVERT + WECTEU + 1 ) / MOREE2 + NBC*2
      EMIN = REAL( DMCN( MNE ) )
      EMAX = REAL( DMCN( MNE ) )
C
C     TEMPS(1) EST LA PREMIERE INFORMATION APRES LE VECTEUR
      MNT  = MNVERT + WECTEU + MOREE2 * NBC * NBV
      TMIN = RMCN( MNT )
      TMAX = RMCN( MNT )
C
      MNT = MNT - 1
      MNE = MNE - 1
C
C     CALCUL DU TEMPS MIN, TEMPS MAX, ERREUR MIN, ERREUR MAX
      DTMAX = 0
      T0    = RMCN( MNT + 1 )
      DO 8 K = 1, NBC
C
C        TEMPS K
         T = RMCN( MNT + K )
         TMIN = MIN( TMIN, T )
         TMAX = MAX( TMAX, T )
         DTMAX = MAX( DTMAX, T-T0 )
         T0    = T
C
C        ERREUR K
         DO NV=3,NBV
            E = REAL( DMCN( MNE + K + (NV-3) * NBC ) )
C           TEST POUR EVITER UN ECRASEMENT DE L'ECHELLE DES ERREURS
            IF( E .GT. 100 * EMAX ) GOTO 8
            EMIN = MIN( EMIN, E )
            EMAX = MAX( EMAX, E )
         ENDDO
C
 8    CONTINUE
C
      WRITE(IMPRIM,*) 'trnlserr: nbc=',NBC,' nbv=',NBV,
     %                ' Time MIN=',TMIN,' Time MAX=',TMAX,
     %                ' Error MIN=',EMIN,' Error MAX=',EMAX
      WRITE(IMPRIM,*) 'trnlserr: Last Time=',T,' Last Error=',E,' k=',K
C
C     ECRETAGE DE L'ERREUR A 30% POUR VOIR LES OSCILLATIONS FAIBLES
      EMIN = MIN( 0.0, EMIN )
      EMAX = MIN( 0.3, EMAX )
C
      TH = ( TMAX - TMIN ) / 15
      EH = ( EMAX - EMIN ) / 15
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
C     LA COULEUR DU FOND EST IMPOSEE
C     THE BACKGROUND COLOR IS IMPOSED
      NC0    = NCOFON
      NCOFON = NCORAN
      CALL XVFOND( NCOFON )
      CALL EFFACE
      CALL FENETRE( TMIN-TH, TMAX+TH, EMIN-EH, EMAX+2*EH )

C     TRACE NOM PROJET UTILISATEUR DATE ...
      CALL TIT1LG
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
C     TEXTE de l'AXE X
      CALL TEXTE2D( NCORAN, TMAX, EMIN, 'TIME' )
C     TEXTE de l'AXE Y
ccc      CALL TEXTE2D( NCBLEU, TMIN, EMAX+EH/3, 'ERROR' )
C
C     LE TRACE DES AXES 2D
C     --------------------
      NETAXE = 0
      CALL TRAXE2
C
C     L'EPAISSEUR DE LA COURBE
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
C
C     LE 1-ER SOMMET
      DT0 = 0.0
      DT  = 0.0
      HDT0= EMAX
      NCT = 0
      NBSTPT = 0
C
C
      DO 10 NV = 3, NBV
C
C        LE TRACE DE L'ERREUR NV-1
C        -------------------------
         NCOUL = NV-2
C
C        LE 1-ER SOMMET
         T0 = RMCN( MNT+K )
         PT(1,2) = T0
         PT(2,2) = REAL( DMCN( MNE+1 ) )
         CALL ENTIER2D(  NCOUL,  PT(1,2), PT(2,2), NV-2 )
ccc      CALL SYMBOLE2D( NCNOIR, PT(1,2), PT(2,2), CAR(NV-2) )
C
C        LE TRACE DES ( TEMPS(K), ERREUR(K,NV) )
         DO K = 2, NBC
C
C           LE POINT INITIAL EST LE POINT FINAL PRECEDENT
            PT(1,1) = PT(1,2)
            PT(2,1) = PT(2,2)
C
c           LE POINT FINAL DE L'ARETE
            T = RMCN( MNT+K )
            IF( T .LT. T0 ) GOTO 10
            PT(1,2) = T
            PT(2,2) = REAL( DMCN( MNE+K ) )
            CALL TRAIT2D(   NCOUL, PT(1,1), PT(2,1), PT(1,2), PT(2,2) )
ccc         CALL SYMBOLE2D( NCNOIR, PT(1,2), PT(2,2), CAR(NV-2) )
C
C           LE PAS DE TEMPS CHANGE T IL?
            IF( T .NE. T0 ) THEN
               NBSTPT = NBSTPT + 1
               DT = T - T0
               T0 = T
            ENDIF
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

         ENDDO
         CALL ENTIER2D(  NCOUL,  PT(1,2), PT(2,2), NV-2 )
C
C        PASSAGE AU VECTEUR SUIVANT
         MNE = MNE + NBC
C
 10   CONTINUE
      CALL TRAIT2D( NCBLAN, TMIN, EMIN, TMAX, EMIN )

C     L'AXE DES PAS DE TEMPS AU COURS DU TEMPS
      TITRE='TIME'
      NBC = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLAN, TMAX, EMAX, TITRE(1:NBC) )
      CALL TRAIT2D( NCBLAN, TMIN, EMAX, TMAX, EMAX )
      CALL TRAIT2D( NCBLAN, TMIN, EMAX, TMIN, EMAX+EH )
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 25 )

      TITRE = KDT0(1:LDT0) // '  Step Time(time)  STNb='
      KDT   = ' '
      WRITE( KDT(1:8), '(I8)' ) NBSTPT
      CALL TEXTSB( KDT, LDT0 )
      NB = NUDCNB( TITRE )
      TITRE(1:NB+LDT0) = TITRE(1:NB) // KDT(1:LDT0)
      NB = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLEU, TMIN+TH/6, EMAX+EH, TITRE(1:NB) )

      TITRE ='ERROR(time)=max|UExact-UComput|(N)/(max UExact(N)-min UExa
     %ct(N))'
      NBC = NUDCNB( TITRE )
      CALL TEXTE2D( NCBLEU, TMIN+TH/6, EMAX-EH, TITRE(1:NBC) )
C
      TITRE ='REAL PART ERROR Curve'
      NBC = NUDCNB( TITRE )
      CALL TEXTE2D( 1, TMIN+TH/6, EMAX-EH*1.5, TITRE(1:NBC) )
C
      TITRE ='IMAGINARY PART ERROR Curve '
      NBC = NUDCNB( TITRE )
      CALL TEXTE2D( 2, TMIN+TH/6, EMAX-EH*1.75, TITRE(1:NBC) )

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
C     MISE SUR FICHIER erreur_temps.eps du TRACE
C     ATTENTION NBC PASSAGE PAR VARIABLE OBLIGATOIRE
      NBC = NUDCNB( 'erreur_temps' )
      CALL xvsauverps( 'erreur_temps', NBC )
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
