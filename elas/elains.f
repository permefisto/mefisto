      SUBROUTINE ELAINS( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES DEPLACEMENTS ET CONTRAINTES DANS UN DOMAINE
C -----    2D OU 3D OU AXISYMETRIQUE EN ELASTICITE INSTATIONNAIRE
C          SELON LES DIFFERENTES METHODES PROGRAMMEES
C          ACTUELLEMENT: PAS CONSTANT ET SCHEMA DE NEWMARK
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET DE ELASTICITE INSTATIONNAIRE A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       MARS 1999
C23456---------------------------------------------------------------012
      PARAMETER     (MXTYEL=7)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/donthe.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___vecteur.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      CHARACTER*(*)     KNOMOB
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      DOUBLE PRECISION  DINFO,  DCPU, DDEPL, PENALI, DT
      DOUBLE PRECISION  DEXMAX, DECMAX
      DOUBLE PRECISION  AMORTM, AMORTK
      DOUBLE PRECISION  RELMIN, D2PI
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('*')/
     %' RESOLUTION DE L''ELASTICITE INSTATIONNAIRE DE L''OBJET: ',
     %A/1X,80('*')/)
20000 FORMAT(/1X,80('*')/
     %' COMPUTATION of the UNSTEADY ELASTICITY of the OBJECT: ',
     %A/1X,80('*')/)
C
C     2 PI
      D2PI  = ATAN(1D0) * 8D0
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PENALISATION DE LA CONDITION DE DEPLACEMENT IMPOSE
      PENALI = 1D30
C
C     INSTANT ACTUEL DU CALCUL
      TEMPS  = 0.0
      MNTIME = 0
C
C     LES TEMPS CALCUL DES DIFFERENTES ETAPES DU CALCUL
      DDEPL  = 0D0
      DCPU   = 0D0
C
C     LES ADRESSES DANS MCN
      MNNPEF = 0
      DO 5 I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
 5    CONTINUE
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNELAS = 0
      MNFORC = 0
      MNX    = 0
      MNLPLK = 0
      MNLPCK = 0
      MNUG   = 0
      MNTEM0 = 0
      MNIP   = 0
      MNNODL = 0
C
C     LECTURE DES MOTS CLE
C     ====================
C     L'ANCIEN HISTORIQUE EST EFFACE
 10   CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
      CALL LIMTCL( 'resoelin', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9999
C     SEUL LE CAS DES COEFFICIENTS INDEPENDANT du TEMPS etdu DEPLACEMENT
C     EST ACTUELLEMENT PROGRAMME
      NMTCL=1
      GOTO( 20, 30, 30, 10 ),NMTCL
C
C     ELASTICITE_INSTATIONNAIRE_LINEAIRE =>
C     MASSE YOUNG INDEPENDANTS DU TEMPS ET DU DEPLACEMENT
C     SEUL SOURCE PEUT DEPENDRE DU TEMPS (ET NON DU DEPLACEMENT)
C     ----------------------------------------------------------
 20   TESTNL = 0
      GOTO 40
C
C     ELASTICITE_INSTATIONNAIRE_LINEAIRE
C     COEFFICIENTS DEPENDANTS DU TEMPS ET INDEPENDANTS DES DEPLACEMENTS
 30   NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'OPTION NON PROGRAMMEE'
         KERR(2) = 'SEULE OPTION 1 ACTIVE ACTUELLEMENT'
      ELSE
         KERR(1) = 'NOT PROGRAMMED OPTION'
         KERR(2) = 'ONLY OPTION 1 ACTIVE ACTUALLY'
      ENDIF
      CALL LEREUR
      NMTCL = 1
      GOTO 20
C
C     ENTREEE DES DONNEES ELASTIQUES POUR L'INSTANT INITIAL
C     POUR UN PROBLEME ELASTIQUE INSTATIONNAIRE D'ORDRE 2 EN TEMPS
C     ============================================================
C     INITIALISATION DU TEMPS CALCUL INITIAL
 40   DCPU  = DINFO( 'DELTA CPU' )
      CALL ELAD2T( KNOMOB, MOREE2, NOAXIS,
     %             NTLXOB, MNTOPO, NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBCODE, NBTYEL, MNNPEF, NBTTEF,
     %             MNTPOB, NBDLMX, MOAUX,  MNTAUX,
     %             NUMIOB, NUMAOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXTYEL, MXDOEL, MNDOEL,
     %             IEMASS, IEYOUN, IEDILA, IECOED, IECOIN, IEFOIN,
     %             IEFIXA, IEFOCL, IEFOPO,
     %             AMORTM, AMORTK,
     %             MNELAS, MNFORC, MOTAEL, MNTAEL, MNX,
     %             MNIP,   MNNODL,
     %             NORESO, MNLPLK, MNLPCK, NIVEAU, NBRDKG,
     %             NCODMG, NBRDMG,
     %             TPSINI, TPSFIN, DT,     NBDEPL, TESTDE, DTSTDE,
     %             NTVEDE, MNVEDE, NTDLDE, MNUG,
     %             NOTHEL, NTVETE, MNVETE, NBTEMP, NTDLTE, MNTEM0,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 50
C
C     ELASTICITE_INSTATIONNAIRE_LINEAIRE
C     MASSE YOUNG INDEPENDANTS DU TEMPS ET DES DEPLACEMENTS
C     FORCE PEUT DEPENDRE DU TEMPS MAIS PAS DU DEPLACEMENT
C     =====================================================
      CALL ELAILC( KNOMOB, MOREE2, D2PI,
     %             NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBCODE, NBTYEL, MNNPEF,
     %             MNTPOB, MNTAUX,
     %             NUMIOB, NUMAOB,
     %             MNDOEL, PENALI, RELMIN, AMORTM, AMORTK,
     %             MNELAS, MNFORC, MNTAEL, MNX,
     %             MNIP,   MNNODL,
     %             NORESO, MNLPLK, MNLPCK, NIVEAU, NBRDKG,
     %             NCODMG, NBRDMG,
     %             TPSINI, TPSFIN, DT,     NBDEPL, TESTDE, DTSTDE,
     %             NTVEDE, MNVEDE, NTDLDE, MNUG,
     %             NOTHEL, MNVETE, NBTEMP, NTDLTE, MNTEM0,
     %             NPIMAX, IERR )
 50   IF( MNLPLK .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDLDE+1, MNLPLK )
      IF( MNUG   .GT. 0 ) CALL TNMCDS( 'REEL2', 3*NTDLDE,  MNUG   )
      IF( MNTEM0 .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDLTE,   MNTEM0 )
      IF( IERR .NE. 0 ) GOTO 9990
C
C     AFFICHAGE DU DERNIER VECTEUR"DEPLACT CALCULE
C     ============================================
      CALL AFDEPL( 10,     NBDEPL, MNXYZN,
     %             NBCODE, MCN(MNXYZN+WNBNOE),
     %             NBDEPL, MCN(MNVEDE+WECTEU),
     %             DECMAX, NOFOTI, DEXMAX )
C
C     COUT CALCUL DES DEPLACEMENTS
C     ============================
      DDEPL = DINFO( 'DELTA CPU' )
C
C     CALCUL DES CONTRAINTES EN CHAQUE POINT D INTEGRATION DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     =======================================================
      IF( IECOIN .GT. 0 ) THEN
         IECOIN = 1
      ELSE
         IECOIN = 0
      ENDIF
      CALL ELASTR( NTLXOB, MNTOPO, NOAXIS, NDIM,   MOREE2,
     %             NPIMAX, NBDEPL, IECOIN, NTDLDE,
     %             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %             MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNELAS, MOTAEL, MNTAEL, MNX,
     %             MNVEDE, NOTHEL, MNVETE,
     %             DECMAX, DEXMAX )
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
C     =========================================
 9990 IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF )
      DO 9991 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
 9991 CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA, MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2' , MOAUX ,   MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2' , MOTAEL,   MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX,   MNNODL )
      IF( MNIP   .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX,   MNIP   )
      IF( MNELAS .GT. 0 ) CALL TNMCDS( 'REEL2' , 128   ,   MNELAS )
      IF( MNFORC .GT. 0 ) CALL TNMCDS( 'REEL2' , 6     ,   MNFORC )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX*3, MNX  )
      IF( MNTHDL .GT. 0 ) CALL TNMCDS( 'REEL2',  NBDLMX,   MNTHDL )
C
C     BILAN SUR LE TEMPS CALCUL UTILISE
C     =================================
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19999) DDEPL, DCPU, DDEPL+DCPU
      ELSE
         WRITE(IMPRIM,29999) DDEPL, DCPU, DDEPL+DCPU
      ENDIF
19999 FORMAT(/' TEMPS CALCUL DES DEPLACEMENTS =',F12.2,' SECONDES CPU'
     %       /' TEMPS CALCUL DES CONTRAINTES  =',F12.2,' SECONDES CPU'
     %       /' TEMPS TOTAL  DE LA RESOLUTION =',F12.2,' SECONDES CPU')
29999 FORMAT(/' DISPLACEMENT COMPUTATION TIME =',F12.2,' CPU SECONDS'
     %       /' STRESS       COMPUTATION TIME =',F12.2,' CPU SECONDS'
     %       /' SOLUTION     TOTAL       TIME =',F12.2,' CPU SECONDS')
C
C     SORTIE
 9999 RETURN
      END
