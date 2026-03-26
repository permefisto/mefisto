      SUBROUTINE THEINS( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES TEMPERATURES ET FLUX DANS UN DOMAINE
C -----    2D OU 3D OU AXISYMETRIQUE EN THERMIQUE INSTATIONNAIRE
C          SELON LES DIFFERENTES METHODES PROGRAMMEES
C          ACTUELLEMENT: PAS CONSTANT ET THETA-METHODE
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET DE THERMIQUE INSTATIONNAIRE A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1998
C23456---------------------------------------------------------------012
      PARAMETER     (MXTYEL=7)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/donthe.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyznoeud.inc"
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
      REAL              DT,     DTSTOC, TPSINI, TPSFIN
      DOUBLE PRECISION  DINFO,  DCPU, DTEMP, BETA(0:1), PENALI
      DOUBLE PRECISION  RELMIN, D2PI
      DATA              RELMIN/-1D28/
C
C     INITIALISATION DU TEMPS CALCUL INITIAL
      DCPU  = DINFO( 'CPU' )
      D2PI  = ATAN(1D0) * 8D0
C
C     PENALISATION DE LA CONDITION DE DIRICHLET
      PENALI = 1D20
C
C     INSTANT ACTUEL DU CALCUL
      TEMPS  = 0.0
      MNTIME = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
      DTEMP  = 0D0
      DCPU   = 0D0
C
      MNNPEF = 0
      MNDOEL(1) = 0
      MNDOEL(2) = 0
      MNDOEL(3) = 0
      MNDOEL(4) = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNTHER = 0
      MNNODL = 0
      MONODL = 0
      MNX    = 0
      NBCOOR = 0
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('*')/
     %' RESOLUTION THERMIQUE INSTATIONNAIRE DE L''OBJET ',A/1X,80('*'))
20000 FORMAT(/1X,80('*')/
     %' UNSTEADY HEAT TRANSFER of the OBJECT ',A/1X,80('*'))
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
C
      CALL LIMTCL( 'resothin', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9999
      GOTO( 33, 34, 35, 30 ),NMTCL
C
C     THERMIQUE_INSTATIONNAIRE_LINEAIRE =>
C     CAPACITE CONDUCTIVITE COEFFICIENT D'ECHANGE INDEPENDANTS DU TEMPS ET DE LA
C     SEULS CONTACT ET SOURCE PEUVENT DEPENDRE DU TEMPS (ET NON DE LA TEMPERATUR
C     --------------------------------------------------------------------------
 33   TESTNL = 0
      GOTO 38
C
C     THERMIQUE_INSTATIONNAIRE_LINEAIRE =>
C     CAPACITE CONDUCTIVITE COEFFICIENT D'ECHANGE INDEPENDANTS DU TEMPS ET DE LA
C     SEULS CONTACT ET SOURCE PEUVENT DEPENDRE DU TEMPS ET DE LA TEMPERATURE
C     --------------------------------------------------------------------------
 34   TESTNL = 0
      GOTO 38
C
C     THERMIQUE_INSTATIONNAIRE NON_LINEAIRE
C     CAPACITE CONDUCTIVITE COEFFICIENT D'ECHANGE CONTACT SOURCE
C     PEUVENT DEPENDRE DU TEMPS ET DE LA TEMPERATURE
C     ----------------------------------------------------------
 35   TESTNL = 1
C
C     ENTREEE DES DONNEES THERMIQUES POUR L'INSTANT INITIAL
C     POUR UN PROBLEME THERMIQUE INSTATIONNAIRE D'ODRE 1 EN TEMPS
C     ===========================================================
 38   CALL THED1T( KNOMOB, MOREE2, NOAXIS,
     %             NTLXOB, MNTOPO, NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBTYEL, MNNPEF, NBTTEF,
     %             MNTPOB, NBDLMX, MOAUX,  MNTAUX,
     %             NUMIOB, NUMAOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXTYEL, MXDOEL, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             MNTHER, MOTAEL, MNTAEL, MNX,
     %             NORESO, MNLPLI, MNLPCO, NIVEAU, NBRDAG,
     %             NCODSK, NCODSM, MNUG0,
     %             DT,     DTSTOC, TPSINI, TPSFIN, NBVECT, MXVECT,
     %             NTDL,   NTVECT, MNVECT,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9990
      NBCOOR = MCN(MNXYZN+WBCOON)
C
C     ***************************************************************
C     INTEGRATION EN TEMPS PAR LA THETA(=BETA ICI) METHODE IMPLICITE
C     u(tn+1) - u(tn) = dt * (Somme i=0,1 de
C     ( betai [MG]**-1 ( {fg(tn+i)} - [KG] {u(tn+i)} )
C     avec [MG] la matrice de CAPACITE CALORIFIQUE
C          [KG] la matrice de CONDUCTIVITE
C     ***************************************************************
C
C     ENTREE DU PARAMETRE BETA COMPRIS ENTRE 1/2 ET 1 (COMPROMIS+2/3)
C     ===============================================================
      CALL INVITE( 1 )
      NCVALS  = 6
      BETA(1) = 2.D0/3.D0
      CALL LIRRDP( NCVALS, BETA(1) )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        IERR = 2
        RETURN
      ENDIF
      IF( BETA(1) .LT. 0.5D0 .OR. BETA(1) .GT. 1D0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'THETA<1/2 ou THETA>1 => THETA=1 IMPOSE'
         ELSE
            KERR(1) = 'THETA<1/2 or THETA>1 => THETA=1 IMPOSED'
         ENDIF
         CALL LEREUR
C        1 EST IMPOSE
         BETA(1) = 1D0
      ENDIF
      BETA(0) = 1D0 - BETA(1)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) BETA(1)
      ELSE
         WRITE(IMPRIM,20100) BETA(1)
      ENDIF
10100 FORMAT(' COEFFICIENT THETA DU SCHEMA EN TEMPS =',G15.7)
20100 FORMAT(' THETA COEFFICIENT of TIME SCHEMA =',G15.7)
C
C     TRAITEMENT SELON LA NATURE DES CARACTERISTIQUES THERMIQUES
C     ==========================================================
      GOTO( 100, 200, 300, 30 ), NMTCL
C
C     THERMIQUE_INSTATIONNAIRE_LINEAIRE
C     COEFFICIENTS INDEPENDANTS DU TEMPS ET DE LA TEMPERATURE
C     POUR LA CAPACITE ET LA CONDUCTIVITE ET COEFFICIENT D'ECHANGE
C     LE SECOND MEMBRE PEUT DEPENDRE DU TEMPS ET DE LA TEMPERATURE
C     ------------------------------------------------------------
C     LE SECOND MEMBRE PEUT DEPENDRE DU TEMPS SEULEMENT
 100  TESTNL = 0
      GOTO 210
C
C     LE SECOND MEMBRE PEUT DEPENDRE DU TEMPS ET DE LA TEMPERATURE
 200  TESTNL = 1
 210  CALL THEILC( KNOMOB, MOREE2, D2PI,
     %             NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %             NUMIOB, NUMAOB, MNDOEL, IEVIFL, IESOPO,
     %             PENALI, RELMIN,
     %             MNTHER, MNTAEL, MNX,
     %             NORESO, MNLPLI, MNLPCO, NIVEAU, NBRDAG,
     %             NCODSM, MNUG0,
     %             BETA,   DT,     DTSTOC, TPSINI, TPSFIN,
     %             NBVECT, MXVECT,
     %             NTDL,   NTVECT, MNVECT,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9990
      GOTO 900
C
C     THERMIQUE_INSTATIONNAIRE_LINEAIRE
C     COEFFICIENTS DEPENDANTS DU TEMPS ET DE LA TEMPERATURE
C     LINEARISATION PAR POINT FIXE
C     -----------------------------------------------------
 300  TESTNL = 1
      CALL THEINL( KNOMOB, MOREE2, D2PI,
     %             NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBTYEL, MNNPEF,
     %             MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %             PENALI, RELMIN, MNTHER, MNTAEL, MNX,
     %             NORESO, MNLPLI, MNLPCO, NIVEAU, NBRDAG,
     %             NCODSM, MNUG0,
     %             BETA,   DT,     DTSTOC, TPSINI, TPSFIN,
     %             NBVECT, MXVECT, NTDL,   NTVECT, MNVECT,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9990
C
C     COUT CALCUL DES TEMPERATURES
C     ============================
 900  DTEMP = DINFO( 'DELTA CPU' )
C
C     CALCUL DES FLUX EN CHAQUE POINT D INTEGRATION DES FACES DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     ==========================================================
      CALL THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     %             NDIM,   MOREE2, MXVECT, NTDL,
     %             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %             MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MOTAEL, MNTAEL, MNX,    MNVECT )
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
C     =========================================
 9990 IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF )
      DO 9991 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
 9991 CONTINUE
      IF( MNTHDL .GT. 0 ) CALL NLDATADS
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX,  MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL, MNTAEL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2',  128,    MNTHER )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', MONODL, MNNODL )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL' ,  NBDLMX*NBCOOR, MNX )
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA,MNTPOB )
C
C     BILAN SUR LE TEMPS CALCUL UTILISE
C     =================================
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19000) DTEMP, DCPU, DTEMP+DCPU
      ELSE
         WRITE(IMPRIM,29000) DTEMP, DCPU, DTEMP+DCPU
      ENDIF
19000 FORMAT(/
     %' TEMPS CALCUL DES TEMPERATURES =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES FLUX         =',F12.2,' SECONDES CPU'/
     %' TEMPS TOTAL  DE LA RESOLUTION =',F12.2,' SECONDES CPU'/)
29000 FORMAT(/
     % ' COMPUTATION TEMPERATURE TIME =',F12.2,' CPU SECONDS'/
     % ' COMPUTATION HEAT FLUX   TIME =',F12.2,' CPU SECONDS'/
     % ' SOLUTION TOTAL          TIME =',F12.2,' CPU SECONDS'/)
C
C     SORTIE
 9999 RETURN
      END
