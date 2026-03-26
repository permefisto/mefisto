      SUBROUTINE SDTHER( KNOMST , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES TEMPERATURES ET LES FLUX DANS UN DOMAINE
C ----- 2D OU 3D OU AXISYMETRIQUE EN THERMIQUE LINEAIRE STATIONNAIRE
C       STATIONNAIRE PAR LA METHODE DES SOUS-DOMAINES
C
C ENTREES :
C ---------
C KNOMST : NOM DE L'OBJET
C
C SORTIES :
C ---------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  FEVRIER 1990
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / MSIMTA / NOIMPR
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___interface.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION  DMCN(1),PROSCD,D0,DN,DNP1,RN,RHON,RAPPORT,EPS0
      DOUBLE PRECISION  R0
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*(*)     KNOMST
      CHARACTER*24      KNOMOB,KNOMSD,KNOM
      CHARACTER*12      KNOMAT
      EQUIVALENCE       (MCN(1),DMCN(1))
      DATA              KNOMAT/'CONDUCTIVITE'/
C
      IERR   = 0
C
C     L'OBJET
C     =======
C
C     NOM DE L'OBJET "SOUS-DOMAINES"
      L = INDEX( KNOMST , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMST )
      ENDIF
      KNOMSD = KNOMST(1:L) // '_SD'
C     RECHERCHE DE L'OBJET KNOMSD DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE , KNOMSD , NTLXSD , MNLXSD )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOMSD
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU DEFINITION DE L'OBJET KNOMSD
      CALL LXTSOU( NTLXSD , 'DEFINITION' , NTDFSD , MNDFSD )
C     IL N'EXISTE PAS
      IF( NTDFSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : DEFINITION INCONNUE OBJET ' // KNOMSD
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C
C     L'INTERFACE
C     ===========
C
C     LE TABLEAU INTERFACE DE L'OBJET KNOMSD
      CALL LXTSOU( NTLXSD , 'INTERFACE' , NTINSD , MNINSS )
C     IL N'EXISTE PAS
      IF( NTINSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : INTERFACE INCONNU OBJET ' // KNOMSD
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C     LE NOMBRE DE NOEUDS DE L'INTERFACE
      NBNOIN = MCN(MNINSS+WBNOIN)
      MNOBNO = MNINSS+WBOBNO
      MNNUNO = MNINSS+WBOBNO+NBNOIN
C
C     CHOIX DE LA METHODE DE RESOLUTION
C     =================================
C
C     LA METHODE DE CHOLESKY EST IMPOSEE ( STOCKAGE PROFIL )
      NORESO = 1
C     LE NOMBRE DE D'INCONNUES PAR NOEUDS
      NBINCO = 1
C
C     1) PREPARATION DU CALCUL PAR SOUS-DOMAINE
C     =========================================
C
C     LE NOMBRE DE SOUS-DOMAINES
      NBDOBJ = MCN(MNDFSD+WBDOBJ)
      MNDF = MNDFSD + WTYOBJ
C     LE TABLEAU DES POIDS SUR L'INTERFACE
      MNPOID = 0
      LOPOID = NBNOIN * NBINCO
      CALL TNMCDC( 'REEL2' , LOPOID , MNPOID )
      CALL AZEROD( LOPOID , MCN(MNPOID) )
C     ADRESSE ET TAILLE DES TABLEAUX UTILITAIRES
      MNADR = 0
      LOADR = NBDOBJ * 4
      CALL TNMCDC( 'ENTIER' , LOADR , MNADR )
      CALL AZEROI( LOADR , MCN(MNADR) )
      MNVAUX = MNADR
      MNLAUX = MNADR  + NBDOBJ
      MNINTE = MNLAUX + NBDOBJ
      MNMNPO = MNINTE + NBDOBJ
      MOREE2 = MOTVAR(6)
      NTDLO  = 0
      NBVECT = 8
C
      DO 100 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C
C           PREPARATION : ASSEMBLAGE ET FACTORISATION
C           -----------------------------------------
C           L'ADRESSE DE L'INTERFACE DU SOUS-DOMAINE
            MNINSD = 0
            LOINSD =  NBNOIN * 2 + 1
            CALL TNMCDC( 'ENTIER' , LOINSD , MNINSD )
            CALL AZEROI( LOINSD , MCN(MNINSD) )
C           STOCKAGE DE L'ADRESSE
            MCN(MNINTE+NO) = MNINSD
C           L'ADRESSE DES POIDS DU SOUS-DOMAINE
            MNPOIL = 0
            LOPOIL =  NBNOIN * NBINCO
            CALL TNMCDC( 'REEL2' , LOPOIL , MNPOIL )
            CALL AZEROD( LOPOIL , MCN(MNPOIL) )
C           STOCKAGE DE L'ADRESSE
            MCN(MNMNPO+NO) = MNPOIL
            IF (NOIMPR.GE.10) WRITE(IMPRIM,4000) KNOM
            CALL SDTHE1( NTLXOB , NUTYSD , NUOBSD , MNINSD ,
     &                   MNPOID , MNPOIL , MNDLTR , NBNOIN ,
     &                   MNOBNO , MNNUNO , IERR )
C           LES ADRESSES DES VECTEURS AUXILIAIRES DU SOUS-DOMAINE
            MNVX = 0
            LOVX =  NTDL * NDSM * NBVECT
            CALL TNMCDC( 'REEL2' , LOVX , MNVX )
            CALL AZEROD( LOVX , MCN(MNVX) )
C           STOCKAGE DE L'ADRESSE ET DE LA LONGUEUR
            MCN(MNVAUX+NO) = MNVX
            MCN(MNLAUX+NO) = NTDL
C           STOCKAGE TEMPORAIRE DE L'ADRESSE DES C.L. SUR L'INTERFACE
            MCN(MNVX)      = MNDLTR
            NTDLO  = NTDLO + NTDL
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 100  CONTINUE
C
C     LES TABLEAUX SUR L'INTERFACE
C     ============================
C
C     LE TABLEAU RESIDU SUR L'INTERFACE
      MNRSIN = 0
      LORSIN = NBNOIN * NBINCO * NDSM
      CALL TNMCDC( 'REEL2' , LORSIN , MNRSIN )
C     LE TABLEAU DE LA TRACE SUR L'INTERFACE
      MNTRIN = 0
      LOTRIN = NBNOIN * NBINCO * NDSM
      CALL TNMCDC( 'REEL2' , LOTRIN , MNTRIN )
C
C     INITIALISATION DE LA TRACE SUR L'INTERFACE :
C     MISE A ZERO PUIS PRISE EN COMPTE DES CONDITIONS AUX LIMITES
C     -----------------------------------------------------------
      CALL AZEROD( LOTRIN , MCN(MNTRIN) )
      MNDF   = MNDFSD + WTYOBJ
      DO 115 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
C           CONTRIBUTION DES C.L. DU SOUS-DOMAINE
            MNDLTR = MCN(MNVECT)
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES5( MCN(MNDLTR) , MCN(MNTRIN) , NBNOIN , NBNOSD ,
     &              MCN(MNINS1) , MCN(MNINS2) , MCN(MNOBNO) , NBINCO )
CC            CALL SDRES7( MCN(MNDLTR) , MCN(MNTRIN) , NBNOIN , NBNOSD ,
CC     &              MCN(MNINS1) , MCN(MNINS2) ,
CC     &              MCN(MNPOIL) , MCN(MNPOID) , NBINCO )
C           DESTRUCTION DU TABLEAU DES C.L. SUR L'INTERFACE
            CALL TNMCDS( 'REEL2' , NTDL*NDSM , MNDLTR )
         END IF
 115  CONTINUE
C=========================================================================
C
C     AFFICHAGE DE LA TRACE INITIALE
C     ------------------------------
      IF (NOIMPR.GE.10) THEN
         NTDLI = NBNOIN*NBINCO
         IA = ( MNTRIN - 1 ) / 2
         WRITE(IMPRIM,10116)
10116    FORMAT(/' LA TRACE INITIALE :',/1X,80(1H=))
         DO 116 I=1,NTDLI
         WRITE(IMPRIM,10117) (I,DMCN(IA+I+(J-1)*NBNOIN*NBINCO),J=1,NDSM)
 116     CONTINUE
10117    FORMAT( ' DL',I6,' : TRACE =',G16.7)
      END IF
C=========================================================================
C
C     2) CALCUL DES TEMPERATURES
C     ==========================
C
C     2.1) INITIALISATION
C     -------------------
C
C     2.1.1) RESOLUTION DES PROBLEMES DE DIRICHLET
C
      MNDF = MNDFSD + WTYOBJ
      DO 101 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C           LE VECTEUR TEMPERATURE DE L'OBJET
            CALL LXTSOU( NTLXOB , 'VECTEUR"TEMPERATURE' , NTU , MNU )
            MNVTNO = MNU + WECTEU
C           LE VECTEUR SECONDE MEMBRE DE L'OBJET
            CALL LXTSOU( NTLXOB , 'VECTEUR"B' , NTB , MNB )
            MNVBNO = MNB + WECTEU
C
C           INITIALISATION DE LA TRACE LAMBDA
C           ---------------------------------
            MNLANO = MNVECT
            CALL AZEROD( NTDL * NDSM , MCN(MNLANO) )
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES4( MCN(MNTRIN) , MCN(MNLANO) , NBNOIN ,
     &                   NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C
C           RESOLUTION DU PROBLEME DE DIRICHLET
C           SECOND MEMBRE B + C.L. DIRICHLET SUR L'INTERFACE
C           ------------------------------------------------
            IF (NOIMPR.GE.10) WRITE(IMPRIM,5000) KNOM
C           LE SECOND MEMBRE
            MNLINO = MNVECT + LAUX
            CALL SDRES1( NTLXOB , MNLANO , MNVBNO , MNLINO ,
     &                   KNOMAT , IERR )
C           LES CONDITIONS AUX LIMITES  SUR L'INTERFACE
            CALL SDRES4( MCN(MNTRIN) , MCN(MNLINO) , NBNOIN ,
     &                   NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C           RESOLUTION DU SYSTEME LINEAIRE
            CALL SDTHE2( NTLXOB , MNVTNO , MNLINO , IERR )
            IF (IERR.NE.0) GO TO 9900
C
C=========================================================================
C           AFFICHAGE DES TEMPERATURES INITIALES
C           ------------------------------------
            IF (NOIMPR.GE.10) THEN
                NTDLI = NTDL
                IA = ( MNVTNO - 1 ) / 2
                WRITE(IMPRIM,10799)
10799           FORMAT(/' LES TEMPERATURES INITIALES :',/1X,80(1H=))
                DO 799 I=1,NTDLI
                WRITE(IMPRIM,10800) (I,DMCN(IA+I+(J-1)*NTDL),J=1,NDSM)
 799            CONTINUE
10800           FORMAT( ' DL',I6,' : TEMPERATURE =',G16.7)
            ENDIF
C=========================================================================
C
C           CALCUL DU RESIDU
C           ----------------
            CALL SDRES1( NTLXOB , MNVTNO , MNVBNO , MNLINO ,
     &                   KNOMAT , IERR )
C           CHANGEMENT DE SIGNE !
            CALL CHSGND( NTDL * NDSM , MCN(MNLINO) )
C
C=========================================================================
C           AFFICHAGE DU RESIDU INITIAL
C           ---------------------------
            IF (NOIMPR.GE.10) THEN
                NTDLI = NTDL
                IA = ( MNLINO - 1 ) / 2
                WRITE(IMPRIM,11799)
11799           FORMAT(/' LES RESIDUS INITIAUX :',/1X,80(1H=))
                DO 1799 I=1,NTDLI
                WRITE(IMPRIM,11800) (I,DMCN(IA+I+(J-1)*NTDL),J=1,NDSM)
1799            CONTINUE
11800           FORMAT( ' DL',I6,' : RESIDU =',G16.7)
            ENDIF
C=========================================================================
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 101  CONTINUE
C
C     2.1.2) MOYENNE DU RESIDU SUR L'INTERFACE
C
C     MISE A ZERO AVANT SOMMATION
      CALL AZEROD( LORSIN , MCN(MNRSIN) )
      MNDF   = MNDFSD + WTYOBJ
      DO 102 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           CONTRIBUTION DU RESIDU DU SOUS-DOMAINE
C           --------------------------------------
            MNLINO = MNVECT + LAUX
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES7( MCN(MNLINO) , MCN(MNRSIN) , NBNOIN , NBNOSD ,
     &              MCN(MNINS1) , MCN(MNINS2) ,
     &              MCN(MNPOIL) , MCN(MNPOID) , NBINCO )
         END IF
 102  CONTINUE
C
C=========================================================================
C     VERIFICATION DU RACCORD
C     -----------------------
      IF (NOIMPR.GE.10) THEN
         RACCO = REAL( PROSCD( MCN(MNRSIN) , MCN(MNRSIN) , LORSIN ) )
         RACCO = SQRT(RACCO) / NBNOIN
         WRITE (IMPRIM,2345) RACCO
 2345    FORMAT(/,1X,80(1H*)/,' INITIALISATION : RACCORD ',
     S         E15.6/1X,80(1H*)/)
      ENDIF
C=========================================================================
C
C     2.1.3) RESOLUTION DES PROBLEMES DE NEUMANN
C
      DN   = 0.D0
      R0   = 0.D0
      MNDF = MNDFSD + WTYOBJ
      DO 103 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           TRANSFERT DU RESIDU
C           -------------------
            MNLINO = MNVECT + LAUX
            CALL AZEROD( NTDL * NDSM , MCN(MNLINO) )
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES4( MCN(MNRSIN) , MCN(MNLINO) , NBNOIN ,
     &                    NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C
C           RESOLUTION DU PROBLEME DE NEUMANN
C           ---------------------------------
            MNFINO = MNVECT + LAUX * 2
            CALL SDTHE3( NTLXOB , MNFINO , MNLINO , IERR )
            IF (IERR.NE.0) GO TO 9900
C
C           CALCUL DE DN
C           ------------
            DN = DN + PROSCD( MCN(MNFINO) , MCN(MNLINO) , NTDL*NDSM )
            R0 = R0 + PROSCD( MCN(MNLINO) , MCN(MNLINO) , NTDL*NDSM )
C
C           SAUVEGARDE DU RESULTAT
C           ----------------------
            MNWINO = MNVECT + LAUX * 3
            CALL TRTATA( MCN(MNFINO) , MCN(MNWINO) , LAUX )
            MNRINO = MNVECT + LAUX * 4
            CALL TRTATA( MCN(MNLINO) , MCN(MNRINO) , LAUX )
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 103  CONTINUE
      D0 = DN
C
C     2.2) ITERATIONS
C     ---------------
C
C     LE NOMBRE MAXIMUM D'ITERATIONS
      ITEMAX = NTDLO
C     LA PRECISION DEMANDEE
      EPS  = EPZERO * EPZERO
      EPS0 = EPS * D0
      IF(EPS0.LT.EPS) EPS0  = EPS
      ITER = 0
 200  ITER = ITER + 1
C
C     2.2.1) MOYENNE DE LA TRACE SUR L'INTERFACE
C
C     MISE A ZERO AVANT SOMMATION
      CALL AZEROD( LOTRIN , MCN(MNTRIN) )
      MNDF   = MNDFSD + WTYOBJ
      DO 104 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           CONTRIBUTION DE LA TEMPERATURE DU SOUS-DOMAINE
C           ----------------------------------------------
            MNWINO = MNVECT + LAUX * 3
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES5( MCN(MNWINO) , MCN(MNTRIN) , NBNOIN , NBNOSD ,
     &              MCN(MNINS1) , MCN(MNINS2) , MCN(MNOBNO) , NBINCO )
CC            CALL SDRES7( MCN(MNWINO) , MCN(MNTRIN) , NBNOIN , NBNOSD ,
CC     &              MCN(MNINS1) , MCN(MNINS2) ,
CC     &              MCN(MNPOIL) , MCN(MNPOID) , NBINCO )
C
         END IF
 104  CONTINUE
C
C     2.2.1) RESOLUTION DES PROBLEMES DE DIRICHLET
C
      MNDF = MNDFSD + WTYOBJ
      DO 105 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           TRANSFERT DE LA TRACE
C           ---------------------
            MNLANO = MNVECT
C           MISE A ZERO DE LA TRACE
            CALL AZEROD( NTDL * NDSM , MCN(MNLANO) )
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            CALL SDRES4( MCN(MNTRIN) , MCN(MNLANO) , NBNOIN ,
     &                   NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C
C=========================================================================
C            AFFICHAGE DE LA TRACE
C            ---------------------
             IF (NOIMPR.GE.10) THEN
                NTDLI = NBNOIN*NBINCO
                IA = ( MNTRIN - 1 ) / 2
                WRITE(IMPRIM,14126)
14126           FORMAT(/' LA TRACE SUR L''INTERFACE :',/1X,80(1H=))
                DO 4216 I=1,NTDLI
                WRITE(IMPRIM,10117) (I,DMCN(IA+I+(J-1)*NTDL),J=1,NDSM)
4216            CONTINUE
             END IF
C=========================================================================
C
C           RESOLUTION DU PROBLEME DE DIRICHLET :
C           SECOND MEMBRE NUL + C.L. DIRICHLET SUR L'INTERFACE
C           --------------------------------------------------
C           LE SECOND MEMBRE
            MNLINO = MNVECT + LAUX
            MNLJNO = MNVECT + LAUX * 7
            CALL AZEROD( NTDL * NDSM , MCN(MNLINO) )
            CALL SDRES1( NTLXOB , MNLANO , MNLINO , MNLJNO ,
     &                   KNOMAT ,  IERR )
C           LES CONDITIONS AUX LIMITES SUR L'INTERFACE
            CALL SDRES4( MCN(MNTRIN) , MCN(MNLJNO) , NBNOIN ,
     &                   NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C           RESOLUTION DU SYSTEME LINEAIRE
            MNZINO = MNVECT + LAUX * 5
            CALL SDTHE2( NTLXOB , MNZINO , MNLJNO , IERR )
            IF (IERR.NE.0) GO TO 9900
C
C           CALCUL DU RESIDU
C           ----------------
            CALL SDRES1( NTLXOB , MNZINO , MNLJNO , MNLINO ,
     &                   KNOMAT ,  IERR )
C           CHANGEMENT DE SIGNE !
            CALL CHSGND( NTDL * NDSM , MCN(MNLINO) )
C
C=========================================================================
C           AFFICHAGE DU RESIDU
C           -------------------
             IF (NOIMPR.GE.10) THEN
                NTDLI = NTDL
                IA = ( MNLJNO - 1 ) / 2
                WRITE(IMPRIM,12799)
12799           FORMAT(/' LES RESIDUS :',/1X,80(1H=))
                DO 2799 I=1,NTDLI
                WRITE(IMPRIM,11800) (I,DMCN(IA+I+(J-1)*NTDL),J=1,NDSM)
2799            CONTINUE
            ENDIF
C=========================================================================
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 105  CONTINUE
C
C     2.2.3) MOYENNE DU RESIDU SUR L'INTERFACE
C
C     MISE A ZERO AVANT SOMMATION
      CALL AZEROD( LORSIN , MCN(MNRSIN) )
      MNDF   = MNDFSD + WTYOBJ
      DO 106 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           CONTRIBUTION DU RESIDU DU SOUS-DOMAINE
C           --------------------------------------
            MNLINO = MNVECT + LAUX
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            MNPOIL = MCN(MNMNPO+NO)
            CALL SDRES7( MCN(MNLINO) , MCN(MNRSIN) , NBNOIN , NBNOSD ,
     &              MCN(MNINS1) , MCN(MNINS2) ,
     &              MCN(MNPOIL) , MCN(MNPOID) , NBINCO )
C
         END IF
 106  CONTINUE
C
C=========================================================================
C     VERIFICATION DU RACCORD
C     -----------------------
      IF (NOIMPR.GE.10) THEN
         RACCO = REAL( PROSCD( MCN(MNRSIN) , MCN(MNRSIN) , LORSIN ) )
         RACCO = SQRT(RACCO)
         WRITE (IMPRIM,2346) ITER,RACCO
 2346    FORMAT(/,1X,80(1H*)/,' ITERATION',I6,' : RACCORD',
     S      E15.6/1X,80(1H*)/)
      ENDIF
C=========================================================================
C
C     2.2.3) RESOLUTION DES PROBLEMES DE NEUMANN
C
      RN   = 0.D0
      MNDF = MNDFSD + WTYOBJ
      DO 107 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
C
C           TRANSFERT DU RESIDU
C           -------------------
            MNLINO = MNVECT + LAUX
            CALL AZEROD( NTDL * NDSM , MCN(MNLINO) )
            MNINSD = MCN(MNINTE+NO)
            NBNOSD = MCN(MNINSD)
            MNINS1 = MNINSD+1
            MNINS2 = MNINS1+NBNOIN
            CALL SDRES4( MCN(MNRSIN) , MCN(MNLINO) , NBNOIN ,
     &                   NBNOSD , MCN(MNINS1) , MCN(MNINS2) , NBINCO )
C
C           RESOLUTION DU PROBLEME DE NEUMANN
C           ---------------------------------
            MNPINO = MNVECT + LAUX * 6
            CALL SDTHE3( NTLXOB , MNPINO , MNLINO , IERR )
            IF (IERR.NE.0) GO TO 9900
C
C           CALCUL DE RN
C           ------------
            MNWINO = MNVECT + LAUX * 3
            RN = RN + PROSCD( MCN(MNWINO) , MCN(MNLINO) , NTDL*NDSM )
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 107  CONTINUE
C
C     2.2.4) CALCUL DE RHO
C
      IF (ABS(RN) .GT. EPZERO) THEN
         RHON = - DN / RN
      ELSE
         GO TO 201
      ENDIF
C
C     2.2.5) MISE A JOUR DES VECTEURS
C
      MNDF = MNDFSD + WTYOBJ
      DNP1 = 0.D0
      DO 108 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C
C           LE TABLEAU TEMPERATURE
C           ----------------------
            CALL LXTSOU( NTLXOB , 'VECTEUR"TEMPERATURE' , NTU , MNU )
            MNVTNO = MNU + WECTEU
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
            MNLINO = MNVECT + LAUX
            MNFINO = MNVECT + LAUX * 2
            MNRINO = MNVECT + LAUX * 4
            MNZINO = MNVECT + LAUX * 5
            MNPINO = MNVECT + LAUX * 6
            CALL SDRES2( RHON , MCN(MNVTNO) , MCN(MNZINO) ,
     &                          MCN(MNFINO) , MCN(MNPINO) ,
     &                          MCN(MNRINO) , MCN(MNLINO) )
C
C=========================================================================
C           AFFICHAGE DES TEMPERATURES
C           --------------------------
            IF (NOIMPR.GE.10) THEN
                NTDLI = NTDL
                IA = ( MNVTNO - 1 ) / 2
                WRITE(IMPRIM,10789)
10789           FORMAT(/' LES TEMPERATURES ACTUELLES :',/1X,80(1H=))
                DO 789 I=1,NTDLI
                WRITE(IMPRIM,10800) (I,DMCN(IA+I+(J-1)*NTDL),J=1,NDSM)
 789            CONTINUE
            ENDIF
C=========================================================================
C
C           CALCUL DE DN+1
C           --------------
            DNP1 = DNP1 + PROSCD( MCN(MNRINO) , MCN(MNFINO) ,
     &             NTDL*NDSM )
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 108  CONTINUE
C
C     2.2.6) TEST DE CONVERGENCE
C
      IF (ABS(DNP1) .LT. EPS0 ) THEN
         GO TO 201
      ELSE
         IF (ITER.EQ.ITEMAX) THEN
            WRITE(IMPRIM,2000) ITER
            WRITE(IMPRIM,3000) DNP1
            GO TO 202
         END IF
      END IF
C
C     2.2.7) CALCUL DE LA NOUVELLE DIRECTION
C
      RAPPORT = DNP1 / DN
      DN      = DNP1
      MNDF = MNDFSD + WTYOBJ
      DO 109 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            MNVECT = MCN(MNVAUX+NO)
            NTDL   = MCN(MNLAUX+NO)
            LAUX   = NTDL * NDSM * MOREE2
            MNFINO = MNVECT + LAUX * 2
            MNWINO = MNVECT + LAUX * 3
            CALL SDRES3( RAPPORT , MCN(MNWINO) , MCN(MNFINO) )
         END IF
 109  CONTINUE
      GO TO 200
C
C     L'ALGORITHME A CONVERGE !
C
 201  WRITE(IMPRIM,1000) ITER
      DNP1=SQRT(DNP1)
      WRITE(IMPRIM,3000) DNP1
C
C     2.2.8) IMPRESSION DE LA SOLUTION PAR SOUS-DOMAINE
C
 202  MNDF = MNDFSD + WTYOBJ
      DO 120 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C
C           LE TABLEAU TEMPERATURE
C           ----------------------
            NTDL   = MCN(MNLAUX+NO)
            CALL LXTSOU( NTLXOB , 'VECTEUR"TEMPERATURE' , NTU , MNU )
            MCN( MNU + WBCOVE ) = NTDL
            MCN( MNU + WBVECT ) = NDSM
            MCN( MNU + WBCPIN ) = 0
C           LA DATE
            CALL ECDATE( MCN(MNU) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNU + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
            MNVTNO = MNU + WECTEU
C=========================================================================
            WRITE(IMPRIM,7000) KNOM
            NTDLI = MIN(NTDL,10)
            IAU0 = ( MNVTNO - 1 ) / 2
            DO 121 I=1,NTDLI
               WRITE(IMPRIM,8000) (I,DMCN(IAU0+I+(J-1)*NTDL),J=1,NDSM)
 121        CONTINUE
C=========================================================================
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 120  CONTINUE
C
C     3) CALCUL DES FLUX PAR SOUS-DOMAINE
C     ===================================
C
C     RECENSEMENT DES OBJETS AUX LIMITES
      MNDF = MNDFSD + WTYOBJ
      NBLIGF = 0
      CALL TNMCDC ( 'ENTIER' , NBDOBJ , MNLIGN )
      CALL AZEROI ( NBDOBJ , MCN(MNLIGN) )
      CALL TNMCDC ( 'REEL2' , NBDOBJ*NDSM , MNFLUX )
      CALL AZEROD ( NBDOBJ*NDSM , MCN(MNFLUX) )
      DO 112 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
C        SEULES LES LIGNES SONT PRISES EN COMPTE
         IF (NUTYSD .EQ. 2 ) THEN
            MCN(MNLIGN+NBLIGF) = NUOBSD
            NBLIGF = NBLIGF + 1
         END IF
 112  CONTINUE
      MNDF = MNDFSD + WTYOBJ
      DO 110 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            I = INDEX( KNOMOB , '_' )
            KNOM = KNOMOB(1:I-1)
C
C           CALCUL DES FLUX
C           ---------------
            WRITE(IMPRIM,6000) KNOM
            CALL SDTHE4( KNOMOB, NTLXOB, NBLIGF, MNLIGN, MNFLUX, IERR )
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 110  CONTINUE
C
C     AFFICHAGE DES FLUX
C     ==================
C
C     NOM DE L'OBJET
      L = INDEX( KNOMST , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMST )
      ENDIF
      WRITE(IMPRIM,6000) KNOMST(1:L)
      WRITE(IMPRIM,10000)
10000 FORMAT(1X,62('-'))
      IAF = ( MNFLUX - 1 ) / 2
      DO 22 NL = 1 , NBLIGF
         NUTYLI = 2
         NUOBLI = MCN(MNLIGN+NL-1)
         KNOMTY = NMTYOB(NUTYLI)
         CALL NMOBNU( KNOMTY , NUOBLI , KNOMOB )
         L = INDEX( KNOMOB , '_' )
         IF( L .GT. 0 ) THEN
            L = L - 1
         ELSE
            L = LEN( KNOMST )
         ENDIF
         KNOM = KNOMOB(1:L-1)
C        AFFICHAGE DU FLUX A TRAVERS CETTE SURFACE
         WRITE(IMPRIM,10020) KNOMTY,KNOM,
     %   (DMCN(IAF+NL+(K-1)*NBLIGF),K=1,NDSM)
10020    FORMAT(1X,'| ',A,': ',A,' | FLUX =',5(G13.5,' |') )
 22   CONTINUE
      WRITE(IMPRIM,10000)
C
C     DESTRUCTION DES TABLEAUX UTILITAIRES
C     ====================================
C
      CALL TNMCDS( 'REEL2'  , LOPOID , MNPOID )
      CALL TNMCDS( 'REEL2'  , LOTRIN , MNTRIN )
      CALL TNMCDS( 'REEL2'  , LORSIN , MNRSIN )
      CALL TNMCDS( 'ENTIER' , LOADR  , MNADR  )
      CALL TNMCDS( 'ENTIER' , LOINSD , MNINSD )
      CALL TNMCDS ( 'ENTIER' , NBDOBJ , MNLIGN )
      CALL TNMCDS ( 'REEL2' , NBDOBJ*NDSM , MNFLUX )
C
C     DESTRUCTION DES TABLEAUX ASSOCIES DANS LES LEXIQUES
C     ===================================================
C
      MNDF = MNDFSD + WTYOBJ
      DO 111 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C           DESTRUCTION DES TABLEAUX INUTILES
            NTDL =  MCN(MNLAUX+NO)
            LOVX =  NTDL * NDSM * NBVECT
            MNVX =  MCN(MNVAUX+NO)
            CALL TNMCDS( 'REEL2' , LOVX , MNVX )
            MNPOIL = MCN(MNMNPO+NO)
            LOPOIL =  NBNOIN * NBINCO
            CALL TNMCDS( 'REEL2' , LOPOIL , MNPOIL )
            IF (NORESO .EQ. 1) THEN
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE' )
               CALL LXTSDS( NTLXOB , 'PROFIL"CONDUCTIVITE_D' )
               CALL LXTSDS( NTLXOB , 'PROFIL"CONDUCTIVITE_N' )
               CALL LXTSDS( NTLXOB , 'VECTEUR"B')
            ELSE IF (NORESO .EQ. 2 ) THEN
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE' )
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE_D' )
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE_N' )
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE_DC' )
               CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE_NC' )
               CALL LXTSDS( NTLXOB , 'VECTEUR"B')
            END IF
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 111  CONTINUE
C
C     FERMETURE DU LEXIQUE DE L'OBJET 'SOUS-DOMAINES'
C     ===============================================
      CALL LXLXFE( NTOBJE , KNOMSD )
C
      RETURN
C
C     GESTION DES ERREURS
C     ===================
 9900 IF( INTERA .GE. 3 ) RETURN
      CALL ARRET(100)
      RETURN
C
 1000 FORMAT(//,' ** CONVERGENCE EN       ',I6,' ITERATIONS **',/)
 2000 FORMAT(//,' ** NON CONVERGENCE APRES',I6,' ITERATIONS **',/)
 3000 FORMAT(' NORME DU RESIDU',T40,G16.7)
 4000 FORMAT(/1X,80('=')/
     %' PREPARATION DES PROBLEMES DE DIRICHLET ET NEUMANN'
     % ,' SUR L''OBJET: ',A/1X,80('=')/)
 5000 FORMAT(/1X,80('=')/,' RESOLUTION DU PROBLEME DE DIRICHLET'
     % ,' SUR L''OBJET: ',A/1X,80('=')/)
 6000 FORMAT(/1X,80('=')/
     %' CALCUL DES FLUX SUR L''OBJET: ',A,/1X,80('=')/)
 7000 FORMAT(/1X,80('=')/,' LES TEMPERATURES'
     % ,' SUR L''OBJET: ',A/1X,80('=')/)
 8000 FORMAT(' DL',I6,' : TEMPERATURE =',G16.7)
C
      END
