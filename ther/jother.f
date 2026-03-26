      SUBROUTINE JOTHER( KNOMST , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES TEMPERATURES ET LES FLUX DANS UN DOMAINE
C ----- 2D OU 3D OU AXISYMETRIQUE EN THERMIQUE LINEAIRE STATIONNAIRE
C       PAR LA METHODE DES JOINTS
C
C ENTREES :
C ---------
C KNOMST : NOM DE L'OBJET
C
C SORTIES :
C ---------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MARS 1993
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / MSIMTA / NOIMPR
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___joint.inc"
      include"./incl/a___vecteur.inc"
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1),PROSCD,RN,RNP1,QDAQD,ALPHA,BETA,UN
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*(*)     KNOMST
      CHARACTER*24      KNOMOB,KNOMSD,KNOM
      CHARACTER*12      KNOMAT
      CHARACTER*14      KNOMAC
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))
      DATA              KNOMAT/'CONDUCTIVITE'/
      DATA              KNOMAC/'CONDUCTIVITE_C'/
C
      IERR   = 0
      IMPRE  = 1
      IRGP   = 0
      IADIAG = 0
      IADG   = 0
      NBNOJ  = 0
C
C     L'OBJET
C     =======
C
C     NOM DE L'OBJET
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
C     LES JOINTS
C     ==========
C
C     LE TABLEAU JOINT DE L'OBJET KNOMSD
      CALL LXTSOU( NTLXSD , 'JOINT' , NTJOIN , MNJOIN )
C     IL N'EXISTE PAS
      IF( NTJOIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = ' ERREUR : OBJET ' // KNOMSD
         KERR(1) = ' SANS TABLEAU JOINT'
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C     LE NOMBRE DE COMPOSANTS DE L'OBJET
      NBDOBJ = MCN(MNDFSD+WBDOBJ)
C     LE TABLEAU JOINT
      NBJOIN = MCN(MNJOIN+WBJOIN)
      MXNOJO = MCN(MNJOIN+WXNOJO)
      MNOBJE = MNJOIN + WBOBJE
      MNNUJO = MNOBJE + NBJOIN * 4
      MNNOJO = MNNUJO + NBJOIN * 2
      MNCOJO = MNNOJO + NBJOIN * MXNOJO * 2
C=========================================================================
C                                  IMPRESSIONS
C                                  ===========
      IF (NOIMPR.GE.10)  THEN
      DO NO = 0 , NBJOIN - 1
         NO1 = NO + 1
         WRITE (IMPRIM,5000) NO1,MCN(MNOBJE+NO),MCN(MNOBJE+NO+NBJOIN),
     S                  MCN(MNOBJE+NO+NBJOIN*2),MCN(MNOBJE+NO+NBJOIN*3),
     S                           MCN(MNNUJO+NO),MCN(MNNUJO+NO+NBJOIN)
         NUT1 = MCN(MNOBJE+NO)
         NOM1 = MCN(MNOBJE+NO+NBJOIN)
         KNOMTY = NMTYOB(NUT1)
         CALL NMOBNU( KNOMTY , NOM1 , KNOMOB )
         I = INDEX( KNOMOB , '_' )
         KNOM = KNOMOB(1:I-1)
         WRITE (IMPRIM,5300) KNOM
         NUT2 = MCN(MNOBJE+NO+NBJOIN*2)
         NOM2 = MCN(MNOBJE+NO+NBJOIN*3)
         KNOMTY = NMTYOB(NUT2)
         CALL NMOBNU( KNOMTY , NOM2 , KNOMOB )
         I = INDEX( KNOMOB , '_' )
         KNOM = KNOMOB(1:I-1)
         WRITE (IMPRIM,5400) KNOM
         NB1 = MCN(MNNUJO+NO)
         NB2 = MCN(MNNUJO+NO+NBJOIN)
         WRITE (IMPRIM,5100)
         DO  NB = 0 , NB1 - 1
            WRITE (IMPRIM,5500)  MCN(MNNOJO+NO+NB*NBJOIN*2),
     S     (RMCN(MNCOJO+NO+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
         WRITE (IMPRIM,5200)
         DO  NB = 0 , NB2 - 1
            WRITE (IMPRIM,5500) MCN(MNNOJO+NO+NBJOIN+NB*NBJOIN*2),
     S   (RMCN(MNCOJO+NO+NBJOIN+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
      ENDDO
      ENDIF
C=========================================================================
      LOJOSD =  NBDOBJ * 8
      MNJOSD = 0
      CALL TNMCDC( 'ENTIER' , LOJOSD , MNJOSD )
      CALL AZEROI( LOJOSD , MCN(MNJOSD) )
      MNNOSD = 0
      CALL TNMCDC( 'ENTIER' , NBDOBJ , MNNOSD )
      CALL AZEROI( NBDOBJ , MCN(MNNOSD) )
      MNDF = MNDFSD + WTYOBJ
      DO 1 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
C           SELECTION DES JOINTS MAITRES ET ESCLAVES
C           ----------------------------------------
            MCN(MNNOSD+NO) = NUOBSD
            MNJOIM = 0
            CALL TNMCDC( 'ENTIER' , NBJOIN, MNJOIM )
            CALL AZEROI( NBJOIN , MCN(MNJOIM) )
            MNJOIE = 0
            CALL TNMCDC( 'ENTIER' , NBJOIN, MNJOIE )
            CALL AZEROI( NBJOIN , MCN(MNJOIE) )
C
C           1) SELECTION SUIVANT LE NOMBRE DE NOEUDS
C
            CALL JOTHE2( NUOBSD , MCN(MNOBJE) , NBJOIN ,
     &                            MCN(MNNUJO) ,
     &                            MCN(MNJOIM) , NBJOIM ,
     &                            MCN(MNJOIE) , NBJOIE )
C
C           2) SELECTION SUIVANT LE NUMERO DE SOUS-DOMAINE
C
C            CALL JOTHE3( NUOBSD , MCN(MNOBJE) , NBJOIN ,
C     &                            MCN(MNJOIM) , NBJOIM ,
C     &                            MCN(MNJOIE) , NBJOIE )
C           STOCKAGE
            NO8 = NO*8
            MCN(MNJOSD+NO8)   = MNJOIM
            MCN(MNJOSD+NO8+1) = NBJOIM
            MCN(MNJOSD+NO8+2) = MNJOIE
            MCN(MNJOSD+NO8+3) = NBJOIE
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 1    CONTINUE
      MNNUSD = 0
      CALL TNMCDC( 'ENTIER' , NBJOIN*2 , MNNUSD )
      CALL AZEROI( NBJOIN*2 , MCN(MNNUSD) )
      DO 2 NJ = 0 , NBJOIN - 1
         NUSD1 = MCN(MNOBJE+NJ+NBJOIN)
         MNDF = MNDFSD + WTYOBJ
         DO 3 NO = 0 , NBDOBJ - 1
            NUTYSD = MCN(MNDF)
            MNDF   = MNDF + 2
            IF (NUTYSD .EQ. 5 ) THEN
               IF (MCN(MNNOSD+NO).EQ.NUSD1) THEN
                  MCN(MNNUSD+NJ*2) = NO+1
                  GOTO 4
               ENDIF
            ENDIF
 3       CONTINUE
 4       NUSD2 = MCN(MNOBJE+NJ+NBJOIN*3)
         MNDF = MNDFSD + WTYOBJ
         DO 5 NO = 0 , NBDOBJ - 1
            NUTYSD = MCN(MNDF)
            MNDF   = MNDF + 2
            IF (NUTYSD .EQ. 5 ) THEN
               IF (MCN(MNNOSD+NO).EQ.NUSD2) THEN
                  MCN(MNNUSD+NJ*2+1) = NO+1
                  GOTO 2
               ENDIF
            ENDIF
 5       CONTINUE
 2    CONTINUE
C
C     CHOIX DE LA METHODE DE RESOLUTION
C     =================================
C
C     RESOLUTION PAR GRADIENT CONJUGUE GLOBAL ( STOCKAGE MORSE )
      NORESO = 2
C     LE NOMBRE DE D'INCONNUES PAR NOEUDS
C     NBINCO = 1
C     LE PRECONDITIONNEMENT :
C        0 = PAS DE PRECONDITIONNEMENT
C        1 = PRECONDITIONNEMENT   NEUMANN / JOINT DANS CHAQUE SOUS-DOMAINE
C        2 = 1 + PRECONDITIONNEMENT DIAGONAL DU SYSTEME DU JOINT
C        3 = PRECONDITIONNEMENT DIRICHLET / JOINT DANS CHAQUE SOUS-DOMAINE
C        4 = 3 + PRECONDITIONNEMENT DIAGONAL DU SYSTEME DU JOINT
      IPREC = 0
      IF (IPREC.EQ.0) THEN
         PRINT 9001
      ELSE IF (IPREC.EQ.1) THEN
         PRINT 9002
      ELSE IF (IPREC.EQ.2) THEN
         PRINT 9003
      ELSE IF (IPREC.EQ.3) THEN
         PRINT 9004
      ELSE IF (IPREC.EQ.4) THEN
         PRINT 9005
      ENDIF
C
C     1) PREPARATION DU CALCUL
C     ========================
C
      MNDF = MNDFSD + WTYOBJ
C     ADRESSE ET TAILLE DES TABLEAUX UTILITAIRES
      MNADR = 0
      LOADR = NBDOBJ * 3
      CALL TNMCDC( 'ENTIER' , LOADR , MNADR )
      CALL AZEROI( LOADR , MCN(MNADR) )
      MNLAUX = MNADR
      MNINTE = MNLAUX + NBDOBJ
      MNDCAL = MNINTE + NBDOBJ
      MOREE2 = MOTVAR(6)
      NTDLO  = 0
      NTOTAL = 0
C     LES TABLEAUX ASSOCIEES AUX JOINTS
      LOQ = MXNOJO*3
      LOS = MXNOJO
      LOQQ = MXNOJO*2
      MNMAJO = 0
      LOMAJO = NBJOIN * 10
      CALL TNMCDC( 'ENTIER' , LOMAJO , MNMAJO )
      CALL AZEROI( LOMAJO , MCN(MNMAJO) )
      DO 13 NJ = 0 , NBJOIN - 1
         MNQM=0
         CALL TNMCDC( 'REEL2' , LOQ , MNQM )
         CALL AZEROD( LOQ , MCN(MNQM) )
         MNQE=0
         CALL TNMCDC( 'REEL2' , LOQ , MNQE )
         CALL AZEROD( LOQ , MCN(MNQE) )
         MNLUM=0
         CALL TNMCDC( 'REEL2' , LOQ , MNLUM )
         CALL AZEROD( LOQ , MCN(MNLUM) )
         MNLUE=0
         CALL TNMCDC( 'REEL2' , LOQ , MNLUE )
         CALL AZEROD( LOQ , MCN(MNLUE) )
         MNSMM=0
         CALL TNMCDC( 'REEL2' , LOS , MNSMM )
         CALL AZEROD( LOS , MCN(MNSMM) )
         MNSME=0
         CALL TNMCDC( 'REEL2' , LOS , MNSME )
         CALL AZEROD( LOS , MCN(MNSME) )
         MNQQM=0
         CALL TNMCDC( 'REEL2' , LOQQ , MNQQM )
         CALL AZEROD( LOQQ , MCN(MNQQM) )
         MNQQE=0
         CALL TNMCDC( 'REEL2' , LOQQ , MNQQE )
         CALL AZEROD( LOQQ , MCN(MNQQE) )
         MNNUM=0
         CALL TNMCDC( 'ENTIER' , LOQQ , MNNUM )
         CALL AZEROI( LOQQ , MCN(MNNUM) )
         MNNUE=0
         CALL TNMCDC( 'ENTIER' , LOQQ , MNNUE )
         CALL AZEROD( LOQQ , MCN(MNNUE) )
         NJ10 = NJ*10
         MCN(MNMAJO+NJ10)   = MNQM
         MCN(MNMAJO+NJ10+1) = MNLUM
         MCN(MNMAJO+NJ10+2) = MNSMM
         MCN(MNMAJO+NJ10+3) = MNQQM
         MCN(MNMAJO+NJ10+4) = MNNUM
         MCN(MNMAJO+NJ10+5) = MNQE
         MCN(MNMAJO+NJ10+6) = MNLUE
         MCN(MNMAJO+NJ10+7) = MNSME
         MCN(MNMAJO+NJ10+8) = MNQQE
         MCN(MNMAJO+NJ10+9) = MNNUE
 13   CONTINUE
C
      DO 6 NO = 0 , NBDOBJ - 1
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
C           PREPARATION : ASSEMBLAGE DE LA MATRICE LOCALE
C           ---------------------------------------------
            IF (NOIMPR.GE.10) WRITE(IMPRIM,4000) KNOM
            CALL JOTHE1( NTLXOB , NUTYSD , NUOBSD , NBJOIN ,
     &                   MXNOJO ,  MCN(MNOBJE) , MCN(MNNUJO) ,
     &                             MCN(MNNOJO) , IPREC , IERR )
C           STOCKAGE DU NOMBRE LOCAL D'INCONNUES
            MCN(MNLAUX+NO) = NTDL
C
C           LE NOMBRE D'INCONNUES DU SYSTEME GLOBAL
C           ---------------------------------------
            NO8 = NO*8
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            NTDLE = 0
            DO 7 NBJ = 0 , NBJOIE - 1
               NJ = MCN(MNJOIE+NBJ) - 1
C              LE NOMBRE DE NOEUDS DE CE JOINT SUR LE SD
               NOBJE1 = MCN(MNOBJE+NJ+NBJOIN)
               NOBJE2 = MCN(MNOBJE+NJ+NBJOIN*3)
               IF (NUOBSD.EQ.NOBJE1) THEN
                 NBNOJ = MCN(MNNUJO+NJ)
               ELSE IF (NUOBSD.EQ.NOBJE2) THEN
                 NBNOJ = MCN(MNNUJO+NJ+NBJOIN)
               ENDIF
               NTDLE = NTDLE + NBNOJ - 2
 7          CONTINUE
            NTDLO = NTDLO + NTDL - NTDLE
            MCN(MNDCAL+NO) = NTOTAL
            NTOTAL = NTOTAL + NTDL
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 6    CONTINUE
      MNDF = MNDFSD + WTYOBJ
      DO 8 NO = 0 , NBDOBJ - 1
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
C           CONSTRUCTION ET FACTORISATION DE LA MATRICE Q
C           ---------------------------------------------
            IF (NOIMPR.GE.10) WRITE(IMPRIM,4200) KNOM
            NO8 = NO*8
            MNJOIM = MCN(MNJOSD+NO8)
            NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            CALL JOTHE4( NUOBSD ,  NBJOIN , MXNOJO ,
     S        MCN(MNOBJE) , MCN(MNNOJO) , MCN(MNNUJO) , MCN(MNCOJO) ,
     S        MCN(MNJOIM) , NBJOIM , MCN(MNJOIE) , NBJOIE , MNMAJO )
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 8    CONTINUE
      CALL TNMCDS( 'ENTIER' , NBDOBJ , MNNOSD )
C
C     2) CALCUL DES TEMPERATURES
C     ==========================
C
C     2.1) INITIALISATION
C     -------------------
C
      LOXG =  NTOTAL
      MNXG = 0
      CALL TNMCDC( 'REEL2' , LOXG , MNXG )
      CALL AZEROD( LOXG , MCN(MNXG) )
      IXG = ( MNXG  - 1 ) / 2
      MNRG = 0
      CALL TNMCDC( 'REEL2' , LOXG , MNRG )
      CALL AZEROD( LOXG , MCN(MNRG) )
      IRG = ( MNRG  - 1 ) / 2
C
C     2.1.1) LE PRECONDITIONNEMENT GLOBAL (DIAGONAL)
C
      IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
         MNADG = 0
         CALL TNMCDC( 'REEL2' , LOXG , MNADG )
         CALL AZEROD( LOXG , MCN(MNADG) )
         IADG = ( MNADG  - 1 ) / 2
         MNDIAG = 0
         CALL TNMCDC( 'REEL2' , LOXG , MNDIAG )
         CALL AZEROD( LOXG , MCN(MNDIAG) )
         IADIAG = ( MNDIAG - 1 ) / 2 + 1
         MNLIDI = 0
         CALL TNMCDC( 'ENTIER' , NTOTAL , MNLIDI )
         CALL AZEROI( NTOTAL , MCN(MNLIDI) )
         NPJOIN = 0
         MNRGP = 0
         CALL TNMCDC( 'REEL2' , LOXG , MNRGP )
         CALL AZEROD( LOXG , MCN(MNRGP) )
         IRGP = ( MNRGP  - 1 ) / 2
      ENDIF
C
C     2.1.2) LA SOLUTION ET LE RESIDU
C
      MNDF = MNDFSD + WTYOBJ
      NBITEM = 0
      DO 9 NO = 0 , NBDOBJ - 1
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
            IF (NOIMPR.GE.10) WRITE(IMPRIM,4100) KNOM
            NBITEM = NBITEM + 1
            NTDL   = MCN(MNLAUX+NO)
C           LE VECTEUR TEMPERATURE DE L'OBJET (INITIALISEE A ZERO DANS JOTHE1)
            CALL LXTSOU( NTLXOB , 'VECTEUR"TEMPERATURE' , NTU , MNU )
            MNVTNO = MNU + WECTEU
C           LE VECTEUR SECOND MEMBRE (CALCULE DANS JOTHE1)
            CALL LXTSOU( NTLXOB , 'VECTEUR"B' , NTB , MNB )
            MNVBNO = MNB + WECTEU
C
C           TRANSFERT DE LA SOLUTION LOCALE
C           -------------------------------
            IXL = ( MNVTNO - 1 ) / 2
            NDCAL = MCN(MNDCAL+NO)
            DO K = 1 , NTDL
               DMCN(IXG+NDCAL+K) = DMCN(IXL+K)
            ENDDO
C           PAS DE CORRECTION DE LA SOLUTION SUR LE JOINT (=0.)
C
C           LE RESIDU LOCAL
C           ---------------
            MNRLO  = 0
            CALL TNMCDC( 'REEL2' , NTDL , MNRLO )
            CALL AZEROD( NTDL , MCN(MNRLO) )
            NO8 = NO*8
            MCN(MNJOSD+NO8+4) = MNRLO
C           ADRESSAGE ET STOCKAGE DES DIRECTIONS (AVANT ITERATIONS)
            MNDLO  = 0
            CALL TNMCDC( 'REEL2' , NTDL , MNDLO )
            CALL AZEROD( NTDL , MCN(MNDLO) )
            MCN(MNJOSD+NO8+5) = MNDLO
            MNADLO  = 0
            CALL TNMCDC( 'REEL2' , NTDL , MNADLO )
            CALL AZEROD( NTDL , MCN(MNADLO) )
            MCN(MNJOSD+NO8+6) = MNADLO
            MNZLO  = 0
            CALL TNMCDC( 'REEL2' , NTDL , MNZLO )
            CALL AZEROD( NTDL , MCN(MNZLO) )
            MCN(MNJOSD+NO8+7) = MNZLO
C           CALCUL DU RESIDU LOCAL
            CALL SDRES1( NTLXOB , MNVTNO, MNVBNO , MNRLO ,
     S                   KNOMAT , IERR )
C
C           PRECONDITIONNEMENT DU RESIDU LOCAL
C           ----------------------------------
C           ANCIENNE FORME DE PRECONDITIONNEMENT (AVEC IPREC=1)
C           CALL SDRES14( NTLXOB , MNRLO , MNZLO , KNOMAC )
            IF (IPREC.EQ.1 .OR. IPREC.EQ.3) THEN
               CALL SDRES15( NTLXOB , MNRLO , MNZLO , KNOMAC )
               CALL TRTATA( MCN(MNZLO) , MCN(MNRLO) , NTDL*MOREE2 )
            ELSE IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
               CALL SDRES15( NTLXOB , MNRLO , MNZLO , KNOMAC )
               CALL TRTATA( MCN(MNZLO) , MCN(MNRLO) , NTDL*MOREE2 )
C              LA DIAGONALE DU PRECONDITIONNEMENT GLOBAL
               CALL SDRES20( NTLXOB , MNADLO , KNOMAT )
            ENDIF
C
C           LE RESIDU GLOBAL
C           ----------------
            IRL = ( MNRLO - 1 ) / 2
C           TRANSFERT DU RESIDU LOCAL
            NDCAL = MCN(MNDCAL+NO)
            DO K = 1 , NTDL
               DMCN(IRG+NDCAL+K) = DMCN(IRL+K)
            ENDDO
C           LA DIAGONALE DU PRECONDITIONNEMENT GLOBAL
            IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
               IADL = ( MNADLO - 1 ) / 2
               DO K = 1 , NTDL
                  DMCN(IADG+NDCAL+K) =  DMCN(IADL+K)
               ENDDO
            ENDIF
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 9    CONTINUE
C
C     2.1.2) CORRECTION DU RESIDU SUR LE JOINT
C
      MNDF = MNDFSD + WTYOBJ
      DO 10 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C           CORRECTION DU RESIDU AUX NOEUDS ESCLAVES
            NO8 = NO*8
            MNJOIM = MCN(MNJOSD+NO8)
            NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            MNRLO  = MCN(MNJOSD+NO8+4)
            CALL JOTHE6( NUOBSD , NBJOIN , MXNOJO ,
     S                   MCN(MNOBJE) , MCN(MNNUJO) ,
     S                   MCN(MNNOJO) , MCN(MNNUSD) , MCN(MNDCAL) ,
     S                   MCN(MNRG)   ,
     S                   MCN(MNJOIM) , NBJOIM , MNMAJO )
C           LA DIAGONALE DU PRECONDITIONNEMENT GLOBAL
            IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
            CALL JOTHE7( NUOBSD , NBJOIN , MXNOJO ,
     S                   MCN(MNOBJE) , MCN(MNNUJO) ,
     S                   MCN(MNNOJO) , MCN(MNNUSD) , MCN(MNDCAL) ,
     S                   MCN(MNADG)  ,
     S                   MCN(MNJOIM) , NBJOIM ,
     S                   MCN(MNDIAG) , MCN(MNLIDI) , NPJOIN )
            ENDIF
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 10   CONTINUE
C
C     LA NORME DU RESIDU INITIAL
      RN = PROSCD( MCN(MNRG) , MCN(MNRG) , NTOTAL )
C
C     2.1.3) PRECONDITIONNEMENT DU RESIDU GLOBAL
C
      IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
         CALL TRTATA( MCN(MNRG) , MCN(MNRGP) , NTOTAL*MOREE2 )
         DO K = 0 , NPJOIN - 1
            L = MCN(MNLIDI+K)
            DMCN(IRGP+L) = DMCN(IRGP+L) / DMCN(IADIAG+K)
         ENDDO
         CALL TNMCDS( 'REEL2' , LOXG , MNADG )
         RN = PROSCD( MCN(MNRG) , MCN(MNRGP) , NTOTAL )
      ENDIF
C
C     2.1.4) LA DIRECTION DE DESCENTE
C
      MNDG = 0
      CALL TNMCDC( 'REEL2' , LOXG , MNDG )
      CALL TRTATA( MCN(MNRG) , MCN(MNDG) , LOXG*MOREE2 )
      IDG = ( MNDG - 1 ) / 2
      MNADG = 0
      CALL TNMCDC( 'REEL2' , LOXG , MNADG )
      IADG = ( MNADG - 1 ) / 2
C=========================================================================
C                                  IMPRESSIONS
C                                  ===========
      IF (NOIMPR.GE.5)  THEN
      DO NO = 0 , NBJOIN - 1
         NJ = NO + 1
         WRITE (IMPRIM,5000) NJ,MCN(MNOBJE+NO),MCN(MNOBJE+NO+NBJOIN),
     S                  MCN(MNOBJE+NO+NBJOIN*2),MCN(MNOBJE+NO+NBJOIN*3),
     S                           MCN(MNNUJO+NO),MCN(MNNUJO+NO+NBJOIN)
         NB1   = MCN(MNNUJO+NO)
         NB2   = MCN(MNNUJO+NO+NBJOIN)
         NOB1  = MCN(MNNUSD+NO*2)
         NOB2  = MCN(MNNUSD+NO*2+1)
         NCAL1 = MCN(MNDCAL+NOB1-1)
         NCAL2 = MCN(MNDCAL+NOB2-1)
         WRITE (IMPRIM,5100)
         DO  NB = 0 , NB1 - 1
            NUM1  = MCN(MNNOJO+NO+NB*NBJOIN*2)
            NUMG1 = NUM1 + NCAL1
            WRITE (IMPRIM,5600)  NUM1,DMCN(IXG+NUMG1),DMCN(IRG+NUMG1),
     S          (RMCN(MNCOJO+NO+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
         WRITE (IMPRIM,5200)
         DO  NB = 0 , NB2 - 1
            NUM2  = MCN(MNNOJO+NO+NBJOIN+NB*NBJOIN*2)
            NUMG2 = NUM2 + NCAL2
            WRITE (IMPRIM,5600)  NUM2,DMCN(IXG+NUMG2),DMCN(IRG+NUMG2),
     S   (RMCN(MNCOJO+NO+NBJOIN+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
      ENDDO
      ENDIF
C=========================================================================
C
C     2.2) ITERATIONS
C     ---------------
C
C     LE NOMBRE MAXIMUM D'ITERATIONS
      ITEMAX = NTDLO
C     LA PRECISION DEMANDEE
      EPS  = EPZERO
      EPS0 = REAL( EPS * EPS * RN )
      IF (RN.LT.1.D0) EPS0  = EPS * EPS
      ITER = 0
      KITER = 25
      UN = 1.D0
 200  ITER = ITER + 1
C
C     2.2.1) CALCUL DU PRODUIT  A * Q * D
C
      CALL AZEROD( LOXG , MCN(MNADG) )
      MNDF   = MNDFSD + WTYOBJ
      DO 11 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C           LES ADRESSES DES TABLEAUX TEMPORAIRES
            NTDL   = MCN(MNLAUX+NO)
            NO8 = NO*8
            NDCAL = MCN(MNDCAL+NO)
C
C           LA DIRECTION LOCALE DK
C           ----------------------
            MNDLO  = MCN(MNJOSD+NO8+5)
            CALL AZEROD( NTDL , MCN(MNDLO) )
            IDL = ( MNDLO - 1 ) / 2
C           TRANSFERT DE LA DIRECTION GLOBALE
            DO K = 1 , NTDL
               DMCN(IDL+K) = DMCN(IDG+NDCAL+K)
            ENDDO
C           CORRECTION DE LA DIRECTION AUX NOEUDS ESCLAVES
            MNJOIM = MCN(MNJOSD+NO8)
            NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            CALL JOTHE5( NUOBSD , NBJOIN , MXNOJO ,
     S                   MCN(MNOBJE) , MCN(MNNUJO) ,
     S                   MCN(MNNOJO) , MCN(MNNUSD) , MCN(MNDCAL) ,
     S                   MCN(MNDG)   , MCN(MNDLO)  ,
     S                   MCN(MNJOIE) , NBJOIE , MNMAJO )
C
C           PRECONDITIONNEMENT DE LA DIRECTION LOCALE
C           -----------------------------------------
C
C           ANCIENNE FORME DE PRECONDITIONNEMENT (AVEC IPREC=1)
C           IF (IPREC.GE.3) THEN
C
            IF (IPREC.GE.1) THEN
               MNZLO  = MCN(MNJOSD+NO8+7)
               CALL AZEROD( NTDL , MCN(MNZLO) )
               CALL SDRES16( NTLXOB , MNDLO , MNZLO , KNOMAC )
               CALL TRTATA( MCN(MNZLO) , MCN(MNDLO) , NTDL*MOREE2 )
            END IF
C
C           LE PRODUIT A * Q * DK LOCAL
C           ---------------------------
            MNADLO = MCN(MNJOSD+NO8+6)
            CALL AZEROD( NTDL , MCN(MNADLO) )
            MNZLO  = MCN(MNJOSD+NO8+7)
            CALL AZEROD( NTDL , MCN(MNZLO) )
            CALL SDRES1( NTLXOB , MNDLO , MNZLO , MNADLO ,
     &                   KNOMAT , IERR )
C           CE SP FOURNIT  AD = Z - A * D , SOIT - AD !
C           ON PEUT DONC SE PASSER DE CE CHANGEMENT DE SIGNE (VOIR 1 ET 2)
C           CALL CHSGND( NTDL , MCN(MNADLO) )
C
C           PRECONDITIONNEMENT DU PRODUIT LOCAL
C           -----------------------------------
C
C           ANCIENNE FORME DE PRECONDITIONNEMENT (AVEC IPREC=1)
C            IF (IPREC.EQ.1 .OR. IPREC.EQ.2) THEN
C               CALL AZEROD( NTDL , MCN(MNZLO) )
C               CALL SDRES14( NTLXOB , MNADLO , MNZLO , KNOMAC )
C               CALL TRTATA( MCN(MNZLO) , MCN(MNADLO) , NTDL*MOREE2 )
C            ELSE IF (IPREC.EQ.3 .OR. IPREC.EQ.4) THEN
C
            IF (IPREC.GE.1) THEN
               CALL AZEROD( NTDL , MCN(MNZLO) )
               CALL SDRES15( NTLXOB , MNADLO , MNZLO , KNOMAC )
               CALL TRTATA( MCN(MNZLO) , MCN(MNADLO) , NTDL*MOREE2 )
            END IF
C
C           LE PRODUIT  A * Q * DK GLOBAL
C           -----------------------------
            IADL = ( MNADLO - 1 ) / 2
C           TRANSFERT AU SYSTEME GLOBAL
            DO K = 1 , NTDL
               DMCN(IADG+NDCAL+K) = DMCN(IADL+K)
            ENDDO
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 11   CONTINUE
C
C     2.2.2) CORRECTION DU PRODUIT SUR LE JOINT
C
      MNDF   = MNDFSD + WTYOBJ
      DO 12 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
C           CORRECTION DU PRODUIT SUR LE JOINT
C           ----------------------------------
            NO8 = NO*8
            MNJOIM = MCN(MNJOSD+NO8)
            NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            MNADLO = MCN(MNJOSD+NO8+6)
            CALL JOTHE6( NUOBSD , NBJOIN , MXNOJO ,
     S                   MCN(MNOBJE) , MCN(MNNUJO) ,
     S                   MCN(MNNOJO) , MCN(MNNUSD) , MCN(MNDCAL) ,
     S                   MCN(MNADG)  ,
     S                   MCN(MNJOIM) , NBJOIM , MNMAJO )
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 12   CONTINUE
C
C     2.2.3) CALCUL DE ALPHA
C
C     ON PEUT SE PASSER DE CE CHANGEMENT (1) MCN(MNADG) CONTIENT - AD
C     QDAQD =  PROSCD( MCN(MNDG) , MCN(MNADG) , NTOTAL )
      QDAQD =  - PROSCD( MCN(MNDG) , MCN(MNADG) , NTOTAL )
      IF (ABS(QDAQD) .GT. EPZERO*RN) THEN
         ALPHA = RN / QDAQD
      ELSE
         GO TO 201
      ENDIF
C
C     2.2.4) MISE A JOUR DES VECTEURS
C
C     LES TEMPERATURES GLOBALES
      CALL SDRES9(MCN(MNXG),UN,MCN(MNXG),ALPHA,MCN(MNDG),NTOTAL)
C     LES RESIDUS GLOBAUX
C     ON PEUT SE PASSER DE CE CHANGEMENT (2) MCN(MNADG) CONTIENT - AD
C     CALL SDRES9(MCN(MNRG),UN,MCN(MNRG),-ALPHA,MCN(MNADG),NTOTAL)
      CALL SDRES9(MCN(MNRG),UN,MCN(MNRG),ALPHA,MCN(MNADG),NTOTAL)
C     LA NORME DU RESIDU GLOBAL
      RNP1 = PROSCD( MCN(MNRG) , MCN(MNRG) , NTOTAL )
C
C     2.2.5) PRECONDITIONNEMENT DU RESIDU GLOBAL (DIAGONAL)
C
      IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
         CALL TRTATA( MCN(MNRG) , MCN(MNRGP) , NTOTAL*MOREE2 )
         DO K = 0 , NPJOIN - 1
            L = MCN(MNLIDI+K)
            DMCN(IRGP+L) = DMCN(IRGP+L) / DMCN(IADIAG+K)
         ENDDO
         RNP1 = PROSCD( MCN(MNRG) , MCN(MNRGP) , NTOTAL )
      ENDIF
C
C=========================================================================
      IF (NOIMPR.GE.10)  THEN
      NTDLI = NTOTAL
         DO I=1,NTDLI
            WRITE(IMPRIM,12811)
     &      (I,DMCN(IXG+I+(J-1)*NTDL),DMCN(IRG+I+(J-1)*NTDL),
     &         DMCN(IDG+I+(J-1)*NTDL),DMCN(IADG+I+(J-1)*NTDL),J=1,NDSM)
         ENDDO
12811    FORMAT( ' DL',I6,' : SOL =',D15.6,' RES =',D15.6,
     &           ' DG =',D15.6,' ADG =',D15.6 )
      ENDIF
C=========================================================================
C
C     2.2.6) TEST DE CONVERGENCE
C
      IF (MOD(ITER,KITER).EQ.0) WRITE(IMPRIM,9000)  ITER,RNP1
      IF (RNP1 .LT. EPS0 ) THEN
         GO TO 201
      ELSE
         IF (ITER.EQ.ITEMAX) THEN
            WRITE(IMPRIM,2000) ITER
            WRITE(IMPRIM,3000) RNP1
            GO TO 202
         END IF
      END IF
C
C     2.2.7) CALCUL DE LA NOUVELLE DIRECTION
C
      BETA = RNP1 / RN
      RN   = RNP1
      IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
         CALL SDRES9(MCN(MNDG),UN,MCN(MNRGP),BETA,MCN(MNDG),NTOTAL)
      ELSE
         CALL SDRES9(MCN(MNDG),UN,MCN(MNRG),BETA,MCN(MNDG),NTOTAL)
      ENDIF
C=========================================================================
C                                  IMPRESSIONS
C                                  ===========
      IF (NOIMPR.GE.5)  THEN
      DO NO = 0 , NBJOIN - 1
         NJ = NO + 1
         WRITE (IMPRIM,5000) NJ,MCN(MNOBJE+NO),MCN(MNOBJE+NO+NBJOIN),
     S                  MCN(MNOBJE+NO+NBJOIN*2),MCN(MNOBJE+NO+NBJOIN*3),
     S                           MCN(MNNUJO+NO),MCN(MNNUJO+NO+NBJOIN)
         NB1   = MCN(MNNUJO+NO)
         NB2   = MCN(MNNUJO+NO+NBJOIN)
         NOB1  = MCN(MNNUSD+NO*2)
         NOB2  = MCN(MNNUSD+NO*2+1)
         NCAL1 = MCN(MNDCAL+NOB1-1)
         NCAL2 = MCN(MNDCAL+NOB2-1)
         WRITE (IMPRIM,5100)
         DO  NB = 0 , NB1 - 1
            NUM1  = MCN(MNNOJO+NO+NB*NBJOIN*2)
            NUMG1 = NUM1 + NCAL1
            WRITE (IMPRIM,5600)  NUM1,DMCN(IXG+NUMG1),DMCN(IRG+NUMG1),
     S          (RMCN(MNCOJO+NO+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
         WRITE (IMPRIM,5200)
         DO  NB = 0 , NB2 - 1
            NUM2  = MCN(MNNOJO+NO+NBJOIN+NB*NBJOIN*2)
            NUMG2 = NUM2 + NCAL2
            WRITE (IMPRIM,5600)  NUM2,DMCN(IXG+NUMG2),DMCN(IRG+NUMG2),
     S   (RMCN(MNCOJO+NO+NBJOIN+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
      ENDDO
      ENDIF
C=========================================================================
      GO TO 200
C
C     L'ALGORITHME A CONVERGE !
C
 201  WRITE(IMPRIM,1000) ITER
      WRITE(IMPRIM,3000) RNP1
C
C     2.2.8) MISE A JOUR DE LA SOLUTION PAR SOUS-DOMAINE
C
 202  MNDF = MNDFSD + WTYOBJ
      DO 14 NO = 0 , NBDOBJ - 1
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
            IXL = ( MNVTNO - 1 ) / 2
C           TRANSFERT DE LA SOLUTION GLOBALE
            NDCAL = MCN(MNDCAL+NO)
            DO K = 1 , NTDL
               DMCN(IXL+K) = DMCN(IXG+NDCAL+K)
            ENDDO
C           CORRECTION DE LA SOLUTION AUX NOEUDS ESCLAVES
            NO8 = NO*8
            MNJOIM = MCN(MNJOSD+NO8)
            NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
            NBJOIE = MCN(MNJOSD+NO8+3)
            CALL JOTHE5( NUOBSD , NBJOIN , MXNOJO ,
     S                   MCN(MNOBJE) , MCN(MNNUJO) ,
     S                   MCN(MNNOJO) , MCN(MNNUSD) , MCN(MNDCAL) ,
     S                   MCN(MNXG)   , MCN(MNVTNO) ,
     S                   MCN(MNJOIE) , NBJOIE , MNMAJO )

C
C           PRECONDITIONNEMENT DE LA SOLUTION LOCALE
C           ----------------------------------------
C            IF (IPREC.GE.3) THEN
            IF (IPREC.GE.1) THEN
               MNZLO  = MCN(MNJOSD+NO8+7)
               CALL AZEROD( NTDL , MCN(MNZLO) )
               CALL SDRES16( NTLXOB , MNVTNO , MNZLO , KNOMAC )
               CALL TRTATA( MCN(MNZLO) , MCN(MNVTNO) , NTDL*MOREE2 )
            END IF
C
C              IMPRESSION DE LA SOLUTION
               WRITE(IMPRIM,7000) KNOM
               NTDLI = MIN(NTDL,10)
               IXL = ( MNVTNO - 1 ) / 2
               DO I=1,NTDLI
                  WRITE(IMPRIM,8000) (I,DMCN(IXL+I+(J-1)*NTDL),J=1,NDSM)
               ENDDO
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 14   CONTINUE
C=========================================================================
C                                  IMPRESSIONS
C                                  ===========
      IF (IMPRE.GE.10) THEN
      DO NO = 0 , NBJOIN - 1
         NJ = NO + 1
         WRITE (IMPRIM,5000) NJ,MCN(MNOBJE+NO),MCN(MNOBJE+NO+NBJOIN),
     S                  MCN(MNOBJE+NO+NBJOIN*2),MCN(MNOBJE+NO+NBJOIN*3),
     S                           MCN(MNNUJO+NO),MCN(MNNUJO+NO+NBJOIN)
         NB1   = MCN(MNNUJO+NO)
         NB2   = MCN(MNNUJO+NO+NBJOIN)
         NOB1  = MCN(MNNUSD+NO*2)
         NOB2  = MCN(MNNUSD+NO*2+1)
         NCAL1 = MCN(MNDCAL+NOB1-1)
         NCAL2 = MCN(MNDCAL+NOB2-1)
         WRITE (IMPRIM,5100)
         DO  NB = 0 , NB1 - 1
            NUM1  = MCN(MNNOJO+NO+NB*NBJOIN*2)
            NUMG1 = NUM1 + NCAL1
            WRITE (IMPRIM,5600)  NUM1,DMCN(IXG+NUMG1),DMCN(IRG+NUMG1),
     S          (RMCN(MNCOJO+NO+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
         WRITE (IMPRIM,5200)
         DO  NB = 0 , NB2 - 1
            NUM2  = MCN(MNNOJO+NO+NBJOIN+NB*NBJOIN*2)
            NUMG2 = NUM2 + NCAL2
            WRITE (IMPRIM,5600)  NUM2,DMCN(IXG+NUMG2),DMCN(IRG+NUMG2),
     S   (RMCN(MNCOJO+NO+NBJOIN+NB*NBJOIN*2+K*NBJOIN*2*MXNOJO),K=0,2)
         ENDDO
      ENDDO
      ENDIF
C=========================================================================
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
      DO 21 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
C        SEULES LES LIGNES SONT PRISES EN COMPTE
         IF (NUTYSD .EQ. 2 ) THEN
            MCN(MNLIGN+NBLIGF) = NUOBSD
            NBLIGF = NBLIGF + 1
         END IF
 21   CONTINUE
      MNDF = MNDFSD + WTYOBJ
      DO 17 NO = 0 , NBDOBJ - 1
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
 17   CONTINUE
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
      MNDF = MNDFSD + WTYOBJ
      DO 19 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
            NO8 = NO*8
            MNJOIM = MCN(MNJOSD+NO8)
C           NBJOIM = MCN(MNJOSD+NO8+1)
            MNJOIE = MCN(MNJOSD+NO8+2)
C           NBJOIE = MCN(MNJOSD+NO8+3)
            MNRLO  = MCN(MNJOSD+NO8+4)
            MNDLO  = MCN(MNJOSD+NO8+5)
            MNADLO = MCN(MNJOSD+NO8+6)
            MNZLO  = MCN(MNJOSD+NO8+7)
            NTDL = MCN(MNLAUX+NO)
            CALL TNMCDS( 'ENTIER' , NBJOIN , MNJOIM )
            CALL TNMCDS( 'ENTIER' , NBJOIN , MNJOIE )
            CALL TNMCDS( 'REEL2'  , NTDL   , MNRLO )
            CALL TNMCDS( 'REEL2'  , NTDL   , MNDLO )
            CALL TNMCDS( 'REEL2'  , NTDL   , MNADLO )
            CALL TNMCDS( 'REEL2'  , NTDL   , MNZLO )
         END IF
 19   CONTINUE
C
      DO 20 NJ = 0 , NBJOIN - 1
         NJ10  = NJ*10
         MNQM  = MCN(MNMAJO+NJ10)
         MNLUM = MCN(MNMAJO+NJ10+1)
         MNSMM = MCN(MNMAJO+NJ10+2)
         MNQQM = MCN(MNMAJO+NJ10+3)
         MNNUM = MCN(MNMAJO+NJ10+4)
         MNQE  = MCN(MNMAJO+NJ10+5)
         MNLUE = MCN(MNMAJO+NJ10+6)
         MNSME = MCN(MNMAJO+NJ10+7)
         MNQQE = MCN(MNMAJO+NJ10+8)
         MNNUE = MCN(MNMAJO+NJ10+9)
         CALL TNMCDS( 'REEL2'  , LOQ  , MNQM )
         CALL TNMCDS( 'REEL2'  , LOQ  , MNLUM )
         CALL TNMCDS( 'REEL2'  , LOS  , MNSMM )
         CALL TNMCDS( 'REEL2'  , LOQQ , MNQQM )
         CALL TNMCDS( 'ENTIER' , LOQQ , MNNUM )
         CALL TNMCDS( 'REEL2'  , LOQ  , MNQE )
         CALL TNMCDS( 'REEL2'  , LOQ  , MNLUE )
         CALL TNMCDS( 'REEL2'  , LOS  , MNSME )
         CALL TNMCDS( 'REEL2'  , LOQQ , MNQQE )
         CALL TNMCDS( 'ENTIER' , LOQQ , MNNUE )
 20   CONTINUE
C
      CALL TNMCDS( 'ENTIER' , LOADR  , MNADR  )
      CALL TNMCDS( 'ENTIER' , LOJOSD , MNJOSD )
      CALL TNMCDS( 'REEL2'  , LOXG   , MNXG )
      CALL TNMCDS( 'REEL2'  , LOXG   , MNRG )
      CALL TNMCDS( 'REEL2'  , LOXG   , MNDG )
      CALL TNMCDS( 'REEL2'  , LOXG   , MNADG )
      IF (IPREC.EQ.2 .OR. IPREC.EQ.4) THEN
         CALL TNMCDS( 'REEL2'  , LOXG   , MNDIAG )
         CALL TNMCDS( 'ENTIER' , NTOTAL , MNLIDI )
      ENDIF
      CALL TNMCDS( 'ENTIER' , NBJOIN*2 , MNNUSD )
      CALL TNMCDS( 'ENTIER' , LOMAJO , MNMAJO )
      CALL TNMCDS ( 'ENTIER' , NBDOBJ , MNLIGN )
      CALL TNMCDS ( 'REEL2' , NBDOBJ*NDSM , MNFLUX )
C
C     DESTRUCTION DES TABLEAUX ASSOCIES DANS LES LEXIQUES
C     ===================================================
C
      MNDF = MNDFSD + WTYOBJ
      DO 18 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C           DESTRUCTION DES TABLEAUX INUTILES
            CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE' )
            IF (IPREC.GT.0)
     S      CALL LXTSDS( NTLXOB , 'MORSE"CONDUCTIVITE_C' )
            CALL LXTSDS( NTLXOB , 'VECTEUR"B')
C
C           FERMETURE DU LEXIQUE DE L'OBJET
            CALL LXLXFE( NTOBJE , KNOMOB )
         END IF
 18   CONTINUE
C
C     FERMETURE DU LEXIQUE DE L'OBJET 'SOUS-DOMAINES'
C     ===============================================
C
      CALL LXLXFE( NTOBJE , KNOMSD )
C
      RETURN
C
 1000 FORMAT(//,' ** CONVERGENCE EN       ',I6,' ITERATIONS **',/)
 2000 FORMAT(//,' ** NON CONVERGENCE APRES',I6,' ITERATIONS **',/)
 3000 FORMAT(' NORME DU RESIDU',T40,D15.6)
 4000 FORMAT(/1X,80('=')/
     %' PREPARATION SUR L''OBJET: ',A/1X,80('=')/)
 4100 FORMAT(/1X,80('=')/
     %' INITIALISATION SUR L''OBJET: ',A/1X,80('=')/)
 4200    FORMAT(//,1X,80(1H=)/,
     % ' CONSTRUCTION DES MATRICES JOINTS SUR L''OBJET: ',A/1X,80('=')/)
ccc 4300 FORMAT(/1X,80('=')/
ccc     %' CALCUL SUR L''OBJET: ',A/1X,80('=')/)
 5000 FORMAT(//,1X,80(1H-)/,' JOINT',I3,
     S ' OBJETS T/N ',4I5,' NOEUDS',2I5/,1X,80(1H-))
 5100 FORMAT(/,' PREMIER COTE :',/)
 5200 FORMAT(/,' SECOND COTE  :',/)
 5300 FORMAT(' PREMIER OBJET : ',A24)
 5400 FORMAT(' SECOND  OBJET : ',A24)
 5500 FORMAT(' NUMERO',I5,' COOR',3D15.6)
 5600 FORMAT(' NUMERO',I5,' SOLUTION',D15.6,' RESIDU',D15.6,' COOR',
     S       3D16.6)
 6000 FORMAT(/1X,80('=')/
     %' CALCUL DES FLUX SUR L''OBJET: ',A,/1X,80('=')/)
 7000 FORMAT(/1X,80('=')/,' LES TEMPERATURES'
     % ,' SUR L''OBJET: ',A/1X,80('=')/)
 8000 FORMAT(' DL',I6,' : TEMPERATURE =',D15.6)
 9000 FORMAT( ' + ITERATION',I4,' RESIDU ',D15.6)
 9001 FORMAT(' PAS DE PRECONDITIONNEMENT')
 9002 FORMAT(' PRECONDITIONNEMENT NEUMANN DANS CHAQUE SOUS-DOMAINE' )
 9003 FORMAT(' PRECONDITIONNEMENT NEUMANN DANS CHAQUE SOUS-DOMAINE +'
     &    ,/,' PRECONDITIONNEMENT DIAGONAL DU SYSTEME DU JOINT' )
 9004 FORMAT(' PRECONDITIONNEMENT DIRICHLET DANS CHAQUE SOUS-DOMAINE' )
 9005 FORMAT(' PRECONDITIONNEMENT DIRICHLET DANS CHAQUE SOUS-DOMAINE +'
     &   ,/,' PRECONDITIONNEMENT DIAGONAL DU SYSTEME DU JOINT' )
C
      END
