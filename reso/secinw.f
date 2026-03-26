      SUBROUTINE SECINW( NBROOT, NC,   N,  TETA, 
     +                   MUA,    A,    NCODSA, 
     +                   MUB,    IAAB, B,  NCODSB, 
     +                   V0, V1, V, W, 
     +                   VV, WW, EIGV, 
     +                   NITE,   NINV, 
     +                   MORAID, RAID, NFMUB1, 
     +                   NOV,    IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES NBROOT PLUS PETITES VALEURS ET VECTEURS PROPRES DU
C -----   PROBLEME  (A-L*B) * X =0  (L VALEUR X VECTEUR PROPRE)
C         A MATRICE DE RAIDEUR  B MATRICE DE MASSE
C         METHODE DU DETERMINANT ET ITERATIONS INVERSES
C         CALCUL DES ZEROS DU POLYNOME CARACTERISTIQUE PAR LA METHODE
C         DE LA SECANTE AVEC ACCELERATION
C
C ENTREES:
C --------
C NBROOT: NOMBRE DES PLUS PETITES VALEURS ET VECTEURS PROPRES A CALCULER
C NC    : NOMBRE MAJORANT DES VALEURS PROPRES >NBROOT AFIN DE PRENDRE EN
C         CONSIDERATION LES VALEURS MULTIPLES OU LES CALCULS DES VALEURS
C         SUPERIEURES A CELLES DEMANDEES
C N     : NOMBRE DE LIGNES ET COLONNES  DES MATRICES DE A ET B
C TETA  : TRANSLATION DE A => A + TETA * B DEJA EFFECTUEE
C
C MUA,IAAA,A,NCODSA : POINTEUR  , ADRESSE , MATRICE , CODE ( RAIDEUR )
C MUB,IAAB,B,NCODSB : POINTEUR  , ADRESSE , MATRICE , CODE ( MASSE )
C
C V0,V1,V,W: TABLEAUX AUXILIAIRES DE TAILLE N*(NOMBRE DE MOTS D UNE VARIABLE)
C VV    : VV(N,NC) TABLEAU DES VECTEURS PROPRES
C WW    : WW(N,NC) TABLEAU DE B FOIS LES VECTEURS PROPRES
C EIGV  : EIGV(N)  VALEURS PROPRES
C
C NITE  : NOMBRE D ITERATIONS TOTALES  PAR VALEUR PROPRE
C NINV  : NOMBRE D ITERATIONS INVERSES PAR VALEUR PROPRE
C MORAID: NOMBRE DE MOTS DE LA MATRICE DE SAUVEGARDE DE A
C RAID  : MATRICE DE SAUVEGARDE DE A
C
C SORTIES:
C --------
C NOV   : NOMBRE DE VALEURS PROPRES SAUTEES
C IERR  : 0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: BROCHARD-PERRONNET UPMC ANALYSE NUMERIQUE PARIS  AOUT 1976-98
C2345X7..............................................................012
      include"./incl/epsvvp.inc"
      COMMON            M(1)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (M(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      INTEGER          MUA(*),MUB(*),NITE(NC),NINV(NC)
      DOUBLE PRECISION RAID(*),A(*),B(*),V0(N),V1(N),V(N),W(N),
     +                 VV(N,NC),WW(N,NC),EIGV(NC),
     +                 TETA,DBLE
      DOUBLE PRECISION RQ0,RQ,XX,ERR,RC1,TOL,AL,RA,RR,ETA,
     +                 RQB,BS,RB,RT,RQT,EI,RC,
     +                 DETA,DETB,DETC,DETR
C
 100  FORMAT(' %% ERREUR SECINW : MATRICE RAIDEUR K NON INVERSIBLE')
 110  FORMAT('  ITERATION ',I5,' : VALEUR DE RB : ',G15.7)
 120  FORMAT(' MASSE INVERSIBLE (FICHIER ',I3,'). AMELIORATION DE LA '
     +       ,'VALEUR DE  RB : ',G15.7)
 130  FORMAT(' %% ATTENTION SECINW : MASSE NON INVERSIBLE , FICHIER',I6)
 140  FORMAT(' NOMBRE DE REAJUSTEMENT DE RB AU DEMARRAGE : ',I6)
 150  FORMAT(' DEBUT ITERATION INVERSE   RC           : '
     +       ,G15.7,' (',I6,' VALEURS SAUTEES ) ')
 160  FORMAT(' %% ERREUR SECINW : RELANCER EN JOUANT SUR TETA '/
     +       '    PRENDRE UNE VALEUR ENTRE 0. ET ',G15.7)
 170  FORMAT(' NOMBRE DE CONVERGENCES LENTES DE L ITERATION INVERSE : '
     +      ,I5/)
 180  FORMAT(1X,10I8)
 190  FORMAT(' RECHERCHE D''UNE BORNE INFERIEURE DE LA PLUS PETITE ',
     +        'VALEUR PROPRE ')
 200  FORMAT(/' **  SECANTE ABANDONNEE . NOMBRE D ITERATIONS : ',I6,
     +        ' POUR LA VALEUR PROPRE',I6/)
 210  FORMAT(' VALEUR DE DEPART DE LA SECANTE  (RB) : ',G15.7,
     +       ' ( ',I6,' VALEURS SAUTEES ) '/)
 220  FORMAT(' %% ERREUR SECINW : VERIFIEZ LES MATRICES ' )
 230  FORMAT(/'%% ERREUR: LE POLYNOME DEFLATIONNE N A PLUS DE NOEUDS'/)
 240  FORMAT(' %% ERREUR SECINW : PLUS DE PLACE POUR LES VECTEURS',
     +       ' PROPRES'/)
 250  FORMAT(' NOMBRE D''ITERATIONS SECANTE POUR CHAQUE VALEUR PROPRE :'
     +      )
 260  FORMAT(' NOMBRE D''ITERATIONS INVERSES POUR CHAQUE VALEUR PROPRE '
     +       ,':')
 270  FORMAT('  ITERATION INVERSE ',I6,' MAUVAISE CONVERGENCE ',
     +       ', REPRISE AVEC RC : ',G15.7)
 280  FORMAT(' CALCUL DE LA VALEUR DE DEPART DE LA SECANTE (RB) ')
 290  FORMAT(1X,19('*')/' VALEUR PROPRE',I6/1X,19('*'))
 300  FORMAT(' DEBUT SECANTE RA RB ETA :',3(G15.7,1X))
 310  FORMAT('  ITERATION SECANTE',I4,
     +       ' RA RB RC ETA          : ',4(G15.7,1X),'(',i3,
     +       ' VALEURS SAUTEES)')
 320  FORMAT(' %% ERREUR SECINW : ',I6,' ITERATIONS INVERSES.ABANDON ',
     +       ' SANS LA VALEUR ACTUELLE ',G15.2,
     +       ' DE CETTE VALEUR PROPRE',I5/)
 340  FORMAT(' VALEUR PROPRE ',I5,' DEFINITIVE         : ',G15.7/
     +       ' NOMBRE DE CONVERGENCES LENTES ITERATION INVERSE : ',I5/)
 350  FORMAT(I6,' VALEURS PROPRES SAUTEES AVANT     ',G15.7)
 360  FORMAT('  ITERATION INVERSE',I4,' RT RC RAYLEIGHT : ',
     +3(G15.7,1X))
 370  FORMAT(' VALEURS PROPRES CALCULEES:'/5(I5,':',G15.7))
C
C     INITIALISATIONS DES SEUILS ET AUTRES CONSTANTES
C     -----------------------------------------------
      IERR   = 0
C
C     SCLS SEUIL DE CONVERGENCE LENTE DE LA SECANTE POUR REAJUSTER ETA
      SCLS   = 1.E-01
C     SCSA SEUIL DE CONVERGENCE DE L ITERATION SECANTE ET DE LA BORNE INF
      SCSA   = 1.E-03
C     SCINV SEUIL DE CONVERGENCE DE L ITERATION INVERSE
      SCINV  = 1.E-5
C     SCLINV SEUIL DE CONVERGENCE LENTE DE L ITERATION INVERSE
      SCLINV = 1.E-02
C     SSDEF SEUIL DE STABILITE DE LA DEFLATION
      SSDEF  = 1.E-05
C     FAST FACTEUR PAR LEQUEL ON TRANSLATE LEGEREMENT MU POUR PRESERVER LA
C     FAST1 IDEM FAST STABILITE DE LA FACTORISATION MAIS POUR LE CALCUL DES BORN
      FAST1  = 1.E-04
C     IITEM NOMBRE D ITERATIONS MAX POUR CALCULER LA DEUXIEME VALEUR
C     DE DEPART DE LA SECANTE
      IITEM  = 10
C     NITEM NOMBRE D ITERATIONS MAX POUR CALCULER VALEURS ET VECTEURS PROPRES
      NITEM  = 30
C     NITIN NOMBRE D ITERATIONS INVERSES MAX PAR LE CALCUL D'UNE VALEUR PROPRE
      NITIN  = 100
C     NINVT NOMBRE D ITERATIONS INVERSES SANS TEST DE CONVERGENCE
      NINVT  = 7
C     NMAXII NOMBRE MAX D ITERATIONS INVERSES APRES LESQUELLES LA SECANTE REPART
      NMAXII = 2*NINVT
C
      KNOV = 0
      ETA  = 2
C     NOV EST LA DIFFENCE ENTRE LE NOMBRE DE PIVOTS NEG ET LE NOMBRE DE VALEURS
C     DEJA CALCULEES DONC LE NOMBRE DE VALEURS PROPRES SAUTEES
      NOV = 0
C     JR NUMERO DU NOEUD QUE L ON CALCULE
      JR  = 1
C     NSK VAUT 0 SAUF SI L ON A SAUTE PLUSIEURS VALEURS PROPRES LORS DE L
C     ITERATION SECANTE OU S IL ON A CONVERGE VERS LE VECTEUR PROPRE SUIVANT
C     LORS DE L ITERATION INVERSE ET ALORS NSK VAUT 1
      NSK = 0
      NIJ = 0
C     NIF NOMBRE MAX DE FOIS OU L ON REAJUSTE LA DEUX VALEUR DE DEPART DE LA SEC
C     AU CAS OU ELLE SE SITUE AU DELA DE LA PREMIERE VALEUR PROPRE
      NIF = 5
C
      RA      = 0D0
      RR      = 0D0
      NDSM    = 1
      NBROOT1 = NBROOT
C
C     RA ET RB SONT LES DEUX VALEURS NECESSAIRES A LA SECANTE
C     RC EST LE RESULTAT DE L'ITERATION DE LA SECANTE
C     RR SERT A STOCKER RA LORSQUE L ON PASSE D UN NOEUD DEJA CALCULE AU SUIVANT
C
      DO 1 I=1,NC
         NITE(I) = 0
         NINV(I) = 0
 1    CONTINUE
C
C     RECHERCHE DE LA BORNE INFERIEURE DE LA PLUS PETITE VALEUR PROPRE
C     ----------------------------------------------------------------
      NENTRE = 0
      CALL TRTATA( RAID, A, MORAID )
C     FACTORISATION A = L D TL
      CALL CRMC1D( MUA,A,N,EPS,NENTRE, A,NRETOU )
      IF( NRETOU .NE. 0 ) THEN
         WRITE(IMPRIM,100)
         IERR = 1
         RETURN
      ENDIF
C     DETERMINANT DE A = DETERMINANT DE D
      CALL DETETD( A,MUA,N, NBVI,NBHA,DETA )
      DETR = DETA
      NBHR = NBHA
      WRITE(IMPRIM,190)
C
C     CALCUL DE V0=DIAG(A)  V1=DIAG(B)
C     --------------------------------
      NIVO = 1
      DO 2 I=1,N
         XX = A(MUA(I+1))
         IF( XX .LT. 1D-20 ) XX=1D0
         V0(I) = XX
         V1(I) = B( MUB(I+1) )
 2    CONTINUE
C
C     IITE NOMBRE D ITERATIONS POUR CALCULER LA DEUXIEME VALEUR DE
C     DEMARRAGE DE LA SECANTE
      IITE   = 0
      NENTRE = 1
      WRITE(IMPRIM,280)
C
C     INITIALISATION DU QUOTIENT DE RAYLEIGH
C     --------------------------------------
      RT = 0
      CALL MAPRVE( 0, 1D0, N, NCODSB, MUB, B, V1,  W )
      NIVEAU = 3
C
 1000 IITE = IITE+1
      CALL DRCRPR(N,NCODSA,MUA,A,W,NIVEAU, V,IERR)
      RQT = 0D0
      DO 6 I=1,N
         RQT = RQT+W(I)*V(I)
 6    CONTINUE
C
C     CALCUL DE (V,BE) E VECTEUR UNITE  V=L(-T)BE
C     CALCUL DE (V,BV)
      CALL MAPRVE( 0, 1D0, N, NCODSB, MUB, B, V,  W )
      RQB = 0D0
      DO 7 I=1,N
         RQB = RQB+W(I)*V(I)
 7    CONTINUE
      RQ = RQT / RQB
      WRITE (IMPRIM,110) IITE,RQ
C
      BS  = SQRT( ABS(RQB) )
      TOL = ABS(RQ-RT) / RQ
      IF( TOL .GE. SCSA ) THEN
C
C        B-NORMALISATION
         DO 8 I=1,N
            W(I) = W(I) / BS
 8       CONTINUE
         RT = RQ
         IF( IITE .LT. IITEM ) GOTO 1000
C
C        SI AU BOUT DU NOMBRE D ITERATIONS MAX LE QUOTIENT DE RAYLEIGH
C        N'EST PAS STABILISE ON PREND RB ARBITRAIREMENT
         RB = RQ * ( 1.D0 - MIN(0.1D0,TOL) )
C
      ELSE
C
C        B EST-ELLE INVERSIBLE ?
C        -----------------------
         CALL TRTATA( RAID, A, MORAID )
         CALL MUA2PD(N,1D0,NCODSA,MUA,A,-RQ,NCODSB,MUB,B,MUA,A)
         CALL CRMC1D(MUB,B,N,EPS,NENTRE,B,NRETOU)
         IF( NRETOU .NE. 0 ) THEN
C
C           B N EST PAS INVERSIBLE
            RB = RQ*(1D0-0.01D0)
            WRITE(IMPRIM,130) NFMUB1
         ELSE
C
C           B EST INVERSIBLE. AMELIORATION DE LA VALEUR DE DEPART RB
            DO 9 I=1,N
               V(I) = V(I)/BS
 9          CONTINUE
            CALL MAPRVE( 0, 1D0, N, NCODSA, MUA, A, V,  W )
            CALL DRCRPR(N,NCODSB,MUB,B,W,NIVEAU, V,IERR)
C
C           CALCUL DE (V,B(-1)V)
            ERR = 0
            DO 10 I=1,N
               ERR = ERR + V(I) * W(I)
 10         CONTINUE
            RB = RQ - SQRT( ABS(ERR) )
            WRITE(IMPRIM,120) NFMUB1,RB
         ENDIF
      ENDIF
C
C     VERIFICATION RB < LAMDA1
C     ------------------------
C     IS COMPTE LE NOMBRE DE FOIS OU L ON A DU REAJUSTER RB AU DEMARRAGE
      IS  = 0
      RQB = 0D0
C
      BACKSPACE NFMUB1
      READ (NFMUB1) LE,IS,(M(IAAB-1+I),I=1,LE)
 1010 CALL TRTATA( RAID, A, MORAID )
      CALL MUA2PD(N,1.D00,NCODSA,MUA,A,-RB,NCODSB,MUB,B,MUA,A)
      CALL CRMC1D(MUA,A,N,EPS,NENTRE,A,NRETOU)
      IF( NRETOU .NE. 0 ) THEN
         RB = RB*(1D0-FAST)
         GOTO 1010
      ENDIF
      CALL DETETD( A,MUA,N, NBVI,NBHB,DETB )
      WRITE (IMPRIM,210) RB,NBVI
      IF( NBVI .NE. 0 ) THEN
         IS = IS+1
         WRITE (IMPRIM,140) IS
         IF( IS .LE. NIF ) THEN
            WRITE (IMPRIM,220)
            STOP
         ENDIF
         RB = RB / (NBVI+1)
         GOTO 1010
      ENDIF
C
C     -----------------------------------------------------------------
C     LES ITERATIONS DE LA SECANTE POUR CALCULER LE ZERO DU DETERMINANT
C     -----------------------------------------------------------------
C
 1111 WRITE (IMPRIM,290) JR
      WRITE (IMPRIM,370) (KK,EIGV(KK),KK=1,JR-1)
      IF(  JR .GT. NC ) GOTO 2222
      IF( (JR .GT. NBROOT1 .AND. NOV .EQ. 0) .OR.
     +    (JR .GT. NBROOT1+NOV)            ) GOTO 2222
      NJ     = 0
      NINVTE = NINVT
      NMAXI  = NMAXII
      WRITE (IMPRIM,300) RA,RB,ETA
C
 1040 IF( NBHA .NE. NBHB ) THEN
C         DETA ET DETB N'ONT PAS MEME NOMBRE D'HOMOTHETIES
          IF( NBHA .GT. NBHB ) THEN
              DETA = DETA / SCALE
              NBHA = NBHA - 1
              GOTO 1040
          ELSE
              DETA = DETA * SCALE
              NBHA = NBHA + 1
              GOTO 1040
          ENDIF
      ENDIF
C
      IF( DETA .EQ. DETB ) THEN
         WRITE (IMPRIM,230)
         RC = RA
         GOTO 1030
      ENDIF
C
C     LA NOUVELLE VALEUR DE LA SECANTE A ESSAYER
      RC   = RB - ETA * ( RB - RA ) / ( DETB - DETA ) * DETB
C     TEMOIN NOV N'EST PAS A JOUR
      KNOV = 1
      IF( ABS(RC-RB) .LE. SCSA * ABS(RC) ) GO TO 1030
C
C     PAS DE CONVERGENCE => UNE ITERATION SECANTE DE PLUS
 1050 CALL TRTATA( RAID, A, MORAID )
      CALL MUA2PD( N,
     +               1.D00,NCODSA,MUA,A,
     +               -RC,  NCODSB,MUB,B,  MUA,A )
      CALL CRMC1D( MUA,A,N,EPS,NENTRE,A,NRETOU )
      IF( NRETOU .NE. 0 ) THEN
         RC = RC * (1.D0-FAST)
         GOTO 1030
      ENDIF
      CALL DETETD( A,MUA,N, NBVI,NBHC,DETC )
      NITE(JR) = NITE(JR)+1
      JR1 = JR-1
      IF( JR1 .NE. 0 ) THEN
C        DEFLATION DU PROBLEME
         DO 11 K=1,JR1
            DETC = DETC / (RC-EIGV(K))
 11      CONTINUE
      ENDIF
C
C     SI L ON A SAUTE UNE VALEUR PROPRE AU MOINS L ITERATION INVERSE COMMENCE
C     NES COMPTE LE NOMBRE DE VALEURS PROPRES INFERIEURES A RC DEJA CALCULEES
      NES = 0
      IF( JR1 .NE. 0 ) THEN
         DO 12 I=1,JR1
            IF( EIGV(I) .LT. RC ) NES = NES+1
 12      CONTINUE
      ENDIF
      NOV  = NBVI - NES
C     TEMOIN NOV EST A JOUR
      KNOV = 0
      WRITE (IMPRIM,310) NITE(JR),RA,RB,RC,ETA,NOV
      IF( NOV .GT. 0 ) THEN
         EIGV(JR) = RC
         IF( NOV .GT. 1 ) NSK = 1
         GOTO 1030
      ENDIF
C
C     NSK VAUT 1 SI L ON A SAUTE PLUS D UNE VALEUR PROPRE
      RR   = RA
      DETR = DETA
      NBHR = NBHA
C
      RA   = RB
      DETA = DETB
      NBHA = NBHB
C
      RB   = RC
      DETB = DETC
      NBHB = NBHC
C
C     CONVERGENCE LENTE DE LA SECANTE ETA=ETA*2
C     -----------------------------------------
      TOL = RB*SCLS
      IF( ABS(RA-RB) .LT. TOL   ) ETA = ETA * 2
      IF( NITE(JR)   .LE. NITEM ) GOTO 1040
      WRITE (IMPRIM,200) NITE(JR),JR
      GOTO 2222
C
C     --------------------------------------------------------------------------
C     ITERATION INVERSE POUR AFFINER LA VALEUR PROPRE ET CALCULER LE VECTEUR PRO
C     --------------------------------------------------------------------------
C
C     PREPARATION DE L ITERATION INVERSE
 1030 IF( JR .GT. NC ) THEN
C        CONTROLE DE STOCKAGE
         WRITE (IMPRIM,240)
         GOTO 2222
      ENDIF
 1070 IF( KNOV .NE. 0 ) THEN
C
C        FACTORISATION DE  A - RC * B = L * D * TL
C        -----------------------------------------
         CALL TRTATA( RAID, A, MORAID )
         CALL MUA2PD( N,1.D00,NCODSA,MUA,A,-RC,NCODSB,MUB,B, MUA,A )
         CALL CRMC1D( MUA,A,N,EPS,NENTRE,A,NRETOU )
         IF(NRETOU.NE.0) THEN
            RC = RC * (1.D0-FAST)
            GOTO 1070
         ENDIF
C
C        CALCUL DE NOV NOMBRE DE VALEURS PROPRES CALCULEES INFERIEURES A RC
C        ------------------------------------------------------------------
         CALL DETETD( A,MUA,N, NBVI,NBHC,DETC )
         NES = 0
         JR1 = JR-1
         IF( JR1 .GT. 0 ) THEN
            DO 13 I=1,JR1
               IF( EIGV(I) .LT. RC ) NES = NES+1
 13         CONTINUE
         ENDIF
         NOV  = NBVI - NES
C        TEMOIN NOV EST A JOUR
         KNOV = 0
         WRITE(IMPRIM,350) NOV,RC
         IF( NOV .GT. 1 ) NSK = 1
         EIGV(JR) = RC
      ENDIF
C
      WRITE (IMPRIM,150) RC,NOV
      IF( NOV .LT. 0 ) THEN
         WRITE (IMPRIM,160) (EIGV(1)-TETA)/2D0
         STOP
      ENDIF
C     RT EST LA VALEUR PROPRE CALCULEE PAR ITERATION INVERSE
C     SA VALEUR INITIALE EST CELLE DE LA VALEUR FINALE DE LA SECANTE
      RT  = RC
      RQ0 = 0D0
C
C     INITIALISATION DU VECTEUR PROPRE
C     --------------------------------
      DO 14 I=1,N
         XX = ABS(V0(I)-RC*V1(I))
         IF( XX .LT. 1.D-20 ) XX = 1D0
         V(I) = V1(I) / XX
 14   CONTINUE
      CALL MAPRVE( 0, 1D0, N, NCODSB, MUB, B, V,  W )
      CALL DRCRPR( N,NCODSA,MUA,A,W,NIVEAU, V,IERR )
      GOTO 1080
C
C     ITERATION INVERSE PROPREMENT DITE
C     ---------------------------------
 1090 NINV(JR) = NINV(JR)+1
      CALL DRCRPR( N,NCODSA,MUA,A,W,NIVEAU, V,IERR )
C
C     RQ QUOTIENT DE RAYLEIGH
      RQ0 = RQ
      RQT = 0D0
      DO 15 I=1,N
         RQT = RQT + W(I) * V(I)
 15   CONTINUE
      CALL MAPRVE( 0, 1D0, N, NCODSB, MUB, B, V,  W )
      RQB = 0D0
      DO 16 I=1,N
         RQB = RQB + W(I) * V(I)
 16   CONTINUE
      RQ = RQT / RQB
      RT = RC  + RQ
C
      WRITE(IMPRIM,360) NINV(JR),RT,RC,RQ
      IF( ABS(RQ-RQ0) .LE. RQ * SCINV ) GOTO 1100
C
C     TEST DE CONVERGENCE LENTE DE L ITERATION INVERSE
C     ------------------------------------------------
 1080 IF( NINV(JR) .GT. NMAXI )  THEN
          NMAXI  = NMAXI    + NMAXII
          NINVTE = NINV(JR) + NINVT
          RC     = 0.6D0 * RB + 0.4D0 * RC
          ETA    = ETA * 0.5D0
          GOTO 1050
      ENDIF
      IF( NINV(JR) .LT. NINVTE ) GOTO 1110
      IF( ABS(RQ-RQ0) .LE. RQ * SCLINV ) GOTO 1110

C     REPRISE DES CONVERGENCES LENTES DE L ITERATION INVERSE
C     ------------------------------------------------------
      NJ     = NJ+1
      NMAXI  = NMAXI+NINVT
      NINVTE = NINV(JR)+NINVT

      GOTO (1120,1130,1130,1130,1130,1250),NJ
 1120 IF( JR.NE.1 .AND. NOV.NE.0 ) THEN
         EI = EIGV(JR-1)

C        3 CAS : EI>RC ET RT ,  EI<RC ET RT , RC OU RT<EI< RT OU RC
C        ----------------------------------------------------------
         IF(EI.GE.RC .OR. EI.GE.RT) THEN
            IF(EI.LE.RC .OR. EI.LE.RT) THEN

C              RC OU RT < EI < RT OU RC
C              ------------------------
               IF(2.*ABS(RT-EI) .LT. ABS(RC-EI)) GOTO 1140
ccc            IF(RC-RT) 1130,1130,1140
               IF( RC .GT. RT ) THEN
                  GOTO 1140
               ELSE
                  GOTO 1130
               ENDIF
            ENDIF

C           EI>RC ET EI>RT
C           --------------
ccc         IF(2.*(EI-RT)-(EI-RC)) 1140,1140,1130
            IF(2.*(EI-RT) .GT. (EI-RC)) THEN
               GOTO 1130
            ELSE
               GOTO 1140
            ENDIF
         ENDIF

 1140    RC = EI*(1.D0+SIGN(DBLE(SCLINV),RT-EI))
         GOTO 1160
      ENDIF
C
 1130 RC1 = RT*(1D0-SCLINV)
 1170 IF( ABS(RC1-RC) .LE. RC*SCLINV ) THEN
         RC1 = RC1*(1D0-SCLINV)
         GOTO 1170
      ENDIF
      RC = RC1
 1160 CALL TRTATA( RAID, A, MORAID )
      CALL MUA2PD(N,1.D00,NCODSA,MUA,A,-RC,NCODSB,MUB,B,MUA,A)
      CALL CRMC1D(MUA,A,N,EPS,NENTRE,A,NRETOU)
      IF(NRETOU.NE.0) THEN
         RC = RC * (1D0-FAST)
         GOTO 1160
      ENDIF
      CALL DETETD( A,MUA,N, NBVI,NBHC,DETC )
      NES = 0
      JR1 = JR-1
      IF( JR1 .NE. 0 ) THEN
         DO 17 I=1,JR1
            IF( EIGV(I) .LT. RC ) NES = NES + 1
 17      CONTINUE
      ENDIF
C     NOMBRE DE VALEURS PROPRES SAUTEES
      NOV  = NBVI - NES
C     TEMOIN NOV EST A JOUR
      KNOV = 0
      WRITE(IMPRIM,350) NOV,RC
      IF( NOV .GT. 1 ) NSK = 1
      WRITE (IMPRIM,270) NINV(JR),RC
C     NIJ = NOMBRE DE CONVERGENCES LENTES DE L ITERATION INVERSE
      NIJ = NIJ + 1
C
C     B-ORTHONORMALISATION DE GRAMM-SCHMIDT
C     -------------------------------------
 1110 IF( RQB .NE. 0D0 ) THEN
         BS = SQRT( ABS(RQB) )
         DO 18 I=1,N
            W(I) = W(I) / BS
 18      CONTINUE
      ENDIF
      JR1 = JR - 1
      IF(JR1 .NE. 0 ) THEN
         DO 19 K=1,JR1
            AL=0D0
            DO 20 I=1,N
               AL = AL + VV(I,K) * W(I)
 20         CONTINUE
            DO 21 I=1,N
               W(I) = W(I) - AL * WW(I,K)
 21         CONTINUE
 19      CONTINUE
      ENDIF
      IF( NINV(JR) .LE. NITIN ) GOTO 1090
 1250 WRITE (IMPRIM,320) NINV(JR),RT,JR
      GOTO 2222
C
C     **************************************
C     CALCUL DES VECTEURS ET VALEURS PROPRES
C     **************************************
 1100 CALL MAPRVE( 0, 1D0, N, NCODSB, MUB, B, V,  W )
      RQT = 0
      DO 22 I=1,N
         RQT = RQT + V(I) * W(I)
 22   CONTINUE
      RQB = RQT
C
C     LA VALEUR PROPRE JR DEFINITIVE
      EIGV(JR) = RT
      WRITE (IMPRIM,340) JR,RT,NJ
C
C     PASSAGE A LA VALEUR PROPRE SUIVANTE
      NJ = 0
      BS = SQRT( ABS(RQB) )
      DO 24 I=1,N
         W(I) = W(I) / BS
         V(I) = V(I) / BS
 24   CONTINUE
      DO 25 I=1,N
         WW(I,JR) = W(I)
         VV(I,JR) = V(I)
 25   CONTINUE
C
C     **************************************************************
C     DEFLATION DU POLYNOME CARACTERISTIQUE ET STRATEGIE DU CHOIX DE
C     RB,RR POUR PRESERVER LA STABILITE
C     **************************************************************
      TOL = SSDEF * EIGV(JR)
      IF( RA .LE. 0 ) THEN
         RA = RB / 2
C        TEMOIN NOV N'EST PAS A JOUR
         KNOV = 1
 1180    CALL TRTATA( RAID, A, MORAID )
         CALL MUA2PD( N,1.D00,NCODSA,MUA,A,-RA,NCODSB,MUB,B,MUA,A )
         CALL CRMC1D( MUA,A,N,EPS,NENTRE,A,NRETOU )
         IF( NRETOU .NE. 0 ) THEN
            RA = RA * (1D0-FAST)
            GOTO 1180
         ENDIF
         CALL DETETD( A,MUA,N, NBVI,NBHA,DETA )
      ENDIF
      IF( NOV .LE. 0 ) THEN
C
C        ON N A PAS SAUTE DE NOEUD
C        -------------------------
         IF( ABS(EIGV(JR)-RB) .GT. TOL ) GOTO 1200
C        LE NOEUD EST DONC PROCHE DE RB
         RB   = RA
         DETB = DETA
         NBHB = NBHA
C
         RA   = RR
         DETA = DETR
         NBHA = NBHR
         GOTO 1200
      ENDIF
C
C     ON A SAUTE AU MOINS UN NOEUD
C     ----------------------------
      IF( EIGV(JR) .GT. RC ) NSK = 1
C
C     DONC NSK VAUT 1 SOIT SI L ON A SAUTE PLUSIEURS NOEUDS
C     SOIT SI L ITERATION INVERSE A CONVERGE VERS LE VECTEUR PROPRE SUIVANT
C
      IF( NSK .EQ. 1 ) GOTO 1225
      IF( ABS(RC-EIGV(JR)) .LT. TOL ) GOTO 1220
C
C     LE NOEUD N EST PAS PROCHE DE RC
      IF( ABS(EIGV(JR)-RB) .GE. TOL ) THEN
C
C        LE NOEUD N EST PAS PROCHE DE RB
         RA   = RB
         DETA = DETB
         NBHA = NBHB
      ENDIF
C
      RB   = RC
      DETB = DETC
      GOTO 1200
C
 1220 IF( ABS(EIGV(JR)-RB) .LE. TOL ) THEN
         RB   = RA
         DETB = DETA
         NBHB = NBHA
C
         RA   = RR
         DETA = DETR
         NBHA = NBHR
      ENDIF
C
 1200 DETA = DETA / (RA-EIGV(JR))
      DETB = DETB / (RB-EIGV(JR))
      JR   = JR + 1
      ETA  = 2
      GOTO 1111
C
 1225 IF( ABS(EIGV(JR)-RB) .LE. TOL ) THEN
         RB   = RA
         DETB = DETA
         NBHB = NBHA
C
         RA   = RR
         DETA = DETR
         NBHA = NBHR
C
C        TEMOIN NOV N'EST PAS A JOUR
         KNOV = 1
      ENDIF
      DETA = DETA / (RA-EIGV(JR))
      DETB = DETB / (RB-EIGV(JR))
      DETR = DETR / (RR-EIGV(JR))
      IF( EIGV(JR) .LE. RC ) NOV = NOV - 1
      JR = JR + 1
      EIGV(JR) = RC
      NMAXI = NMAXII
      IF( NOV .GT. 0 ) GOTO 1030
      NSK = 0
      ETA = 2
      GOTO 1111
C
C     *****************************************************************
C     LES NBROOT1 VALEURS ET VECTEURS PROPRES ONT ETE CALCULES
C     CLASSEMENT DES VALEURS ET VECTEURS PROPRES DANS L ORDRE CROISSANT
C     *****************************************************************
 2222 WRITE (IMPRIM,350) NOV,RC
      NBROOT = JR - 1
      IF( NBROOT .EQ. 0 ) RETURN
      IF( JR .NE. 2 ) THEN
         JR = JR - 2
 2240    IS = 0
         DO 2260 I=1,JR
            IF( EIGV(I+1) .LT. EIGV(I) ) THEN
               IS        = IS + 1
               RT        = EIGV(I+1)
               EIGV(I+1) = EIGV(I)
               EIGV(I)   = RT
               K         = NITE(I+1)
               NITE(I+1) = NITE(I)
               NITE(I)   = K
               K         = NINV(I+1)
               NINV(I+1) = NINV(I)
               NINV(I  ) = K
               DO 2250 K=1,N
                  RT        = VV(K,I+1)
                  VV(K,I+1) = VV(K,I)
                  VV(K,I)   = RT
 2250          CONTINUE
            ENDIF
 2260    CONTINUE
         IF( IS .GT. 0 ) GOTO 2240
      ENDIF
      NBROOT = MIN( NBROOT, NBROOT1 )
C
      WRITE (IMPRIM,250)
      WRITE (IMPRIM,180) (NITE(J),J=1,NBROOT)
      WRITE (IMPRIM,260)
      WRITE (IMPRIM,180) (NINV(J),J=1,NBROOT)
      WRITE (IMPRIM,170) NIJ
      END
