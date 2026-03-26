        SUBROUTINE SUEX01( NTLXSU, LADEFI, RADEFI,
     %                     NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN QUADRANGLE DEFINI PAR LES
C -----    LIGNES STRUCTUREES DE CHACUN DES 4 COTES AVEC OU SANS TG
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU QUADRANGLE
C LADEFI : TABLEAU DE DEFINITION DE LA SURFACE PARTITIONNEE
C          CF ~/td/d/a_surface__definition
C RADEFI : TABLEAU REEL DE DEFINITION DE LA SURFACE
C          CES DEUX TABLEAUX LADEFI ET RADEFI ONT MEME ADRESSE A L'APPEL
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF ~/td/d/a___nsef
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       JUIN      1988
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS       NOVEMBRE  1988
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       OCTOBRE   1996
C....................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      INTEGER           NBSOCT(4),MNSOCT(4),MNSTLI(4),NUCOTE(4)
      INTEGER           NBARLI(4),MNARLI(4),
     %                  MNXYTG(4),MNNTGL(4),MNCGEF(4)
      REAL              XYZI(3,4),XYZF(3,4),FADIST(4)
C
      CHARACTER*24      KNOMLG
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      IERR   = 0
      NBCOOR = 3
      CALL AZEROI( 4, NBARLI )
      CALL AZEROI( 4, MNNTGL )
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE'
      NOFOTI = NOFOTIEL()
C
C     VERIFICATION DES TYPES DE CHAQUE COTE
C     =====================================
      NUTYSU = LADEFI(WUTYSU)
      NBTGS  = 0
      DO 20 N=1,4
C
C        LE NUMERO DE LA LIGNE N
         IF (NUTYSU.EQ.1) THEN
           NOLI = LADEFI(WU4COT-1+N)
         ENDIF
         IF (NUTYSU.EQ.2) THEN
           NOLI = LADEFI(WU4COR-1+N)
         ENDIF
CCC         IF (NUTYSU.EQ.7) THEN
CCC           NOLI = LADEFI(WU4CAT-1+N)
CCC         ENDIF
         IF( NOLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') NOLI
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='LIGNE DE NUMERO' // KERR(MXLGER)(1:10)
     %               // ' est INCORRECT'
            ELSE
               KERR(1)='NUMBER of LINE ' // KERR(MXLGER)(1:10)
     %               // ' is INCORRECT'
            ENDIF
            CALL LEREUR
            IERR = 9
            GOTO 10
         ENDIF
C
C        LE TABLEAU LEXIQUE DE CETTE LIGNE
         CALL LXNLOU( NTLIGN, NOLI, NTLXLI, MN )
         IF( NTLXLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') N
            IF( LANGAG .EQ. 0 ) THEN
              KERR(1) ='LIGNE INCONNUE SUR LE COTE '//KERR(MXLGER)(1:10)
            ELSE
              KERR(1) ='UNKNOWN LINE on the EDGE '//KERR(MXLGER)(1:10)
            ENDIF
            CALL LEREUR
            IERR = 5
            GOTO 10
         ENDIF
C
C        LE NOM DE LA LIGNE
         CALL NMOBNU( 'LIGNE', NOLI, KNOMLG )
C
C        LE TABLEAU 'NSEF' DE CETTE LIGNE
         CALL LXTSOU( NTLXLI, 'NSEF', NTARLI, MNARLI(N) )
         IF( NTARLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='LIGNE SANS ARETES: ' // KNOMLG
            ELSE
               KERR(1)='LINE without EDGES: ' // KNOMLG
            ENDIF
            CALL LEREUR
            IERR = 6
            GOTO 10
         ENDIF
C
C        LE TYPE DE LA LIGNE
         NUTYLI = MCN( MNARLI(N) + WUTYMA )
C
C        LA LIGNE EST ELLE FERMEE ?
         NUTFMA = MCN( MNARLI(N) + WUTFMA )
C
C        LE NOMBRE DE SOMMETS DE LA LIGNE
         NBSOCT(N) = MCN( MNARLI(N) + WBEFOB ) + 1
C
C        LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
         CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
         IF( NTSOLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='LIGNE SANS XYZSOMMET: ' // KNOMLG
            ELSE
               KERR(1)='LINE without XYZSOMMET: ' // KNOMLG
            ENDIF
            CALL LEREUR
            IERR = 8
            GOTO 10
         ENDIF
C        L'ADRESSE DU TABLEAU DES COORDONNEES DES SOMMETS DU COTE N
         MNSOCT(N) = MNSOLI + WYZSOM
         MNSTLI(N) = MNSOLI
C
         IF( NUTYLI .LE. 0 ) THEN
C           LIGNE NON STRUCTUREE
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LIGNE NON CONTINUE: ' // KNOMLG
            ELSE
               KERR(1) = 'LINE NOT CONTINUE: ' // KNOMLG
            ENDIF
            CALL LEREUR
            IERR = 1
         ENDIF
C
C        TEST SUR LE NOMBRE DE SOMMETS
         IF( NBSOCT(N) .LE. 1 .OR. NBSOCT(N) .GT. 100000 ) THEN
            NBLGRC(NRERR) = 3
            WRITE(KERR(MXLGER)(1:10),'(I10)') NBSOCT(N)-1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LIGNE : ' // KNOMLG
               KERR(2) = 'PAS ASSEZ OU TROP D''ARETES >100000'
               KERR(3) = KERR(MXLGER)(1:10) // ' ARETES'
            ELSE
               KERR(1) = 'LINE : ' // KNOMLG
               KERR(2) = 'NOT ENOUGH or TOO MORE EDGES >100000'
               KERR(3) = KERR(MXLGER)(1:10) // ' EDGES'
            ENDIF
            CALL LEREUR
            IERR = 4
         ENDIF
C
C        CALCUL DES 2 NUMEROS DES 2 TANGENTES DES ARETES DE LA LIGNE
         CALL TGARLI( MNARLI(N), MNSOLI,
     %                NBARLI(N), MNXYTG(N), MNNTGL(N), MNCGEF(N) )
C        LE NOMBRE DE COTES AVEC DES TANGENTES
         IF( NBARLI(N) .GT. 0 ) NBTGS = NBTGS + 1
C
 10      IF ( IERR .NE. 0 ) GOTO 800
 20   CONTINUE
C
C     LE QUADRANGLE ELLIPTIQUE EST SUPPOSE SANS TG DANS UN PREMIER TEMPS
      IF( NUTYSU .EQ. 2 ) NBTGS=0
C
C     VERIFICATION DE LA GEOMETRIE
C     ============================
C
C     1) PAS DE LIGNE FERMEE
C        RANGEMENT DES COTES SUIVANT L'ORDRE ET
C        POUR LE SENS:  C1:S1S2 C2:S2S3 C3:S3S4 C4:S4S1
C
C                       C3
C                S4----<-------S3
C                |             |
C                |             |
C         C4    \/             /\  C2
C                |             |
C                |             |
C                S1---->-------S2
C                       C1
C
      DO 30 N=1,4
C
C        LES COORDONNEES DU POINT INITIAL DU COTE N
         IAI = MNSOCT(N)
         XYZI(1,N)=RMCN(IAI  )
         XYZI(2,N)=RMCN(IAI+1)
         XYZI(3,N)=RMCN(IAI+2)
C
C        LES COORDONNEES DU POINT FINAL DU COTE N
         IAF = MNSOCT(N) + 3 * (NBSOCT(N)-1)
         XYZF(1,N)=RMCN(IAF  )
         XYZF(2,N)=RMCN(IAF+1)
         XYZF(3,N)=RMCN(IAF+2)
  30  CONTINUE
      CALL SUEXQ1( NUCOTE, XYZI, XYZF, IERR )
      IF( IERR .NE. 0 ) GOTO 800
C     NUCOTE(3) NUMERO DE 1 A 4 DU COTE 3 PARMI LES 4 LIGNES...
C               SI NUCOTE(3)=-4 LE COTE 3 EST LA LIGNE 4 A PARCOURIR
C                               EN SENS INVERSE DE SON RANGEMENT...
C               POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S4 C4:S4S1
C
C     VERIFICATION DU NOMBRE DE SOMMETS => CAS A TRAITER
C     2 COTES OPPOSES ONT-ILS MEME NOMBRE DE SOMMETS ?
      DO 55 I=1,3
         DO 42 J=I+1,4
            IF( ABS(NUCOTE(I)) .EQ. ABS(NUCOTE(J)) ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='ERREUR: QUADRANGLE NON FERME'
               ELSE
                  KERR(1) ='ERROR: QUADRANGLE NOT CLOSED'
               ENDIF
               CALL LEREUR
               IERR = 4
               GOTO 60
            ENDIF
 42     CONTINUE
 55   CONTINUE
C
C     SI LA FONCTION TAILLE_IDEALE(x,y,z) EXISTE ALORS TRIANGULATION FORCEE
      IF( NOFOTI .GT. 0 ) GOTO 58
      IF( .NOT.( (NUTYSU .EQ. 1) .AND.
     %    ( NBSOCT(ABS(NUCOTE(1))) .NE. NBSOCT(ABS(NUCOTE(3)))  .AND.
     %      NBSOCT(ABS(NUCOTE(2))) .NE. NBSOCT(ABS(NUCOTE(4))))))GOTO 61
C
C     QUADRANGLE ALGEBRIQUE A TRIANGULER
C     ----------------------------------
 58   DO 59 I=1,4
C        L'ADRESSE MCN DU TABLEAU XYZSOMMET DES 4 LIGNES COTES DU QUADRANGLE
C        DANS L'ORDRE DE LA DONNEE UTILISATEUR
         MNSOCT(I) = MNSOCT(I) - WYZSOM
 59   CONTINUE
C     NUCOTE(I) = +-NO 1 A 4 DE LA LIGNE COTE I DU QUADRANGLE POUR LE SENS DIREC
      CALL QUNSAL( NTLXSU, NUCOTE, NBTGS,  NBSOCT, MNSOCT,
     %             NBARLI, MNXYTG, MNNTGL,
     %             NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
      GOTO 9000
C
 61   IF( (NUTYSU .NE. 1) .AND.
     %    ( NBSOCT(ABS(NUCOTE(1))) .NE. NBSOCT(ABS(NUCOTE(3)))  .OR.
     %      NBSOCT(ABS(NUCOTE(2))) .NE. NBSOCT(ABS(NUCOTE(4))) ) )THEN
         NBLGRC(NRERR) = 3
         WRITE(KERR(MXLGER)(1:40),'(4(I10))')
     %                     (NBSOCT(ABS(NUCOTE(K)))-1,K=1,4)
         IF( LANGAG .EQ. 0 ) THEN
        KERR(1)='LE NOMBRE DES ARETES DES COTES 1-3 DOIT ETRE IDENTIQUE'
        KERR(2)='LE NOMBRE DES ARETES DES COTES 2-4 DOIT ETRE IDENTIQUE'
        KERR(3)='NOMBRE DES ARETES DES 4 COTES :'//KERR(MXLGER)(1:40)
         ELSE
            KERR(1)='The NUMBER of SUB-EDGES of EDGES 1-3 MUST BE EQUAL'
            KERR(2)='The NUMBER of SUB-EDGES of EDGES 2-4 MUST BE EQUAL'
      KERR(3)='NUMBER of SUB-EDGES of the 4 EDGES :'//KERR(MXLGER)(1:40)
         ENDIF
         CALL LEREUR
         IERR = 4
      ENDIF
C
 60   IF( IERR .GT. 0 ) GOTO 9999
C
C     PERMUTATIONS POUR AMENER LES COTES 1-3 A AVOIR LE MEME NOMBRE DE SOMMETS
C     LE COTE 2 A AVOIR PLUS DE SOMMETS QUE LE COTE 4
C     ========================================================================
      IF( NBSOCT(ABS(NUCOTE(1))) .NE. NBSOCT(ABS(NUCOTE(3))) ) THEN
         CALL PECI4R( NUCOTE )
      ENDIF
      IF( NBSOCT(ABS(NUCOTE(2))) .LT. NBSOCT(ABS(NUCOTE(4))) ) THEN
         DO 68 J=1,2
            CALL PECI4R( NUCOTE )
68       CONTINUE
      ENDIF
      IF( NBSOCT(ABS(NUCOTE(2)))-NBSOCT(ABS(NUCOTE(4))) .GT.
     %    NBSOCT(ABS(NUCOTE(1)))-1 ) THEN
         GOTO 58
      ENDIF
C
C     GENERATION DU TABLEAU XYZSOMMET (COORDONNEES DES SOMMETS ET TANGENTES)
C     ======================================================================
C     LE NOMBRE DE NSEF DES COTES 1-3 ET 2-4
      NBS1 = NBSOCT(ABS(NUCOTE(1)))
      NBA1 = NBS1 - 1
      NBS2 = NBSOCT(ABS(NUCOTE(2)))
      NBA2 = NBS2 - 1
      NBS4 = NBSOCT(ABS(NUCOTE(4)))
      N24  = NBS2 - NBS4
      NDSQ = NBS1 - N24
      IF( N24 .EQ. 0 ) THEN
C        LE NOMBRE DE SOMMETS DE LA QUADRANGULATION
         NBS = NBS1 * NBS2
C        LE NOMBRE DE TANGENTES STOCKEES
         IF( NBTGS .GT. 0 ) NBTGS = 8 * NBA1 * NBA2
      ELSE
C        LE NOMBRE DE SOMMETS DE LA TRIANGULATION-QUADRANGULATION
         NBS  = NBS1 * NBS4 + N24 * (N24+1) / 2
         NBA4 = NBS4 - 1
         NDSQ = NBS1 - N24
C        LE NOMBRE DE QUADRANGLES
         NBQ  = (NDSQ-1) * NBA4
C        LE NOMBRE DE TRIANGLES
         I = 1
         DO 69 J=2,N24
            I = I + 2 * J - 1
 69      CONTINUE
         NBT = (NBS1-NDSQ)*2*NBA4 + I
C        LE NOMBRE D'EF DU MAILLAGE
         NBEF = NBT + NBQ
C        LE NOMBRE DE TANGENTES STOCKEES
         IF( NBTGS .GT. 0 ) NBTGS = 6*NBT + 8*NBQ
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %             WYZSOM + NBCOOR * (NBS+NBTGS) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )
C     LE NOMBRE DE COORDONNEES DES SOMMETS
      MCN( MNSOFA + WBCOOR ) = NBCOOR
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM ) = NBS
C     LE NOMBRE DE TANGENTES
      MCN( MNSOFA + WNBTGS ) = NBTGS
C     L'ADRESSE DU DEBUT DES COORDONNEES DES SOMMETS
      MNCOSO = MNSOFA + WYZSOM
C     L'ADRESSE DU DEBUT DES COMPOSANTES DES EVENTUELLES TANGENTES
      MNCOTG = MNCOSO + NBCOOR * NBS
C
C     REMPLISSAGE DU TABLEAU 'XYZSOMMET' : LES SOMMETS DU BORD
C     ----------------------------------- --------------------
C     LE COTE 1
      N = NUCOTE(1)
      DO 75 NP=1,NBSOCT(ABS(N))
         IF( N .GE. 0 ) THEN
            IA = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA = MNSOCT(-N) - 1 + 3 * (NBSOCT(-N)-NP)
         ENDIF
         IN = MNCOSO - 1 + 3 * (NP-1)
         DO 70 I=1,3
            RMCN(IN+I) = RMCN(IA+I)
  70     CONTINUE
  75  CONTINUE
C
C     LE COTE 2
      N = NUCOTE(2)
      DO 85 NP=2,NBSOCT(ABS(N))-1
         IF( N .GE. 0 ) THEN
            IA = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA = MNSOCT(-N) - 1 + 3 * (NBSOCT(-N)-NP)
         ENDIF
         NUM = NUSOTQ( NBS1, NBS2, NBS4, NBS1, NP )
         IN  = MNCOSO - 1 + 3 * (NUM-1)
         DO 80 I=1,3
            RMCN(IN+I) = RMCN(IA+I)
  80     CONTINUE
  85  CONTINUE
C
C     LE COTE 3
      N = NUCOTE(3)
      DO 95 NP=1,NBSOCT(ABS(N))
         IF( N .LE. 0 ) THEN
            IA = MNSOCT(-N) - 1 + 3 * (NP-1)
         ELSE
            IA = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         IF( NP .LE. NDSQ ) THEN
            I = NBS4
         ELSE
            I = NBS4 + NP - NDSQ
         ENDIF
         NUM = NUSOTQ( NBS1, NBS2, NBS4, NP, I )
         IN  = MNCOSO - 1 + 3 * (NUM-1)
         DO 90 I=1,3
            RMCN(IN+I) = RMCN(IA+I)
  90     CONTINUE
  95  CONTINUE
C
C     LE COTE 4
      N = NUCOTE(4)
      DO 102 NP=2,NBSOCT(ABS(N))-1
         IF( N .LE. 0 ) THEN
            IA = MNSOCT(-N) - 1 + 3 * (NP-1)
         ELSE
            IA = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         NUM = NBS1 * (NP-1) + 1
         IN  = MNCOSO - 1 + 3 * (NUM-1)
         DO 100 I=1,3
            RMCN(IN+I) = RMCN(IA+I)
 100     CONTINUE
 102  CONTINUE
C
C     CONSTRUCTION DU TABLEAU 'NSEF' TOPOLOGIE DE CETTE SURFACE
C     =========================================================
      IF( N24 .GT. 0 )  THEN
C
C        TRIANGULATION-QUADRANGULATION ALGEBRIQUE (COOK GORDON HALL)
C        -----------------------------------------------------------
C        GENERATION DU TABLEAU NSEF
         CALL TQUADR( NTLXSU, NBS1, NBS2, NBS4, NBTGS,
     %                NTFASU, MNFASU )
C        GENERATION DES SOMMETS INTERNES ET DES TANGENTES DU TABLEAU XYZSOMMET
         CALL TNMCDC( 'REEL', 2*NBS, MNCARR )
         CALL TNMCDC( 'REEL', 6*NBS, MN2DER )
         DO 105 I=1,4
C           L'ADRESSE MCN DU TABLEAU XYZSOMMET DES 4 LIGNES COTES DU QUADRANGLE
            MNSOCT(I) = MNSOCT(I) - WYZSOM
 105     CONTINUE
         CALL SQUADR( NUCOTE, NBTGS,  NBSOCT, MNSOCT,
     %                NBARLI, MNXYTG, MNNTGL, MNSOFA, MNFASU,
     %                MNCARR, RMCN(MNCARR), MN2DER, RMCN(MN2DER) )
         CALL TNMCDS( 'REEL', 6*NBS, MN2DER )
         CALL TNMCDS( 'REEL', 2*NBS, MNCARR )
         GOTO 800
      ENDIF
C
C     CAS D'UNE QUADRANGULATION STRUCTUREE ALGEBRIQUE (COOK GORDON HALL)
C     ------------------------------------------------------------------
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0 OU C1
      IF( NBTGS .LE. 0 ) THEN
C        EF SANS TG
         NBTGEF = 0
         N      = 0
      ELSE
C        EF AVEC TG
         NBTGEF = 8
         N      = NBA1 * NBA2
      ENDIF
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER',
     %             1 + WBARYQ + N * ( 2 + NBTGEF ))
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU, MNFASU )
C     TYPE DE L'OBJET : SURFACE
      MCN ( MNFASU + WUTYOB ) = 3
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN ( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0 OU C1
      MCN ( MNFASU + WBTGEF ) = NBTGEF
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBA1 * NBA2
C     LE NOMBRE D'EF AVEC TANGENTES DE LA SURFACE
      MCN ( MNFASU + WBEFTG ) = N
C     LE NOMBRE D'EF AVEC POINTEUR SUR LES EF  A TG DE LA SURFACE
      MCN ( MNFASU + WBEFAP ) = N
C     NUMERO DU TYPE DU MAILLAGE : QUADRANGLE STRUCTURE
      MCN ( MNFASU + WUTYMA ) = 4
C     VARIABLE NBARXQ : NOMBRE DE SEGMENTS SUIVANT X
      MCN ( MNFASU + WBARXQ ) = NBA1
C     VARIABLE NBARYQ : NOMBRE DE SEGMENTS SUIVANT Y
      MCN ( MNFASU + WBARYQ ) = NBA2
C
      IF( NBTGS .GT. 0 ) THEN
C        LE TABLEAU DES POINTEURS ET DES TANGENTES
         MN = MNFASU + WBARYQ
         DO 107 NP=1,N
            MN = MN + 1
            MCN(MN) = NP
C           CODE GEOMETRIQUE 'INTERPOLATION TRANSFINI'
            MCN(MN+N) = 16
 107     CONTINUE
         MN = MN + N
         DO 108 NP=1,8*N
            MCN(MN+NP) = NP
 108     CONTINUE
      ENDIF
C
      IF( NUTYSU .EQ. 1 .OR. NUTYSU .EQ. 2 ) THEN
C
C        LAGRANGE AVEC QUADRANGULATION STRUCTUREE DU CARRE UNITE
C                 METHODE DE COOK-GORDON-HALL
C        -------------------------------------------------------
C        RANGEMENT DES TANGENTES DES COTES SUIVANT L'ORDRE ET
C        POUR LE SENS:  C1:S1S2 C2:S2S3 C3:S4S3 C4:S1S4
C
C                       C3
C                S4---->-------S3
C                |             |
C                |             |
C         C4    /\             /\  C2
C                |             |
C                |             |
C                S1---->-------S2
C                       C1
         DO 111 NP=1,2
C           LE COTE NP EST FAIT LA LIGNE N
            N = NUCOTE(NP)
            IF( N .LT. 0 ) THEN
C              LA LIGNE DU COTE NP EST EN SENS INVERSE
C              LES DERNIERS NUMEROS DEVIENNENT LES PREMIERS AVEC INVERSION DU SI
               CALL RENVER( 2*NBARLI(-N), MCN(MNNTGL(-N)) )
            ENDIF
 111     CONTINUE
         DO 112 NP=3,4
C           LE COTE NP EST FAIT LA LIGNE N
            N = NUCOTE(NP)
            IF( N .GT. 0 ) THEN
C              LA LIGNE DU COTE NP EST EN SENS INVERSE
C              LES DERNIERS NUMEROS DEVIENNENT LES PREMIERS AVEC INVERSION DU SI
               CALL RENVER( 2*NBARLI(N), MCN(MNNTGL(N)) )
            ENDIF
 112     CONTINUE
C
C        LES TABLEAUX STCARR (XY DES SOMMETS DU CARRE UNITE) ET TGPURE
         NBCARR = 8 * NBS1 * NBS2
         CALL TNMCDC( 'REEL', NBCARR, MNSTCA )
         MNTGPU = MNSTCA + 2 * NBS1 * NBS2
C
         CALL SUEXQ2( NBS1, NBS2, NBA1, NBA2, NBTGS,
     S                MNXYTG(ABS(NUCOTE(1))), MNNTGL(ABS(NUCOTE(1))),
     S                MNXYTG(ABS(NUCOTE(2))), MNNTGL(ABS(NUCOTE(2))),
     S                MNXYTG(ABS(NUCOTE(3))), MNNTGL(ABS(NUCOTE(3))),
     S                MNXYTG(ABS(NUCOTE(4))), MNNTGL(ABS(NUCOTE(4))),
     S                RMCN(MNCOSO), RMCN(MNCOTG),
     S                RMCN(MNSTCA), RMCN(MNTGPU) )
C
C        DESTRUCTION DES TABLEAUX AUXILIAIRES
         CALL TNMCDS( 'REEL', NBCARR, MNSTCA )
C
C        FIN DU CAS ALGEBRIQUE
         IF( NUTYSU .EQ. 1 ) GOTO 800
C
C        CAS ELLIPTIQUE: REMISE EN ETAT DES TANGENTES => DANS LE SENS DIRECT
         DO 113 NP=1,2
C           LE COTE NP EST FAIT LA LIGNE N
            N = NUCOTE(NP)
            IF( N .LT. 0 ) THEN
C              LA LIGNE DU COTE NP EST EN SENS INVERSE
C              LES DERNIERS NUMEROS DEVIENNENT LES PREMIERS AVEC INVERSION DU SI
               CALL RENVER( 2*NBARLI(-N), MCN(MNNTGL(-N)) )
            ENDIF
 113     CONTINUE
         DO 114 NP=3,4
C           LE COTE NP EST FAIT LA LIGNE N
            N = NUCOTE(NP)
            IF( N .GT. 0 ) THEN
C              LA LIGNE DU COTE NP EST EN SENS INVERSE
C              LES DERNIERS NUMEROS DEVIENNENT LES PREMIERS AVEC INVERSION DU SI
               CALL RENVER( 2*NBARLI(N), MCN(MNNTGL(N)) )
            ENDIF
 114     CONTINUE
      ENDIF
C
      IF ( NUTYSU .EQ. 2 ) THEN
C
C        MAILLAGE ELLIPTIQUE METHODE DE WINSLOW
C        --------------------------------------
C        L'INITIALISATION DE LA QUADRANGULATION PAR COOK GORDON HALL
C        A ETE FAITE AU DESSUS EN NE PRENANT PAS EN COMPTE LES TANGENTES!
C
C        RESERVATION DE LA PLACE NECESSAIRE AUX TABLEAUX TEMPORAIRES :
C        DANS L'ORDRE LES MATRICES 3 COLONNES MAT ET A
C                     LES TABLEAUX 2 COLONNES B ET TRAV
C                     LES VECTEURS VEC ITRAV FCI ET FCJ
C
        LMAT   = 3*(NBS1*NBS2)
        NADCOS = 0
        CALL TNMCDC( 'REEL', LMAT, NADCOS )
        LMAT   = 3*(NBS1*NBS2)
        NADMAT = 0
        CALL TNMCDC( 'REEL', LMAT, NADMAT )
        NADA   = 0
        CALL TNMCDC( 'REEL', LMAT, NADA )
        LTAB   = 2*(NBS1*NBS2)
        NADB   = 0
        CALL TNMCDC( 'REEL', LTAB, NADB )
        NADTRA = 0
        CALL TNMCDC( 'REEL', LTAB, NADTRA )
        LVEC   = NBS1*NBS2
        NADVEC = 0
        CALL TNMCDC( 'REEL', LVEC, NADVEC )
        NADITR = 0
        CALL TNMCDC( 'REEL', LVEC, NADITR )
        NADFCI = 0
        CALL TNMCDC( 'REEL', LVEC, NADFCI )
        NADFCJ = 0
        CALL TNMCDC( 'REEL', LVEC, NADFCJ )
        LCOB   = 6*(NBS1+NBS2)-9
        NADCOB = 0
        CALL TNMCDC( 'REEL', LCOB, NADCOB )
        LALP   = 2*(NBS1+NBS2)-3
        NADALP = 0
        CALL TNMCDC( 'REEL', LALP, NADALP )
        LDIS   = 2*(NBS1+NBS2-2)
        NADDIS = 0
        CALL TNMCDC( 'REEL', LDIS, NADDIS )
C
C       LE NOMBRE D'ITERATIONS DE CORRECTION
        NITCOR = LADEFI( WITC2D )
        IF( NITCOR .LE. 0 .OR. NITCOR .GT. 100 ) THEN
C          PAS D'ITERATIONS DE CORRECTION
           NITCOR = 0
           DO 119 I=1,4
C             PAS DE CORRECTION SUR LES DISTANCES
              FADIST(I) = 0.
 119       CONTINUE
        ELSE
           DO 120 I=1,4
              FADIST(I) = RADEFI(WADI2D+ABS(NUCOTE(I))-1)
  120      CONTINUE
        ENDIF
C
        CALL SUEXW3( NBS1,  NBS2,  NBS1*NBS2,  RMCN(MNCOSO),
     S               RMCN(NADMAT), RMCN(NADA), RMCN(NADB),
     S               RMCN(NADFCI), RMCN(NADFCJ),
     S               RMCN(NADITR), RMCN(NADTRA),
CCC     %                RMCN(NADVEC),
     S               RMCN(NADCOB), RMCN(NADALP), RMCN(NADDIS),
     S               NITCOR, FADIST, IERR )
C
C       DESTRUCTION DES TABLEAUX DE TRAVAIL
        LMAT   = 3*(NBS1*NBS2)
        CALL TNMCDS( 'REEL', LMAT, NADCOS )
        LMAT   = 3*(NBS1*NBS2)
        CALL TNMCDS( 'REEL', LMAT, NADMAT )
        CALL TNMCDS( 'REEL', LMAT, NADA )
        LTAB   = 2*(NBS1*NBS2)
        CALL TNMCDS( 'REEL', LTAB, NADB )
        CALL TNMCDS( 'REEL', LTAB, NADTRA )
        LVEC   = NBS1*NBS2
        CALL TNMCDS( 'REEL', LVEC, NADVEC )
        CALL TNMCDS( 'REEL', LVEC, NADITR )
        CALL TNMCDS( 'REEL', LVEC, NADFCI )
        CALL TNMCDS( 'REEL', LVEC, NADFCJ )
        LCOB   = 6*(NBS1+NBS2)-9
        CALL TNMCDS( 'REEL', LCOB, NADCOB )
        LALP   = 2*(NBS1+NBS2)-3
        CALL TNMCDS( 'REEL', LALP, NADALP )
        LDIS   = 2*(NBS1+NBS2-2)
        CALL TNMCDS( 'REEL', LDIS, NADDIS )
C
C       AJOUT DES EVENTUELLES TANGENTES SUR LA FRONTIERE DU QUADRANGLE STRUCTURE
C       PRESENTATION DES DONNEES DES LIGNES
C       SELON LE SENS C1:S1S2 C2:S2S3  C3:S4S3 C4:S1S4
        NUCOTE(3) = -NUCOTE(3)
        NUCOTE(4) = -NUCOTE(4)
C
C                       C3
C                S4---->-------S3
C                |             |
C                |             |
C           C4  / \           / \ C2
C                |             |
C                |             |
C                S1---->-------S2
C                       C1
C        NBSOCT(1:4) SERT DE TABLEAU AUXILIAIRE POUR LA REMISE EN ORDRE DE MNSOC
         DO 122 I=1,4
            NBSOCT(I) = MNSOCT( ABS( NUCOTE(I) ) )
 122     CONTINUE
         DO 123 I=1,4
            MNSOCT(I) = NBSOCT(I)
 123     CONTINUE
         DO 124 I=1,4
            NBSOCT(I) = MNNTGL( ABS( NUCOTE(I) ) )
 124     CONTINUE
         DO 127 I=1,4
C           L'ADRESSE MCN DU TABLEAU XYZSOMMET DES 4 LIGNES COTES DU QUADRANGLE
            MNSOCT(I) = MNSOCT(I) - WYZSOM
C           REMISE EN ORDRE DE MNNTGL ET INVERSION EVENTUELLE
            MNNTGL(I) = NBSOCT(I)
            N = ABS( NUCOTE(I) )
            IF( MNNTGL(I) .GT.0 .AND. NUCOTE(I) .LT. 0 ) THEN
               CALL RENVER( 2*NBARLI(N), MCN(MNNTGL(I)) )
            ENDIF
 127     CONTINUE
         CALL TGQUST( MNSOCT, MNNTGL,
     %                NTLXSU, MNSOFA, MNFASU )
C        POUR LA DESTRUCTION DES TABLEAUX MNNTGL LE TABLEAU NBARLI
C        EST REMIS EN ORDRE
         DO 128 I=1,4
            NBSOCT(I) = NBARLI( ABS( NUCOTE(I) ) )
 128     CONTINUE
         DO 129 I=1,4
            NBARLI(I) = NBSOCT(I)
 129     CONTINUE
C
CCC      ELSE IF ( NUTYSU .EQ. 7 ) THEN
CCCC
CCCC        MAILLAGE ALGEBRIQUE : PROJECTIONS SUR CERCLES OU CYLINDRES
CCCC        ----------------------------------------------------------
CCCC        RESERVATION DE LA PLACE NECESSAIRE AUX TABLEAUX TEMPORAIRES :
CCCC        DANS L'ORDRE LES MATRICES 3 COLONNES MAT ET A
CCCC                     LES TABLEAUX 2 COLONNES B ET TRAV
CCCC                     LES VECTEURS VEC ITRAV FCI ET FCJ
CCCC
CCC        NBRPRO = LADEFI(WBRPRO)
CCC        IF (NBRPRO.LT.1) THEN
CCC           NBLGRC(NRERR) = 1
CCC           KERR(1) =  'SUEX01: NOMBRE ERRONE DE PROJECTIONS'
CCC           CALL LEREUR
CCC           IERR = 1
CCC           GOTO 9999
CCC        ENDIF
CCC        MNRAYO = 0
CCC        CALL TNMCDC(  'REEL' , NBRPRO, MNRAYO )
CCC        LONG   = 4*NBRPRO
CCC        MNLIMI = 0
CCC        CALL TNMCDC( 'ENTIER', LONG, MNLIMI )
CCC        LONG   = 3*NBRPRO
CCC        MNCEN  = 0
CCC        MNVECN = 0
CCC        CALL TNMCDC(  'REEL' , LONG, MNCEN )
CCC        CALL TNMCDC(  'REEL' , LONG, MNVECN )
CCCC
CCC        CALL SUEX07( NBS1,NBS2,LADEFI,RADEFI,NBRPRO,
CCC     S               RMCN(MNLIMI),RMCN(MNCEN),RMCN(MNVECN),
CCC     S               RMCN(MNRAYO),RMCN(MNCOSO), IERR )
CCCC
CCC        CALL TNMCDS(  'REEL' , NBRPRO, MNRAYO )
CCC        CALL TNMCDS( 'ENTIER', 4*NBRPRO, MNLIMI )
CCC        LONG = 3*NBRPRO
CCC        CALL TNMCDS(  'REEL' , LONG, MNCEN )
CCC        CALL TNMCDS(  'REEL' , LONG, MNVECN )
C
      ELSE
C
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:2),'(I2)') NUTYSU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'METHODE INCONNUE NUTYSU=' // KERR(MXLGER)(1:2)
         ELSE
            KERR(1) = 'UNKNOWN METHOD NUTYSU=' // KERR(MXLGER)(1:2)
         ENDIF
         CALL LEREUR
         IERR = 9
C
      ENDIF
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
 800  DO 810 N=1,4
C        LES EVENTUELS TABLEAUX CREES DANS L'EXECUTION DU SP TGARLI
         IF( MNNTGL(N) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', 3*NBARLI(N), MNNTGL(N) )
         ENDIF
 810  CONTINUE
      IF( IERR .NE. 0 ) THEN
         IF( NTSOFA .GT. 0 ) CALL LXTSDS( NTLXSU, 'XYZSOMMET' )
         IF( NTFASU .GT. 0 ) CALL LXTSDS( NTLXSU, 'NSEF' )
         GOTO 9999
      ENDIF
C
C     FIN DE LA CONSTRUCTION DE LA SURFACE
C     ====================================
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     SURFACE NON FERMEE
      MCN( MNFASU + WUTFMA ) = 0
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     SUPPRESSION DES TANGENTES DOUBLES
C     =================================
 9000 CALL MOINTG( NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C
C     IDENTIFICATION ET IMPOSITION DES COORDONNEES DES SOMMETS
C     DES LIGNES DU CONTOUR AU SOMMET LE PLUS PROCHE DE LA SURFACE
C     ============================================================
      IF( MNSOFA .GT. 0 ) THEN
         DO 9010 I=1,4
            IF( MNSTLI(I) .GT. 0 ) THEN
               CALL IDLISU( MNSTLI(I), MNSOFA )
            ENDIF
 9010    CONTINUE
      ENDIF
C
 9999 RETURN
      END
