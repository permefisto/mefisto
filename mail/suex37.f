      SUBROUTINE SUEX37( NTLXSU, LADEFI,
     %                   NTFASE, MNFASE, NTSOSE, MNSOSE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   EXTRAIRE UNE SURFACE D'UN VOLUME A PARTIR D'UN CRITERE LOGIQUE
C -----
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI : TABLEAU ENTIER DE DEFINITION DU VOLUME
C          CF $MEFISTO/td/d/a_surface__definition
C
C SORTIES:
C --------
C NTFASE : NUMERO      DU TMS 'NSEF' DE LA SURFACE EXTRAITE
C MNFASE : ADRESSE MCN DU TMS 'NSEF' DE LA SURFACE EXTRAITE
C          CF $MEFISTO/td/d/a___nsef
C NTSOSE : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE EXTRAITE
C MNSOSE : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE EXTRAITE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          1 SI VOLUME INITIAL INCONNU
C          2 SI VOLUME INITIAL SANS NSEF
C          3 SI VOLUME INITIAL SANS XYZSOMMET
C          4 SI FONCTION INCONNUE
C          5 SI TYPE INCONNU DE FACES
C          6 SI IMPOSSIBLE DE CREER LES FACES DU VOLUME
C         10 SI AUCUNE FACE EXTRAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1991
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___face.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      DOUBLE PRECISION  DXYZ(3),DOUI
      LOGICAL           OUI
      CHARACTER*24      NMOBJT
C
C     LE VOLUME INITIAL
C     =================
C     LE NOM DE CE VOLUME
      NUVOSE = LADEFI( WUVOSE )
C     LE TABLEAU LEXIQUE DE CE VOLUME
      CALL LXNLOU( NTVOLU, NUVOSE, NTLXVL, MN )
      IF( NTLXVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INITIAL INCONNU'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
      CALL NMOBNU( 'VOLUME', NUVOSE, NMOBJT )
C
C     LE TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'NSEF', NTCUVL, MNCUVL )
      IF( NTCUVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS NSEF'
         ELSE
            KERR(1) = 'VOLUME WITHOUT NSEF'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'XYZSOMMET', NTSOVL, MNSOVL )
      IF( NTSOVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS XYZSOMMET'
         ELSE
            KERR(1) = 'VOLUME WITHOUT XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DU VOLUME
      NBSOVL = MCN( MNSOVL + WNBSOM )
C     LE NOMBRE DE TANGENTES DU VOLUME
      NBTGVL = MCN( MNSOVL + WNBTGS )
C
C     LE NUMERO DE LA FONCTION
      NUFOCS = LADEFI( WUFOCS )
      IF( NUFOCS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FONCTION CRITERE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN CRITERION FUNCTION'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     LE TYPE DES FACES A EXTRAIRE
      NTYFSE = LADEFI( WTYFSE )
      IF( NTYFSE .LE. 0 .OR. NTYFSE .GT. 3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TYPE INCONNU DE FACES A EXTRAIRE'
         ELSE
            KERR(1) = 'UNKNOWN TYPE FOR FACES TO BE EXTRACTED'
         ENDIF
         CALL LEREUR
         IERR = 5
         RETURN
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNCUVL),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE DE FACES, SOMMETS, ... DE LA SURFACE EXTRAITE
      NBFAEX = 0
      NBFTEX = 0
      NBSOEX = 0
      NBTGEX = 0
      MNOUIS = 0
      MNFAEX = 0
      MNNEW  = 0
      MNNEWT = 0
      MNFAIN = 0
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU CRITERE
      CALL TNMCDC( 'ENTIER', NBSOVL, MNOUIS )
      MNO = MNOUIS - 1
C
C     CALCUL DU CRITERE EN CHACUN DES SOMMETS DU MAILLAGE
C     ===================================================
      MNS  = MNSOVL + WYZSOM - 3
      DO 10 N=1,NBSOVL
C
C        LES 3 COORDONNEES EN DOUBLE PRECISION DU SOMMET N
         MN = MNS + 3 * N
         DXYZ(1) = RMCN( MN   )
         DXYZ(2) = RMCN( MN+1 )
         DXYZ(3) = RMCN( MN+2 )
C
C        LE SOMMET N VERIFIE T IL LE CRITERE ?
         CALL FONVAL( NUFOCS, 3, DXYZ, NCODEV, DOUI )
         IF( NCODEV .EQ. 0 ) DOUI = 0
C
C        1 SI CRITERE VERIFIE, 0 SINON
         I = NINT( DOUI )
         MCN( MNO+N ) = I
         IF( I .EQ. 1 ) NBSOEX = NBSOEX + 1
 10   CONTINUE
C
C     EXISTE T IL DES SOMMETS VERIFIANT LE CRITERE ?
      IF( NBSOEX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AUCUN SOMMET NE VERIFIE LE CRITERE'
         ELSE
            KERR(1) = 'NO VERTEX VERIFIES THE CRITERION'
         ENDIF
         CALL LEREUR
         IERR = 10
         GOTO 9999
      ENDIF
C
C     GENERATION EVENTUELLE PAR HACHAGE DES FACES DES CUBES
C     CHAINAGE DES FACES FRONTALIERES EN POSITION 7
C     AVEC UN LIEN POUR LES FACES FRONTALIERES
C     =====================================================
      CALL HAFAVO( NMOBJT, 3, NTFAVO, MNFAVO, N, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME: ' // NMOBJT
            KERR(2) = 'IMPOSSIBLE DE CREER SES FACES'
         ELSE
            KERR(1) = 'VOLUME: ' // NMOBJT
            KERR(2) = 'IMPOSSIBLE TO CREATE ITS FACES'
         ENDIF
         CALL LEREUR
         IERR    = 6
         GOTO 9999
      ENDIF
C
C     LE NOMBRE D'ENTIERS PAR FACE
      MOFACE = MCN( MNFAVO + WOFACE )
C     LA MAJORATION DU NOMBRE DE FACES
      MXFACE = MCN( MNFAVO + WXFACE )
C     LE NOMBRE DE FACES FRONTALIERES
      NBFAFR = MCN( MNFAVO + WBFAFR )
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      L1FAFR = MCN( MNFAVO + W1FAFR )
C     LE NOMBRE DE FACES AVEC DES TANGENTES
      NBFATG = MCN( MNFAVO + WBFATG )
C     ADRESSE DANS MCN POUR ATTEINDRE LES FACES
      MNF0   = MNFAVO + WFACES - MOFACE - 1
C     ADRESSE DANS MCN POUR ATTEINDRE LES 8 NUMEROS DE TG DES FACES A TG
      MNTGFA = MNFAVO + WFACES + MOFACE * MXFACE
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU TEMPORAIRE
C     DES NUMEROS DES SOMMETS DES FACES EXTRAITES
      MXFAEX = 5 * MXFACE
      CALL TNMCDC( 'ENTIER', MXFAEX, MNFAEX )
      MNFAE  = MNFAEX - 1
      MNFTEX = MNFAEX - 1 + 4 * MXFACE
      CALL AZEROI( MXFACE, MCN(MNFTEX) )
C
C     RECENSEMENT DES FACES A TRAITER
C     ===============================
      MNFAI1 = 0
      IF( NTYFSE .EQ. 1 .OR. NTYFSE .EQ. 2 ) THEN
         CALL TNMCDC( 'ENTIER', MXFACE, MNFAIN )
         MNFAI1 = MNFAIN - 1
         DO 20 N=1,MXFACE
            MCN(MNFAI1+N) = N
 20      CONTINUE
C        MISE A ZERO DES FACES NON A TRAITER
C        LA PREMIERE FACE FRONTALIERE
         NF   = L1FAFR
         DO 30 N = NBFAFR, 1, -1
            MCN(MNFAI1+NF) = -NF
C           LA FACE FRONTALIERE SUIVANTE
            NF = MCN( MNF0 + MOFACE * NF +7 )
 30      CONTINUE
C
C        ICI MCN(MNFAI1+N) >0 SI FACE INTERNE OU NON INITIALISEE
C                          <0 SI FACE FRONTALIERE
C
      ENDIF
C
C     LA BOUCLE SUR LES FACES DU VOLUME
C     =================================
      DO 290 N = 1, MXFACE
C
C        ADRESSE MCN DU TABLEAU DE LA FACE
         MNF = MNF0 + MOFACE * N
C
C        LA FACE EST-ELLE INITIALISEE ?
         IF( MCN( MNF+1 ) .EQ. 0 ) GOTO 290
C
C        OUI LES FACES A EXTRAIRE SONT ELLES QUELCONQUES ?
         IF( NTYFSE .EQ. 3 ) GOTO 200
C
C        TYPE DE LA FACE N : INTERNE OU FRONTALIERE ?
         NTY = MCN( MNFAI1+N )
         IF( NTY .GT. 0 .AND. NTYFSE .EQ. 2 ) GOTO 200
         IF( NTY .LT. 0 .AND. NTYFSE .EQ. 1 ) GOTO 200
         GOTO 290
C
C        LA FACE N EXISTE ET A TRAITER
C        -----------------------------
C        LE NOMBRE DE SOMMETS DE LA FACE
 200     IF( MCN( MNF+4 ) .GT. 0 ) THEN
            NBSF = 4
         ELSE
            NBSF = 3
         ENDIF
C
         OUI  = .TRUE.
         DO 220 J=1,NBSF
C           LE NUMERO DU SOMMET J DE LA FACE
            NSOM = MCN( MNF+J )
C           LE SOMMET VERIFIE T IL LE CRITERE ?
            OUI = OUI .AND. ( MCN(MNO+NSOM).EQ. 1 )
 220     CONTINUE
C
         IF( OUI ) THEN
C           LA FACE VERIFIE LE CRITERE
            NBFAEX = NBFAEX + 1
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO 230 I=1,NBSF
C              LE NUMERO DU SOMMET I DE LA FACE
               MCN( MNFAE + I ) = MCN( MNF+I )
 230        CONTINUE
C           COMPLETION EVENTUELLE AVEC DES ZEROS
            DO 240 I=NBSF+1,4
               MCN( MNFAE + I ) = 0
 240        CONTINUE
            MNFAE = MNFAE + 4
C
            NUFATG = MCN( MNF + 8 )
            IF( NUFATG .GT. 0 ) THEN
C              UNE FACE A TG DE PLUS
               NBFTEX = NBFTEX + 1
C              LE NUMERO DE FACE A EXTRAIRE DU VOLUME
               MCN( MNFTEX + N ) = N
            ENDIF
         ENDIF
 290  CONTINUE
C
C     EXISTE T IL DES FACES VERIFIANT LE CRITERE ?
      IF( NBFAEX .LE. 0 ) THEN
         NBSOEX = 0
         GOTO 460
      ENDIF
C
C     RECENSEMENT DES SOMMETS DE CETTE SURFACE EXTRAITE
C     =================================================
      CALL TNMCDC( 'ENTIER', NBSOVL+1, MNNEW )
      CALL AZEROI( NBSOVL+1, MCN(MNNEW) )
C
C     LE NOUVEAU NUMERO DES SOMMETS RECENSES
      DO 400 I=MNFAEX,MNFAE
C        LE NUMERO DU SOMMET EST RECENSE
         N = MCN(I)
C        ET REMPLACE LE ZERO INITIAL
         MCN(MNNEW+N) = N
 400  CONTINUE
C
      NBSOEX = 0
      DO 450 I=1,NBSOVL
         IF( MCN(MNNEW+I) .NE. 0 ) THEN
            NBSOEX = NBSOEX + 1
C           LE NOUVEAU NUMERO DU SOMMET EXTRAIT
            MCN(MNNEW+I) = NBSOEX
         ENDIF
 450  CONTINUE
C
 460  IF( NBSOEX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SELON CE CRITERE AUCUNE FACE EXTRAITE'
         ELSE
            KERR(1) = 'FOR THIS CRITERION NO EXTRACTED FACE'
         ENDIF
         CALL LEREUR
         IERR = 10
         GOTO 9999
      ENDIF
C
C     CREATION DE LA SURFACE EXTRAITE
C     ===============================
C
C     CONSTRUCTION DU TABLEAU 'NSEF'
C     ------------------------------
      IF( NBFTEX .GT. 0 ) THEN
         NBTGEF = 8
      ELSE
         NBTGEF = 0
      ENDIF
      CALL LXTNDC( NTLXSU, 'NSEF', 'MOTS',
     %             WUSOEF+4*NBFAEX+NBFAEX+NBFTEX*9)
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASE , MNFASE )
C
C     LE TYPE DE L'OBJET : ICI SURFACE
      MCN( MNFASE + WUTYOB ) = 3
C     LE TYPE DU MAILLAGE : ICI SURFACE NON STRUCTUREE
      MCN( MNFASE + WUTYMA ) = 0
C     LE NOMBRE DE FACES DE LA SURFACE
      MCN( MNFASE + WBEFOB ) = NBFAEX
C     LE NOMBRE DE SOMMETS DES FACES
      MCN( MNFASE + WBSOEF ) = 4
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNFASE + WUTFMA ) = -1
C     LES TANGENTES STOCKEES
      MCN( MNFASE + WBTGEF ) = NBTGEF
      MCN( MNFASE + WBEFAP ) = NBFTEX
      MCN( MNFASE + WBEFTG ) = NBFTEX
C
C     LE CHANGEMENT DE NUMERO DES SOMMETS DE LA SURFACE EXTRAITE
      MN1 = MNFASE + WUSOEF
      DO 700 I=MNFAEX,MNFAE
C        L'ANCIEN NUMERO DU SOMMET
         N = MCN(I)
C        LE NOUVEAU NUMERO DU SOMMET
         MCN(MN1) = MCN(MNNEW+N)
         MN1 = MN1 + 1
 700  CONTINUE
C
C     LES EVENTUELLES FACES A TG
C     ==========================
      IF( NBFTEX .GT. 0 ) THEN
         NBFTEX = 0
         MNEFAP = MNFASE + WUSOEF + 4*NBFAEX
         MNCGEF = MNEFAP + NBFAEX
         MNNUT0 = MNCGEF + NBFTEX
         MNNUTG = MNNUT0 - 1
         DO 720 N=1,MXFACE
            IF( MCN(MNFTEX+N) .GT. 0 ) THEN
C              LA FACE EXTRAITE EST ELLE A TG?
               MNF    = MNF0 + MOFACE * N
               NUFATG = MCN( MNF + 8 )
               IF( NUFATG .GT. 0 ) THEN
C                 FACE A TG
                  NBFTEX = NBFTEX + 1
                  MCN( MNEFAP ) = NBFTEX
C                 CODE GEOMETRIQUE
                  MCN( MNCGEF ) = 0
                  MNCGEF = MNCGEF + 1
C                 LES 8 NUMEROS DES TANGENTES DE LA FACE
                  MN = MNTGFA - 9 + NUFATG * 8
                  DO 710 I=1,8
C                    LE NUMERO DE LA TG DANS LE VOLUME
                     MCN( MNNUTG + I ) = MCN( MN + I )
 710              CONTINUE
                  MNNUTG = MNNUTG + 8
               ELSE
C                 FACE SANS TG
                  MCN( MNEFAP ) = 0
               ENDIF
               MNEFAP = MNEFAP + 1
            ENDIF
 720     CONTINUE
C
C        RECENSEMENT DES TANGENTES DE CETTE SURFACE EXTRAITE
C        ===================================================
         CALL TNMCDC( 'ENTIER', NBTGVL+1, MNNEWT )
         CALL AZEROI( NBTGVL+1, MCN(MNNEWT) )
C
C        LE NOUVEAU NUMERO DES SOMMETS RECENSES
         DO 730 I=MNNUT0,MNNUTG
C           LE NUMERO DE LA TG EST RECENSEE
            N = MCN(I)
C           ET REMPLACE LE ZERO INITIAL
            MCN(MNNEWT+N) = N
 730     CONTINUE
C
         DO 750 I=1,NBTGVL
            IF( MCN(MNNEWT+I) .NE. 0 ) THEN
               NBTGEX = NBTGEX + 1
C              LE NOUVEAU NUMERO DE LA TANGENTE EXTRAITE
               MCN(MNNEWT+I) = NBTGEX
            ENDIF
 750     CONTINUE
C
         DO 760 I=MNNUT0,MNNUTG
C           L'ANCIEN NUMERO DE TG
            N = MCN(I)
C           LE NOUVEAU NUMERO DE TG
            MCN(I) = MCN(MNNEWT+N)
 760     CONTINUE
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFASE) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASE + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS'  ,
     %             WYZSOM+3*(NBSOEX+NBTGEX) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOSE , MNSOSE )
C
C     LE NOMBRE DE SOMMETS DE LA SURFACE
      MCN( MNSOSE + WNBSOM ) = NBSOEX
C     LE NOMBRE DE TANGENTES DE LA SURFACE EXTRAITE
      MCN( MNSOSE + WNBTGS ) = NBTGEX
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOSE + WBCOOR ) = 3
C     LES 3 COORDONNEES DES SOMMETS DE LA SURFACE
      MN0 = MNSOVL + WYZSOM - 3
      MN1 = MNSOSE + WYZSOM
      DO 800 I=1,NBSOVL
         N = MCN(MNNEW+I)
         IF( N .NE. 0 ) THEN
            MN  = MN0 + 3 * I
            RMCN(MN1  ) = RMCN(MN  )
            RMCN(MN1+1) = RMCN(MN+1)
            RMCN(MN1+2) = RMCN(MN+2)
            MN1 = MN1 + 3
         ENDIF
 800  CONTINUE
C
C     LES 3 COORDONNEES DES TANGENTES DE LA SURFACE
      IF( NBTGEX .GT. 0 ) THEN
         MN0 = MN0 + 3 * NBSOVL
         DO 810 I=1,NBTGVL
            N = MCN(MNNEWT+I)
            IF( N .NE. 0 ) THEN
               MN  = MN0 + 3 * I
               RMCN(MN1  ) = RMCN(MN  )
               RMCN(MN1+1) = RMCN(MN+1)
               RMCN(MN1+2) = RMCN(MN+2)
               MN1 = MN1 + 3
            ENDIF
 810     CONTINUE
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOSE) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOSE + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
 9999 IF( MNFAEX .GT. 0 ) CALL TNMCDS( 'ENTIER', MXFAEX  , MNFAEX )
      IF( MNNEW  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOVL+1, MNNEW  )
      IF( MNNEWT .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTGVL+1, MNNEWT )
      IF( MNOUIS .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOVL  , MNOUIS )
      IF( MNFAIN .GT. 0 ) CALL TNMCDS( 'ENTIER', MXFACE  , MNFAIN )
      RETURN
      END
