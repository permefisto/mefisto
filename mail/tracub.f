      SUBROUTINE TRACUB( NMVOLU, NUVOLU, MNTSMA, MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES CUBES DU VOLUME NUVOLU
C -----
C ENTREES:
C --------
C NMVOLU : NOM DU VOLUME A TRACER
C NUVOLU : NUMERO DU VOLUME DANS SON LEXIQUE
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET' A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        MARS 1991
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___trace.inc"
      include"./incl/a___face.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      INTEGER, allocatable, dimension(:)   :: NOFAFR, NOFAFF
      REAL,    allocatable, dimension(:)   :: DIOEIL
      REAL,    allocatable, dimension(:,:) :: BARFAC

      CHARACTER*(*)     NMVOLU
      REAL              XYZ(3)

      IERNOFA = 1
      IERNOFF = 1
      IERDIOE = 1
      IERBARF = 1

C     COULEUR ET EPAISSEUR DES ARETES DES FACES
      CALL XVCOULEUR( NCOUAF )
      CALL XVEPAISSEUR( NEPARF )

      IF( IAVFAC .EQ. 0 .OR. MCN(MNSOMM+WBCOOR) .EQ. 6 ) THEN

C        TRACE DES ARETES DES 6-CUBES EN PROJECTION XYZ (OUBLI DE UVW)
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
         GOTO 5

C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
 3       CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) RETURN

C        TRACE DES ARETES DES 3-CUBES DES 6-CUBES
 5       CALL T3ARVO( NMVOLU )
C        TRACE DES AXES DU CUBE ENGLOBANT
         CALL TRAXE3
C        TRACE DU TITRE ET FERMETURE
         CALL TRFINS( NMVOLU )

C        REPRISE DE L'ORBITE
         IF( LORBITE .NE. 0 ) GOTO 3
C        SORTIE DU TRACE
         RETURN
      ENDIF

C     TRACE SI LA DIMENSION DE L'ESPACE VAUT 3
      IF( MCN(MNSOMM+WBCOOR) .NE. 3 ) RETURN

C     GENERATION EVENTUELLE PAR HACHAGE DES FACES DES CUBES
C     CHAINAGE DES FACES FRONTALIERES EN POSITION 7
C     AVEC UN LIEN NEGATIF POUR LES FACES FRONTALIERES
C     =====================================================
      CALL HAFAVO( NMVOLU, 3, NTFAVO, MNFAVO, MNS, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'VOLUME: ' // NMVOLU
         KERR(2) = 'IMPOSSIBLE DE CREER SES FACES'
         CALL LEREUR
         RETURN
      ENDIF

C     LE NOMBRE D'ENTIERS PAR FACE
      MOFACE = MCN( MNFAVO + WOFACE )
C     LA MAJORATION DU NOMBRE DE FACES
      MXFACE = MCN( MNFAVO + WXFACE )
C     LE NOMBRE DE FACES FRONTALIERES
      NBFAFR = MCN( MNFAVO + WBFAFR )
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      L1FAFR = MCN( MNFAVO + W1FAFR )
C     LE NOMBRE DE FACES INTERFACES
      NBFA2M = MCN( MNFAVO + WBFA2M )
C     LE NUMERO DE LA PREMIERE FACE INTERFACE
      L1FA2M = MCN( MNFAVO + W1FA2M )

C     RESERVATION D'UN TABLEAU POUR LE TRI
      NBFAFI = NBFAFR + NBFA2M
      ALLOCATE ( NOFAFR(1:NBFAFI), STAT=IERNOFA )
      ALLOCATE ( NOFAFF(1:NBFAFI), STAT=IERNOFF )
      IF( IERNOFA .NE. 0 .OR. IERNOFF .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tracub: NOFAFR ou NOFAFF de',NBFAFI,
     %             ' entiers NON ALLOUES'
         ELSE
            PRINT*,'tracub: NOFAFR or NOFAFF of',NBFAFI,
     %             ' INTEGERS are NOT ALLOCATED'
         ENDIF
         GOTO 9000
      ENDIF

      ALLOCATE ( DIOEIL(1:NBFAFI), STAT=IERDIOE )
      IF( IERDIOE .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tracub: DIOEIL de',NBFAFI,
     %             ' reels NON ALLOUES'
         ELSE
            PRINT*,'tracub: DIOEIL of',NBFAFI,
     %             ' REAL are NOT ALLOCATED'
         ENDIF
         GOTO 9000
      ENDIF

      ALLOCATE ( BARFAC(1:3,1:NBFAFI), STAT=IERBARF )
      IF( IERBARF .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tracub: BARFAC de',3*NBFAFI,
     %             ' reels NON ALLOUES'
         ELSE
            PRINT*,'tracub: BARFAC of',3*NBFAFI,
     %             ' REAL are NOT ALLOCATED'
         ENDIF
         GOTO 9000
      ENDIF

C     ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS = MNSOMM + WYZSOM

C     BARYCENTRE DES FACES FRONTALIERES
C     =================================
      CALL BAFAFR( MOFACE, MXFACE, MCN(MNFAVO+WFACES), RMCN(MNS),
     %             0,      L1FAFR, NBFAFR, L1FA2M, NBFA2M,
     %             NOFAFR, BARFAC )

C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
C     ===========================================
      IF( LORBITE .EQ. 0 ) GOTO 19
      CALL ORBITE0( NOTYEV )
      GOTO 19

C     ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ====================================
 12   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000

C     TRACE DES AXES
      CALL TRAXE3

C     CALCUL DE LA DISTANCE DU BARYCENTRE DES FACES A L'OEIL
C     ======================================================
 19   DO NF = 1, NBFAFI
C        COTE Z AXONOMETRIQUE DANS LA DIRECTION DE VISEE
         CALL XYZAXO( BARFAC(1,NF), XYZ )
         DIOEIL(NF) = XYZ(3)
         NOFAFF(NF) = NOFAFR(NF)
      ENDDO

C     LE TRI PAR TAS DE CETTE DISTANCE
C     ================================
C     LA FACE LA PLUS PROCHE EST LA PREMIERE
      CALL TRITRP( NBFAFI, DIOEIL, NOFAFF )
C     ICI NOFAFF = NUMERO DE LA FACE LA PLUS PROCHE DE L'OEIL

      IF( NBRCOU .EQ. 0 ) THEN
C
C        PEAU DES FACES EN NOIR ET BLANC DE LA FRONTIERE DU VOLUME
C        ----------------------------)-----------------------------
         CALL T3FPNB( NMVOLU, NUVOLU ,
     %                MOFACE, MXFACE, MCN(MNFAVO+WFACES), RMCN(MNS),
     %                NBFAFI, NOFAFF )

      ELSE IF( IAVFAC .EQ. 0 ) THEN
C
C        LA PEAU DES FACES EN FIL DE FER DE LA FRONTIERE DU VOLUME
C        ---------------------------------------------------------
         CALL T3TOUR( NMVOLU, NUVOLU ,
     %                MOFACE, MXFACE, MCN(MNFAVO+WFACES), RMCN(MNS),
     %                NBFAFI, NOFAFF )

      ELSE
C
C        LES FACES FRONTALIERES ALGORITHME DU PEINTRE
C        --------------------------------------------
         CALL T3FACO( NMVOLU, NUVOLU ,
     %                MOFACE, MXFACE, MCN(MNFAVO+WFACES),
     %                MCN(MNSOMM+WNBSOM), RMCN(MNS),
     %                NBFAFI, NOFAFF, DIOEIL, MNTSMA )

      ENDIF
      CALL TRFINS( NMVOLU )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 12

C     DESTRUCTION DES TABLEAUX SERVANT AU TRI
 9000 IF( IERNOFA .EQ. 0 ) DEALLOCATE( NOFAFR )
      IF( IERNOFF .EQ. 0 ) DEALLOCATE( NOFAFF )
      IF( IERDIOE .EQ. 0 ) DEALLOCATE( DIOEIL )
      IF( IERBARF .EQ. 0 ) DEALLOCATE( BARFAC )

      RETURN
      END
