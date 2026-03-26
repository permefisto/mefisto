      SUBROUTINE VOEX26( NUVOLU, LADEFI, RADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : MISE A JOUR de la LONGUEUR des ARETES des TETRAEDRES SELON
C ----- SOIT la LONGUEUR>0 des ARETES par DEFAUT (code 0 du Menu DEBUT)
C       SOIT la FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z) DONNEE ou 
C            the USER''s Function EDGE_LENGTH(x,y,z) GIVEN

C ENTREES:
C --------
C NUVOLU   : NUMERO DU VOLUME A TRAITER DANS LE LEXIQUE DES VOLUMES
C LADEFI   : TABLEAU DE DEFINITION Vu de TYPE ENTIER DU VOLUME A TRAITER
C RADEFI   : TABLEAU DE DEFINITION Vu de TYPE REEL   DU VOLUME A TRAITER
C            CF '~/td/d/a_volume__definition'

C SORTIES:
C --------
C NTNSEF : NUMERO du TMS DES NUMEROS DES TETRAEDRES DU VOLUME
C MNNSEF : ADRESSE MCN DU TABLEAU DU NO des SOMMETS des EF du VOLUME
C          CF '~/td/d/a___nsef'
C NTXYZS : NUMERO du TABLEAU TMS DES XYZ des SOMMETS de la TETRAEDRISATION
C MNXYZS : ADRESSE MCN du TABLEAU XYZSOMMET de la TETRAEDRISATION
C          CF '~/td/d/a___xyzsommet'
C IERR   : =0 SI PAS D'ERREUR,  >0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE DU PERRAY              Juin 2019
C2345X7..............................................................012
C     QTEAME: QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C             AMELIORATION DE LA QUALITE DES TETRAEDRES EST DEMANDEE
ccc      PARAMETER        (QTEAME=0.001)  9/9/17
ccc      PARAMETER        (QTEAME=0.01 ) 13/10/17
ccc      PARAMETER        (QTEAME=0.02 )  3/12/18
ccc      PARAMETER        (QTEAME=0.05)
      PARAMETER        (QTEAME=0.08)

C     QUAMINEX: QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C               C-A-D AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE
      PARAMETER        (QUAMINEX=0.005)

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ

      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/a___face.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/darete.inc"

      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))

      DOUBLE PRECISION, allocatable, dimension(:,:) :: PTXYZD
      INTEGER, allocatable, dimension(:,:) ::  NOTETR
      INTEGER, allocatable, dimension(:)   ::  NVOLTE
      INTRINSIC         ALLOCATED

C     FORMULES : TETRAEDRES + ARETES + NB_FRONTIERES = FACES + SOMMETS
C     ========== 4*TETRAEDRES = FACES_FRONTALIERES + 2*FACES_INTERNES
C                ASYMPTOTIQUEMENT  FACES = 2 * TETRAEDRES

      CHARACTER*24      KNMVOLU, KNMVOLUIN
      CHARACTER*80      KINFO
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      INTEGER           NUFACE(16), NOSOCU(1:8), NOSOTR(3),
     %                  LEFACO(1), N1FASC(1)
      EQUIVALENCE     ( NUFACE(1), NOSOCU(1) )
      INTEGER           NBIPAV(3)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
      DOUBLE PRECISION  DINFO, VOLUMT, VOLMOY, CPU, CPUT
ccc      REAL              QUAMIN, QUAMOY

C     INITIALISATION DU TEMPS CPU
C     ---------------------------
      CPU  = DINFO( 'DELTA CPU' )
      CPUT = 0D0
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE
      STOPTE = .FALSE.
      TRACTE = .FALSE.
ccc      TRACTE = .TRUE.

      NTNSEF = 0
      MNNSEF = 0
      NTXYZS = 0
      MNXYZS = 0

C     LE TABLEAU LEFACO N'EXISTE PAS
      INFACO = 0
      MXFACO = 0

      MNOPTSUIV=0
      MN1SPAVE=0
      MNLIPO = 0
      MNTETO = 0
      MNSOFR = 0
      MNTETS = 0
      MNNEWS = 0

      IERAOALLOC = 1
      IERATALLOC = 1
      IERZDALLOC = 1
      IERTEALLOC = 1
      IERNVALLOC = 1

C     RECUPERATION DU NOM DU VOLUME DE NUMERO NUVOLU
C     DANS LE LEXIQUE DES VOLUMES
      CALL NMOBNU( 'VOLUME', NUVOLU, KNMVOLU )

C     OUVERTURE DU MMS LEXIQUE DU VOLUME A TRAITER
C     NTLXVOLU : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME A TRAITER
      CALL LXNLOU( NTVOLU, NUVOLU, NTLXVOLU, MNLXVOLU )

C     NUMERO DU VOLUME INITIAL TETRAEDRISE DE LONGUEUR D'ARETE A TRAITER
C     variable NUVOIN 'nom du volume de longueur d''arete a ameliorer' ^~>VOLUME
C     --------------------------------------------------------------------------
      NUVOIN = LADEFI( WUVOIN )
      IF( NUVOIN .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'voex26: VOLUME a TRAITER INCONNU ',KNMVOLU
         ELSE
            PRINT*,'voex26: UNKNOWN VOLUME to TREAT ',KNMVOLU
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

      CALL NMOBNU( 'VOLUME', NUVOIN, KNMVOLUIN )
      IF( LANGAG .EQ. 0 ) THEN
       PRINT*,'voex26: MISE a JOUR de la LONGUEUR deS ARETES du VOLUME '
     %       ,KNMVOLUIN
      ELSE
         PRINT*,'voex26: UPDATE of EDGE LENGTH of VOLUME TETRAHEDRA '
     %         ,KNMVOLUIN
      ENDIF

C     LE TABLEAU LEXIQUE DE CE VOLUME INITIAL
      CALL LXLXOU( NTVOLU, KNMVOLUIN, NTLXVI, MNLXVI )

C     variable NBPTIM 'nombre points internes imposes par l''utilisateur' entier
C     --------------------------------------------------------------------------
      NBPTIM = LADEFI( WBPTIM )

C     tableau XYZDIM(1..4,1..NBPTIM) 'XYZ et distance souhaitee aux voisins' reel
C     RADEFI(WYZDIM:WYZDIM+NBPTIM-1)
C     ---------------------------------------------------------------------------

C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE'
C     -----------------------------------------------
      ARETGR = ABS( RADEFI( WRETGR ) )
      DARETE = ARETGR
      NOFOTI = NOFOTIEL()
      IF( LANGAG .EQ. 0 ) THEN
         IF( NOFOTI .GT. 0 ) THEN
C           LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
            WRITE(IMPRIM,10000) 'AVEC'
         ELSE
C           LA FONCTION TAILLE_IDEALE(X,Y,Z) N'EXISTE PAS
            WRITE(IMPRIM,10000) 'SANS'
            WRITE(IMPRIM,10001) DARETE
         ENDIF
      ELSE
         IF( NOFOTI .GT. 0 ) THEN
C           LA FONCTION EDGE_LENGTH(X,Y,Z) EXISTE
            WRITE(IMPRIM,20000) '  WITH '
         ELSE
C           LA FONCTION EDGE_LENGTH(X,Y,Z) N'EXISTE PAS
            WRITE(IMPRIM,20000) 'WITHOUT'
            WRITE(IMPRIM,20001) DARETE
         ENDIF
      ENDIF
10000 FORMAT(' voex26: TETRAEDRISATION ',A,
     %     ' la FONCTION TAILLE_IDEALE(X,Y,Z) (ou EDGE_LENGTH(X,Y,Z))' )
20000 FORMAT(' voex26: TETRAHEDRISATION ',A,
     %    ' the FUNCTION EDGE_LENGTH(X,Y,Z) (or TAILLE_IDEALE(X,Y,Z))' )
10001 FORMAT(' LONGUEUR de l''ARETE SOUHAITEE =',G14.6)
20001 FORMAT(' WISHED EDGE LENGTH =',G14.6)

C     RECUPERATION DES TABLEAUX DE LA TETRAEDRISATION INITIALE NUVOIN
C     ---------------------------------------------------------------
      CALL LXNLOU( NTVOLU, NUVOIN, NTLXQU, MNLXQU )
      IF( NTLXQU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TETRAEDRISATION INITIALE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL TETRAHEDRIZATION'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

C     LE TABLEAU 'XYZSOMMET'
      CALL LXTSOU( NTLXQU, 'XYZSOMMET', NTXYZSVI, MNXYZSVI )
      IF( NTXYZSVI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS le TMS XYZSOMMET'
         ELSE
            KERR(1) = 'VOLUME WITHOUT the XYZSOMMET TMS'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF

C     LE NOMBRE DE SOMMETS ET DE COORDONNEES DE LA TETRAEDRISATION INITIALE
      NBSOMMVI = MCN( MNXYZSVI + WNBSOM )
      NBCOOR   = MCN( MNXYZSVI + WBCOOR )

C     LE TABLEAU 'NSEF'
      CALL LXTSOU( NTLXQU, 'NSEF', NTNSEFVI, MNNSEFVI )
      IF( NTNSEFVI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS le TMS NSEF'
         ELSE
            KERR(1) = 'VOLUME WITHOUT the NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF

C     LES CARACTERISTIQUES DES NSEF DE CETTE TETRAEDRISATION
      CALL NSEFPA( MCN(MNNSEFVI) ,
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBTETRVI,
     %             NX,     NY,     NZ,
     %             IERR   )
C     NUTYMA : 'NUMERO DE TYPE DU MAILLAGE'    ENTIER
C              0 : 'NON STRUCTURE'       , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TETRAEDRE  STRUCTURE', 4 : 'TETRAEDRE  STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE' , 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
C     NBSOEL : NOMBRE DE SOMMETS DES EF
C              0 SI MAILLAGE NON STRUCTURE
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DES NSEF
C              ( TETRAEDRE NBSOEF=8 )
C     NBTETRVI : NOMBRE DE EF DU MAILLAGE
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTIONS X Y Z
C                CF LE TMS ~td/d/a___nsef

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'TYPE de STRUCTURE du MAILLAGE =',NUTYMA
         PRINT*, 'NOMBRE INITIAL de SOMMETS     =',NBSOMMVI
         PRINT*, 'NOMBRE INITIAL de TETRAEDRES  =',NBTETRVI
      ELSE
         PRINT*, 'MESH STRUCTURE TYPE          =',NUTYMA
         PRINT*, 'INITIAL NUMBER of VERTICES   =',NBSOMMVI
         PRINT*, 'INITIAL NUMBER of TETRAHEDRA =',NBTETRVI
      ENDIF

      IF( NBSOEF .NE. 8 .OR. NUTYMA .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MAILLAGE DIFFERENT d''une TETRAEDRISATION'
         ELSE
            KERR(1) = 'MESH is NOT a TETRAHEDRIZATION'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF

C     VERIFICATION PAS D'EF NON TETRAEDRES DANS LE MAILLAGE
C     -----------------------------------------------------
      MNTE = MNNSEFVI + WUSOEF -1
      DO NUELEM = 1, NBTETRVI

C        L'EF NUELEM A T IL UN 5-EME SOMMET?
         IF( MCN( MNTE + 5 ) .GT. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PRESENCE INTERDITE d''UN EF NON TETRAEDRE'
            ELSE
               KERR(1) = 'FORBIDDEN PRESENCE of a NOT TETRAHEDRON FE'
            ENDIF
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF

         MNTE = MNTE + NBSOEF

      ENDDO

C     VISEE AVEC LES COORDONNEES COOEXT DE L'HEXAEDRE ENGLOBANT
C     ---------------------------------------------------------
      INIEXT = 0
      CALL MAJEXT( MNXYZSVI )

C     AXONOMETRIE
      DISTMX=MAX( COOEXT(1,2) - COOEXT(1,1),
     %            COOEXT(2,2) - COOEXT(2,1),
     %            COOEXT(3,2) - COOEXT(3,1) )
      AXOAVA = 0
      AXOARR = 0
      AXOLAR = DISTMX * 0.6
      AXOHAU = DISTMX * 0.6
      AXOPTV(1) = ( COOEXT(1,1) + COOEXT(1,2) ) / 2
      AXOPTV(2) = ( COOEXT(2,1) + COOEXT(2,2) ) / 2
      AXOPTV(3) = ( COOEXT(3,1) + COOEXT(3,2) ) / 2
      AXOEIL(1) =   COOEXT(1,2)
      AXOEIL(2) =   COOEXT(2,2)
      AXOEIL(3) =   COOEXT(3,2)
      CALL AXONOMETRIE( AXOPTV,AXOEIL, AXOLAR,AXOHAU, AXOARR,AXOAVA )

      IF( TRACTE ) THEN
C        TRACE DES EF DU VOLUME NUVOIN
         CALL EFFACE
         LCRITR = 0
         PREDUF = 25.
         CALL TRACUB( KNMVOLUIN, NUVOIN, MNNSEFVI, MNXYZSVI )

C        TRACE DES ARETES DES TETRAEDRES DEFINIS DANS XYZSOMMET et NSEF
         CALL EFFACE
         MXTET = NBTETRVI
         CALL TRNSTETR( NBSOEF, MXTET, MCN(MNNSEFVI+WUSOEF),
     %                  NBCOOR, RMCN(MNXYZSVI+WYZSOM),
     %                  NBTETR, VOLET0, QUAMIN0, QUAMOY0 )
      ENDIF

C     MAJORATION DU NOMBRE DE SOMMETS DES TETRAEDRES
C     MAJORATION DU NOMBRE DES TETRAEDRES
C     ----------------------------------------------
      print*,'voex26: VOLUME: ',KNMVOLU
C     LE NOMBRE MAXIMAL DE SOMMETS DECLARABLES
      print*,'voex26: Nombre de Sommets     =',NBSOMMVI
      MXSOMM = 64 * NBSOMMVI
C     MAJORATION DU NOMBRE DE TETRAEDRES SELON LE MAXIMUM DE SOMMETS
      print*,'voex26: Nombre de Tetraedres  =',NBTETRVI
      MXTETR = 5 * MXSOMM

C     IMPRESSIONS DE CES MAJORATIONS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10010) MXSOMM, MXTETR
      ELSE
         WRITE(IMPRIM,20010) MXSOMM, MXTETR
      ENDIF
10010 FORMAT(/
     %' NOMBRE MAXIMUM DECLARE de SOMMETS    =',I9/
     %' NOMBRE MAXIMUM DECLARE de TETRAEDRES =',I9)
20010 FORMAT(/
     %' DECLARED VERTEX     MAXIMUM NUMBER =',I9/
     %' DECLARED TETRAHEDRA MAXIMUM NUMBER =',I9)

C     ALLOCATION FORTRAN DYNAMIQUE du tableau PTXYZD(4,MXSOMM)
C     --------------------------------------------------------
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'DEMANDE ALLOCATION de',4*MXSOMM,
     %                   ' DOUBLE PRECISION pour PtXYZd'
         ALLOCATE ( PTXYZD( 1:4, 1:MXSOMM ), STAT=IERZDALLOC )
         IF( IERZDALLOC .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION des',4*MXSOMM,
     %                     ' DOUBLE PRECISION pour PtXYZd'
            IERR = IERZDALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'CORRECTE ALLOCATION de PtXYZD(1:4,1:',MXSOMM,
     %                ') DOUBLE PRECISION'
      ELSE
         PRINT*, 'ALLOCATION DEMAND  of',4*MXSOMM,
     %                ' DOUBLE PRECISION of PtXYZd'
         ALLOCATE ( PTXYZD( 1:4, 1:MXSOMM ), STAT=IERZDALLOC )
         IF( IERZDALLOC .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',4*MXSOMM,
     %                ' DOUBLE PRECISION of PtXYZd'
            IERR = IERZDALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'ALLOCATION CORRECT of PtXYZD(1:4,1:',MXSOMM,
     %                ') DOUBLE PRECISION'
      ENDIF

C     ALLOCATION FORTRAN de NOTETR(8,MXTETR)
C     no des SOMMETS des TETRAEDRES et no des TETRAEDRES VOISINS
C     ----------------------------------------------------------
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'DEMANDE ALLOCATION de',8*MXTETR,
     %                ' ENTIERS No DES SOMMETS des TETRAEDRES'
         ALLOCATE ( NOTETR( 1:8, 1:MXTETR ), STAT=IERTEALLOC )
         IF( IERTEALLOC .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION de',8*MXTETR,
     %                ' ENTIERS No DES SOMMETS des TETRAEDRES'
            IERR = IERTEALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'CORRECTE ALLOCATION de NOTETR(1:8,1:',MXTETR,
     %                ') ENTIERS No DES SOMMETS des TETRAEDRES'
      ELSE
         PRINT*, 'ALLOCATION DEMAND  of',8*MXTETR,
     %                ' INTEGER of TETRAHEDRA VERTICES NUMBER'
         ALLOCATE ( NOTETR( 1:8, 1:MXTETR ), STAT=IERTEALLOC )
         IF( IERTEALLOC .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',8*MXTETR,
     %                ' INTEGER of TETRAHEDRA VERTICES NUMBER'
            IERR = IERTEALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'ALLOCATION CORRECT of NOTETR(1:8,1:',MXTETR,
     %                ') INTEGER of TETRAHEDRA VERTICES NUMBER'
      ENDIF

C     ALLOCATION FORTRAN du tableau NVOLTE
C     LE NUMERO DE VOLUME DE 1 A NBVOPA DES TETRAEDRES DU VOLUME
C     ----------------------------------------------------------
C     LE TABLEAU 'MATERIAUX'
      CALL LXTSOU( NTLXQU, 'MATERIAUX', NTMATE, MNMATE )
      IF( NTMATE .LE. 0 ) THEN

C        UN SEUL MATERIAU = PAS DE NO DE MATERIAU DE CHAQUE EF
C                           LE TABLEAU NVOLTE N'EXISTE PAS
         IVOLTE  = 0
         MXVOLTE = 1
         NBVOPA  = 0
         NBDMEF  = 0
         MNUDMEF = 0

      ELSE

C        EN FAIT PLUSIEURS MATERIAUX = STOCKAGE DU NO DE VOLUME DE CHAQUE EF
C        CHAQUE EF A SON NUMERO DE MATERIAU=VOLUME RANGE A L'ADRESSE MNUDMEF
C        LE TABLEAU NVOLTE EXISTE
         IVOLTE  = 1
         MXVOLTE = MXTETR
C        NOMBRE DE MATERIAUX
         NBVOPA = MCN( MNMATE + WNBDM )
C        NOMBRE D''ELEMENTS FINIS DU MAILLAGE AVEC NO DE MATERIAU
         NBDMEF = MCN( MNMATE + WBDMEF )

         IF( NBVOPA .GT. 0 .AND. NBTETRVI .NE. NBDMEF ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOMBRE DE TETRAEDRES NON EGAL ENTRE VOLUME et 
     %MATERIAUX'
            ELSE
               KERR(1) = 'NUMBER of TETRAHEDRA NOT EQUAL BETWEEN VOLUME 
     %and MATERIALS'
            ENDIF
            CALL LEREUR
            IERR = 5
            RETURN
         ENDIF

C        ADRESSE MCN DU NO DE MATERIAU DE CHAQUE TETRAEDRE
         MNUDMEF= MNMATE + WUDMEF

      ENDIF

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'DEMANDE ALLOCATION de',MXVOLTE,
     %           ' ENTIERS pour NVOLTE'
         ALLOCATE ( NVOLTE( 1:MXVOLTE ), STAT=IERNVALLOC )
         IF( IERNVALLOC .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION de',MXVOLTE,
     %                ' ENTIERS pour NVOLTE'
            IERR = IERNVALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'CORRECTE ALLOCATION de NVOLTE(1:',MXVOLTE,
     %                ') ENTIERS pour NVOLTE'
      ELSE
         PRINT*, 'ALLOCATION DEMAND  of',MXVOLTE,
     %        ' INTEGER of NVOLTE'
         ALLOCATE ( NVOLTE( 1:MXVOLTE ), STAT=IERNVALLOC )
         IF( IERNVALLOC .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',MXVOLTE,
     %           ' INTEGER of NVOLTE'
            IERR = IERNVALLOC
            GOTO 9900
         ENDIF
         PRINT*, 'ALLOCATION CORRECT of NVOLTE(1:',MXVOLTE,
     %              ') INTEGER'
      ENDIF

      IF( IVOLTE .NE. 0 ) THEN
C        COPIE DU TMS MATERIAUX DANS LE TABLEAU NVOLTE
         CALL TRTATA( MCN(MNUDMEF), NVOLTE, NBDMEF )
C        COMPLETION PAR LA VALEUR -1 => VOLUME INCONNU
         DO N = NBDMEF+1, MXVOLTE
            NVOLTE( N ) = -1
         ENDDO
      ENDIF

C     INITIALISATION XYZ DU TABLEAU PTXYZD(1:4,1:NBSOMMVI)
C     ----------------------------------------------------
      MN = MNXYZSVI + WYZSOM - 1
      DO NS = 1, NBSOMMVI

         DO K = 1, 3
            PTXYZD( K, NS ) = RMCN( MN + K )
         ENDDO

C        LONGUEUR SOUHAITEE DE L'ARETE ISSUE DU POINT NS
C        PAR LE CALCUL DE TAILLE_IDEALE( NS ) SI LA FONCTION UTILISATEUR EXISTE
         CALL TAILIDEA( NOFOTI, PTXYZD(1,NS), NCODEV, PTXYZD(4,NS) )

         MN = MN + 3

      ENDDO
      NBSOMM = NBSOMMVI

C     INITIALISATION DU TABLEAU N1TETS(1:MXSOMM)
C     ------------------------------------------
C     MCN(MNTETS) ( MXSOMM )  NO 1 TETRAEDRE DE CHAQUE SOMMET
      MNTETS = -1
      CALL TNMCDC( 'ENTIER', MXSOMM, MNTETS )
      IF( MNTETS .LE. 0 ) GOTO 9995
      CALL AZEROI( MXSOMM, MCN(MNTETS) )

C     INITIALISATION DU TABLEAU NOTETR(1:8,1:NBTETRVI) et N1TETS
C     ----------------------------------------------------------
      MNTE = MNNSEFVI + WUSOEF -1
      MN1T = MNTETS - 1
      DO NUELEM = 1, NBTETRVI

C        REMPLISSAGE DU TABLEAU NSEF DE LA TETRAEDRISATION
         DO K = 1, 4

C           NO DU SOMMET K DU TETRAEDRE NUELEM
            NS = MCN( MNTE + K )
            NOTETR( K, NUELEM ) = NS

C           NO DU TETRAEDRE VOISIN PAR LA FACE K INCONNU
            NOTETR( 4+K, NUELEM ) = -1

C           NO D'UN TETRAEDRE DE SOMMET NS
            MCN( MN1T + NS ) = NUELEM

         ENDDO

         MNTE = MNTE + NBSOEF

      ENDDO

C     CONSTRUCTION DU TMS 'FACES' des TETRAEDRES du MAILLAGE
C     GENERATION DU TABLEAU NPSOFR(1:NBSOMM) a PARTIR de LFACES
C     =0 SI SOMMET INTERNE
C     =1 SI SOMMET FRONTALIER
C     =2 SI SOMMET SUR INTERFACE ENTRE 2 MATERIAUX
C     =3 SI SOMMET IMPOSE PAR L'UTILISATEUR LORS DE LA DEFINITION DU VOLUME
C     ---------------------------------------------------------------------
      CALL PTSUFR( KNMVOLUIN, NBSOMM, MNSOFR )
      IF( MNSOFR .LE. 0 ) GOTO 9995
C     AUGMENTATION DE LA TAILLE DU TABLEAU SOFR
      CALL TNMCAU( 'ENTIER', NBSOMM, MXSOMM, NBSOMM, MNSOFR )
C     PAR DEFAUT LE POINT A AJOUTER EST INTERNE
      DO K = MNSOFR+NBSOMM, MNSOFR-1+MXSOMM
         MCN( K ) = 0
      ENDDO

C     OUVERTURE du TABLEAU 'FACE' du VOLUME INITIAL et CONSTRUIT par PTSUFR
      CALL LXTSOU( NTLXVI, 'FACE', NTFAVO, MNFAVO )
      IF( MNFAVO .LE. 0 ) THEN
         IERR = 6
         GOTO 9999
      ENDIF
C     LE NOMBRE D'ENTIERS PAR FACE
      MOFACE = MCN( MNFAVO + WOFACE )
C     LE NOMBRE MAXIMAL DE FACES
      MXFACE = MCN( MNFAVO + WXFACE )
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      L1FAFR = MCN( MNFAVO + W1FAFR )
C     LE NOMBRE DE FACES FRONTALIERES
      NBFAFR = MCN( MNFAVO + WBFAFR )
C     LE NUMERO DE LA PREMIERE FACE INTERFACE
      L1FA2M = MCN( MNFAVO + W1FA2M )
C     LE NOMBRE DE FACES INTERFACES
      NBFA2M = MCN( MNFAVO + WBFA2M )
C     LE NOMBRE DE FACES AVEC DES TANGENTES
      NBFATG = MCN( MNFAVO + WBFATG )
      print*,'voex26: Nombre Faces Fontiere =',NBFAFR
      print*,'voex26: Nombre Faces Interface=',NBFA2M


C     PARCOURS DES FACES DES TETRAEDRES DU VOLUME
C     POUR LA MISE A JOUR DE NOTETR(5:8,*) TETRAEDRES VOISINS PAR LES FACES
C     ---------------------------------------------------------------------
      NBFACES = 0
      MNFACE  = MNFAVO + WFACES
      DO NOFACE = 1, MXFACE

         IF( MCN( MNFACE ) .GT. 0 ) THEN

C           LA FACE EXISTE: NO DES 3 SOMMETS DE LA FACE NOFACE
            NBFACES = NBFACES + 1
            DO K=1,3
               NOSOTR( K ) = MCN( MNFACE-1 + K )
            ENDDO

C           NO DU TETRAEDRE 1 DE CETTE FACE
            NOTET1 = ABS( MCN( MNFACE + 5 ) )

C           RECHERCHE DU NO DE FACE DANS LE TETRAEDRE 1
            IF( NOTET1 .GT. 0 ) THEN
               CALL NUFATRTE( NOSOTR, NOTETR(1,NOTET1), NUFATE1 )
               IF( NUFATE1 .EQ. 0 ) THEN
C                 LA FACE NOFACE N'EST PAS UNE FACE DU TETRAEDRE NOTET1
                  PRINT*,'voex26: LaFACE',NOFACE,' St',NOSOTR,
     %                   ' NON FACE du TETRAEDRE',NOTET1,
     %                   ':',(NOTETR(KK,NOTET1),KK=1,8),' ????'
                  NOTET1 = 0
               ENDIF
            ENDIF

C           NO DU TETRAEDRE 2 DE CETTE FACE
C           OU CHAINAGE SUR LES FACES DE LA FRONTIERE
            NOTET2 = ABS( MCN( MNFACE + 6 ) )

C           RECHERCHE DU NO DE FACE DANS LE TETRAEDRE 2
            IF( NOTET2 .GT. 0 ) THEN
               CALL NUFATRTE( NOSOTR, NOTETR(1,NOTET2), NUFATE2 )
               IF( NUFATE2 .EQ. 0 ) THEN
C                 C'EST UN CHAINAGE SUR LA FRONTIERE
ccc                  PRINT*,'voex26: LaFACE ',NOFACE,' St',NOSOTR,
ccc     %                   ' EST FRONTALIERE'
                  NOTET2 = 0
               ENDIF
            ENDIF

            IF( NOTET1 .GT. 0 .AND. NOTET2 .GT. 0 ) THEN
C              LA FACE EST COMMUNE A 2 TETRAEDRES
C              MISE A JOUR DES 2 TETRAEDRES VOISINS PAR CETTE FACE
               NOTETR( 4+NUFATE1, NOTET1 ) = NOTET2
               NOTETR( 4+NUFATE2, NOTET2 ) = NOTET1
               GOTO 20
            ENDIF

            IF( NOTET1 .GT. 0 .AND. NOTET2 .EQ. 0 ) THEN
C              NOTET2=0 LA FACE DE NOTET1 EST FRONTALIERE
               NOTETR( 4+NUFATE1, NOTET1 ) = 0
               GOTO 20
            ENDIF

            IF( NOTET2 .GT. 0 .AND. NOTET1 .EQ. 0 ) THEN
C              NOTET1=0 LA FACE DE NOTET2 EST FRONTALIERE
               NOTETR( 4+NUFATE2, NOTET2 ) = 0
               GOTO 20
            ENDIF

         ENDIF

 20      MNFACE = MNFACE + MOFACE

      ENDDO
      print*,'voex26: Nombre Faces Total    =',NBFACES

C     COMPLETION DE L'INITIALISATION DU TABLEAU NOTETR
C     LE CHAINAGE DES TETRAEDRES VIDES DEBUTE EN POSITION N1TEVI
C     PUIS SE FAIT SUR NOTETR(5,.)
      NUDTETR = NBTETRVI
      N1TEVI  = NBTETRVI+1
      DO NUELEM = N1TEVI, MXTETR
         NOTETR( 1, NUELEM ) = 0
         NOTETR( 5, NUELEM ) = NUELEM + 1
      ENDDO
C     FIN DU CHAINAGE DES TETRAEDRES VIDES
      NOTETR( 5, MXTETR ) = 0
      NBTETR = NBTETRVI


cccC     TENTATIVE D'AMELIORER LA QUALITE DE TOUS LES TETRAEDRES
cccC     DE QUALITE MEDIOCRE PAR DIFFERENTES TECHNIQUES SELON
cccC     LE TYPE DU TETRAEDRE PLAT ( TRIANGLE ou QUADRANGLE ) :
cccC     SUPPRESSION DES TETRAEDRES D'UNE ARETE TRES COURTE
cccC     CHANGEMENTS 3 TETRAEDRES -> 2 TETRAEDRES
cccC     DECOMPOSITION DE TETRAEDRES EN 2 TETRAEDRES
cccC     AMELIORER LA QUALITE DES TETRAEDRES AUTOUR 
cccC     DU TETRAEDRE NTE DE MEDIOCRE QUALITE ET DE LUI-MEME
cccC     -------------------------------------------------------
ccc      CALL AMQALLTE( QTEAME,  MXTETR, N1TEVI, NOTETR, MCN(MNTETS),
ccc     %               MXSOMM,  NBSOMM, PTXYZD, MCN(MNSOFR),
ccc     %               INFACO,  MXFACO, LEFACO, IVOLTE, NVOLTE,
ccc     %               MXTE1S,   MCN(MNTE1S), MXETOI, MCN(MNFETO),
ccc     %               MXTRCF,   MCN(MNTRCF),
ccc     %               NUDTETR, VOLMOY, MODIFT, IERR )

ccc      CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
ccc     %               NBTETR, NUDTETR,   QUAMIN, QUAMOY, VOLUMT )
ccc      VOLMOY = VOLUMT / NBTETR

ccc      CPU  = DINFO('DELTA CPU')
ccc      CPUT = CPUT + CPU
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,12001) 'amqallte',CPU
ccc      ELSE
ccc         WRITE(IMPRIM,22001) 'amqallte',CPU
ccc      ENDIF

12001 FORMAT(' ',A,T12,F11.2,'  SECONDES CPU')
22001 FORMAT(' ',A,T12,F11.2,'  CPU SECONDS')


cccC     AJOUT DE POINTS DANS LES TETRAEDRES POUR AJUSTER LES ARETES
cccC     A LA DISTANCE SOUHAITEE ENTRE LEURS SOMMETS
cccC     ===========================================================
cccC     TETRAEDRISATION DE POINTS SUR L'ARETE MAX DES TETRAEDRES
cccC     SI SA TAILLE EST SUPERIEURE A LA TAILLE SOUHAITEE EN AU
cccC     MOINS UN DES SOMMETS PUIS
cccC     AMELIORIATION DES TETRAEDRES AU DELA DES DECOUPES PAR
cccC     2T->3T ou mT->2m-4T
ccc      MXTE1A = 3*MXFETO
ccc      MXSSTE = 2*MXFETO
ccc      CALL AJPTARMX( ARETGR, NOFOTI, NBVOPA, VOLMOY,
ccc     %               NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
ccc     %               MXTETR, N1TEVI, NOTETR, NUDTETR,
ccc     %               MCN(MNTETS), MCN(MNF1VO),
ccc     %               MXFACO, LEFACO, IVOLTE, NVOLTE,
ccc     %               MXPILE, MCN(MNPILE),
ccc     %               MXTE1A, MCN(MNFETO),
ccc     %               MXSSTE, MCN(MNFETO+MXTE1A),
ccc     %               IERR )
cccccc      DEJA EXECUTE DANS ajptarmx
cccccc      CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
cccccc     %               NBTETR, NUDTETR,   QUAMIN, QUAMOY, VOLUMT )
ccc      CPU  = DINFO('DELTA CPU')
ccc      CPUT = CPUT + CPU
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,12001) 'ajptarmx:',CPU
ccc      ELSE
ccc         WRITE(IMPRIM,22001) 'ajptarmx:',CPU
ccc      ENDIF

cccC     MISE A JOUR DU NOMBRE DE SOMMETS ACTIFS ET
cccC                 DU TABLEAU N1TETS(NS)=1 TETRAEDRE DE SOMMET NS
cccC     ----------------------------------------------------------
ccc      CALL AZEROI( MXSOMM, MCN(MNTETS) )
ccc      NBSOMM = 0
ccc      MNT1   = MNTETS - 1
ccc      DO NTE = 1, NUDTETR
ccc         IF( NOTETR(1,NTE) .GT. 0 ) THEN
ccc            DO K=1,4
ccc               NS = NOTETR( K, NTE )
ccc               IF( NS .GT. NBSOMM ) THEN
ccc                  NBSOMM = NS
ccc               ENDIF
ccc               MCN( MNT1 + NS ) = NTE
ccc            ENDDO
ccc         ENDIF
ccc      ENDDO
ccc      PRINT*,'voex26: NOMBRE DE SOMMETS ACTIFS NBSOMM=',NBSOMM


C     CREATION DE L'HEXAEDRE DU PAVAGE DE L'HEXAEDRE ENGLOBANT
C     POUR ACCELERER LA RECHERCHE DU PLUS PROCHE SOMMET D'UN POINT
C     ------------------------------------------------------------
      DO K=1,3
         HEXAPAVE(K,1) = COOEXT(K,1) - ARETGR
         HEXAPAVE(K,2) = COOEXT(K,2) + ARETGR
      ENDDO

C     CREATION DU PAVAGE DE L'HEXAEDRE ENGLOBANT
      CALL CREERPAVAGE( ARETGR, HEXAPAVE, NBIPAV, ECHPAV )

C     NOMBRE DE PAVES DE L'HEXAEDRE
      NBPAVE = (NBIPAV(1)+1) * (NBIPAV(2)+1) * (NBIPAV(3)+1)

C     CREATION ET INITIALISATION DU TABLEAU DU PAVAGE
C     POINTEUR SUR UN SOMMET PTXYZD DE CHAQUE PAVE
      CALL TNMCDC( 'ENTIER', NBPAVE, MN1SPAVE )
      CALL AZEROI(  NBPAVE,  MCN(MN1SPAVE) )

C     NOMBRE DE SOMMETS DES CUBES DE L'HEXAEDRE ENGLOBANT
      NBSTPV = (NBIPAV(1)+2) * (NBIPAV(2)+2) * (NBIPAV(3)+2)

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10002) ((HEXAPAVE(I,J),I=1,3),J=1,2),
     %                       NBSTPV, NBPAVE, ARETGR
      ELSE
         WRITE(IMPRIM,20002) ((HEXAPAVE(I,J),I=1,3),J=1,2),
     %                       NBSTPV, NBPAVE, ARETGR
      ENDIF

10002 FORMAT(/
     % ' PAVAGE de l''HEXAEDRE: COORDONNEES MIN-MAX'/
     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 /
     % I9,' SOMMETS des',I9,' PAVES de LONGUEUR d''ARETE',G15.7)

20002 FORMAT(/
     % ' HEXAHEDRON PAVEMENT: MIN-MAX COORDINATES'/
     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 /
     % I9,' VERTICES of',I9,' CUBES of EDGE LENGTH',G15.7)

C     NO DU SOMMET SUIVANT DE CHAQUE SOMMET DANS SON PAVE
      CALL TNMCDC( 'ENTIER', MXSOMM, MNOPTSUIV )
      CALL AZEROI(  MXSOMM,  MCN(MNOPTSUIV) )

      DO NS = 1, NBSOMM
C        AJOUT DANS LE PAVAGE DU POINT NS
         IF( MCN(MNTETS-1+NS) .GT. 0 ) THEN
            CALL NUPAVEST( NS, PTXYZD, HEXAPAVE, NBIPAV, ECHPAV,
     %                     MCN(MN1SPAVE), MCN(MNOPTSUIV) )
         ENDIF
      ENDDO

C     DECLARATION DES PILES SUR L'ETOILE NFETOI ET NTETOI
C     ---------------------------------------------------
      MXFETO = MAX( 4096, NBFAFR )
      MNFETO = -1
      CALL TNMCDC( 'ENTIER', 5*MXFETO, MNFETO )
      IF( MNFETO .LE. 0 ) GOTO 9995

      MNTETO = -1
      CALL TNMCDC( 'ENTIER',   MXFETO, MNTETO )
      IF( MNTETO .LE. 0 ) GOTO 9995

C     MCN(MNLIPO) (0:MXSOMM) NO DU POINT A TETRAEDRISER AVANT ET APRES TRI
C     --------------------------------------------------------------------
      MNLIPO = -1
      CALL TNMCDC( 'ENTIER', 1+MXSOMM, MNLIPO )
      IF( MNLIPO .LE. 0 ) GOTO 9995

C     VOLUME ET QUALITE DU MAILLAGE INITIAL
C     -------------------------------------
      print*
      PRINT*,'voex26: VOLUME et QUALITE du MAILLAGE INITIAL'
      CALL QUALTETR( PTXYZD, MXTETR,  NOTETR,
     %               NBTETR, NUDTETR, QUAMIN, QUAMOY, VOLUMT )
      VOLMOY = VOLUMT / NBTETR

C     TETRAEDRISATION DES NBPTIM POINTS IMPOSES COMME SOMMETS PAR L'UTILISATEUR
C     =========================================================================
      IF( NBPTIM .GT. 0 ) THEN
         CALL TETRPTIM( KNMVOLU, QTEAME, VOLMOY, NBPTIM, RADEFI(WYZDIM),
     %                  NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
     %                  HEXAPAVE, NBIPAV, ECHPAV,
     %                  MCN(MN1SPAVE), MCN(MNOPTSUIV),
     %                  MXTETR, N1TEVI, NOTETR, NUDTETR, NBTETR,
     %                  MCN(MNTETS),
     %                  INFACO, MXFACO, LEFACO, N1FASC,
     %                  IVOLTE, NVOLTE,
     %                  MXFETO, MCN(MNFETO), MCN(MNTETO), MCN(MNLIPO),
     %                  IERR )
         CPU  = DINFO('DELTA CPU')
         CPUT = CPUT + CPU
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,12001) 'tetrptim:',CPU
         ELSE
            WRITE(IMPRIM,22001) 'tetrptim:',CPU
         ENDIF
         print*
         PRINT*,'voex26: VOLUME et QUALITE du MAILLAGE apres TETRPTIM'
         CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
     %                  NBTETR, NUDTETR,   QUAMIN, QUAMOY, VOLUMT )
         VOLMOY = VOLUMT / NBTETR

      ENDIF


C     TETRAEDRISATION DE POINTS SUR UNE GRILLE D'OCTAEDRES REGULIERS
C     DE TAILLE D'ARETE ARETGR ou FONCTAION TAILLE_IDEALE(x,y,z)
C     DONNEES DE L'UTILISATEUR
C     ==============================================================
      IF( NOFOTI .EQ. 0 ) THEN
         CALL TETRGROC( KNMVOLU, QTEAME, COOEXT, VOLMOY,
     %                  NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
     %                  HEXAPAVE, NBIPAV, ECHPAV,
     %                  MCN(MN1SPAVE), MCN(MNOPTSUIV),
     %                  MXTETR, N1TEVI, NOTETR, NUDTETR, NBTETR,
     %                  MCN(MNTETS),
     %                  INFACO, MXFACO, LEFACO, N1FASC,
     %                  IVOLTE, NVOLTE,
     %                  MXFETO, MCN(MNFETO), MCN(MNTETO), MCN(MNLIPO),
     %                  IERR )
         CPU  = DINFO('DELTA CPU')
         CPUT = CPUT + CPU
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,12001) 'tetrgroc:',CPU
         ELSE
            WRITE(IMPRIM,22001) 'tetrgroc:',CPU
         ENDIF
         print*
         PRINT*,'voex26: VOLUME et QUALITE du MAILLAGE apres TETRGROC'
         CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
     %                  NBTETR, NUDTETR,   QUAMIN, QUAMOY, VOLUMT )
         VOLMOY = VOLUMT / NBTETR
      ENDIF


C     ITERATIONS DE TETRAEDRISATION DE POINTS A L'INTERIEUR DES
C     TETRAEDRES TROP VOLUMINEUX DE LA BOULE DE CENTRE UN SOMMET
C     ET DE RAYON LA TAILLE SOUHAITEE DES ARETES AUTOUR DU SOMMET
C     AINSI QUE LES POINTS SUR LES ARETES TROP LONGUES ISSUES DU SOMMET
C     =================================================================
      CALL TETRARBO( NUVOLU, KNMVOLU, QUAMINEX, QTEAME, VOLMOY,
     %               NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
     %               HEXAPAVE, NBIPAV, ECHPAV,
     %               MCN(MN1SPAVE), MCN(MNOPTSUIV),
     %               MXTETR, N1TEVI, NOTETR, NUDTETR, NBTETR,
     %               MCN(MNTETS),
     %               INFACO, MXFACO, LEFACO, N1FASC, IVOLTE, NVOLTE,
     %               MXFETO, MCN(MNFETO), MCN(MNTETO), MCN(MNLIPO),
     %               IERR )
      CPU  = DINFO('DELTA CPU')
      CPUT = CPUT + CPU
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) 'tetrarbo:',CPU
      ELSE
         WRITE(IMPRIM,22001) 'tetrarbo:',CPU
      ENDIF
      print*
      PRINT*,'voex26: VOLUME et QUALITE du MAILLAGE apres TETRARBO'
      CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
     %               NBTETR, NUDTETR,   QUAMIN, QUAMOY, VOLUMT )
      VOLMOY = VOLUMT / NBTETR

      IF( NBTETR .GT. 0 ) GOTO 8000

cccC     TETRAEDRISATION DE POINTS SUR LES ARETES ET A L'INTERIEUR
cccC     DES TETRAEDRES TROP VOLUMINEUX AU REGARD DE LA TAILLE
cccC     SOUHAITEE DES ARETES DEFINIES EN LEURS 4 SOMMETS
cccC     ===========================================================
ccc      CALL TETROPVO( NUVOLU, KNMVOLU, QUAMINEX, QTEAME, NBVOPA, VOLMOY,
ccc     %               NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
ccc     %               HEXAPAVE, NBIPAV, ECHPAV,
ccc     %               MCN(MN1SPAVE), MCN(MNOPTSUIV),
ccc     %               MXTETR, N1TEVI, NOTETR, NUDTETR, NBTETR,
ccc     %               MCN(MNTETS), MCN(MNF1VO),
ccc     %               INFACO, MXFACO, LEFACO, N1FASC, IVOLTE, NVOLTE,
ccc     %               MXFETO, MCN(MNFETO), MCN(MNTETO), MCN(MNLIPO),
ccc     %               IERR )
ccc      CPU  = DINFO('DELTA CPU')
ccc      CPUT = CPUT + CPU
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,12001) 'tetropvo:',CPU
ccc      ELSE
ccc         WRITE(IMPRIM,22001) 'tetropvo:',CPU
ccc      ENDIF


C     ITERATIONS de TETRAEDRISATION DU BARYCENTRE DES TETRAEDRES
C     TROP VOLUMINEUX AU REGARD DE LA TAILLE SOUHAITEE D'ARETES
C     AUX 4 SOMMETS
C     ==========================================================
      CALL TETROPVL( NUVOLU, KNMVOLU, QUAMINEX, QTEAME,
ccc     %               VOLMOY,
     %               NBSOMM, MXSOMM, PTXYZD, MCN(MNSOFR),
     %               HEXAPAVE, NBIPAV, ECHPAV,
     %               MCN(MN1SPAVE), MCN(MNOPTSUIV),
     %               MXTETR, N1TEVI, NOTETR, NUDTETR, NBTETR,
     %               MCN(MNTETS),
     %               INFACO, MXFACO, LEFACO, N1FASC, IVOLTE, NVOLTE,
     %               MXFETO, MCN(MNFETO), MCN(MNTETO), MCN(MNLIPO),
     %               IERR )
      CPU  = DINFO('DELTA CPU')
      CPUT = CPUT + CPU
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) 'tetropvl:',CPU
      ELSE
         WRITE(IMPRIM,22001) 'tetropvl:',CPU
      ENDIF


C     RENUMEROTATION DES SOMMETS DE LA TETRAEDRISATION AMELIOREE
C     ----------------------------------------------------------
C     CONSTRUCTION DU NOUVEAU NO DES SOMMETS DE LA TETRAEDRISATION
 8000 MXNEWS = 1+NBSOMM
      CALL TNMCDC( 'ENTIER', MXNEWS, MNNEWS )
      IF( MNNEWS .LE. 0 ) GOTO 9995
      CALL AZEROI( MXNEWS, MCN(MNNEWS) )
      DO NUELEM = 1, NUDTETR
         IF( NOTETR(1,NUELEM) .GT. 0 ) THEN
C           LE TETRAEDRE NUELEM EXISTE
            DO K = 1, 4
               NS = NOTETR( K, NUELEM )
C              TEMOIN DE SOMMET ACTIF
               MCN( MNNEWS + NS ) = 1
            ENDDO
         ENDIF
      ENDDO

      NBSOM = 0
      DO NS = 1, NBSOMM
         IF( MCN( MNNEWS + NS ) .GT. 0 ) THEN
C           LE NOUVEAU NUMERO DU SOMMET NS
            NBSOM = NBSOM + 1
            MCN( MNNEWS + NS ) = NBSOM
         ENDIF
      ENDDO

C     CONSTRUCTION DU TMS 'XYZSOMMET' DE LA TETRAEDRISATION AMELIOREE
C     ---------------------------------------------------------------
      CALL LXTSOU( NTLXVOLU, 'XYZSOMMET',  NTXYZS, MNXYZS )
      IF( MNXYZS .GT. 0 ) CALL LXTSDS( NTLXVOLU, 'XYZSOMMET' )
      CALL LXTNDC( NTLXVOLU, 'XYZSOMMET', 'MOTS', WYZSOM + 3 * NBSOM )
      CALL LXTSOU( NTLXVOLU, 'XYZSOMMET',  NTXYZS, MNXYZS )
C     NOMBRE DE SOMMETS DE LA TETRAEDRISATION
      MCN( MNXYZS + WNBSOM ) = NBSOM
C     NOMBRE DE TANGENTES DE LA TETRAEDRISATION
      MCN( MNXYZS + WNBTGS ) = 0
C     NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNXYZS + WBCOOR ) = 3
C     XYZ DES SOMMETS
      DO NS = 1, NBSOMM
C        LE NOUVEAU NUMERO DU SOMMET NS
         NSNEW = MCN( MNNEWS + NS )
         IF( NSNEW .GT. 0 ) THEN
            MN = MNXYZS + WYZSOM + 3*NSNEW - 4
            DO K = 1, 3
               RMCN( MN + K ) = REAL( PTXYZD( K, NS ) )
            ENDDO
         ENDIF
         MN = MN + 3
      ENDDO
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNXYZS) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C     NOMBRE EFFECTIF DE TETRAEDRES NBT APRES AMELIORATION
      NBT = 0
      DO NUELEM = 1, NUDTETR
         IF( NOTETR(1,NUELEM) .GT. 0 ) THEN
C           LE TETRAEDRE EXISTE
            NBT = NBT + 1
         ENDIF
      ENDDO

C     CONSTRUCTION DU TABLEAU 'NSEF' DE LA TETRAEDRISATION AMELIOREE
C     --------------------------------------------------------------
      CALL LXTSOU( NTLXVOLU, 'NSEF',  NTNSEF, MNNSEF )
      IF( MNNSEF .GT. 0 ) CALL LXTSDS( NTLXVOLU, 'NSEF' )
      CALL LXTNDC( NTLXVOLU, 'NSEF', 'MOTS', WUSOEF+NBT*NBSOEF )
      CALL LXTSOU( NTLXVOLU, 'NSEF',  NTNSEF, MNNSEF )
C     TYPE DE L'OBJET : VOLUME
      MCN( MNNSEF + WUTYOB ) = 4
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNNSEF + WUTFMA ) = -1
C     NOMBRE DE SOMMETS PAR EF
      MCN( MNNSEF + WBSOEF ) = 8
C     PAS DE TANGENTES STOCKEES
      MCN( MNNSEF + WBTGEF ) = 0
C     NOMBRE DE TETRAEDRES
      MCN( MNNSEF + WBEFOB ) = NBT
C     NOMBRE DE TETRAEDRES A TG
      MCN( MNNSEF + WBEFTG ) = 0
C     NOMBRE D'EF AVEC POINTEUR SUR LES EF A TG
      MCN( MNNSEF + WBEFAP ) = 0
C     NUMERO DU TYPE DU MAILLAGE : NON STRUCTURE
      MCN( MNNSEF + WUTYMA ) = 0
C     MISE A JOUR DU NO DES SOMMETS DES TETRAEDRES ACTIFS
      NBTETR = 0
      MN     = MNNSEF + WUSOEF -1
      DO NUELEM = 1, NUDTETR
         IF( NOTETR(1,NUELEM) .GT. 0 ) THEN
C           LE TETRAEDRE EXISTE
            NBTETR = NBTETR + 1
            DO K = 1, 4
C              ANCIEN  NO DU SOMMET
               NS    = NOTETR( K, NUELEM )
C              NOUVEAU NO DU SOMMET
               NSNEW = MCN( MNNEWS + NS )
               IF( NSNEW .EQ. 0 ) THEN
                  PRINT*,'voex26: SOMMET',NS,' DE NOUVEAU NO',
     %                    NSNEW,' A PROBLEME'
               ENDIF
               MCN( MN + K ) = NSNEW
C              NO DU SOMMET AU DELA DE 4 EST NUL
               MCN( MN + 4 + K ) = 0
            ENDDO
            MN = MN + 8
         ENDIF
      ENDDO
      IF( NBTETR .NE. NBT ) THEN
         PRINT*,'voex26: PB NBTETR=',NBTETR,' NON EGAL a NBT=',NBT
      ENDIF
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNNSEF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
C     -----------------------------------------
 9900 PRINT*
      IF( LANGAG .EQ. 0 ) THEN

         IF( IERNVALLOC .EQ. 0 ) then
            PRINT*,'NVOLTE est DESALLOUE'
            DEALLOCATE( NVOLTE )
         ENDIF
         IF( IERTEALLOC .EQ. 0 ) then
            PRINT*,'NOTETR est DESALLOUE'
            DEALLOCATE( NOTETR )
         ENDIF
         IF( IERZDALLOC .EQ. 0 ) then
            PRINT*,'PtXYZd est DESALLOUE'
            DEALLOCATE( PTXYZD )
         ENDIF

      ELSE

        IF( IERNVALLOC .EQ. 0 ) then
            PRINT*,'NVOLTE is DESALLOCATED'
            DEALLOCATE( NVOLTE )
         ENDIF
         IF( IERTEALLOC .EQ. 0 ) then
            PRINT*,'NOTETR is DESALLOCATED'
            DEALLOCATE( NOTETR )
         ENDIF
         IF( IERZDALLOC .EQ. 0 ) then
            PRINT*,'PtXYZd is DESALLOCATED'
            DEALLOCATE( PTXYZD )
         ENDIF

      ENDIF

      IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNEWS, MNNEWS )
      IF( MNOPTSUIV.GT.0) CALL TNMCDS( 'ENTIER', MXSOMM, MNOPTSUIV)
      IF( MN1SPAVE .GT.0) CALL TNMCDS( 'ENTIER', NBPAVE, MN1SPAVE )
      IF( MNTETO .GT. 0 ) CALL TNMCDS( 'ENTIER', MXFETO  , MNTETO )
      IF( MNFETO .GT. 0 ) CALL TNMCDS( 'ENTIER', 5*MXFETO, MNFETO )
      IF( MNLIPO .GT. 0 ) CALL TNMCDS( 'ENTIER', MXSOMM+1, MNLIPO )
      IF( MNSOFR .GT. 0 ) CALL TNMCDS( 'ENTIER', MXSOMM,   MNSOFR )
      IF( MNTETS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXSOMM,   MNTETS )


C     AFFICHAGE DU TEMPS CPU TOTAL UTILISE PAR LA TETRAEDRISATION
C     -----------------------------------------------------------
      IF( IERR .NE. 0 ) THEN

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'voex26: SUITE AU PROBLEME VU AVANT sur le VOLUME ',
     %              KNMVOLU
       PRINT*,'voex26: la TETRAEDRISATION AMELIOREE N''EST PAS REALISEE'
         ELSE
            PRINT*,'voex26: From the PROBLEM SEEN BEFORE on VOLUME',
     %              KNMVOLU
           PRINT*,'voex26: the IMPROVED TETRAHEDRIZATION is NOT CREATED'
         ENDIF

         GOTO 9997

      ELSE

C        TEMPS CPU DE LA GENERATION FINALE DES TETRAEDRES
         CPU  = DINFO('DELTA CPU')
         CPUT = CPUT + CPU
         WRITE(IMPRIM,19800) KNMVOLU
         IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,19900) CPU,CPUT,NBSOM,NBT,VOLUMT,KINFO('MACHINE')
         ELSE
       WRITE(IMPRIM,29900) CPU,CPUT,NBSOM,NBT,VOLUMT,KINFO('MACHINE')
         ENDIF

      ENDIF

19800 FORMAT(' voex26 FIN: VOLUME:',A)
19900 FORMAT(' TEMPS CPU CREATION FINALE des TETRAEDRES et STOCKAGE du M
     %AILLAGE',F11.2,' SECONDES'//
     %' voex26: TEMPS CPU TOTAL DE LA TETRAEDRISATION',F11.2,
     %' SECONDES pour',I9,' SOMMETS et',I9,' TETRAEDRES et un VOLUME TOT
     %AL de ',G20.12,' sur l''ordinateur ',A/166('='),/)

29900 FORMAT(' CPU TIME of FINAL TETRAHEDRA CREATION and MESH STORAGE ',
     %F11.2,' SECONDS'//
     %' voex26: TETRAHEDRIZATION TOTAL CPU TIME',F11.2,
     %' SECONDS for',I9,' VERTICES and',I9,' TETRAHEDRA of TOTAL VOLUME'
     %,G20.12,' on COMPUTER ',A/ 166('='),/)

      GOTO 9999


C     SATURATION DE LA PLACE MEMOIRE
C     ------------------------------
 9995 NBLGRC(NRERR) = 4
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='voex26: SATURATION MEMOIRE OCCUPEE par les TABLEAUX'
         KERR(2)='AUGMENTEZ la TAILLE du SUPER-TABLEAU MCN'
         KERR(3)='ou AUGMENTEZ le NOMBRE MAXIMAL de SOMMETS DE LA TETRAE
     %DRISATION'
         KERR(4) ='ou REDUISEZ LE NOMBRE DE TRIANGLES INITIAUX'
      ELSE
         KERR(1)='voex26: ALL MEMORY is USED by the TMS'
         KERR(2)='AUGMENT the LENGTH of SUPER-ARRAY MCN'
       KERR(3)='or AUGMENT the TETRAHEDRIZATION VERTICES MAXIMUM NUMBER'
         KERR(4)='or REDUCE the INITIAL TRIANGLES NUMBER'
      ENDIF
      CALL LEREUR
      GOTO 9900

C     ERREUR
 9997 IERR = 100


C     FIN DU TRAITEMENT
 9999 RETURN
      END
