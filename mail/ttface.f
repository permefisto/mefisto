      SUBROUTINE TTFACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER TOUTES LES FACES TRIANGULAIRES OU QUADRANGULAIRES
C -----    SANS LES TANGENTES DES SURFACES ACTUELLES DU LEXIQUE SURFACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C....................................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/xyzext.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
C
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
C
      REAL          X(4), Y(4)
      REAL          XYZP(3), CNORFA(3), COUL(5), XYZ(3,5), XYZB(3)
      INTEGER       NOSOEL(1:64)
      CHARACTER*24  NMSURF
      CHARACTER*8   NMSOMM
C
      MNTSOM = 0
      MNTEF  = 0
      MNBAEF = 0
      MNDIST = 0
      MNNOFA = 0
      MNQUEF = 0
      NBCOUL = NDCOUL - N1COUL
      ZMIN   = 0
      ZMAX   = 0
C
C     RECHERCHE DE L'ADRESSE MCN DU LEXIQUE DES SURFACES
      CALL TAMSOU( NTSURF, MNSURF )
C
C     LE DEBUT DU CHAINAGE DES SURFACES OCCUPEES DANS LE LEXIQUE
      NUSF0 = MCN( MNSURF + 5 )
      IF( NUSF0 .LE. 0 ) RETURN
C
C     NOMBRE D'ENTIERS POUR STOCKER UN NOM DE SURFACE
      NBENNM = MCN( MNSURF + 2 )
C
C     LA BOUCLE SUR LES SURFACES OCCUPEES POUR CALCULER LES DIMENSIONS
C     ================================================================
      NBTSOM = 0
      NBTEF  = 0
      NUSURF = NUSF0
C
 10   IF( NUSURF .GT. 0 ) THEN
C
C        ADRESSE MCN DU DEBUT DE LA SURFACE DANS LE LEXIQUE
         MNSF = MNSURF + MCN(MNSURF) * NUSURF
C
C        LE LEXIQUE DE CETTE SURFACE EXISTE-T-IL ?
         NTLXSF = MCN( MNSF + NBENNM + 2 )
C
         IF( NTLXSF .GT. 0 ) THEN
C
C           CET OBJET EXISTE : RECHERCHE DE SON NOM
            CALL ENTNOM( NBENNM, MCN(MNSF), NMSURF )
C           RECHERCHE OU GENERATION DU MAILLAGE DE LA SURFACE
            NBLGRC(NRERR) = 0
            CALL MAILEX( 'SURFACE', NMSURF,
     %                    NTNSEF, MNNSEF, NTSOM, MNSOM, IERR )
C           SI SORTIE DE LA DESTRUCTION DU PLSV (IERR=82)
            IF( IERR .EQ. 82 ) RETURN
            IF( IERR .NE. 0  ) THEN
                NBLGRC(NRERR) = NBLGRC(NRERR) + 1
                IF( LANGAG .EQ. 0 ) THEN
                KERR(NBLGRC(NRERR)) = 'SURFACE SANS MAILLAGE: '//NMSURF
                ELSE
                KERR(NBLGRC(NRERR)) = 'SURFACE WITHOUT MESH: '//NMSURF
                ENDIF
                CALL LEREUR
C               DESTRUCTION DU PLSV NON MAILLE
                CALL LXLXDS( NTSURF, NMSURF )
                RETURN
            ENDIF
C
C           UNE SURFACE ACTIVE
C           NOMBRE DE SOMMETS
            NBTSOM = NBTSOM + MCN( MNSOM + WNBSOM )
C
C           NOMBRE DE TRIANGLES-QUADRANGLES
            NBTEF = NBTEF + MCN( MNNSEF + WBEFOB )
C
C           NOMBRE DE SOMMETS PAR EF
            NBSOEF = MCN( MNNSEF + WBSOEF )
            IF( NBSOEF .NE. 4 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1)='TTFACE: SURFACE A EF N''AYANT PAS 4 SOMMETS'
               ELSE
                  KERR(1)='TTFACE: SURFACE WITH FE WITH NOT 4 VERTICES'
               ENDIF
               CALL LEREUR
               RETURN
            ENDIF
C
         ENDIF
C
C        PASSAGE A LA SURFACE SUIVANTE
         NUSURF = MCN( MNSF + NBENNM )
         GOTO 10
C
      ENDIF
C
C     RESERVATION DES TABLEAUX DES XYZ DES SOMMETS ET NO DES SOMMETS
      CALL TNMCDC( 'REEL',   3*NBTSOM,     MNTSOM )
      CALL TNMCDC( 'ENTIER', NBSOEF*NBTEF, MNTEF  )
C     XYZ DU BARYCENTRE DE CHAQUE EF
      CALL TNMCDC( 'REEL', 3*NBTEF, MNBAEF )
      IF( MNTSOM .LE. 0 .OR. MNTEF .LE. 0 .OR. MNBAEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TTFACE: PAS ASSEZ DE MC POUR LES SURFACES'
         ELSE
            KERR(1) = 'TTFACE: NOT ENOUGH MEMORY FOR THE SURFACES'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     CONCATENATION DES TABLEAUX XYZSOMMET ET NSEF
C     ============================================
      MNXYZ  = MNTSOM
      MNEF   = MNTEF
      MNBA   = MNBAEF
      NBS    = 0
      NUSURF = NUSF0
C
 40   IF( NUSURF .GT. 0 ) THEN
C
C        ADRESSE MCN DU DEBUT DE LA SURFACE DANS LE LEXIQUE
         MNSF = MNSURF + MCN(MNSURF) * NUSURF
C
C        LE LEXIQUE DE CETTE SURFACE EXISTE-T-IL ?
         NTLXSF = MCN( MNSF + NBENNM + 2 )
C
         IF( NTLXSF .GT. 0 ) THEN
C
C           OUVERTURE DU TABLEAU 'XYZSOMMET'
            CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTSOM, MNSOM  )
            IF( NTSOM .LE. 0 ) GOTO 9990
C           NOMBRE DE SOMMETS
            NBSOM = MCN( MNSOM + WNBSOM )
            IF( NBSOM .LE. 0 ) GOTO 9990
C
            CALL LXTSOU( NTLXSF, 'NSEF', NTNSEF, MNNSEF )
            IF( NTNSEF .LE. 0 ) GOTO 9990
C
C           LES PARAMETRES DES NO SOMMET DU MAILLAGE
            CALL NSEFPA( MCN(MNNSEF),
     %                   NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF, NBEF,
     %                   NX    , NY    , NZ    , IERR   )
            IF( IERR .NE. 0 ) GOTO 9990
C
C           COPIE DU TABLEAU DES XYZ DES SOMMETS
            CALL TRTATA( MCN(MNSOM+WYZSOM), MCN(MNXYZ), 3*NBSOM )
C
C           COPIE AVEC DECALAGE DES NUMEROS DE SOMMETS DES EF
            DO 60 NEF=1,NBEF
C
C              LE NUMERO DES NBSOEF SOMMETS DE L'EF NEF
               CALL NSEFNS( NEF   , NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNNSEF, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C              LE NOMBRE DE SOMMETS DE CET EF = NCOGEL 3 OU 4
               DO 50 J=1,NCOGEL
                  MCN(MNEF-1+J) = NOSOEL(J) + NBS
 50            CONTINUE
               IF( NCOGEL .EQ. 3 ) THEN
C                 C'EST UN TRIANGLE
                  MCN(MNEF+3) = 0
               ENDIF
C
C              CONSTRUCTION DU BARYCENTRE DE CHAQUE FACE
               XB = 0
               YB = 0
               ZB = 0
               DO 54 J=1,NCOGEL
                  NS = MCN(MNEF-1+J)
                  MN = MNTSOM - 3 + 3 * NS
                  XB = XB + RMCN( MN     )
                  YB = YB + RMCN( MN + 1 )
                  ZB = ZB + RMCN( MN + 2 )
 54            CONTINUE
               RMCN( MNBA     ) = XB / NCOGEL
               RMCN( MNBA + 1 ) = YB / NCOGEL
               RMCN( MNBA + 2 ) = ZB / NCOGEL
C
               MNBA = MNBA + 3
               MNEF = MNEF + 4
 60         CONTINUE
C
C           POINTEURS POUR LA SURFACE SUIVANTE
            MNXYZ = MNXYZ + 3*NBSOM
            NBS   = NBS   +   NBSOM
C
         ENDIF
C
C        PASSAGE A LA SURFACE SUIVANTE
         NUSURF = MCN( MNSF + NBENNM )
         GOTO 40
C
      ENDIF
C
C     SI QUALITE DEMANDEE CALCUL DE LA QUALITE DE CHAQUE EF
      IF( LCRITR .GT. 0 ) THEN
         CALL TNMCDC( 'REEL', NBTEF, MNQUEF )
         QEFMIN = 2
         NEFMIN = 0
         MNEF   = MNTEF
         DO 96 NEF = 1, NBTEF
            IF( MCN(MNEF+3) .EQ. 0 ) THEN
C              TRIANGLE
               NCOGEL = 3
            ELSE
C              QUADRANGLE
               NCOGEL = 4
            ENDIF
            CALL QUALEF( NCOGEL,   MCN(MNEF), NBTSOM, RMCN(MNTSOM),
     %                   SURFVOLU, QUALIT,    IERR )
            IF( IERR .NE. 0 ) GOTO 9990
            IF( QUALIT .LT. QEFMIN ) THEN
C               REPERAGE DE L'EF DE PLUS MAUVAISE QUALITE
                QEFMIN = QUALIT
                NEFMIN = NEF
            ENDIF
            RMCN( MNQUEF-1+NEF ) = QUALIT
            MNEF = MNEF + 4
 96      CONTINUE
      ENDIF
      GOTO 102
C
C     =====================================================
C     SAISIE DE L'OPTION DE TRACE DE LA SURFACE EN 2D OU 3D
C     =====================================================
 100  LORBITE = 0
C
      CALL LEOPSU( LOPTRA )
      IF( LOPTRA .LE. 0 ) GOTO 9990
 102  CALL T3PLAV
C
C     REDUCTION DES FACES TRACEES
      REDUCF = PREDUF * 0.01
      REDUC1 = 1.0 - REDUCF
C
      IF( NDIMLI .EQ. 3 ) GOTO 300
C
C     ==================================================
C     SURFACES EN DIMENSION 2 SANS ALGORITHME DU PEINTRE
C     ==================================================
      IF( LORBITE .EQ. 0 ) GOTO 110
C     INITIALISATION DU ZOOM TRANSLATION
      CALL ZOOM2D0( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 100
      GOTO 110
C
C     ZOOM OU TRANSLATION ACTIFS
 105  CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 100
C
C     LE TRACE DES AXES EN 2D
 110  CALL TRAXE2
C
C     AU DELA LES TRAITS SONT CONTINUS ET DE NEPARF EPAISSEURS
      CALL XVTYPETRAIT( 0 )
      CALL XVEPAISSEUR( NEPARF )
C
C     LA BOUCLE SUR LES TRIANGLES-QUADRANGLES DES SURFACES 2D
      MNEF = MNTEF
C
      DO 200 NEF = 1, NBTEF
C
C        LES 3 COORDONNEES DES NCOGEL SOMMETS DE LA FACE
         IF( MCN(MNEF+3) .EQ. 0 ) THEN
C           TRIANGLE
            NCOGEL = 3
         ELSE
C           QUADRANGLE
            NCOGEL = 4
         ENDIF
C
C        REDUCTION DES FACES DE CET EF SANS TG
         DO 140 J=1,NCOGEL
            NS   = MCN(MNEF-1+J)
            MN   = MNTSOM + 3 * NS - 3
            MNBA = MNBAEF-3+3*NEF
            X(J) = RMCN(MN  ) * REDUC1 + RMCN(MNBA  ) * REDUCF
            Y(J) = RMCN(MN+1) * REDUC1 + RMCN(MNBA+1) * REDUCF
 140     CONTINUE
C
C        TRACE DE LA FACE DE NCOGEL ARETES ET DES ARETES
         IF( IAVFAC .EQ. 0 ) THEN
C           PAS DE TRACE DE L'INTERIEUR DE LA FACE
            NCF = -1
         ELSE
C           LA COULEUR DE LA FACE
            NCF = NCOUFA
            IF( LCRITR .GT. 0 ) THEN
C              TRACE DE LA COULEUR DE LA QUALITE DE L'ELEMENT FINI
               QUALIT = RMCN( MNQUEF-1+NEF )
               IF( QUALIT .LE. 0.1 ) THEN
                  NCF = NCROUG
               ELSE
                  NCF = N1COUL + 9 - INT( 10.0 * ( 1.0 - QUALIT ) )
               ENDIF
            ENDIF
         ENDIF
C
C        LA COULEUR DES ARETES DE L'EF
         IF( IAVARE .EQ. 0 ) THEN
            NCA = -1
         ELSE
            NCA = NCOUAF
         ENDIF
C
C        TRACE DE L'EF SANS SES EVENTUELLES TG
         CALL FACE2D( NCF, NCA, NCOGEL, X, Y )
C
C        TRACE EVENTUEL DU NO DE L'EF EN SON BARYCENTRE
         IF( IAVNEF .NE. 0 ) THEN
            WRITE( NMSOMM , '(I8)' ) NEF
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, RMCN(MNBAEF-3+3*NEF), NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO 190 J=1,NCOGEL
               NS = MCN(MNTEF-5+4*NEF+J)
               MN = MNTSOM - 3 + 3 * NS
               WRITE( NMSOMM , '(I8)' ) NS
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, RMCN(MN), NMSOMM(1:L) )
 190        CONTINUE
         ENDIF
C
C        PASSAGE A L'EF SUIVANT
         MNEF = MNEF + 4
 200  CONTINUE
C
C     LA MEMOIRE PIXELS EST COPIEE DANS LA FENETRE => VISIBLE
      CALL MEMPXFENETRE
      IF( LORBITE .NE. 0 ) GOTO 105
      GOTO 100
C
C     ==================================================
C     SURFACES EN DIMENSION 3 AVEC ALGORITHME DU PEINTRE
C     ==================================================
C     DISTANCE A L'OEIL DU BARYCENTRE DES NBTEF FACES
 300  CALL TNMCDC( 'REEL', NBTEF, MNDIST )
      MNDIS1 = MNDIST - 1
C     NUMERO AVANT TRI DE CHAQUE EF
      CALL TNMCDC( 'ENTIER', NBTEF, MNNOFA )
      MNNOF1 = MNNOFA - 1
C
ccc      IF( LORBITE .NE. 0 ) GOTO 330
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      CALL ORBITE0( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 100
      GOTO 330
C
C     ORBITE OU ZOOM OU TRANSLATION ACTIFS
 325  CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 100
C
C     CALCUL DE LA DISTANCE A L'OEIL DU BARYCENTRE DES FACES
 330  MNBA = MNBAEF
      DO 340 NEF = 1, NBTEF
C        AXONOMETRIE DU BARYCENTRE
         CALL XYZAXO( RMCN(MNBA), XYZB )
         RMCN( MNDIS1+NEF ) = -XYZB(3)
C        IDENTITE INITIALE POUR LE NUMERO DE FACE
         MCN( MNNOF1+NEF ) = NEF
         MNBA = MNBA + 3
 340  CONTINUE
C
C     LE TRI CROISSANT SELON LA DISTANCE (COTE Z ) A L'OEIL
C     =====================================================
      CALL TRITRP( NBTEF , RMCN(MNDIST) , MCN(MNNOFA) )
C
C     DISTANCES MIN ET MAX DES BARYCENTRES A L'OEIL
      DMIN = RMCN(MNDIST)
      DMAX = RMCN(MNDIS1+NBTEF)
C
C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
      EP1   = 0.05
      EP2   = SQRT(EP1)
      EP3   = SQRT(1+EP1) - EP2
      DELTA = DMAX - DMIN
      IF ( DELTA .EQ. 0 ) THEN
         DELTA = 1.0
         POID  = 1.0
      ELSE
         POID  = 0.6
      ENDIF
      POID1 = 1.0 - POID
C
C     TRACE DES FACES SELON L'ELOIGNEMENT ET LA DIRECTION DE VISEE
      DIREVI(1) = AXOEIL(1) - AXOPTV(1)
      DIREVI(2) = AXOEIL(2) - AXOPTV(2)
      DIREVI(3) = AXOEIL(3) - AXOPTV(3)
C     LA DIRECTION DE VISEE EST NORMALISEE A 1.
      CALL NORMER( 3, DIREVI, IERR )
C
      IF( LCRITR .EQ. -1 ) THEN
C        TRACE ARC EN CIEL DU MAILLAGE SELON LA COORDONNEE Z DES SOMMETS
C        CALCUL DES Z MIN ET MAX DES SOMMETS
         ZMIN = 1E28
         ZMAX =-1E28
         MN   = MNTSOM - 1
         DO 350 N=1,NBTSOM
            MN = MN + 3
            IF( ZMIN .GT. RMCN(MN) ) ZMIN=RMCN(MN)
            IF( ZMAX .LT. RMCN(MN) ) ZMAX=RMCN(MN)
 350     CONTINUE
      ENDIF
C
C     LE TRACE DES AXES EN 3D
      CALL TRAXE3
C
C     AU DELA LES TRAITS SONT CONTINUS ET DE NEPARF EPAISSEURS
      CALL XVTYPETRAIT( 0 )
      CALL XVEPAISSEUR( NEPARF )
C
      DO 400 NEF = NBTEF,1,-1
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF = MCN( MNNOF1 + NEF )
C
C        ADRESSE DE L'EF A TRACER
         MNEF = MNTEF - 4 + 4 * NF
C
         IF( MCN(MNEF+3) .EQ. 0 ) THEN
C           TRIANGLE
            NCOGEL = 3
         ELSE
C           QUADRANGLE
            NCOGEL = 4
         ENDIF
C
C        COULEUR D'ARETE PAR DEFAUT
         NCA = NCOUAF
C
         IF( IAVFAC .EQ. 0 ) THEN
C
C           PAS DE TRACE DE LA FACE
C           -----------------------
            NCF = -1
C
         ELSE IF( NBRCOU .EQ. 0 ) THEN
C
C           ECRAN NOIR ET BLANC
C           -------------------
            NCF = NCBLAN
            NCA = NCNOIR
C
         ELSE IF( LCRITR .EQ. 0 ) THEN
C
C           CALCUL DE LA COULEUR SELON DIRECTION DE VISEE ET ELOIGNEMENT
C           ------------------------------------------------------------
C           L'ELOIGNEMENT DE LA FACE
            R  = ( RMCN(MNDIS1+NF)  - DMIN ) / DELTA
            R2 = ( SQRT(ABS(R+EP1)) - EP2  ) / EP3
C
C           L'OTHOGONALITE A LA FACE
C           CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
            CALL NORF34( NCOGEL, MCN(MNEF), RMCN(MNTSOM),
     %                   CNORFA, IERR )
            IF( IERR .NE. 0 ) GOTO 400
C
C           LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
            R = PROSCR( DIREVI, CNORFA, 3 )
            R = 1.0 - ABS(R)
C
            IF( IAVELO .NE. 0 ) THEN
C              LA COULEUR PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
               R =  R * POID + POID1 * R2
            ENDIF
C
            NCF = NINT( NBCOUL * R + N1COUL )
C
         ELSE IF( LCRITR .GT. 0 ) THEN
C
C           COULEUR DE LA FACE SELON LA QUALITE DE CET EF
C           ---------------------------------------------
            QUALIT = RMCN( MNQUEF-1+NF )
            IF( QUALIT .LE. 0.1 ) THEN
               NCF = NCROUG
            ELSE
               NCF = N1COUL + 9 - INT( 10.0 * ( 1.0 - QUALIT ) )
            ENDIF
C
         ELSE
C
C           COULEUR DE LA FACE SELON UN DEGRADE DES 3 COULEURS AUX SOMMETS
C           --------------------------------------------------------------
            DO 370 J=1,NCOGEL
               MN = MNTSOM - 4 + 3 * MCN(MNEF-1+J)
C              LES 3 COORDONNEES DES NCOGEL SOMMETS
               DO 360 K=1,3
                  XYZ(K,J) = RMCN( MN + K )
 360           CONTINUE
C              LA COULEUR DES NCOGEL SOMMETS
               COUL(J) = (RMCN(MN+3)-ZMIN) / (ZMAX-ZMIN)
               IF( COUL(J) .LT. 0.0 ) COUL(J)=0.0
               IF( COUL(J) .GT. 1.0 ) COUL(J)=1.0
               COUL(J) = N1COUL + NBCOUL * COUL(J)
 370        CONTINUE
            IF( NCOGEL .EQ. 4 ) THEN
C              LE 5-EME SOMMET EST EN FAIT LE PREMIER
               XYZ(1,5) = XYZ(1,1)
               XYZ(2,5) = XYZ(2,1)
               XYZ(3,5) = XYZ(3,1)
               COUL(5)  = COUL(1)
            ENDIF
C
         ENDIF
C
C        LA COULEUR FINALE DES ARETES DE LA FACE
C        ---------------------------------------
         IF( IAVARE .EQ. 0 ) NCA = -1
C
C        SI 2 COULEURS NEGATIVES => RIEN A TRACER => RETOUR
         IF( NCF .LT. 0 .AND. NCA .LT. 0 ) GOTO 9990
C
C        TRACE DE LA FACE EN 3D
C        ======================
         IF( LCRITR .GE. 0 ) THEN
C
C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            MNBA = MNBAEF - 3 + 3 * NF
            DO 380 J=1,NCOGEL
               NS = MCN( MNEF - 1 + J)
               MN = MNTSOM - 3 + 3 * NS
               XYZ(1,J) = RMCN(MN  ) * REDUC1 + RMCN(MNBA  ) * REDUCF
               XYZ(2,J) = RMCN(MN+1) * REDUC1 + RMCN(MNBA+1) * REDUCF
               XYZ(3,J) = RMCN(MN+2) * REDUC1 + RMCN(MNBA+2) * REDUCF
 380        CONTINUE
C
            IF( LCRITR .EQ. 0 ) THEN
C              LES 3 COORDONNEES DE LA NORMALE A LA FACE PLANE
               CALL NORFA3( XYZ(1,1), XYZ(1,2), XYZ(1,NCOGEL),
     %                      CNORFA, IERR )
C              CALCUL DE LA COULEUR EN CE POINT A PARTIR
C              DE LA NORMALE ET DES ECLAIRAGES
               CALL ECLAIR( CNORFA, RCOUL )
               NCF = NINT( RCOUL )
            ENDIF
C
C           TRACE DE LA FACE P1 SANS TG
            CALL FACE3D( NCF, NCA, NCOGEL, XYZ )
         ELSE
C
C           LE TRACE EFFECTIF DU TRIANGLE DE SOMMETS 123
            CALL TRIACOUL3DBORD( XYZ(1,1), COUL(1), NCA, NCF )
            IF( NCOGEL .EQ. 4 ) THEN
C              LE TRACE EFFECTIF DU TRIANGLE DE SOMMETS 345=341
               CALL TRIACOUL3DBORD( XYZ(1,3), COUL(3), NCA, NCF )
            ENDIF
C
         ENDIF
C
C        TRACE EVENTUEL DU NO DE L'EF EN SON BARYCENTRE
C        ----------------------------------------------
         IF( IAVNEF .NE. 0 ) THEN
            WRITE( NMSOMM , '(I8)' ) MCN(MNNOF1+NF)
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, RMCN(MNBAEF-3+3*NF), NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
C        --------------------------------
         IF( IAVNSO .NE. 0 ) THEN
            DO 390 J=1,NCOGEL
               NS = MCN(MNTEF-5+4*NF+J)
               MN = MNTSOM - 3 + 3 * NS
               WRITE( NMSOMM , '(I8)' ) NS
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, RMCN(MN), NMSOMM(1:L) )
 390        CONTINUE
         ENDIF
C
C        TRACE EVENTUEL DU VECTEUR NORMAL A LA FACE
C        ------------------------------------------
         IF( IAVNRF .NE. 0 ) THEN
C
C           CALCUL DES COORDONNEES DU BARYCENTRE A LA FACE
            XYZP(1) = 0
            XYZP(2) = 0
            XYZP(3) = 0
            DO J=1,NCOGEL
               NS = MCN(MNEF-1+J)
               MN = MNTSOM - 3 + 3 * NS
C              LES 3 COORDONNEES DES NCOGEL SOMMETS
               DO K=1,3
                  XYZ(K,J) = RMCN( MN + K )
                  XYZP(K)  = RMCN( MN + K ) + XYZP(K)
               ENDDO
            ENDDO
            XYZP(1) = XYZP(1) / NCOGEL
            XYZP(2) = XYZP(2) / NCOGEL
            XYZP(3) = XYZP(3) / NCOGEL
C
C           CALCUL DU VECTEUR NORMAL A LA FACE
            CALL NORMTQ( NCOGEL, XYZ, CNORFA, IERR )
            IF( IERR .NE. 0 ) GOTO 400
C
C           TRACE DE LA FLECHE DE LONGUEUR 1CM
            DISTCM = 1.
            CALL T3FLEC( NCONRF, XYZP, DISTCM, CNORFA )
C
         ENDIF
C
 400  CONTINUE
C
cccC     LA MEMOIRE PIXELS EST COPIEE DANS LA FENETRE => VISIBLE
ccc      CALL MEMPXFENETRE
cccC
cccC     POUR VIDER LE BUFFER DE X11
ccc      CALL XVVOIR
      CALL TRFINS('SURFACES')
C
      IF( LORBITE .NE. 0 ) GOTO 325
      GOTO 100
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
C     =========================================
 9990 IF( NDIMLI .EQ. 3 ) THEN
         IF( MNDIST .GT. 0 ) CALL TNMCDS( 'REEL',   NBTEF, MNDIST )
         IF( MNNOFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTEF, MNNOFA )
      ENDIF
      IF( LCRITR .GT. 0 .AND. MNQUEF .GT. 0 ) THEN
         IF( MNQUEF .GT. 0 ) CALL TNMCDS( 'REEL', NBTEF, MNQUEF )
      ENDIF
      IF( MNBAEF .GT. 0 ) CALL TNMCDS( 'REEL', 3*NBTEF, MNBAEF )
      IF( MNTSOM .GT. 0 ) CALL TNMCDS( 'REEL',3*NBTSOM, MNTSOM )
      IF( MNTEF  .GT. 0 ) CALL TNMCDS( 'ENTIER',NBSOEF*NBTEF,MNTEF )
      END
