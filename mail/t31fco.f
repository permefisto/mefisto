      SUBROUTINE T31FCO( NMSURF, NUSURF, MNNSEF, MNXYZS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FACES D'UNE SURFACE MAILLEE C0 OU C1
C -----    SELON LA TECHNIQUE DU PEINTRE
C          ET LES COULEURS SELON LA VALEUR DE LCRITR (cf incl/mecoit.inc)
C          LCRITR=-2 => ARC EN CIEL SELON LA VALEUR DE Z
C                       SANS CALCUL DE ZMIN ET ZMAX POUR LA COULEUR
C                =-1 => ARC EN CIEL SELON LA VALEUR DE Z
C                       AVEC CALCUL DE ZMIN ET ZMAX POUR LA COULEUR
C                = 0 => ELOIGNEMENT ET DIRECTION DE VISEE
C                =+1 => QUALITE GEOMETRIQUE DES ELEMENTS FINIS

C ENTREES:
C --------
C NMSURF : NOM DE LA SURFACE A TRACER
C NUSURF : NUMERO DE LA SURFACE DANS SON LEXIQUE
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF'       NO DES SOMMETS DES EF
C MNXYZS : ADRESSE MCN DU TABLEAU 'XYZSOMMET'  COORDONNEES DES SOMMETS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS    JUIN 1996
C MODIF : PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS OCTOBRE 2003
C MODIF : PERRONNET ALAIN Saint PIERRE du PERRAY              Avril 2020
C2345X...............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"

C     DECLARATION DU SUPER-TABLEAU NUMERIQUE MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)

C     LES PARAMETRES
      CHARACTER*(*)     NMSURF

C     LES VARIABLES LOCALES
      CHARACTER*8       NMSOMM
      INTEGER           NOSOEL(1:64)
      REAL              XYZP(3), CNORFA(3), COUL(5), XYZ(3,5)
      REAL              X(12), Y(12), Z(12), COORNO(3)

      IF( MNNSEF .LE. 0 .OR. MNXYZS .LE. 0 ) RETURN

C     SAUVEGARDE
      NCOUAF0 = NCOUAF
      NCOUFA0 = NCOUFA

C     LES PARAMETRES DES NO SOMMET DES EF DU MAILLAGE
C     -----------------------------------------------
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX,     NY,     NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNST  = MNXYZS + WYZSOM - 3
      MNST1 = MNST - 1

C     LE NOMBRE DE SOMMETS DU MAILLAGE
      NBSOM = MCN(MNXYZS+WNBSOM)

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DE LA 1-ERE TANGENTE DU TMS 'XYZSOMMET'
      MNTG = MNST + 3 * NBSOM

C     LE NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL

C     INITIALISATION DE LA DIRECTION DE VISEE PTV-OEIL NORMALISEE A 1
C     ---------------------------------------------------------------
      DIREVI(1) = AXOEIL(1) - AXOPTV(1)
      DIREVI(2) = AXOEIL(2) - AXOPTV(2)
      DIREVI(3) = AXOEIL(3) - AXOPTV(3)
      CALL NORMER( 3, DIREVI, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'OEIL et POINT VU sont IDENTIQUES. VISEE IMPOSSIBLE'
         CALL LEREUR
         RETURN
      ENDIF

C     CONSTRUIRE LES TABLEAUX DES DISTANCES A L'OEIL ET
C     NUMERO DES FACES TRIEES SELON LEUR ELOIGNEMENT A L'OEIL DE VISEE
C     et NUMERO DE LA COULEUR DE TRACE DES FACES
C     ================================================================
      CALL TNMCDC( 'REEL'  , NBEFOB, MNDIST )
      MNDIST1 = MNDIST - 1
      CALL TNMCDC( 'ENTIER', NBEFOB, MNNOFT )
      MNNOFT1 = MNNOFT - 1

C     LE CALCUL DE LA DISTANCE DU BARYCENTRE DES FACES A L'OEIL
C     ---------------------------------------------------------
      CALL DISTFACE( NBSOM, RMCN(MNXYZS+WYZSOM),
     %               NBEFOB, MCN(MNNSEF+WUSOEF),
     %               RMCN(MNDIST), MCN(MNNOFT) )

C     DISTANCES MIN ET MAX des BARYCENTRES des FACES a l'OEIL
      DMIN = RMCN( MNDIST )
      DMAX = RMCN( MNDIST-1+NBEFOB )

C     CALCUL DES COTES EXTREMES DES NBSOM SOMMETS DU MAILLAGE
C     -------------------------------------------------------
      ZMIN = 0
      ZMAX = 1

      IF( LCRITR .EQ. -1 ) THEN

C        TRACE ARC EN CIEL DU MAILLAGE SELON LA COORDONNEE Z DES SOMMETS
C        CALCUL DES Z MIN ET MAX DES SOMMETS
         ZMIN = 1E28
         ZMAX =-1E28
         MN = MNXYZS+WYZSOM
         DO N=1,NBSOM
            ZN = RMCN(MN+2)
            IF( ZN .LT. ZMIN ) ZMIN = ZN
            IF( ZN .GT. ZMAX ) ZMAX = ZN
            MN = MN + 3
         ENDDO

      ELSE IF( LCRITR .EQ. -2 ) THEN

C        TRACE ARC EN CIEL DU MAILLAGE SELON LA COORDONNEE Z DES SOMMETS
C        AVEC ZMIN ET ZMAX GLOBAUX POUR PLUSIEURS TRACES
         ZMIN = ZZZMIN
         ZMAX = ZZZMAX

      ENDIF

C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
C     ----------------------------------------------------
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
      POID1 = 1. - POID

C     LE TRACE DES FACES SELON L'ELOIGNEMENT ET LA DIRECTION DE VISEE
C     ===============================================================
      MNSTS = MNXYZS + WYZSOM
      DO 100 N=NBEFOB,1,-1

C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF = MCN( MNNOFT1 + N )

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( NF    , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( NOSOEL(1) .EQ. 0 ) GOTO 100

C        LE NOMBRE DE SOMMETS DE CET ELEMENT FINI EST NCOGEL=3 OU 4
         IF( NBRCOU .EQ. 0 ) THEN

C           ECRAN NOIR ET BLANC
C           -------------------
            NCOUAF = NCBLAN
            NCOUFA = NCNOIR
C
         ELSE IF( LCRITR .EQ. 0 ) THEN

C           CALCUL DE LA COULEUR SELON DIRECTION DE VISEE ET ELOIGNEMENT
C           ------------------------------------------------------------
C           L'ELOIGNEMENT DE LA FACE
            R  = ( RMCN(MNDIST1+N) - DMIN ) / DELTA
            R2 = ( SQRT(R+EP1)     - EP2  ) / EP3

C           L'OTHOGONALITE A LA FACE
C           CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
            CALL NORF34( NCOGEL, NOSOEL, RMCN(MNSTS),
     %                   CNORFA, IERR )
            IF( IERR .NE. 0 ) GOTO 100
C
C           LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
            R = PROSCR( DIREVI, CNORFA, 3 )
            R = 1. - ABS(R)

            IF( IAVELO .NE. 0 ) THEN
C              LA COULEUR PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
               R =  R * POID + POID1 * R2
            ENDIF

            NCOUFA = NINT( NBCOUL * R + N1COUL )

         ELSE IF( LCRITR .GT. 0 ) THEN
C
C           COULEUR DE LA FACE SELON LA QUALITE DE CET EF
C           ---------------------------------------------
C           CALCUL DE LA QUALITE DE L'ELEMENT FINI
            CALL QUALEF( NCOGEL,   NOSOEL, NBSOM, RMCN(MNSTS),
     %                   SURFVOLU, QUALIT, IERR )
            IF ( IERR .NE. 0 ) GOTO 9900
C           LA COULEUR VISUALISE LA QUALITE
            NCOUFA = NCOQUA( QUALIT )

         ELSE

C           COULEUR DE LA FACE SELON UN DEGRADE DES 3 COULEURS AUX SOMMETS
C           --------------------------------------------------------------
            DO J=1,NCOGEL
               MN = MNST1 + 3 * NOSOEL(J)
C              LES 3 COORDONNEES DES NCOGEL SOMMETS
               DO K=1,3
                  XYZ(K,J) = RMCN( MN + K )
               ENDDO
C              LA COULEUR DES NCOGEL SOMMETS
               COUL(J) = (RMCN(MN+3)-ZMIN) / (ZMAX-ZMIN)
               IF( COUL(J) .LT. 0.0 ) COUL(J)=0.0
               IF( COUL(J) .GT. 1.0 ) COUL(J)=1.0
               COUL(J) = N1COUL + NBCOUL * COUL(J)
            ENDDO
            IF( NCOGEL .EQ. 4 ) THEN
C              LE 5-EME SOMMET EST EN FAIT LE PREMIER
               XYZ(1,5) = XYZ(1,1)
               XYZ(2,5) = XYZ(2,1)
               XYZ(3,5) = XYZ(3,1)
               COUL(5)  = COUL(1)
            ENDIF

         ENDIF

C        LA COULEUR FINALE DE LA FACE
C        ----------------------------
         IF( IAVFAC .EQ. 0 ) THEN
            NCF = -1
         ELSE
            NCF = NCOUFA
         ENDIF

C        LA COULEUR FINALE DES ARETES DE LA FACE
C        ---------------------------------------
         IF( IAVARE .EQ. 0 ) THEN
            NCA = -1
         ELSE
            NCA = NCOUAF
         ENDIF

C        TRACE DE LA FACE EN 3D
C        ======================
         IF( LCRITR .GE. 0 ) THEN
            CALL TRFA3D( NCF, NCA, PREDUF, MNST, MNTG,
     %                   NCOGEL, NOSOEL, NUEFTG,
     %                   X, Y, Z )
         ELSE
C           LE TRACE EFFECTIF DU TRIANGLE DE SOMMETS 123
            CALL TRIACOUL3DBORD( XYZ(1,1), COUL(1), NCA, NCF )
            IF( NCOGEL .EQ. 4 ) THEN
C              LE TRACE EFFECTIF DU TRIANGLE DE SOMMETS 345=341
               CALL TRIACOUL3DBORD( XYZ(1,3), COUL(3), NCA, NCF )
            ENDIF
         ENDIF

C        TRACE EVENTUEL DU NO DE L'EF
C        ----------------------------
         IF( IAVNEF .NE. 0 ) THEN
C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAST( NCOGEL, X, Y, Z, XYZP )
            WRITE( NMSOMM, '(I8)' ) MCN(MNNOFT1+N)
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
         ENDIF

C        TRACE EVENTUEL DU NO DES SOMMETS DE L'EF
C        ----------------------------------------
         IF( IAVNSO .NE. 0 ) THEN
            DO J=1,NCOGEL
               CALL TRST3D( NCONSO, NOSOEL(J), RMCN(MNSTS) )
ccc               MN = MNST + 3 * NOSOEL(J)
ccc               WRITE( NMSOMM, '(I8)' ) NOSOEL(J)
ccc               CALL SANSBL( NMSOMM, L )
ccc               CALL TEXTE3D( NCONSO, RMCN(MN), NMSOMM(1:L) )
            ENDDO
         ENDIF

C        TRACE EVENTUEL DU VECTEUR NORMAL A LA FACE
C        ------------------------------------------
         IF( IAVNRF .NE. 0 ) THEN

C           CALCUL DES COORDONNEES DU BARYCENTRE A LA FACE
            XYZP(1) = 0
            XYZP(2) = 0
            XYZP(3) = 0
            DO J=1,NCOGEL
               MN = MNST1 + 3 * NOSOEL(J)
C              LES 3 COORDONNEES DES NCOGEL SOMMETS
               DO K=1,3
                  XYZ(K,J) = RMCN( MN + K )
                  XYZP(K)  = RMCN( MN + K ) + XYZP(K)
               ENDDO
            ENDDO
            XYZP(1) = XYZP(1) / NCOGEL
            XYZP(2) = XYZP(2) / NCOGEL
            XYZP(3) = XYZP(3) / NCOGEL

C           CALCUL DU VECTEUR NORMAL A LA FACE
            CALL NORMTQ( NCOGEL, XYZ, COORNO, IERR )
            IF( IERR .NE. 0 ) GOTO 100

C           TRACE DE LA FLECHE DE LONGUEUR 1CM
            DISTCM = 1.
            CALL T3FLEC( NCONRF, XYZP, DISTCM, COORNO )

         ENDIF

 100  ENDDO

C     TRACE DE LA POIGNEE ET DU NOM DE LA SURFACE
C     -------------------------------------------
C     LES 3 COORDONNEES DE LA POIGNEE
C     CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA DERNIERE FACE
      CALL COBAST( NCOGEL, X, Y, Z, XYZP )
      CALL ITEMS3( XYZP, NMSURF, NUSURF )

C     DESTRUCTION DU TABLEAU DES DISTANCES et NUMERO DES FACES
 9900 CALL TNMCDS( 'REEL'  , NBEFOB, MNDIST )
      CALL TNMCDS( 'ENTIER', NBEFOB, MNNOFT )

C     RESTAURATION
      NCOUAF = NCOUAF0
      NCOUFA = NCOUFA0

      RETURN
      END
