      SUBROUTINE TRETCE( KTITRE, PTXYZD, NOTETR, NUPOIN,
     %                   NBTETR, NUTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DE LA LISTE NUTETR DES TETRAEDRES DE NOTETR
C -----    DU CENTRE ET RAYON DE LA BOULE CIRCONSCRITE DES TETRAEDRES

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUPOIN : NUMERO DE PTXYZD DU POINT A TRACER
C NBTETR : NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NUTETR : NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS St PIERRE du PERRAYJanvier 2017
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*), CENTRE(4), V
      INTEGER           NOTETR(8,*)
      INTEGER           NUTETR(NBTETR)
      REAL              XYZ(3), XYZ2(3)

      IF( .NOT. TRACTE .OR. NBTETR .LE. 0 ) RETURN

C     CADRE COOEXT RESTREINT AUX NBTETR ETRAEDRES
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      DO N=1,NBTETR
         NT = NUTETR(N)
         DO K=1,4
C           NUMERO DU SOMMET K DU TETRAEDRE NT
            NS = NOTETR(K,NT)
            DO L=1,3
               XYZ(L) = REAL( PTXYZD(L,NS) )
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), XYZ(L) )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), XYZ(L) )
            ENDDO
         ENDDO
      ENDDO

      AXOLAR = 0
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO
      AXOLAR = AXOLAR * 0.5
      AXOHAU = AXOLAR * 0.75
C     PAS DE PLAN ARRIERE ET AVANT
      AXOARR = 0
      AXOAVA = 0
      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )

      DISMOY = AXOLAR/30

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 5.0
      LORBITE = 1

      IF( LORBITE .EQ. 0 ) GOTO 20

C     INITIALISATION DE L'ORBITE
C     ==========================
      CALL ORBITE0( NOTYEV )
      GOTO 20

C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
 10   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000

C     TRACE DES AXES 3D
 20   CALL TRAXE3

      DO N=1,NBTETR

C        NUMERO DU TETRAEDRE
         NT = NUTETR( N )

         IF( NT .GT. 0 ) THEN

C           TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
            CALL TRTETRA( NCBLEU, NOTETR(1,NT), PTXYZD )

C           CENTRE ET RAYON DE LA BOULE CIRCONSCRITE
            CALL CEBOQU( 0., NOTETR(1,NT), PTXYZD, CENTRE, V, Q )

            XYZ(1) = REAL( CENTRE(1) )
            XYZ(2) = REAL( CENTRE(2) )
            XYZ(3) = REAL( CENTRE(3) )
ccc            RAYON  = REAL( SQRT( CENTRE(4) ) )

            DO I = 1, 4
               NS = NOTETR(I,NT)
               XYZ2(1) = REAL( PTXYZD(1,NS) )
               XYZ2(2) = REAL( PTXYZD(2,NS) )
               XYZ2(3) = REAL( PTXYZD(3,NS) )
C              TRACE DU SEGMENT CENTRE-SOMMET I DE NT
               CALL TRAIT3D(  NCORAN, XYZ, XYZ2 )
C              TRACE DU NUMERO DU SOMMET I
               CALL ENTIER3D( NCMAGE, XYZ2, NS )
            ENDDO

C           TRACE DU CENTRE DE LA BOULE CIRCONSCRITE
            CALL ENTIER3D(  NCVERT, XYZ, NT )
            CALL SYMBOLE3D( NCNOIR, XYZ, '@')

         ELSE

            PRINT*,'tretce: NUPOIN=',NUPOIN,' NT=',NT,
     %             ' NBTETR=', NBTETR
            PRINT*,'NUTETR=',(NUTETR(I),I=1,NBTETR)

         ENDIF

      ENDDO

C     LE TRACE DU POINT NUPOIN DE PXYZD
      XYZ(1) = REAL( PTXYZD(1,NUPOIN) )
      XYZ(2) = REAL( PTXYZD(2,NUPOIN) )
      XYZ(3) = REAL( PTXYZD(3,NUPOIN) )
      CALL SYMBOLE3D( NCNOIR, XYZ, '#' )
      CALL SYMBOLE3D( NCNOIR, XYZ, '*' )
      CALL ENTIER3D(  NCROUG, XYZ, NUPOIN )

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
