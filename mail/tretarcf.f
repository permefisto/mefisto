      SUBROUTINE TRETARCF( KTITRE, PTXYZD, NOTETR,
     %                     NBTETR, NUTETR, NBSTCF, NOSTCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DE LA LISTE NUTETR DES TETRAEDRES DE NOTETR
C -----
C ENTREES:
C --------
C KTITRE : TITRE DU TRACE
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NBTETR : NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NUTETR : NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C NBSTCF : NOMBRE DE SOMMETS DU CF DES FACES PERDUES
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS DU CF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray          Novembre 2018
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
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*), NUTETR(NBTETR), NOSTCF(NBSTCF)
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
         IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN
            DO K=1,4
C              NUMERO DU SOMMET K DU TETRAEDRE NT
               NS = NOTETR(K,NT)
               DO L=1,3
                  XYZ(L) = REAL( PTXYZD(L,NS) )
C                 LE MINIMUM
                  COOEXT(L,1) = MIN( COOEXT(L,1), XYZ(L) )
C                 LE MAXIMUM
                  COOEXT(L,2) = MAX( COOEXT(L,2), XYZ(L) )
               ENDDO
            ENDDO
         ENDIF
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
      PREDUF  = 10.0
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
      CALL XVEPAISSEUR( 0 )
      DO N=1,NBTETR
         NT = NUTETR(N)
         IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN

cccC           LE TRACE DES 6 ARETES DU TETRAEDRE OPPOSE A UNE FACE DE NT
ccc            NC     = NCVERT
ccc            PREDUF = 20.0
ccc            DO K=1,4
ccc               NTEOP = NOTETR( 4+K, NT )
ccc               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP).GT.0 ) THEN
ccc                  CALL TRTETRA( NC, NOTETR(1,NTEOP), PTXYZD )
ccc               ENDIF
ccc            ENDDO

C           LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
            PREDUF = 10.0
            NC     = NCNOIR
            CALL TRTETRA( NC, NOTETR(1,NT), PTXYZD )

C           LE TRACE DU NUMERO DES SOMMETS DU TETRAEDRE NT
            DO I = 1, 4
               NS = NOTETR(I,NT)
               XYZ(1) = REAL( PTXYZD(1,NS) )
               XYZ(2) = REAL( PTXYZD(2,NS) )
               XYZ(3) = REAL( PTXYZD(3,NS) )
               CALL ENTIER3D( NCBLEU, XYZ, NS )
            ENDDO

         ENDIF
      ENDDO

C     LE TRACE DES NBSTCF ARETES DU CF DES FACES PERDUES
      IF( NBSTCF .GT. 0 ) THEN
         CALL XVEPAISSEUR( 3 )
         NS1 = NOSTCF( NBSTCF )
         DO N = 1, NBSTCF
            XYZ(1) = REAL( PTXYZD(1,NS1) )
            XYZ(2) = REAL( PTXYZD(2,NS1) )
            XYZ(3) = REAL( PTXYZD(3,NS1) )

            NS2 = NOSTCF( N )
            XYZ2(1) = REAL( PTXYZD(1,NS2) )
            XYZ2(2) = REAL( PTXYZD(2,NS2) )
            XYZ2(3) = REAL( PTXYZD(3,NS2) )

C           TRACE DE L'ARETE NS1-NS2
            CALL TRAIT3D(  NCROUG, XYZ,  XYZ2 )
            CALL ENTIER3D( NCMAGE, XYZ2, NS2  )
            NS1 = NS2
         ENDDO
      ENDIF

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0
      CALL XVEPAISSEUR( 0 )
      RETURN
      END
