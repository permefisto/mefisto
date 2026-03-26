      SUBROUTINE TRANBTET( NOPASS, KTITRE, XYZPOI, NBEF, NBSOTE, NUSOTE,
     %                     NBTETRA, NOTETRA, XYZPT1, XYZPT2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DES ARETES ET NO DES TETRAEDRES D'UNE LISTE DE TETRAEDRES
C -----  ET LES 2 POINTS XYZPT1 et XYZPT2

C ENTREES:
C --------
C NOPASS : =0        CALCUL DE L'ISOFENETRE de VISEE
C          =1 PAS DE CALCUL DE L'ISOFENETRE de VISEE
C KTITRE : TITRE DU TRACE
C XYZPOI : X  Y  Z DES SOMMETS DE LA TETRAEDRISATION
C NBEF   : NOMBRE TOTAL DE TETRAEDRES DU MAILLAGE
C NBSOTE : NOMBRE DE POINTS DES TETRAEDRES (=4)
C NUSOTE : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,

C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NUSOTE DES TETRAEDRES A TRACER
C XYZPT1 2: LES 2 POINTS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2020
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES

      INTEGER           NUSOTE(NBEF,NBSOTE), NOTETRA(NBTETRA), NOSOTE(4)
      REAL              XYZPOI(3,*), XYZPT1(3), XYZPT2(3), R,
     %                  XYZBAR1(3), XYZBAR2(3)
      CHARACTER*(*)     KTITRE

      IF( .NOT. TRACTE ) RETURN

ccc      print *
ccc      DO N=1,NBTETRA
ccc         NT = NOTETRA(N)
ccc         IF( NT .GT. 0 ) THEN
ccc            print *,'trntetra: nusote(',NT,')=',(NUSOTE(NT,M),M=1,4)
ccc         ENDIF
ccc      ENDDO

      IF( NOPASS .NE. 0 ) GOTO 8

C     CADRE COOEXT RESTREINT AUX NBTETRA TETRAEDRES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      DO 5 N = 1, NBTETRA
         NT = ABS( NOTETRA(N) )
         DO K = 1, 4
            NS = NUSOTE( NT, K )
            IF( NS .LE. 0 ) GOTO 5
            DO L=1,3
               R = XYZPOI(L,NS)
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
         ENDDO
 5    ENDDO

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

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
 8    PREDU0 = PREDUF
      PREDUF = 20.0
      LORBITE= 1

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

C     TRACE DES TETRAEDRES
      DO 30 N = 1, NBTETRA

C        NUMERO DU TETRAEDRE DANS NUSOTE
         NT = ABS( NOTETRA(N) )
         IF( NT .GT. 0 ) THEN
            IF( NUSOTE(NT,1) .LE. 0 ) GOTO 30

C           LE TRACE DES 6 ARETES DU TETRAEDRE NUSOTE(*,NT)
            DO I=1,4
               NOSOTE( I ) = NUSOTE( NT, I )
            ENDDO

            IF( N .EQ. 1 ) THEN
C              POUR LE 1-ER TETRAEDRE TRACE EN VERT
               NCL = NCVERT
               CALL XVEPAISSEUR( 5 )
            ELSE
               IF( N .EQ. NBTETRA ) THEN
C                 POUR LE DERNIER TETRAEDRE TRACE EN CYAN
                  CALL XVEPAISSEUR( 3 )
                  NCL = NCCYAN
               ELSE
C                 POUR LES TETRAEDRES INTERMEDIAIRES TRACE EN GRIS
                  NCL = NCGRIS
                  CALL XVEPAISSEUR( 1 )
               ENDIF
            ENDIF
            CALL TRATETRA( NCL, NOSOTE, XYZPOI )

C           LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
            DO I = 1, 4
               NS = NOSOTE( I )
               CALL ENTIER3D( NCORAN, XYZPOI(1,NS), NS )
            ENDDO

C           LE TRACE DU NUMERO DU TETRAEDRE NT EN SON BARYCENTRE XYZBAR2
            DO K=1,3
               XYZBAR2( K ) = 0
            ENDDO
            DO I=1,4
               NS = NOSOTE( I )
               DO K=1,3
                  XYZBAR2(K) = XYZBAR2(K) + XYZPOI(K,NS)
               ENDDO
            ENDDO
            DO K=1,3
               XYZBAR2( K ) = XYZBAR2( K ) / 4
            ENDDO
            CALL XVEPAISSEUR( 1 )
            IF( N .GT. 1 ) THEN
C              TRACE DU PARCOURS XYZBAR1->XYZBAR2
               CALL TRAIT3D( NCTURQ, XYZBAR1, XYZBAR2 )
            ENDIF
            CALL ENTIER3D( NCL, XYZBAR2, NT )

            DO K=1,3
               XYZBAR1( K ) = XYZBAR2( K )
            ENDDO

C           TRACE DU SEGMENT XYZPT1-XYZPT2
            CALL TRAIT3D( NCROSE, XYZPT1, XYZPT2 )

C           TRACE DU POINT XYZPT1
            CALL SYMBOLE3D( NCMAGE, XYZPT1, '.xyz1' )

C           TRACE DU POINT XYZPT2
            CALL SYMBOLE3D( NCROUG, XYZPT2, '+XYZ2' )

         ENDIF

 30   ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0
      CALL XVEPAISSEUR( 0 )
      RETURN
      END
