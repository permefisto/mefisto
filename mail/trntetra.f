      SUBROUTINE TRNTETRA( NOPASS,  KTITRE, XYZSOM, NBSOTE, NBTETRA,
     %                     NOTETRA, NSTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACE DES ARETES ET NO DES 4 SOMMETS D'UNE LISTE DE TETRAEDRES
C -----
C ENTREES:
C --------
C NOPASS : =0        CALCUL DE L'ISOFENETRE
C          =1 PAS DE CALCUL DE L'ISOFENETRE
C KTITRE : TITRE DU TRACE
C XYZSOM : X  Y  Z DES SOMMETS DE LA TETRAEDRISATION
C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NSTETR DES TETRAEDRES A TRACER
C NSTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  Janvier 2016
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES

      INTEGER           NSTETR(NBSOTE,*), NOTETRA(NBTETRA)
      REAL              XYZSOM(3,*), R
      CHARACTER*(*)     KTITRE

      IF( .NOT. TRACTE ) RETURN

ccc      print *
ccc      DO N=1,NBTETRA
ccc         NT = NOTETRA(N)
ccc         IF( NT .GT. 0 ) THEN
ccc            print *,'trntetra: nstetr(',NT,')=',(NSTETR(M,NT),M=1,8)
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
            NS = NSTETR( K, NT )
            IF( NS .LE. 0 ) GOTO 5
            DO L=1,3
               R = XYZSOM(L,NS)
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
C
C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
 8    PREDU0 = PREDUF
      PREDUF = 20.0

      LORBITE = 1
      IF( LORBITE .EQ. 0 ) GOTO 20
C
C     INITIALISATION DE L'ORBITE
C     ==========================
      CALL ORBITE0( NOTYEV )
      GOTO 20
C
C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
 10   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C     TRACE DES AXES 3D
 20   CALL TRAXE3
C
C     TRACE DES TETRAEDRES
      DO 30 N = NBTETRA, 1, -1

C        NUMERO DU TETRAEDRE DANS NSTETR
         NT = ABS( NOTETRA(N) )
         IF( NT .GT. 0 ) THEN
            IF( NSTETR(1,NT) .LE. 0 ) GOTO 30

C           LE TRACE DES 6 ARETES DU TETRAEDRE NSTETR(*,NT)
            CALL TRATETRA( NCORAN, NSTETR(1,NT), XYZSOM )
C
C           LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
            DO I = 1, 4
               NS1 = NSTETR(I,NT)
               CALL ENTIER3D( NCNOIR, XYZSOM(1,NS1), NS1 )
            ENDDO

         ENDIF

 30   ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )
C
C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10
C
 9000 PREDUF = PREDU0
      RETURN
      END
