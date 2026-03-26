      SUBROUTINE TR1TEVV( NOPASS, KTITRE, XYZSOM,  NT0, NBSOTE, NSTETR,
     %                    NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLUMET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACE DES ARETES ET NO DES 4 SOMMETS DU TETRAEDRE NT0 DE NSTETR
C -----   DES TETRAEDRES OPPOSES ET DES TETRAEDRES OPPOSES OPPOSES

C ENTREES:
C --------
C NOPASS : 0 FORMER ET TRACER LA LISTE DES TETRAEDRES
C          1 TRACER LA LISTE DES TETRAEDRES
C KTITRE : TITRE A TRACER
C XYZSOM : X  Y  Z DES SOMMETS DE LA TETRAEDRISATION
C NT0    : NUMERO DANS NSTETR DU TETRAEDRE A TRAITER
C NBSOTE : NOMBRE DE SOMMETS DECLARABLES PAR TETRAEDRE DANS NSTETR(>3)
C NBSOTE : NOMBRE DE SOMMETS DECLARABLES PAR TETRAEDRE DANS NSTETR(>3)
C NSTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DE CHAQUE SOMMET
C NOTESO : NUMERO DU TETRAEDRE ET SUIVANT DES SOMMETS

C SORTIES:
C --------
C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NSTETR DES TETRAEDRES A TRACER
C VOLUMET: VOLUME DES NBTETRA TETRAEDRES
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
      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NTOP(4),NSOP(4),
     %                  NOTETRA(*)
      REAL              XYZSOM(3,*), XYZ(3), VOLUMET
      DOUBLE PRECISION  VOLTER, V, VOLUME
      CHARACTER*(*)     KTITRE
      SAVE              NBTETRA0

      VOLUMET = 0
      VOLUME  = 0D0
      IF( .NOT. TRACTE ) RETURN

      IF( NOPASS .NE. 0 ) GOTO 8

      CALL EFFACEMEMPX

C     RECHERCHE DES TETRAEDRES OPPOSES AU TETRAEDRE NT0
      CALL NOTSOPTE( NT0,  NO1TSO, NOTESO, NBSOTE, NSTETR,
     %               NFR, NTOP,   NSOP )

C     LE PREMIER TETRAEDRE A TRACER
      NBTETRA = 1
      NOTETRA( 1 ) = NT0

C     LES TETRAEDRES OPPOSES A NT0
      DO 3 N = 1, 4
         NTO = NTOP(N)
         IF( NTO .GT. 0 ) THEN
            DO K=1,NBTETRA
               IF( NOTETRA(K) .EQ. NTO ) GOTO 3
            ENDDO
            NBTETRA = NBTETRA + 1
            NOTETRA( NBTETRA ) = NTO
         ENDIF
 3    ENDDO
      NBTETRA0 = NBTETRA

C     LES TETRAEDRES OPPOSES AUX TETRAEDRES OPPOSES A NT0
      DO NN = 1, NBTETRA0
         NT1 = NOTETRA(NN)
C        RECHERCHE DES TETRAEDRES OPPOSES AU TETRAEDRE NT1
         CALL NOTSOPTE( NT1, NO1TSO, NOTESO, NBSOTE, NSTETR,
     %                  NFR, NTOP,   NSOP )
         DO 6 N = 1, 4
            NTO = NTOP(N)
            IF( NTO .GT. 0 ) THEN
               DO K=1,NBTETRA
                  IF( NOTETRA(K) .EQ. NTO ) GOTO 6
               ENDDO
               NBTETRA = NBTETRA + 1
               NOTETRA( NBTETRA ) = NTO
            ENDIF
 6       ENDDO
      ENDDO

ccc      print *
ccc      DO N=1,NBTETRA
ccc         NT = NOTETRA(N)
ccc         IF( NT .GT. 0 ) THEN
ccc            print *,'tr1tevv: nstetr(',NT,')=',(NSTETR(M,NT),M=1,4)
ccc         ENDIF
ccc      ENDDO

C     CADRE COOEXT RESTREINT AUX NBTETRA TETRAEDRES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      DO 5 N = 1, NBTETRA

         NT1 = ABS( NOTETRA(N) )
         IF( NT1           .EQ. 0 ) GOTO 5
         IF( NSTETR(1,NT1) .LE. 0 ) GOTO 5

C        VOLUME DU TETRAEDRE NT
         V = VOLTER( XYZSOM(1,NSTETR(1,NT1)),
     %               XYZSOM(1,NSTETR(2,NT1)),
     %               XYZSOM(1,NSTETR(3,NT1)),
     %               XYZSOM(1,NSTETR(4,NT1)) )
         VOLUME = VOLUME + ABS( V )

         DO K = 1, 4
            NS = NSTETR( K, NT1 )
            DO L=1,3
               XYZ(L) = XYZSOM(L,NS)
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), XYZ(L) )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), XYZ(L) )
            ENDDO
         ENDDO

 5    ENDDO

      AXOLAR = 0
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO
      AXOLAR = AXOLAR * 0.5  * 0.6
      AXOHAU = AXOLAR * 0.75 * 0.6
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
         NT1 = ABS( NOTETRA(N) )

         IF( NT1 .GT. 0 ) THEN

            IF( NSTETR(1,NT1) .LE. 0 ) GOTO 30

C           LE TRACE DES 6 ARETES DU TETRAEDRE NSTETR(*,NT)
            IF( N .LE. 1 ) THEN
               NCA = NCROUG
               CALL XVEPAISSEUR( 4 )
            ELSE
               IF( N .LE. NBTETRA0 ) THEN
                  NCA = NCBLEU
                  CALL XVEPAISSEUR( 2 )
               ELSE
                  NCA = NCVERT
                  CALL XVEPAISSEUR( 0 )
               ENDIF
            ENDIF
            CALL TRATETRA( NCA, NSTETR(1,NT1), XYZSOM )
C
C           LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
            DO I = 1, 4
               NS1 = NSTETR(I,NT1)
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

C     PASSAGE DOUBLE PRECISION EN SIMPLE PRECISION
 9000 VOLUMET = REAL( VOLUME )

      PREDUF = PREDU0
      RETURN
      END
