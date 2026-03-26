      SUBROUTINE TRETOPT( KTITRE, XYZSOM, POINT, N1FEOC, NFETOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACE DES ARETES DES TETRAEDRES FORMES PAR UNE FACE DE L'ETOILE
C -----   ET UN POINT

C ENTREES:
C --------
C KTITRE : TITRE A TRACER
C XYZSOM : X  Y  Z DES SOMMETS DE LA TETRAEDRISATION
C POINT  : X Y Z DU POINT A RELIER AUX FACES DE L'ETOILE

C N1FEOC : NUMERO DANS NFETOI DE LA PREMIERE FACE OCCUPEE
C NFETOI : NUMERO DES 3 SOMMETS, TETRAEDRE ET SUIVANTE DES FACES
C          VUES UNE FOIS DE L'ETOILE CHAINEES N1FEOC PUIS NFETOI(5,*)
C NBSOTE : NOMBRE DE SOMMETS DECLARABLES PAR TETRAEDRE DANS NSTETR(>3)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  Fevrier 2016
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

      INTEGER           NFETOI(5,*), NOFEVN(256)
      REAL              XYZSOM(3,*), POINT(3), R
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  VOLTER, V

      IF( .NOT. TRACTE ) RETURN

      CALL EFFACEMEMPX

C     CADRE COOEXT RESTREINT AUX FACES DE L'ETOILE
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NF = N1FEOC
 5    IF( NF .GT. 0 ) THEN
         DO K = 1, 3
            NS = NFETOI( K, NF )
            DO L=1,3
               R = XYZSOM(L,NS)
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
         ENDDO
C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI(5,NF)
         GOTO 5
      ENDIF

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
      PREDU0 = PREDUF
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

C     TRACE DES TETRAEDRES FACE ETOILE + POINT
      NBTEVN = 0
      NF = N1FEOC
 30   IF( NF .GT. 0 ) THEN

C        VOLUME DU TETRAEDRE FACE + POINT
         V = VOLTER( XYZSOM( 1, NFETOI(1,NF) ),
     %               XYZSOM( 1, NFETOI(2,NF) ),
     %               XYZSOM( 1, NFETOI(3,NF) ),
     %               POINT )
         IF( V .LE. 0D0 ) THEN
C           SI VOLUME<0 STOCKAGE DU NO DE FACE
            NBTEVN = NBTEVN + 1
            NOFEVN( NBTEVN ) = NF
         ENDIF

         NS2 = NFETOI( 3, NF )
         DO K = 1, 3
            NS1 = NFETOI( K, NF )

C           TRACE SOMMET1-SOMMET2  ARETE K FACE NF
            CALL TRAIT3D(  NCBLEU, XYZSOM(1,NS2), XYZSOM(1,NS1) )
            CALL ENTIER3D( NCNOIR, XYZSOM(1,NS1), NS1 )

C           TRACE SOMMET K FACE NF -> POINT
            CALL TRAIT3D( NCVERT, XYZSOM(1,NS1), POINT )

            NS2 = NS1
         ENDDO

C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI(5,NF)
         GOTO 30

      ENDIF

C     TRACE EN ROUGE DES TETRAEDRES DE VOLUME<0
      DO L = 1, NBTEVN

C        NO NFETOI DE LA FACE
         NF = NOFEVN( L )

         NS2 = NFETOI( 3, NF )
         DO K = 1, 3
            NS1 = NFETOI( K, NF )

C           TRACE SOMMET1-SOMMET2  ARETE K FACE NF
            CALL TRAIT3D(  NCROUG, XYZSOM(1,NS2), XYZSOM(1,NS1) )
            CALL ENTIER3D( NCNOIR, XYZSOM(1,NS1), NS1 )

C           TRACE SOMMET K FACE NF -> POINT
            CALL TRAIT3D( NCORAN, XYZSOM(1,NS1), POINT )

            NS2 = NS1
         ENDDO

      ENDDO

C     TRACE DU POINT
      CALL SYMBOLE3D( NCNOIR, POINT, 'Pt' )

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
