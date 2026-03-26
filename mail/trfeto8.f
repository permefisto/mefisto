      SUBROUTINE TRFETO8( KTITRE,  PTXYZD,
     %                    NBTRCF,  NOTRCF,  LEFACO, NO0FAR,
     %                    NBTETR0, NOTETR0,
     %                    NBTETR1, NOTETR1, NOTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES TRIANGLES DU CF + TETRAEDRES PRIMITIFS et AJOUTES
C -----
C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NBTRCF : NOMBRE DE TRIANGLES DU TABLEAU NOTRCF
C          C-A-D DU POLYGONE ENCORE DIT ENSUITE ETOILE
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DE L'ETOILE
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C          NORMALE VERS L'INTERIEUR DU TETRAEDRE LA CONTENANT
C NBTETR0: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES INITIAUX
C NOTETR0: NUMERO DANS NOTETR DES TETRAEDRES INITIAUX A TRACER
C NBTETR1: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES FINAUX
C NOTETR1: NUMERO DANS NOTETR DES TETRAEDRES FINAUX A TRACER

C NOTETR : LISTE DES TETRAEDRES DU MAILLAGE
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY    Juin 2016
C2345X7..............................................................012
      PARAMETER         (MXITEM=8192)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE
C
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           LEFACO(1:11,0:*), NOTRCF(NBTRCF), NO0FAR(3,*),
     %                  NOTETR0(NBTETR0), NOTETR1(NBTETR1), NOTETR(8,*)
      INTEGER           NOSOTR(3)
      CHARACTER*(*)     KTITRE

      REAL              R, XYZ(3), XYZBAR(3,MXITEM), DISTOE(MXITEM)
      INTEGER           NOSTTR(MXITEM), ITEMTR(MXITEM)

C
      IF( .NOT. TRACTE ) RETURN

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NBITEM = 0

C     CALCUL DU BARYCENTRE DES FACES TRIANGULAIRES DU CF
      DO K = 1, NBTRCF

         IF( NBITEM .GE. MXITEM ) GOTO 2
         NBITEM = NBITEM + 1
C        NUMERO DE LA FACE K DU CF  <0 ou >0
         NF = NOTRCF( K )
         NOSTTR( NBITEM ) = NF

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO M=1,3
            IF( NF .GT. 0 ) THEN
C              SOMMET M DU TRIANGLE NF DE LEFACO
               NS = LEFACO( M, NF )
            ELSE
C              SOMMET M DU TRIANGLE NF DE NO0FAR
               NS = NO0FAR( M, -NF )
            ENDIF
            DO L=1,3
C              LE BARYCENTRE DU TRIANGLE NBITEM
               R = REAL( PTXYZD(L,NS) )
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
         ENDDO
         DO L=1,3
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 3
         ENDDO

      ENDDO

C     CALCUL DU BARYCENTRE DES TETRAEDRES INITIAUX
 2    NBITEM1 = NBITEM
      DO 3 N = 1, NBTETR0

         IF( NBITEM .GE. MXITEM ) GOTO 4
         NBITEM = NBITEM + 1
         NT = ABS( NOTETR0(N) )

         IF( NOTETR(1,NT) .EQ. 0 ) THEN
ccc            PRINT*,'trfeto8: Tetraedre',NT,' St:?',(NOTETR(L,NT),L=1,8)
            GOTO 3
         ENDIF

C        NUMERO NOTETR DU TETRAEDRE
         NOSTTR(NBITEM) = NT

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO K=1,4
C           SOMMET L DU TETRAEDRE
            NS = NOTETR(K,NT)
            DO L=1,3
C              LE BARYCENTRE DU TETRAEDRE NBITEM
               R = REAL( PTXYZD(L,NS) )
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
         ENDDO
         DO L=1,3
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 4
         ENDDO
 3    ENDDO

C     CALCUL DU BARYCENTRE DES TETRAEDRES AJOUTES
 4    NBITEM2 = NBITEM
      DO 5 N = 1, NBTETR1

         NT = ABS( NOTETR1(N) )

         IF( NOTETR(1,NT) .EQ. 0 ) THEN
ccc            PRINT*,'trfeto8: Tetraedre',NT,' St?:',(NOTETR(L,NT),L=1,8)
            GOTO 5
         ENDIF

         IF( NBITEM .GE. MXITEM ) GOTO 6
         NBITEM = NBITEM + 1

C        NUMERO NOTETR DU TETRAEDRE
         NOSTTR(NBITEM) = NT

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO K=1,4
C           SOMMET L DU TETRAEDRE
            NS = NOTETR(K,NT)
            DO L=1,3
C              LE BARYCENTRE DU TETRAEDRE NBITEM
               R = REAL( PTXYZD(L,NS) )
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
         ENDDO
         DO L=1,3
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 4
         ENDDO

 5    ENDDO

 6    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto8: AUGMENTER MXITEM=',MXITEM
      ENDIF

C     PARAMETRES DE VISEE
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
      PREDU0  = PREDUF
      PREDUF  = 0.0
      LORBITE = 1

C     LE TITRE EST AUGMENTE
      CALL SANSDBL( KTITRE, LK )
      KTITRE = KTITRE(1:LK) //' Tetra INITIAUX BLEUS et CREES en VERT'
      CALL SANSDBL( KTITRE, LK )

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
C     CALCUL DE LA DISTANCE A L'OEIL DU BARYCENTRE DES NBITEM FACES
      DO K = 1, NBITEM
C        AXONOMETRIE DU BARYCENTRE
         CALL XYZAXO( XYZBAR(1,K), XYZ )
C        DISTANCE A L'OEIL
         DISTOE( K ) = -XYZ(3)
C        IDENTITE INITIALE POUR LE NUMERO DE FACE DANS NOSTTR
         ITEMTR( K ) = K
      ENDDO

C     LE TRI CROISSANT SELON LA DISTANCE (COTE Z ) A L'OEIL
      CALL TRITRP( NBITEM, DISTOE, ITEMTR )
C
C     LE TRACE DES TETRAEDRES, DES TRIANGLES DU CF
      DO 50 K = NBITEM, 1, -1

C        NUMERO DE L'ITEM LE PLUS ELOIGNE NON TRACE
         NF = ITEMTR( K )
         NT = NOSTTR( NF )

         IF( NF .GT. NBITEM2 ) THEN

C           TRACE DU TETRAEDRE NOTETR( NT ) FINAL
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
               PREDUF = 6.0
               CALL TRTETRA( NCVERT, NOTETR(1,NT), PTXYZD )
C
C              LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETR1
               DO I = 1, 4
                  NS = NOTETR(I,NT)
                  XYZ(1) = REAL( PTXYZD(1,NS) )
                  XYZ(2) = REAL( PTXYZD(2,NS) )
                  XYZ(3) = REAL( PTXYZD(3,NS) )
                  CALL ENTIER3D( NCTURQ, XYZ, NS )
               ENDDO

            ENDIF
            GOTO 50

         ELSE IF( NF .GT. NBITEM1 ) THEN

C           TRACE DU TETRAEDRE NOTETR( NT ) INITIAL
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
               PREDUF = 3.0
               CALL TRTETRA( NCBLEU, NOTETR(1,NT), PTXYZD )
C
C              LE TRACE DU NO DES SOMMETS DES TETRAEDRES DE NOTETR0
               DO I = 1, 4
                  NS = NOTETR(I,NT)
                  XYZ(1) = REAL( PTXYZD(1,NS) )
                  XYZ(2) = REAL( PTXYZD(2,NS) )
                  XYZ(3) = REAL( PTXYZD(3,NS) )
                  CALL ENTIER3D( NCCYAN, XYZ, NS )
               ENDDO

            ENDIF
            GOTO 50

         ELSE

C           TRIANGLE NT DU CF LEFACO ou NO0FAR
            IF( NT .GT. 0 ) THEN
               NCF = NCROUG
               NCA = NCORAN
               NOSOTR(1) = LEFACO( 1, NT )
               NOSOTR(2) = LEFACO( 2, NT )
               NOSOTR(3) = LEFACO( 3, NT )
            ELSE
               NCF = NCROSE
               NCA = NCORAN
               NT  = -NT
C              SOMMETS DU TRIANGLE -NT DE NO0FAR
               NOSOTR(1) = NO0FAR( 1, NT )
               NOSOTR(2) = NO0FAR( 2, NT )
               NOSOTR(3) = NO0FAR( 2, NT )
            ENDIF

         ENDIF

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATR( NCF, NCA, NOSOTR, PTXYZD )
C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
            CALL ENTIER3D( NCGRIS, XYZ, NS )
         ENDDO

 50   ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE(1:LK) )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
