      SUBROUTINE TRFETO6( KTITRE,  PTXYZD, NBTRCF, NOTRCF, LEFACO,
     %                    NBTETRA, NOTETRA, NOTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DES TRIANGLES DU CF et de leur NORMALE en leur BARYCENTRE
C -----  et DE NBTETRA TETRAEDRES

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
C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY Septembre 2015
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
      INTEGER           LEFACO(1:11,0:*), NOTRCF(NBTRCF),
     %                  NOTETRA(NBTETRA), NOTETR(8,*)
      INTEGER           NOSOTR(3)
      CHARACTER*(*)     KTITRE
      INTRINSIC         MIN, MAX

      REAL              R, XYZ(3), XYZBAR(3,MXITEM), DISTOE(MXITEM)
      INTEGER           NOSTTR(MXITEM), ITEMTR(MXITEM)
      DOUBLE PRECISION  VECNOR(3), D

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
      DO 5 K = 1, NBTRCF

         IF( NBITEM .GE. MXITEM ) GOTO 2
         NBITEM = NBITEM + 1
C        NUMERO DE LA FACE K DU CF
         NF = NOTRCF( K )
         IF( NF .LE. 0 ) GOTO 5
         NOSTTR( NBITEM ) = NF

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO M=1,3
C           SOMMET M DU TRIANGLE NF DE LEFACO
            NS = LEFACO( M, NF )
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

 5    ENDDO

C     CALCUL DU BARYCENTRE DES TETRAEDRES
 2    NBITEM1 = NBITEM
      DO N = 1, NBTETRA

         NT = NOTETRA(N)
         IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN

C           NUMERO NOTETR DU TETRAEDRE
            IF( NBITEM .GE. MXITEM ) GOTO 4
            NBITEM = NBITEM + 1
            NOSTTR(NBITEM) = NT

            DO L=1,3
               XYZBAR(L,NBITEM) = 0
            ENDDO
            DO K=1,4
C              SOMMET L DU TETRAEDRE
               NS = NOTETR(K,NT)
               DO L=1,3
C                 LE BARYCENTRE DU TETRAEDRE NBITEM
                  R = REAL( PTXYZD(L,NS) )
                  XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
C                 LE MINIMUM
                  COOEXT(L,1) = MIN( COOEXT(L,1), R )
C                 LE MAXIMUM
                  COOEXT(L,2) = MAX( COOEXT(L,2), R )
               ENDDO
            ENDDO
            DO L=1,3
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 4
            ENDDO
         ENDIF

      ENDDO

ccc      print*,'trfeto6: nbitem=',NBITEM,' cooext=',cooext

 4    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto6: AUGMENTER MXITEM=',MXITEM
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

      DISMOY = AXOLAR/30

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 10.0
      LORBITE = 1

ccc      print*,'trfeto6: AXOPTV=',AXOPTV,' AXOEIL=',AXOEIL,
ccc     %       ' AXOLAR=',AXOLAR,' AXOHAU=',AXOHAU

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
C     LE TRACE DES TETRAEDRES, DES TRIANGLES DU CF ET DE SON ETOILE
      DO 50 K = NBITEM, 1, -1

         CALL XVEPAISSEUR( 1 )
C        NUMERO DE L'ITEM LE PLUS ELOIGNE NON TRACE
         NF = ITEMTR( K )
         NT = NOSTTR( NF )

         IF( NF .GT. NBITEM1 ) THEN

C           TRACE DU TETRAEDRE NOTETR( NT )
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)

ccc               IF( NF .EQ. NBITEM ) THEN
cccC                 LE DERNIER TETRAEDRE EST TRACE EN NOIR
ccc                  NC = NCNOIR
ccc                  CALL XVEPAISSEUR( 3 )
ccc               ELSE
C                 LES AUTRES TETRAEDRES SONT TRACES EN MAGENTA

               NC = NCNOIR
               CALL XVEPAISSEUR( 1 )

ccc               ENDIF

               CALL TRTETRA( NC, NOTETR(1,NT), PTXYZD )
C
C              LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
               DO I = 1, 4
                  NS = NOTETR(I,NT)
                  XYZ(1) = REAL( PTXYZD(1,NS) )
                  XYZ(2) = REAL( PTXYZD(2,NS) )
                  XYZ(3) = REAL( PTXYZD(3,NS) )
                  CALL ENTIER3D( NCBLEU, XYZ, NS )
               ENDDO

            ENDIF
            GOTO 50

         ELSE

C           TRIANGLE NT DU CF LEFACO
            NCF = NCROUG
            NCA = NCROSE
            NOSOTR(1) = LEFACO(1,NT)
            NOSOTR(2) = LEFACO(2,NT)
            NOSOTR(3) = LEFACO(3,NT)

         ENDIF

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATR( NCF, NCA, NOSOTR, PTXYZD )

C        TRACE DE LA NORMALE A LA FACE EN SON BARYCENTRE
         CALL VECNOR3D( PTXYZD(1,NOSOTR(1)), PTXYZD(1,NOSOTR(2)),
     %                  PTXYZD(1,NOSOTR(3)), VECNOR )

C        NORME DU VECTEUR NORMAL
         D = SQRT( VECNOR(1)**2 + VECNOR(2)**2 + VECNOR(3)**2 )
         DO L=1,3
            XYZ(L) = XYZBAR(L,NF) + REAL( VECNOR(L) / D * DISMOY )
         ENDDO
         CALL SYMBOLE3D( NCJAUN, XYZBAR(1,NF), '*' )
         CALL TRAIT3D(   NCJAUN, XYZBAR(1,NF), XYZ )

C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
            CALL ENTIER3D( NCROSE, XYZ, NS )
         ENDDO

 50   ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )
C
C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
