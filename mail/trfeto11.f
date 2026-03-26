      SUBROUTINE TRFETO11( KTITRE,  PTXYZD,  N1FEOC, NFETOI,
     %                     NBTRCF,  NOTRCF,  LEFACO, NO0FAR,
     %                     NBTETRA, NOTETRA, NOTETR,
     %                     NBFAN2,  NOFAN2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES SIMPLES D'UNE ETOILE
C -----    + TRIANGLES DU CF + TETRAEDRES
C          + ARETES DU CF APPARTENANT PAS SEULEMENT A 2 FACES
C          AVEC LA VERSION 1 DU TABLEAU NFETOI
C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 1
C          1: NUMERO DU TETRAEDRE DANS NOTETR
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3  =0
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C          EN VERSION 2 (NFETOI(3,.)>0)
C          1: NUMERO NOTETR DU TETRAEDRE AU DELA DE LA FACE
C             =0 SI INCONNU
C          2: NUMERO PTXYZD DU 1-ER SOMMET DE LA FACE
C          3: NUMERO PTXYZD DU 2-ME SOMMET DE LA FACE
C          4: NUMERO PTXYZD DU 3-ME SOMMET DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE
C             L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C NBTRCF : NOMBRE DE TRIANGLES DU TABLEAU NOTRCF
C          C-A-D DU POLYGONE ENCORE DIT ENSUITE ETOILE
C NOTRCF : SI NOTRCF(*)>0 NUMERO LEFACO DU TRIANGLE
C                      <0 NUMERO NO0FAR DU TRIANGLE AJOUTE AU CF
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C          NORMALE VERS L'INTERIEUR DU TETRAEDRE LA CONTENANT

C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C NBFAN2: NOMBRE D'ARETES DU CF A PROBLEME DANS L'ETOILE
C NOFAN2: NOFAN2(NDFAN2+1)=NOMBRE NBF DE FACES NFETOI DE L'ARETE
C         NOFAN2(NDFAN2+2)=NUMERO PTXYZD DU SOMMET 1  DE L'ARETE
C         NOFAN2(NDFAN2+3)=NUMERO PTXYZD DU SOMMET 2  DE L'ARETE
C         NOFAN2(NDFAN2+4)=NUMERO DANS N1ARCF DU CIRCUIT DU CF
C         NOFAN2(NDFAN2+4+1...NBF)=NUMERO NF DANS NFETOI DE LA FACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& ST PIERRE DU PERRAY Septembre 2015
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
     %                  NOTETRA(NBTETRA), NOTETR(8,*),
     %                  NOFAN2(*)
      INTEGER           NFETOI(5,*), NOSOTR(3)
      CHARACTER*(*)     KTITRE

      DOUBLE PRECISION  VECNOR(3), D
      REAL              R, XYZ(3), XYZ2(3), XYZBAR(3,MXITEM),
     %                  DISTOE(MXITEM)
      INTEGER           NOSTTR(MXITEM), ITEMTR(MXITEM), NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /
C                       NORMALE VERS L'EXTERIEUR DU TETRAEDRE

      IF( .NOT. TRACTE ) RETURN

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NBITEM = 0
      NF1    = N1FEOC
 2    IF( NF1 .GT. 0 ) THEN

         IF( NBITEM .GE. MXITEM ) GOTO 4

         NBITEM = NBITEM + 1
C        STOKAGE DU NUMERO DE LA FACE DE L'ETOILE
         NOSTTR(NBITEM) = NF1

C        CALCUL DU BARYCENTRE DE LA FACE NBITEM DE L'ETOILE
         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO

         DO K=1,3

            IF( NFETOI(3,NF1) .EQ. 0 ) THEN
C              NFETOI VERSION 1
C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               NTE = NFETOI(1,NF1)
C              NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
               I = ABS( NFETOI(2,NF1) )
C              NUMERO DU SOMMET K DE LA FACE NF1 DE NFETOI
               NS = NOTETR( NOSOFATE(K,I) , NTE )
            ELSE
C              NFETOI VERSION 2
               NS = ABS( NFETOI(1+K,NF1) )
            ENDIF

            DO L=1,3
               R = REAL( PTXYZD(L,NS) )
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
C              LE BARYCENTRE DE LA FACE
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
            ENDDO
         ENDDO
         DO L=1,3
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 3
         ENDDO

C        LA FACE SUIVANTE
         NF1 = NFETOI(5,NF1)
         GOTO 2

      ENDIF
 4    NBITEM0 = NBITEM

cccC     AJOUT DANS LE TITRE
ccc      M = NUDCNB( KTITRE )
ccc      WRITE(KTITRE(M+2:M+5),'(I4)') NBITEM0
ccc      KTITRE = KTITRE(1:M+5) // ' FACES DE L''ETOILE'

C     CALCUL DU BARYCENTRE DES FACES TRIANGULAIRES DU CF
      DO K = 1, NBTRCF

         IF( NBITEM .GE. MXITEM ) GOTO 6
         NBITEM = NBITEM + 1

C        NUMERO DE LA FACE K DU CF
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

C     CALCUL DU BARYCENTRE DES TETRAEDRES
 6    NBITEM1 = NBITEM
      DO N = 1, NBTETRA

         IF( NBITEM .GE. MXITEM ) GOTO 8
         NBITEM = NBITEM + 1
         NT = ABS( NOTETRA(N) )

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
      ENDDO

C     CALCUL DU BARYCENTRE DES ARETES A NBF FACES NON 2
 8    NBITEM2 = NBITEM
      NDFAN2  = 0
      DO N = 1, NBFAN2

C        NOMBRE DE FACES DE L'ARETE DU CF
         NBF = NOFAN2( NDFAN2+1 )

C        LE NO DES 2 SOMMETS DE L'ARETE DU CF
         NS1 = NOFAN2( NDFAN2+2 )
         NS2 = NOFAN2( NDFAN2+3 )

cccC        NUMERO DU CIRCUIT DU CF (DANS N1ARCF)
ccc         NCF = NOFAN2( NDFAN2 + 4 )

         IF( NBITEM .GE. MXITEM ) GOTO 9
         NBITEM = NBITEM + 1

C        NUMERO NOFAN2 DE L'ARETE
         NOSTTR(NBITEM) = NDFAN2

C        POINTEUR SUR LA DERNIERE INFO DE L'ARETE
         NDFAN2 = NDFAN2 + 4 + NBF

         DO L=1,3
C           LE BARYCENTRE DE L'ARETE NBITEM
            R = REAL( (PTXYZD(L,NS1) +  PTXYZD(L,NS2))/2 )
            XYZBAR(L,NBITEM) = R
C           LE MINIMUM
            COOEXT(L,1) = MIN( COOEXT(L,1), R )
C           LE MAXIMUM
            COOEXT(L,2) = MAX( COOEXT(L,2), R )
         ENDDO

      ENDDO       

C     PARAMETRES DE VISEE
 9    AXOLAR = 0
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
      PREDUF  = 30.0
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

         IF( NF .GT. NBITEM2 ) THEN

C           NOMBRE DE FACES DE L'ARETE NT
            NBF = NOFAN2(NT+1)

C           LE NO DES 2 SOMMETS DE L'ARETE DU CF
            NS1 = NOFAN2(NT+2)
            XYZ(1) = REAL( PTXYZD(1,NS1) )
            XYZ(2) = REAL( PTXYZD(2,NS1) )
            XYZ(3) = REAL( PTXYZD(3,NS1) )

            NS2 = NOFAN2(NT+3)
            XYZ2(1) = REAL( PTXYZD(1,NS2) )
            XYZ2(2) = REAL( PTXYZD(2,NS2) )
            XYZ2(3) = REAL( PTXYZD(3,NS2) )

            CALL XVEPAISSEUR( 7 )
C           COULEUR DE TRACE DE L'ARETE NS1-NS2
ccc            NC = 4 + NBF
            NC = NCORAN
            CALL TRAIT3D( NC, XYZ, XYZ2 )
            CALL ENTIER3D( NCBLEU, XYZ,  NS1 )
            CALL ENTIER3D( NCBLEU, XYZ2, NS2 )
            GOTO 50

         ELSE IF( NF .GT. NBITEM1 ) THEN

C           TRACE DU TETRAEDRE NOTETR( NT )
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
ccc               IF( NF .EQ. NBITEM2 ) THEN
cccC                 LE DERNIER TETRAEDRE EST TRACE EN NOIR
ccc                  NC = NCNOIR
ccc                  CALL XVEPAISSEUR( 3 )
ccc               ELSE
C                 LES AUTRES TETRAEDRES SONT TRACES EN CYAN
                  NC = NCCYAN
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
                  CALL ENTIER3D( NCGRIS, XYZ, NS )
               ENDDO

            ENDIF
            GOTO 50

         ELSE IF( NF .GT. NBITEM0 ) THEN

C           TRIANGLE NT DU CF LEFACO
            IF( NT .GT. 0 ) THEN
               NCF = NCROUG
               NCA = NCROSE
               NOSOTR(1) = LEFACO(1,NT)
               NOSOTR(2) = LEFACO(2,NT)
               NOSOTR(3) = LEFACO(3,NT)
            ELSE
               NCF = NCGRIM
               NCA = NCJAUN
               NT = -NT
               NOSOTR(1) = NO0FAR(1,NT)
               NOSOTR(2) = NO0FAR(2,NT)
               NOSOTR(3) = NO0FAR(3,NT)
            ENDIF

         ELSE IF( NBITEM0 .GT. 0 . AND. NF .LE. NBITEM0 ) THEN

C           TRIANGLE SIMPLE NT DE L'ETOILE
            NCF = NCGRIM
            NCA = NCBLAN

            IF( NFETOI(3,NT) .EQ. 0 ) THEN
C              NFETOI VERSION 1
C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               NTE = NFETOI(1,NT)
C              NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
               I = ABS( NFETOI(2,NT) )
C              NUMERO DES 3 SOMMETS
               NOSOTR(1) = NOTETR( NOSOFATE(1,I) , NTE )
               NOSOTR(2) = NOTETR( NOSOFATE(3,I) , NTE )
               NOSOTR(3) = NOTETR( NOSOFATE(2,I) , NTE )
            ELSE
               NOSOTR(1) = ABS( NFETOI(2,NT) )
               NOSOTR(2) = ABS( NFETOI(3,NT) )
               NOSOTR(3) = ABS( NFETOI(4,NT) )
            ENDIF

         ELSE

            GOTO 50

         ENDIF

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATR( NCF, NCA, NOSOTR, PTXYZD )

C        TRACE DE LA NORMALE A LA FACE EN SON BARYCENTRE
         CALL VECNOR3D( PTXYZD(1,NOSOTR(1)),
     %                  PTXYZD(1,NOSOTR(2)),
     %                  PTXYZD(1,NOSOTR(3)), VECNOR )

C        NORME DU VECTEUR NORMAL
         D = SQRT( VECNOR(1)**2 + VECNOR(2)**2 + VECNOR(3)**2 )
         DO L=1,3
            XYZ(L) = XYZBAR(L,NF) + REAL( VECNOR(L) / D * DISMOY )
         ENDDO

C        TRACE DU VECTEUR NORMAL
         CALL SYMBOLE3D( NCNOIR, XYZBAR(1,NF), '*' )
         CALL TRAIT3D(   NCNOIR, XYZBAR(1,NF), XYZ )

C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
            CALL ENTIER3D( NCNOIR, XYZ, NS )
         ENDDO

 50   ENDDO

C     TITRE ET TRACE EFFECTIF
ccc      KTITRE='      FACES DE L''ETOILE,        TRIANGLES du CF,        T
ccc     %ETRAEDRES CREES'
ccc      WRITE(KTITRE(1:5)  ,'(I5)') NBITEM0
ccc      WRITE(KTITRE(27:31),'(I5)') NBTRCF
ccc      WRITE(KTITRE(51:55),'(I5)') NBTETRA
      CALL TRFINS( KTITRE )
C
C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
