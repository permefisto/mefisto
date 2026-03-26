      SUBROUTINE TRFETO3( KTITRE, XYZPT1, XYZPT2, PTXYZD,
     %                    N1FEOC, NFETOI, NBVPSI, NFVPSI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES SIMPLES D'UNE ETOILE
C -----    VERSION 2 DU TABLEAU NFETOI ET NBVPSI FACES TRIANGULAIRES

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE
C XYZPT1 : POINT 1 A TRACER
C XYZPT2 : POINT 2 A TRACER
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 2
C          1: NUMERO NOTETR DU TETRAEDRE AU DELA DE LA FACE
C             =0 SI INCONNU
C          2: NUMERO PTXYZD DU 1-ER SOMMET DE LA FACE
C          3: NUMERO PTXYZD DU 2-ME SOMMET DE LA FACE
C          4: NUMERO PTXYZD DU 3-ME SOMMET DE LA FACE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C          S1S2 x S1S3 VERS L'INTERIEUR DE L'ETOILE
C NBVPSI : NOMBRE DE FACES NFVPSI A TRACER
C NFVPSI : NUMERO NFETOI DES NBVPSI FACES A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY Novembre 2014
C2345X7..............................................................012
      PARAMETER        (MXITEM=8192)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      DOUBLE PRECISION  XYZPT1(3), XYZPT2(3), PTXYZD(4,*)
      INTEGER           NFETOI(5,*), NFVPSI(*)
      CHARACTER*(*)     KTITRE
      INTRINSIC         MIN, MAX, REAL

      DOUBLE PRECISION  VECNOR(3), D
      REAL              XYZBAR(3,MXITEM), DISTOE(MXITEM), XYZB(3),
     %                  XYZ0(3), XYZ(3), R
      INTEGER           NOSTTR(3,MXITEM), ITEMTR(MXITEM)

      IF( .NOT. TRACTE ) RETURN
      IF(  NBVPSI.LE.0 .AND. N1FEOC.LE.0 ) RETURN
      NBITEM = 0

      DO L=1,3
         R = 1E27
C        CADRE COOEXT INITIAL
C        LE MINIMUM
         COOEXT(L,1) = R
C        LE MAXIMUM
         COOEXT(L,2) = -R
      ENDDO

      IF( XYZPT1(1).EQ.1D10 .AND. N1FEOC.GT.0 ) THEN
C        PREMIER SOMMET DE LA PREMIERE FACE DE NFETOI
         NS = ABS( NFETOI(2,N1FEOC) )
         DO L=1,3
C           COORDONNEE L DU SOMMET 1 DE LA FACE N1FEOC DE NFETOI
            R = REAL( PTXYZD(L,NS) )
C           CADRE COOEXT RESTREINT AUX FACES A TRACER
C           LE MINIMUM
            COOEXT(L,1) = MIN( COOEXT(L,1), R )
C           LE MAXIMUM
            COOEXT(L,2) = MAX( COOEXT(L,2), R )
         ENDDO
         GOTO 2
      ENDIF

C     LE POINT XYZPT1 EST LE PREMIER ITEM
      NBITEM = 1
      DO L=1,3
C        COORDONNEE L DE XYZPT1
         R = REAL( XYZPT1(L) )
         XYZBAR(L,NBITEM) = R
C        CADRE COOEXT RESTREINT AUX FACES A TRACER
C        LE MINIMUM
         COOEXT(L,1) = MIN( COOEXT(L,1), R )
C        LE MAXIMUM
         COOEXT(L,2) = MAX( COOEXT(L,2), R )
      ENDDO

C     LE POINT XYZPT2 EST LE PREMIER ITEM
      NBITEM = 2
      DO L=1,3
C        COORDONNEE L DE XYZPT2
         R = REAL( XYZPT2(L) )
         XYZBAR(L,NBITEM) = R
C        CADRE COOEXT RESTREINT AUX FACES A TRACER
C        LE MINIMUM
         COOEXT(L,1) = MIN( COOEXT(L,1), R )
C        LE MAXIMUM
         COOEXT(L,2) = MAX( COOEXT(L,2), R )
      ENDDO

C     LES NBVPSI FACES NFVPSI
      DO K=1, NBVPSI
         IF( NBITEM .GE. MXITEM ) GOTO 2
         NBITEM = NBITEM + 1
         NF = ABS( NFVPSI( K ) )
C        COORDONNEES DU BARYCENTRE DE LA FACE NBITEM
         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO M=1,3
C           SOMMET M DE LA FACE NF DE NFETOI
            NS = ABS( NFETOI(1+M,NF) )
            NOSTTR(M,NBITEM) = NS
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
      ENDDO

C     LES FACES DE L'ETOILE NFETOI
 2    N1FETOI = NBITEM
      NF1 = N1FEOC
 5    IF( NF1 .GT. 0 ) THEN

         DO L=1,NBVPSI
C           ON RETIRE LES FACES VPSI
            IF( NF1 .EQ. ABS(NFVPSI(L)) ) GOTO 8
         ENDDO

C        COORDONNEES DU BARYCENTRE DE LA FACE NBITEM
         IF( NBITEM .GE. MXITEM ) GOTO 9
         NBITEM = NBITEM + 1
         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO M=1,3
C           NUMERO DU SOMMET M DE LA FACE NBITEM DE NFETOI
            NS = ABS( NFETOI(1+M,NF1) )
            NOSTTR(M,NBITEM) = NS
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
 8       NF1 = NFETOI(5,NF1)
         GOTO 5

      ENDIF

ccc      print*,'trfeto3: nbitem=',nbitem,' cooext=',cooext

 9    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto3: AUGMENTER MXITEM=',MXITEM
      ENDIF

C     PARAMETRES DE TRACE
      AXOLAR = 0
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO
      AXOLAR = AXOLAR * 0.25
      AXOHAU = AXOLAR * 0.75
C     PAS DE PLAN ARRIERE ET AVANT
      AXOARR = 0
      AXOAVA = 0
      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )

      DISMOY = AXOLAR/30

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 20.0
      LORBITE = 1

ccc      print*,'trfeto3: AXOPTV=',AXOPTV,' AXOEIL=',AXOEIL,
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
         CALL XYZAXO( XYZBAR(1,K), XYZB )
C        DISTANCE A L'OEIL
         DISTOE( K ) = -XYZB(3)
C        IDENTITE INITIALE POUR LE NUMERO DE FACE DANS NOSTTR
         ITEMTR( K ) = K
      ENDDO

C     LE TRI CROISSANT SELON LA DISTANCE (COTE Z ) A L'OEIL
      CALL TRITRP( NBITEM, DISTOE, ITEMTR )
C
C     LE TRACE DES TRIANGLES DU CF ET DE SON ETOILE
      DO 50 K = NBITEM, 1, -1

C        NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF = ITEMTR( K )

         IF( NF .EQ. 1 ) THEN

C           TRACE DU POINT XYZPT1
            XYZ(1) = REAL( XYZPT1(1) )
            XYZ(2) = REAL( XYZPT1(2) )
            XYZ(3) = REAL( XYZPT1(3) )
            CALL SYMBOLE3D( NCVERT, XYZ, '[' )
            CALL SYMBOLE3D( NCBLEU, XYZ, '1' )

            IF( XYZPT1(1) .NE. XYZPT2(1) .OR.
     %          XYZPT1(2) .NE. XYZPT2(2) .OR.
     %          XYZPT1(3) .NE. XYZPT2(3) ) THEN
               XYZ0(1) = REAL( XYZPT2(1) )
               XYZ0(2) = REAL( XYZPT2(2) )
               XYZ0(3) = REAL( XYZPT2(3) )
               CALL TRAIT3D( NCCYAN, XYZ, XYZ0 )
            ENDIF

            GOTO 50

         ENDIF

         IF( NF .EQ. 2 ) THEN

C           TRACE DU POINT XYZPT2
            XYZ(1) = REAL( XYZPT2(1) )
            XYZ(2) = REAL( XYZPT2(2) )
            XYZ(3) = REAL( XYZPT2(3) )
            CALL SYMBOLE3D( NCVERT, XYZ, ']' )
            CALL SYMBOLE3D( NCBLEU, XYZ, '2' )

            IF( XYZPT1(1) .NE. XYZPT2(1) .OR.
     %          XYZPT1(2) .NE. XYZPT2(2) .OR.
     %          XYZPT1(3) .NE. XYZPT2(3) ) THEN
               XYZ0(1) = REAL( XYZPT1(1) )
               XYZ0(2) = REAL( XYZPT1(2) )
               XYZ0(3) = REAL( XYZPT1(3) )
               CALL TRAIT3D( NCCYAN, XYZ0, XYZ )
            ENDIF

            GOTO 50

         ENDIF

         IF( NF .GT. N1FETOI ) THEN
C           TRIANGLE DE L'ETOILE NFETOI
            NCF = NCGRIC
            NCA = NCNOIR
         ELSE
C           TRIANGLE VPSI
            NCF = NCORAN
            NCA = NCBLEU
         ENDIF

C        TRACE DE LA FACE ET DU NO DES 3 SOMMETS
         CALL TRFATR( NCF, NCA, NOSTTR(1,NF), PTXYZD )

C        TRACE DE LA NORMALE A LA FACE EN SON BARYCENTRE
         CALL VECNOR3D( PTXYZD(1,NOSTTR(1,NF)), PTXYZD(1,NOSTTR(2,NF)),
     %                  PTXYZD(1,NOSTTR(3,NF)), VECNOR )

C        NORME DU VECTEUR NORMAL
         D = SQRT( VECNOR(1)**2 + VECNOR(2)**2 + VECNOR(3)**2 )
         DO L=1,3
            XYZ(L) = XYZBAR(L,NF) + REAL( VECNOR(L) / D * DISMOY )
         ENDDO
         CALL SYMBOLE3D( NCMAGE, XYZBAR(1,NF), '*' )
         CALL TRAIT3D(   NCMAGE, XYZBAR(1,NF), XYZ )

         DO L=1,3
            NS = NOSTTR(L,NF)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
            CALL ENTIER3D( NCGRIS, XYZ, NS )
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
