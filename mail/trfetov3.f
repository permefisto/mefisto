      SUBROUTINE TRFETOV3( KTITRE,  XYZSOM,  NBSSET, N1SSET, NFETOI,
     %                     NBTETRA, NOTETRA, NBSOTE, NSTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES SIMPLES D'UNE ETOILE
C -----    + TETRAEDRES AJOUTES + NUMERO DU SOMMET NS1
C          AVEC LA VERSION 3 DU TABLEAU NFETOI
C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C XYZSOM : PAR POINT : X  Y  Z
C NBSSET : NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET SOUS ETOILES
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 3
C          1: NUMERO XYZSOM DU 1-ER SOMMET DE LA FACE
C          2: NUMERO XYZSOM DU 2-ME SOMMET DE LA FACE
C          3: NUMERO XYZSOM DU 3-ME SOMMET DE LA FACE
C          4: NUMERO NSTETR DU TETRAEDRE
C             =0 SI INCONNU
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NSTETR DES TETRAEDRES A TRACER
C NBSOTE : MAXIMUM DU 1-ER INDICE DU TABLEAU NSTETR
C NSTETR : LISTE DES TETRAEDRES
C          SOMMET1, SOMMET2, SOMMET3, SOMMET4,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY Fevrier 2016
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

      CHARACTER*(*)     KTITRE
      REAL              XYZSOM(3,*)
      INTEGER           N1SSET(NBSSET), NOTETRA(NBTETRA),
     %                  NSTETR(NBSOTE,*), NFETOI(5,*)

      INTRINSIC         MIN, MAX
      REAL              R, XYZ(3), XYZBAR(3,MXITEM), DISTOE(MXITEM)
      INTEGER           NOSTTR(MXITEM), ITEMTR(MXITEM), NOSOTR(3)

      IF( .NOT. TRACTE ) RETURN

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NBITEM = 0
      DO NSS = 1, NBSSET

C        SOUS-ETOILE NSS
         NF1 = N1SSET( NSS )
 2       IF( NF1 .GT. 0 ) THEN

            IF( NBITEM .GE. MXITEM ) GOTO 4
            NBITEM = NBITEM + 1
C           STOKAGE DU NUMERO DE LA FACE DE L'ETOILE
            NOSTTR(NBITEM) = NF1

C           CALCUL DU BARYCENTRE DE LA FACE NBITEM DE L'ETOILE
            DO L=1,3
               XYZBAR(L,NBITEM) = 0
            ENDDO

            DO K=1,3
C              NUMERO DU SOMMET K DE LA FACE NF1 DE NFETOI
               NS = ABS( NFETOI(K,NF1) )
               DO L=1,3
                  R = XYZSOM(L,NS)
C                 LE MINIMUM
                  COOEXT(L,1) = MIN( COOEXT(L,1), R )
C                 LE MAXIMUM
                  COOEXT(L,2) = MAX( COOEXT(L,2), R )
C                 LE BARYCENTRE DE LA FACE
                  XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R
               ENDDO
            ENDDO
            DO L=1,3
               XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 3
            ENDDO

C           LA FACE SUIVANTE
            NF1 = NFETOI(5,NF1)
            GOTO 2

         ENDIF
      ENDDO

C     CALCUL DU BARYCENTRE DES TETRAEDRES
 4    NBITEM0 = NBITEM
      DO N = 1, NBTETRA

         IF( NBITEM .GE. MXITEM ) GOTO 6
         NBITEM = NBITEM + 1
         NT = ABS( NOTETRA(N) )

C        NUMERO NSTETR DU TETRAEDRE
         NOSTTR(NBITEM) = NT

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO
         DO K=1,4
C           SOMMET L DU TETRAEDRE
            NS = NSTETR(K,NT)
            DO L=1,3
C              LE BARYCENTRE DU TETRAEDRE NBITEM
               R = XYZSOM(L,NS)
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

ccc      print*,'trfetov3: nbitem=',NBITEM,' cooext=',cooext

 6    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfetov3: AUGMENTER MXITEM=',MXITEM
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
      PREDUF  = 30.0
      LORBITE = 1

ccc      print*,'trfetov3: AXOPTV=',AXOPTV,' AXOEIL=',AXOEIL,
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

         IF( NF .GT. NBITEM0 ) THEN

C           TRACE DU TETRAEDRE NSTETR( NT )
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NSTETR(*,NT)
               IF( NF .EQ. NBITEM ) THEN
C                 LE DERNIER TETRAEDRE EST TRACE EN NOIR
                  NC = NCNOIR
                  CALL XVEPAISSEUR( 3 )
               ELSE
C                 LES AUTRES TETRAEDRES SONT TRACES EN MAGENTA
                  NC = NCMAGE
                  CALL XVEPAISSEUR( 1 )
               ENDIF
               CALL TRTETRAR( NC, NSTETR(1,NT), XYZSOM )
C
C              LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
               DO I = 1, 4
                  NS = NSTETR(I,NT)
                  CALL ENTIER3D( NCBLEU, XYZSOM(1,NS), NS )
               ENDDO

            ENDIF
            GOTO 50

         ELSE

C           TRIANGLE SIMPLE NT DE L'ETOILE
            NCF = NCGRIM
            NCA = NCBLAN
            NOSOTR(1) = ABS( NFETOI(1,NT) )
            NOSOTR(2) = NFETOI(2,NT)
            NOSOTR(3) = NFETOI(3,NT)
 
         ENDIF

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATRR( NCF, NCA, NOSOTR, XYZSOM )
C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            CALL ENTIER3D( NCNOIR, XYZSOM(1,NS), NS )
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
