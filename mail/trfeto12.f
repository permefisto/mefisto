      SUBROUTINE TRFETO12( KTITRE,  PTXYZD,
     %                     NBTETRA, NOTETRA, NOTETR,
     %                     NBFANR,  NOFANR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES DES NBTETRA TETRAEDRES
C -----    + FACES TRIANGULAIRES DE TETRAEDRES

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE

C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NOTETRA: NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C NBFANR: NOMBRE DE FACES TRIANGULAIRES DES TETRAEDRES
C NOFANR: NOFANR(1,i)=NUMERO NOTETR DU TETRAEDRE
C         NOFANR(2,i)=NUMERO DE LA FACE DU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& ST PIERRE DU PERRAY  Decembre 2016
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
      INTEGER           NOTETRA(NBTETRA), NOTETR(8,*), NOFANR(2,NBFANR)
      INTEGER           NOSOTR(3)
      CHARACTER*(*)     KTITRE
      REAL              R, XYZ(3), XYZBAR(3,MXITEM),
     %                  DISTOE(MXITEM)
      INTEGER           NOSTTR(MXITEM), ITEMTR(MXITEM), NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /

      IF( .NOT. TRACTE ) RETURN

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NBITEM = 0

C     CALCUL DU BARYCENTRE DES TETRAEDRES
      DO N = 1, NBTETRA

         IF( NBITEM .GE. MXITEM ) GOTO 2

         NTE = ABS( NOTETRA(N) )

         IF( NTE.GT.0 .AND. NOTETR(1,NTE).GT.0 ) THEN
            NBITEM = NBITEM + 1

C           NUMERO NOTETR DU TETRAEDRE
            NOSTTR( NBITEM ) = NTE

            DO L=1,3
               XYZBAR(L,NBITEM) = 0
            ENDDO
            DO K=1,4

C              SOMMET L DU TETRAEDRE
               NS = NOTETR(K,NTE)

               IF( NS .LE. 0 ) THEN
                  PRINT*,'trfeto12: PROBLEME NOTETR(',NTE,')=',
     %                (NOTETR(L,NTE),L=1,8),' A deboguer!...'
                  RETURN
               ENDIF

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

C     CALCUL DU BARYCENTRE DES FACES FANR
 2    NBITEM0 = NBITEM
      DO N = 1, NBFANR

C        NUMERO DU TETRAEDRE
         NTE = NOFANR(1,N)

C        NUMERO DE LA FACE
         NFE = NOFANR(2,N)

         IF( NBITEM .GE. MXITEM ) GOTO 4
         NBITEM = NBITEM + 1

C        NUMERO NOFANR DE LA FACE TRIANGULAIRE
         NOSTTR( NBITEM ) = N

         DO L=1,3
            XYZBAR(L,NBITEM) = 0
         ENDDO

         DO L=1,3

C           NUMERO DU SOMMET L DE LA FACE NFE DU TETRAEDRE NTE
            NS = NOTETR( NOSOFATE(L,NFE), NTE )

C           LE BARYCENTRE DU TRIANGLE NBITEM
            R = REAL( PTXYZD(L,NS) )
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) + R

C           LE MINIMUM
            COOEXT(L,1) = MIN( COOEXT(L,1), R )
C           LE MAXIMUM
            COOEXT(L,2) = MAX( COOEXT(L,2), R )

         ENDDO

         DO L=1,3
            XYZBAR(L,NBITEM) = XYZBAR(L,NBITEM) / 3
         ENDDO

      ENDDO       

 4    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto12: AUGMENTER MXITEM=',MXITEM
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
      PREDUF  = 20.0
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

C        NUMERO DE L'ITEM LE PLUS ELOIGNE NON TRACE
         NF = ITEMTR( K )
         NT = NOSTTR( NF )

         IF( NF .GT. NBITEM0 ) THEN

C           NUMERO DU TETRAEDRE
            NTE = NOFANR(1,NT)

C           TRACE DES TETRAEDRES OPPOSES AU TETRAEDRE NTE
            CALL XVEPAISSEUR( 0 )
            DO KK=1,4
               NTEOP = NOTETR( 4+KK, NTE )
               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
                  CALL TRTETRA( NCGRIM, NOTETR(1,NTEOP), PTXYZD )
               ENDIF
            ENDDO

C           NUMERO DE LA FACE
            NFE = NOFANR(2,NT)

C           NUMERO DES 3 SOMMETS DE LA FACE NFE DU TETRAEDRE NTE
            NOSOTR(1) = NOTETR( NOSOFATE(1,NFE) , NTE )
            NOSOTR(2) = NOTETR( NOSOFATE(3,NFE) , NTE )
            NOSOTR(3) = NOTETR( NOSOFATE(2,NFE) , NTE )

C           TRACE DU TRIANGLE NOSOTR
            CALL XVEPAISSEUR( 5 )
            NCF = NCROUG
            NCA = NCORAN
            CALL TRFATR( NCF, NCA, NOSOTR, PTXYZD )

         ELSE

C           TRACE DU TETRAEDRE NOTETR( NT )
            IF( NT .GT. 0 ) THEN

C              LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
               CALL XVEPAISSEUR( 1 )
               CALL TRTETRA( NCNOIR, NOTETR(1,NT), PTXYZD )
C
C              LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETRA
               DO I = 1, 4
                  NS = NOTETR(I,NT)
                  XYZ(1) = REAL( PTXYZD(1,NS) )
                  XYZ(2) = REAL( PTXYZD(2,NS) )
                  XYZ(3) = REAL( PTXYZD(3,NS) )
                  CALL ENTIER3D( NCCYAN, XYZ, NS )
               ENDDO

            ENDIF

         ENDIF

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
