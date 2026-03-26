      SUBROUTINE TRFETO7R( XYZSOM, N1FEOC, NFETOI, NBSSET, N1SSET,
     %                     NBAR3F, NSAR3F )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES SIMPLES D'UNE ETOILE
C -----    FORME DES SOUS-ETOILES N1SSET
C          AVEC LA VERSION 3 DU TABLEAU NFETOI
C          VERSION XYZ EN SIMPLE PRECISION

C ENTREES:
C --------
C XYZSOM : PAR POINT : X  Y  Z
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : 1: NUMERO XYZSOM DU 1-ER SOMMET DE LA FACE
C          2: NUMERO XYZSOM DU 2-ME SOMMET DE LA FACE
C          3: NUMERO XYZSOM DU 3-ME SOMMET DE LA FACE
C          4: NUMERO NOTETR DU TETRAEDRE AU DELA DE LA FACE
C             =0 SI INCONNU
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C NBSSET : NOMBRE DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET SOUS ETOILES
C NBAR3F : NOMBRE D'ARETES APPARTENANT A AU MOINS 3 FACES
C NSAR3F : NUMERO DES 2 SOMMETS DES ARETES 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& St PIERRE DU PERRAY   Fevrier 2016
C2345X7..............................................................012
      PARAMETER         (MXSSET=1024, MXITEM=8192)
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
      REAL              XYZSOM(3,*)
      INTEGER           NFETOI(5,*), N1SSET(NBSSET), NSAR3F(2,NBAR3F)

      CHARACTER*80      KTITRE
      INTRINSIC         MIN, MAX
      REAL              R, XYZ(3), XYZ2(3),
     %                  XYZBAR(3,MXITEM), DISTOE(MXITEM)
      INTEGER           NDITEM(0:MXSSET), NOSTTR(MXITEM), ITEMTR(MXITEM)
      INTEGER           NOSOTR(3)
      DOUBLE PRECISION  VECNOR(3), D

      IF( .NOT. TRACTE ) RETURN

      IF( NBSSET .GT. MXSSET ) THEN
         print *,'trfeto7r: AUGMENTER MXSSET=',MXSSET,' nbsset=',nbsset,
     %           ' n1sset=',(n1sset(K),k=1,nbsset)
      ENDIF

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      NBITEM = 0
      DO M = 0, NBSSET
         IF( M .EQ. 0 ) THEN
C           L'ETOILE
            NF1 = N1FEOC
         ELSE
C           LA SOUS-ETOILE M
            NF1 = N1SSET( M )
         ENDIF

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

C        NUMERO DU DERNIER ITEM DE LA SOUS-ETOILE M
         NDITEM(M) = NBITEM

      ENDDO
      IF( NBITEM .LE. 0 ) RETURN

ccc      print*,'trfeto7r: nbitem=',NBITEM,' cooext=',cooext

 4    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto7r: AUGMENTER MXITEM=',MXITEM
      ENDIF

C     PARAMETRES DE VISEE
      AXOLAR = 0 
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO
      AXOLAR = AXOLAR * 0.7
      AXOHAU = AXOLAR * 0.75
C     PAS DE PLAN ARRIERE ET AVANT
      AXOARR = 0
      AXOAVA = 0
      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )

      DISMOY = AXOLAR/30

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 25.0
      LORBITE = 1

ccc      print*,'trfeto7r: AXOPTV=',AXOPTV,' AXOEIL=',AXOEIL,
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

         DO M=0,NBSSET
            IF( NF .LE. NDITEM(M) ) GOTO 30
         ENDDO

C        TRIANGLE NT D'UNE SOUS ETOILE
C        COULEUR DE LA FACE SELON L'ORDRE
C         0:NCNOIR,  1:NCROUG,  2:NCVERT,  3:NCBLEU,
C         4:NCCYAN,  5:NCJAUN,  6:NCMAGE,  7:NCBLAN,
C         8:NCGRIS,  9:NCGRIM, 10:NCGRIC, 11:NCBEIG,
C        12:NCORAN, 13:NCSAUM, 14:NCROSE, 15:NCTURQ
 30      IF( M .EQ. 0 ) THEN
C           FACES ISSUES DE N1FEOC EN GRIS MOYEN
            NCF = NCGRIM
         ELSE
C           AUTRES FACES AVEC UNE SEULE COULEUR POUR CHAQUE SOUS-ETOILE
            NCF = 3+M
            IF( NCF .GT. 15 ) THEN
               NCF = MOD( NCF, 15) + 3
            ENDIF
         ENDIF

C        3 SOMMETS DE LA FACE NT DE NFETOI
         NOSOTR(1) = ABS( NFETOI(1,NT) )
         NOSOTR(2) = NFETOI(2,NT)
         NOSOTR(3) = NFETOI(3,NT)

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATRR( NCF, NCNOIR, NOSOTR, XYZSOM )

C        TRACE DE LA NORMALE A LA FACE EN SON BARYCENTRE
         CALL VECNOR3( XYZSOM(1,NOSOTR(1)), XYZSOM(1,NOSOTR(2)),
     %                 XYZSOM(1,NOSOTR(3)), VECNOR )

C        NORME DU VECTEUR NORMAL
         D = SQRT( VECNOR(1)**2 + VECNOR(2)**2 + VECNOR(3)**2 )
         DO L=1,3
            XYZ(L) = XYZBAR(L,NF) + REAL( VECNOR(L) / D * DISMOY )
         ENDDO
         CALL SYMBOLE3D( NCROUG, XYZBAR(1,NF), '*' )
         CALL TRAIT3D(   NCROUG, XYZBAR(1,NF), XYZ )

C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            CALL ENTIER3D( NCNOIR, XYZSOM(1,NS), NS )
         ENDDO

 50   ENDDO

C     TRACE DES ARETES D'AU MOINS 3 FACES
      DO K = 1, NBAR3F

         NS = NSAR3F(1,K)
         XYZ(1) = XYZSOM(1,NS)
         XYZ(2) = XYZSOM(2,NS)
         XYZ(3) = XYZSOM(3,NS)

         NS2 = NSAR3F(2,K)
         XYZ2(1) = XYZSOM(1,NS2)
         XYZ2(2) = XYZSOM(2,NS2)
         XYZ2(3) = XYZSOM(3,NS2)

         CALL XVEPAISSEUR( 5 )
         CALL TRAIT3D(  NCROUG, XYZ,  XYZ2 )
         CALL ENTIER3D( NCBLEU, XYZ,  NS   )
         CALL ENTIER3D( NCVERT, XYZ2, NS2  )

      ENDDO

C     TITRE ET TRACE EFFECTIF
      KTITRE = '      FACES de l''ETOILE,       SOUS-ETOILES, N1FEOC=   
     %               '
      WRITE( KTITRE(1:5),  '(I5)') NBITEM
      WRITE( KTITRE(26:28),'(I3)') NBSSET
      WRITE( KTITRE(54:57),'(I4)') N1FEOC
      CALL SANSDBL( KTITRE, L )
      KTITRE(L+1:L+9) = ' NBAR3F= '
      CALL SANSDBL( KTITRE, L )
      WRITE( KTITRE(L+1:L+2),'(I2)') NBAR3F
      CALL SANSDBL( KTITRE, L )
      CALL TRFINS( KTITRE(1:L) )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
