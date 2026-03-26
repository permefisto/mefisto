      SUBROUTINE TRFETO7( PTXYZD,  N1FEOC, NFETOI, NBSSET, N1SSET,
     %                    NBARP2F, NSARP2F )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES FACES TRIANGULAIRES SIMPLES D'UNE ETOILE
C -----    ET DES SOUS-ETOILES N1SSET
C          AVEC LA VERSION 2 DU TABLEAU NFETOI

C ENTREES:
C --------
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
C NBSSET : NOMBRE DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET SOUS ETOILES
C          CHACUNE D'AU MOINS 4 FACES
C NBARP2F: NOMBRE D'ARETES APPARTENANT A PLUS DE 2 FACES
C NSARP2F: NUMERO DES 2 SOMMETS DES ARETES A PLUS DE 2 FACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& St PIERRE DU PERRAY Septembre 2015
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
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           N1SSET(NBSSET), NFETOI(5,*), NSARP2F(3,NBARP2F)

      CHARACTER*96      KTITRE
      INTRINSIC         MIN, MAX
      REAL              R, XYZ(3), XYZ2(3),
     %                  XYZBAR(3,MXITEM), DISTOE(MXITEM)
      INTEGER           NDITEM(0:MXSSET), NOSTTR(MXITEM), ITEMTR(MXITEM)
      INTEGER           NOSOTR(3)
      DOUBLE PRECISION  VECNOR(3), D

      IF( .NOT. TRACTE  ) RETURN
      IF( NBSSET .LE. 0 ) RETURN

ccc      print*,'trfeto7: nbsset=',nbsset,' n1sset=',(n1sset(K),k=1,nbsset)

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
               NS = ABS( NFETOI(1+K,NF1) )
               DO L=1,3
                  R = REAL( PTXYZD(L,NS) )
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

ccc      print*,'trfeto7: nbitem=',NBITEM,' cooext=',cooext

 4    IF( NBITEM .GE. MXITEM ) THEN
         PRINT*,'trfeto7: AUGMENTER MXITEM=',MXITEM
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
      PREDUF  = 16.0
      LORBITE = 1

ccc      print*,'trfeto7: AXOPTV=',AXOPTV,' AXOEIL=',AXOEIL,
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
            IF( NF .LE. NDITEM(M) ) THEN
C              POURCENTAGE DE REDUCTION DE LA FACE
C              SELON SON NO DE SOUS-ETOILE
               PREDUF = 16 - 2*M
               GOTO 30
            ENDIF
         ENDDO

C        TRIANGLE NT D'UNE SOUS ETOILE
C        COULEUR DE LA FACE SELON L'ORDRE
C         0:NCNOIR,  1:NCROUG,  2:NCVERT,  3:NCBLEU,
C         4:NCCYAN,  5:NCJAUN,  6:NCMAGE,  7:NCBLAN,
C         8:NCGRIS,  9:NCGRIM, 10:NCGRIC, 11:NCBEIG,
C        12:NCORAN, 13:NCSAUM, 14:NCROSE, 15:NCTURQ

 30      IF( M .EQ. 0 ) THEN
C           FACES ISSUES DE N1FEOC EN GRIS CLAIR
            NCF = NCGRIC
         ELSE
C           AUTRES FACES AVEC UNE SEULE COULEUR POUR CHAQUE SOUS-ETOILE
            NCF = 3+M
            IF( NCF .GT. 15 ) THEN
               NCF = MOD( NCF, 15) + 3
            ENDIF
         ENDIF

C        3 SOMMETS DE LA FACE NT DE NFETOI
         NOSOTR(1) = ABS( NFETOI(2,NT) )
         NOSOTR(2) = NFETOI(3,NT)
         NOSOTR(3) = NFETOI(4,NT)

C        TRACE DU TRIANGLE NOSOTR
         CALL TRFATR( NCF, NCGRIS, NOSOTR, PTXYZD )

C        TRACE DE LA NORMALE A LA FACE EN SON BARYCENTRE
         CALL VECNOR3D( PTXYZD(1,NOSOTR(1)), PTXYZD(1,NOSOTR(2)),
     %                  PTXYZD(1,NOSOTR(3)), VECNOR )

C        NORME DU VECTEUR NORMAL
         D = SQRT( VECNOR(1)**2 + VECNOR(2)**2 + VECNOR(3)**2 )
         DO L=1,3
            XYZ(L) = XYZBAR(L,NF) + REAL( VECNOR(L) / D * DISMOY )
         ENDDO
         CALL SYMBOLE3D( NCMAGE, XYZBAR(1,NF), '*' )
         CALL TRAIT3D(   NCMAGE, XYZBAR(1,NF), XYZ )

C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO L=1,3
            NS = NOSOTR(L)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
            CALL ENTIER3D( NCGRIS, XYZ, NS )
         ENDDO

 50   ENDDO

C     TRACE DES ARETES D'AU MOINS 3 FACES
      DO K = 1, NBARP2F

         NS = NSARP2F(1,K)
         XYZ(1) = REAL( PTXYZD(1,NS) )
         XYZ(2) = REAL( PTXYZD(2,NS) )
         XYZ(3) = REAL( PTXYZD(3,NS) )

         NS2 = NSARP2F(2,K)
         XYZ2(1) = REAL( PTXYZD(1,NS2) )
         XYZ2(2) = REAL( PTXYZD(2,NS2) )
         XYZ2(3) = REAL( PTXYZD(3,NS2) )

         CALL XVEPAISSEUR( 5 )
         CALL TRAIT3D(  NCROUG, XYZ,  XYZ2 )
         CALL ENTIER3D( NCNOIR, XYZ,  NS   )
         CALL ENTIER3D( NCNOIR, XYZ2, NS2  )

      ENDDO

C     TITRE ET TRACE EFFECTIF
      KTITRE = '      FACES de l''ETOILE,       SOUS-ETOILES, N1FEOC=   
     %               dans souseto.'
      WRITE( KTITRE(1:5),  '(I5)') NBITEM
      WRITE( KTITRE(26:28),'(I3)') NBSSET
      WRITE( KTITRE(54:57),'(I4)') N1FEOC
      CALL SANSDBL( KTITRE, L )
      KTITRE(L+1:L+10) = ' NBARP2F= '
      CALL SANSDBL( KTITRE, L )
      WRITE( KTITRE(L+1:L+2),'(I2)') NBARP2F
      CALL SANSDBL( KTITRE, L )
      CALL TRFINS( KTITRE(1:L) )
C
C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
