      SUBROUTINE TRIA1TRIA( NOFOTI, MXSOM,  NBSOM,  XYZSOM,
     %                      NTR,    MXTRIA, NBTRIA, NOTRIA,
     %                      NOPTAR, NDPTNTR,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRIANGULER UN TRIANGLE DEFINI PAR LA LISTE DES NUMEROS XYZSOM
C -----    DES POINTS DE SES 3 ARETES

C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION 'TAILLE_IDEALE' DES ARETES SINON 0
C MXSOM  : NOMBRE MAXIMAL DE SOMMETS DE LA TRIANGULATION
C NTR    : NUMERO NOTRIA DU TRIANGLE A TRIANGULER
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOTRIA
C NOPTAR : NUMERO XYZSOM DES POINTS DES 3 ARETES DU TRIANGLE NTR
C          NOPTAR( 1 ) = NUMERO XYZSOM DU 1-ER SOMMET DU TRIANGLE NTR
C          NOPTAR( NDPTNTR(3)+1 ) = NOPTAR(1)
C NDPTNTR: NUMERO NOPTAR DU DERNIER POINT SUR L'ARETE I DU TRIANGLE NTR
C          REMARQUE: NDPTNTR(3)=NBS NOMBRE DE POINTS DU CONTOUR FERME
C                    DU TRIANGLE NTR
C MODIFIES:
C ---------
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : X Y Z LES 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION
C NOTRIA : NUMERO XYZSOM DES 3 SOMMETS ET 0 POUR CHACUN DES NBTRIA TRIANGLES

C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR DETECTEE
C          1 SI TAILLE_IDEALE NON CALCULABLE EN UN SOMMET
C          2 SI TABLEAU XYZSOM SATURE (AUGMENTER MXSOM)
C          3 SI TABLEAU NOTRIA SATURE (AUGMENTER MXTRIA)
C          4 SI TABLEAU LAPILE SATURE (AUGMENTER MXPILE?)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC PARIS & VEULETTES SUR MER Avril 2017
C....................................................................012
      PARAMETER        (MXPILE=16384)
      PARAMETER        (MXTRIT=16384)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      REAL              XYZSOM(3,MXSOM)
      INTEGER           NOTRIA(4,MXTRIA),  NDPTNTR(3), NBPTAR(3),
     %                  NOPTAR(0: NDPTNTR(3) )

      INTEGER           LAPILE(2,MXPILE)
      INTEGER           LITRIT(MXTRIT)

C     COSMAX: COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C             LES 2 TRIANGLES SONT JUGES NON COPLANAIRES
C     ICI TOUS LES POINTS SONT DANS LE MEME PLAN => COSMAX=-1
      COSMAX = -1.0

      TRACTE0  = TRACTE
      NBSOM0   = NBSOM
      NBSOM00  = NBSOM
      NBTRIA0  = NBTRIA
      NBTRIA00 = NBTRIA

C     LE BARYCENTRE DU TRIANGLE NTR EST JOINT A TOUS LES POINTS DU CONTOUR
C     --------------------------------------------------------------------
C     NOMBRE DE POINTS DU CONTOUR FERME DU TRIANGLE NTR
      NBS = NDPTNTR(3)

C     NOMBRE DE POINTS SUR LES ARETES  SOMMETS EXTREMITES NON COMPRISES
      NBPTAR(1) = NDPTNTR(1) - 1
      NBPTAR(2) = NDPTNTR(2) - NDPTNTR(1) - 1
      NBPTAR(3) = NDPTNTR(3) - NDPTNTR(2) - 1

C     CONSTRUCTION DU BARYCENTRE NBSOM DU TRIANGLE NTR
      IF( NBSOM .GE. MXSOM ) GOTO 9992
      NBSOM = NBSOM + 1

C     CONSTRUCTION DES NBS TRIANGLES DE SOMMET NBSOM
      NS1 = NOTRIA( 1, NTR )
      NS2 = NOTRIA( 2, NTR )
      NS3 = NOTRIA( 3, NTR )

      DO K = 2, NBS

C        AJOUT DU NOUVEAU TRIANGLE NBTRIA
         IF( NBTRIA .GE. MXTRIA ) GOTO 9993
         NBTRIA = NBTRIA + 1
         NT = NBTRIA

         NOTRIA( 1, NT ) = NOPTAR( K )
         NOTRIA( 2, NT ) = NOPTAR( K+1 )
         NOTRIA( 3, NT ) = NBSOM

      ENDDO

C     NBSOM EST LE BARYCENTRE DES NBS POINTS NOPTAR
      DO K=1,3
         XYZSOM(K,NBSOM) = 0
      ENDDO
      DO M = 1, NBS
         NS = NOPTAR( M )
         DO K=1,3
            XYZSOM(K,NBSOM)= XYZSOM(K,NBSOM) + XYZSOM(K,NS)
         ENDDO
      ENDDO
      DO K=1,3
         XYZSOM(K,NBSOM) = XYZSOM(K,NBSOM) / NBS
      ENDDO

C     COIN DU SOMMET 1 DE NTR INITIAL
C     PREMIER NTR ET DERNIER NBTRIA DE LA TRIANGULATION DE NTR
      NOTRIA( 1, NTR ) = NS1
      NOTRIA( 2, NTR ) = NOPTAR(2)
      NOTRIA( 1, NBTRIA ) = NOPTAR(NBS)
      IF( NBPTAR(1) .GT. 0 ) THEN
         NOTRIA( 3, NTR    ) = NOPTAR(NBS)
         NOTRIA( 2, NBTRIA ) = NOPTAR(2)
      ELSE
         NOTRIA( 3, NTR    ) = NBSOM
         NOTRIA( 2, NBTRIA ) = NS1
      ENDIF
      NOTRIA( 3, NBTRIA ) = NBSOM

C     RECTIFICATION DES 2 SOUS-TRIANGLES DES 2 DERNIERS COINS DU TRIANGLE
C     CE QUI REVIENT A EVENTUELLEMENT ECHANGER LES 2 DIAGONALES
C     COIN DU SOMMET 2 DE NTR INITIAL
      NBPAJ = NBPTAR(1)
      IF( NBPAJ .GT. 0 ) THEN
         NT1 = NBTRIA00 + NBPAJ
         NT2 = NT1 + 1
      ELSE
         NT1 = NTR
         NT2 = NBTRIA00+1
      ENDIF
      CALL EC2D2TRI( COSMAX, XYZSOM, NT1, NT2, NOTRIA, MODIF )

C     COIN DU SOMMET 3 DE NTR INITIAL
      NBPAJ = NBPTAR(2)
      NT1 = NT2 + NBPAJ
      NT2 = NT1 + 1
      CALL EC2D2TRI( COSMAX, XYZSOM, NT1, NT2, NOTRIA, MODIF )

C     TRACE DE LA SOUS-TRIANGULATION AVEC CE MILIEU
      IF( NBTRIA .GT. NBTRIA00 ) THEN
         TRACTE = .TRUE.
         LITRIT( 1 ) = NTR
         NBT = MIN( NBTRIA-NBTRIA00, MXTRIT )
         DO K = 1, NBT
            LITRIT( 1+K ) = NBTRIA00 + K
         ENDDO
         NBT = 1 + NBT
         CALL TRTRIAN( 'tria1tria', XYZSOM, 4, MXTRIT, NBT, LITRIT,
     %                  NOTRIA )
         TRACTE = TRACTE0
      ENDIF

C     FAUT IL AJOUTER UN MILIEU A L'ARETE COMMUNE DE 2 TRIANGLES?
C     -----------------------------------------------------------
 5    NBECDG  = 0
      NBTRIA0 = NBTRIA
      NBSOM0  = NBSOM

      DO NTR1 = NBTRIA00+1, NBTRIA0

         DO 15 NAR1 = 1, 3

            LHPILE = 0
C           RECHERCHE DU TRIANGLE NTR2 D'ARETE NAR2 IDENTIQUE
C           A L'ARETE NAR1 DE NTR1
C           ECHANGER CETTE DIAGONALE POUR L'AUTRE SI ELLE MAXIMISE
C           LE MINIMUM DES QUALITES DES 2 COUPLES DE TRIANGLES
            CALL EC2D1AT( XYZSOM, NTR1, NAR1, 4, NBTRIA00+1, NBTRIA0,
     %                    NOTRIA, NTR2, NAR2, MODIF,
     %                    MXPILE, LHPILE, LAPILE )

            IF( MODIF .LT. 0 ) THEN
C              SATURATION DU TABLEAU LAPILE
               IERR = 4
               GOTO 9994
            ENDIF

            IF( MODIF .GT. 0 ) THEN

C              ECHANGE DES DIAGONALES DE NTR1+NTR2 EFFECTUE
               NBECDG = NBECDG + 1
               GOTO 10

            ENDIF

            IF( NTR2 .GT. 0 .AND. NAR2 .GT. 0 ) THEN

C              FAUT IL AJOUTER UN POINT SUR CETTE ARETE NAR1 DE NTR1?
               CALL AJMIAR( NOFOTI, NBSOM, MXSOM, XYZSOM,
     %                      NAR1,NTR1, NAR2,NTR2, NBTRIA, MXTRIA,NOTRIA,
     %                      MXPILE, LHPILE, LAPILE, IERR )
C              IERR=0 SI AJOUT D'UN MILIEU et PAS D'ERREUR DETECTEE
C                  -1 SI PAS D'AJOUT DE MILIEU ET SANS ERREUR DETECTEE
C                   1 SI TAILLE_IDEALE NON CALCULABLE EN UN SOMMET
C                   2 SI TABLEAU XYZSOM SATURE (AUGMENTER MXSOM)
C                   3 SI TABLEAU NOTRIA SATURE (AUGMENTER MXTRIA)

               IF( IERR .GT. 1 ) THEN
                  GOTO 9999
               ELSE IF( IERR .EQ. -1 .OR. IERR .EQ. 1 ) THEN
                  IERR = 0
                  GOTO 15
               ENDIF

C              UN MILIEU NBSOM A ETE AJOUTE
C              NTR1+NTR2 ONT DONNE 4 SOUS-TRIANGLES

C              TRACE DE LA SOUS-TRIANGULATION AVEC CE MILIEU
               IF( NBTRIA .GT. NBTRIA00 ) THEN
                  TRACTE = .TRUE.
                  NBT = MIN( NBTRIA-NBTRIA00, MXTRIT )
                  DO K = 1, NBT
                     LITRIT( K ) = NBTRIA00 + K
                  ENDDO
                  CALL TRTRIAN( 'tria1tria', XYZSOM, 4,
     %                           MXTRIT, NBT, LITRIT, NOTRIA )
                  TRACTE = TRACTE0
               ENDIF

            ENDIF

C           EXISTE T IL DES ARETES A ECHANGER COMME DIAGONALES DE 2 TRIANGLES?
 10         IF( LHPILE .GT. 0 ) THEN

               NT1 = LAPILE( 1, LHPILE )
               NA1 = LAPILE( 2, LHPILE )

C              L'ARETE NA1 DE NT1 EST DEPILEE
               LHPILE = LHPILE - 1

C              RECHERCHE DU TRIANGLE NT2 D'ARETE NA2 IDENTIQUE
C              A L'ARETE NA1 DE NT1
C              ECHANGER CETTE DIAGONALE POUR L'AUTRE SI ELLE MAXIMISE
C              LE MINIMUM DES QUALITES DES 2 COUPLES DE TRIANGLES
               CALL EC2D1AT( XYZSOM, NT1, NA1, 4, NBTRIA00+1, NBTRIA,
     %                       NOTRIA, NT2, NA2, MODIF,
     %                       MXPILE, LHPILE, LAPILE )

               IF( MODIF .NE. 0 ) THEN

C                 ECHANGE DES DIAGONALES DE NTR1+NTR2 EFFECTUE
                  NBECDG = NBECDG + 1

C                 TRACE DE LA SOUS-TRIANGULATION APRES ECHANGE
                  IF( NBTRIA .GT. NBTRIA00 ) THEN
                     TRACTE = .TRUE.
                     NBT = MIN( NBTRIA-NBTRIA00, MXTRIT )
                     DO K = 1, NBT
                        LITRIT( K ) = NBTRIA00 + K
                     ENDDO
                     CALL TRTRIAN( 'tria1tria', XYZSOM, 4,
     %                              MXTRIT, NBT, LITRIT, NOTRIA )
                     TRACTE = TRACTE0
                  ENDIF

               ENDIF

               GOTO 10

            ENDIF

 15      ENDDO

      ENDDO

      PRINT*
      PRINT*,'tria1tria:',NBSOM, ' SOMMETS   et',
     %        NBSOM -NBSOM0,' SOMMETS   AJOUTES'
      PRINT*,'tria1tria:',NBTRIA,' TRIANGLES et',
     %        NBTRIA-NBTRIA0,' TRIANGLES AJOUTES'
      PRINT*,'tria1tria:',NBECDG,' ECHANGES de DIAGONALES'

      IF( NBTRIA .GT. NBTRIA0 ) GOTO 5


C        CI_DESSOUS INUTILE CAR 0 ECHANGE
cccC     ESSAI FINAL D'ECHANGER LES DIAGONALES DES COUPLES DE SOUS-TRIANGLES
cccC     -------------------------------------------------------------------
ccc      DO ITER = 1, 4
ccc         NBECDG = 0
ccc         DO NTR1 = NBTRIA00+1, NBTRIA
ccc            DO NAR1 = 1, 3

cccC              RECHERCHE DU TRIANGLE NTR2 D'ARETE NAR2 IDENTIQUE
cccC              A L'ARETE NAR1 DE NTR1
cccC              ECHANGER CETTE DIAGONALE POUR L'AUTRE SI ELLE MAXIMISE
cccC              LE MINIMUM DES QUALITES DES 2 COUPLES DE TRIANGLES
ccc               CALL EC2D1AT( XYZSOM, NTR1, NAR1, 4, NBTRIA00+1, NBTRIA,
ccc     %                       NOTRIA, NTR2, NAR2, MODIF,
ccc     %                       MXPILE, LHPILE, LAPILE )
ccc               IF( MODIF .NE. 0 ) THEN
cccC                 ECHANGE DES DIAGONALES DE NTR1+NTR2 EFFECTUE
ccc                  NBECDG = NBECDG + 1
cccC                 TRACE DE LA SOUS-TRIANGULATION AVEC CE MILIEU
ccc                  IF( NBTRIA .GT. NBTRIA00 ) THEN
ccc                     TRACTE = .TRUE.
ccc                     NBT = MIN( NBTRIA-NBTRIA00, MXTRIT )
ccc                     DO K = 1, NBT
ccc                        LITRIT( K ) = NBTRIA00 + K
ccc                     ENDDO
ccc                     CALL TRTRIAN( 'tria1tria', XYZSOM, 4,
ccc     %                              MXTRIT, NBT, LITRIT, NOTRIA )
ccc                     TRACTE = TRACTE0
ccc                  ENDIF
ccc               ENDIF
ccc            ENDDO          
ccc         ENDDO
ccc         PRINT*,'tria1tria: ITERATION',ITER,':',NBECDG,
ccc     %          ' ECHANGES DE DIAGONALES *****************************'
ccc         IF( NBECDG .EQ. 0 ) GOTO 20
ccc      ENDDO
ccc 20   CONTINUE


C     FIN DU TRAITEMENT DU TRIANGLE NTR
C     ---------------------------------
 30   PRINT*
      PRINT*,'tria1tria: TRIANGLE',NTR,' :',NBSOM, ' SOMMETS   dont',
     %NBSOM-NBSOM00,  ' SOMMETS   AJOUTES'
      PRINT*,'tria1tria: TRIANGLE',NTR,' :',NBTRIA,' TRIANGLES dont',
     %NBTRIA-NBTRIA00,' TRIANGLES AJOUTES'

C     TRACE DE LA SOUS-TRIANGULATION FINALE
      IF( NBTRIA .GT. NBTRIA00 ) THEN
         TRACTE = .TRUE.
C        AJOUT DU TRIANGLE INITIAL PEUT ETRE MODIFIE
         LITRIT( 1 ) = NTR
         NBT = MIN( NBTRIA-NBTRIA00, MXTRIT )
         DO K = 1, NBT
            LITRIT( 1+K ) = NBTRIA00 + K
         ENDDO
         NBT = 1 + NBT
         CALL TRTRIAN( 'tria1tria', XYZSOM, 4,
     %                  MXTRIT, NBT, LITRIT, NOTRIA )
         TRACTE = TRACTE0
      ENDIF

      IERR = 0
      GOTO 9999

C     TABLEAU XYZSOM SATURE
 9992 PRINT *,'tria1tria: TABLEAU XYZSOM SATURE MXSOM=',MXSOM
      IERR = 2
      GOTO 9999

C     TABLEAU NOTRIA SATURE
 9993 PRINT*,'tria1tria: TABLEAU NOTRIA SATURE MXTRIA=',MXTRIA
      IERR = 3

C     TABLEAU LAPILE SATURE
 9994 PRINT*,'tria1tria: TABLEAU LAPILE SATURE MXPILE=',MXPILE
      IERR = 4
      GOTO 30

 9999 RETURN
      END
