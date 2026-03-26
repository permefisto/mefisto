      SUBROUTINE CIASFAET( NBVPSI,  NFVPSI,  NFETOI,
     %                     NBCIAS,  MXCIAS,  N1CIAS,
     %                     MXASFVP, MIARSICI, NBASFVP, NSASFVP, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DECTECTER ET CONSTRUIRE LES CYCLES D'ARETES SIMPLES
C -----    D'UNE SELECTION DE FACES SIMPLES D'UNE ETOILE

C ENTREES:
C --------
C NBVPSI : NOMBRE DE FACES SELECTIONNEES DANS NFETOI
C NFVPSI : NUMERO NFETOI DES NBVPSI FACES
C NFETOI : LES FACES TRIANGULAIRES DE L'ETOILE EN VERSION 2
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C MXCIAS : NOMBRE MAXIMAL DE CYCLES
C MXFASFVP: NOMBRE MAXIMAL DE FACES SIMPLES DES CYCLES
C MIARSICI: NOMBRE MINIMAL ACTUEL D'ARETES SIMPLES DES CYCLES
C           DES FACES NFVPSI POUR EVITER DES CALCULS INUTILES
C           SUR DES CYCLES TROP LONGS

C SORTIES:
C --------
C NBCIAS : NOMBRE DE CYCLES DES ARETES SIMPLES DES FACES NFVPSI DANS NFETOI
C N1CIAS : NO DANS NSASFVP DE LA PREMIERE ARETE DES NBCIAS CYCLES
C NBASFVP: NOMBRE DE FACES DES ARETES SIMPLES DES NBCIAS CYCLES
C NSASFVP: 1: NUMERO DU SOMMET 1  DE L'ARETE SIMPLE
C          2: NUMERO DU SOMMET 2  DE L'ARETE SIMPLE
C          3: NUMERO DANS NSASFVP DE L'ARETE SIMPLE SUIVANTE
C          4: NUMERO DANS NFVPSI DE LA FACE DE CETTE ARETE SIMPLE
C IERR   : =9 MXASFVP TROP PETIT. A AUGMENTER
C          =8 MXCIAS  TROP PETIT. A AUGMENTER
C          =1 CYCLE MAL FERME
C          >0 NUMERO DU CYCLE DE MOINS DE 3 ARETES
C          =2 TROP D'ARETES SUR LE BORD (>MIARSICI)
C          =0 PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  LJLL UPMC & St PIERRE du PERRAY Octobre 2017
C MODIFS : ALAIN PERRONNET  Veulettes sur mer               Janvier 2020
C2345X7..............................................................012
      INTEGER  NFVPSI(NBVPSI), NFETOI(5,*),
     %         N1CIAS(MXCIAS), NSASFVP(4,MXASFVP)

      IERR   = 0
      NBCIAS = 0

C     RECENSEMENT DES ARETES SIMPLES DES NBVPSI FACES NFVPSI DE NFETOI
C     ----------------------------------------------------------------
      NBASFVP = 0
      DO K=1,NBVPSI

C        NUMERO DANS NFETOI DE LA FACE K DE NFVPSI
         NF0 = NFVPSI( K )

C        LE 2 SOMMETS DE L'ARETE L DE NF0 SONT NS1-NS2
C        DANS LE SENS DES ARETES DU TRIANGLE NF0 DEFINISSANT LE
C        VECTEUR NORMAL A LA FACE
         NS1 = NFETOI(4,NF0)
         DO L=1,3
            NS2 = NFETOI(1+L,NF0)

C           RECHERCHE DE L'ARETE NS1-NS2 PARMI LES AUTRES ARETES
C           DES NBVPSI FACES NFVPSI
            DO M=1,NBVPSI
               NF1 = NFVPSI( M )
               IF( NF1 .NE. NF0 ) THEN
                  NS3 = NFETOI(4,NF1)
C                 SENS DE L'ARETE COMMUNE
                  DO N=1,3
                     NS4 = NFETOI(1+N,NF1)
                     IF( (NS4 .EQ. NS1 .AND. NS3 .EQ. NS2) .OR.
     %                   (NS4 .EQ. NS2 .AND. NS3 .EQ. NS1) ) THEN
C                        L'ARETE NS1-NS2 EST DOUBLE PARMI
C                        LES NBVPSI TRIANGLES
                         GOTO 10
                     ENDIF
                     NS3 = NS4
                  ENDDO
               ENDIF
            ENDDO

C           L'ARETE NS1-NS2 DE NF0 EST SIMPLE PARMI LES NBVPSI TRIANGLES
            IF( NBASFVP .GE. MXASFVP ) THEN
               PRINT*,'ciasfaet: AUGMENTER MXASFVP=',MXASFVP
               IERR = 9
               GOTO 9999
            ENDIF
            NBASFVP = NBASFVP + 1
C           LE NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE
C           DANS LE SENS INVERSE DE L'ARETE DE NF0 et
C           POUR LE FUTUR TRIANGLE NFETOI A AJOUTER
C           L'ARETE COMMUNE DE 2 FACES SIMPLES EST PARCOURUE DANS
C           LES 2 SENS POUR DESIGNER LA NORMALE VERS LA DIRECTION
C           VERS L'INTERIEUR DE L'ETOILE
            NSASFVP(1,NBASFVP) = NS2
            NSASFVP(2,NBASFVP) = NS1
C           ARETE SIMPLE SUIVANTE DANS UN CYCLE -1 NON ACTIVE
            NSASFVP(3,NBASFVP) = -1
C           NUMERO DE SA FACE DANS NFSVPI
            NSASFVP(4,NBASFVP) = K

            PRINT *,'ciasfaet: ARETE SIMPLE',NBASFVP,' :',
     %              (NSASFVP(N,NBASFVP),N=1,4),' des',
     %               NBVPSI,' FACES+XYZ a V>0'

C           PASSAGE A L'ARETE SUIVANTE DE NF0
 10         NS1 = NS2
         ENDDO

      ENDDO

C     LE BORD D'ARETES A T IL TROP D'ARETES?
      IF( NBASFVP .GT. MIARSICI ) THEN
C        OUI: ABANDON DES CALCULS SUR LES CYCLES
         NBCIAS = 0
         GOTO 9999
      ENDIF


C     RECHERCHE DES CYCLES DE CES NBASFVP ARETES SIMPLES
C     --------------------------------------------------
 20   DO NAC0=1,NBASFVP
         IF( NSASFVP(3,NAC0) .LT. 0 ) GOTO 24
      ENDDO
C     TOUTES LES ARETES SIMPLES ONT ETE RECENSEES DANS NBCIAS CYCLES
      GOTO 30

C     NAC0 1-ERE ARETE NS1 NS2 NON RECENSEE -> 1-ERE DU NOUVEAU CYCLE
 24   NS1 = NSASFVP(1,NAC0)
      NS2 = NSASFVP(2,NAC0)
      IF( NBCIAS .GE. MXCIAS ) THEN
         PRINT*,'ciasfaet: AUGMENTER MXCIAS=',MXCIAS
         IERR = 10
         GOTO 9999
      ENDIF
C     UN NOUVEAU CYCLE NBCIAS EST CREE
      NBCIAS = NBCIAS + 1
C     D'ARETE INITIALE NAC0
      N1CIAS( NBCIAS ) = NAC0
C     DE SOMMET INITIAL NS0
      NS0 = NS1

C     RECHERCHE DE L'ARETE SUIVANTE DU CYCLE
 26   DO NAC1=1,NBASFVP

         IF( NAC1 .NE. NAC0 ) THEN
            IF( NSASFVP(3,NAC1) .LT. 0 ) THEN

C              PARCOURS DES 2 SOMMETS DE L'ARETE NAC1
               DO M=1,2
                  NS3 = NSASFVP(M,NAC1)
                  IF( NS3 .EQ. NS2 ) THEN

C                    L'ARETE NAC0 EST SUIVIE PAR L'ARETE NAC1
                     NSASFVP(3,NAC0) = NAC1
C                    L'ARETE NAC1 EST LA DERNIERE ACTUELLE DU CYCLE
                     NSASFVP(3,NAC1) = 0

C                    NS4 L'AUTRE SOMMET DE L'ARETE NAC1
                     IF( M .EQ. 1 ) THEN
                        NS4 = NSASFVP(2,NAC1)
                     ELSE
                        NS4 = NSASFVP(1,NAC1)
                     ENDIF

C                    NS3 DEVIENT LE PREMIER SOMMET DE L'ARETE NAC1
                     NSASFVP(1,NAC1) = NS3
C                    NS4 DEVIENT LE SECOND  SOMMET DE L'ARETE NAC1
                     NSASFVP(2,NAC1) = NS4

C                    LA NOUVELLE ARETE DE RECHERCHE DU CYCLE
                     NAC0 = NAC1
C                    ELLE SE TERMINE AU SOMMET NS2
                     NS2 = NS4

                     IF( NS2 .EQ. NS0 ) THEN
C                       LE CYCLE NBCIAS EST FERME AU SOMMET NS0
C                       LE CYCLE DEVIENT UN CYCLE:
C                       L'ARETE SUVANTE DE LA DERNIERE NAC1
C                       EST LA PREMIERE
                        NSASFVP(3,NAC1) = N1CIAS( NBCIAS )
                        GOTO 20
                     ENDIF

                     GOTO 26

                  ENDIF
               ENDDO

            ENDIF

         ENDIF

      ENDDO

C     RECHERCHE D'UN AUTRE CYCLE
      GOTO 20


C     TOUTES LES ARETES SIMPLES ONT ETE RECENSEES DANS NBCIAS CYCLES
C     AFFICHAGE DES NBCIAS CYCLES AVANT DETECTION DE LEURS CYCLES INTERNES
C     ---------------------------------------------------------------------
 30   CALL AFCYCLE( NBCIAS, N1CIAS, NSASFVP )


C     DETECTION D'EVENTUELS CYCLES DANS CHACUN DES NBCIAS CYCLES
C     ----------------------------------------------------------
      IF( NBASFVP .LE. 5 ) GOTO 85
C     AU MINIMUM 6 ARETES PEUVENT PRESENTER 2 CYCLES DE 3 ARETES
C     MAIS C'EST IMPOSSIBLE POUR AU PLUS 5 ARETES
      NAC2 = 0
      NAC3 = 0
      NOC  = 0

 50   NOC = NOC + 1
      IF( NOC .LE. NBCIAS ) THEN

C        LA PREMIERE ARETE N1AC DU CYCLE NOC
         N1AC  = N1CIAS( NOC )
         NA1NS = N1AC
         NS0   = 0

C        NS 1-ER SOMMET A CONTROLER DANS LE CYCLE NOC
         NS = NSASFVP( 1, NA1NS )

C        SI NS APPARTIENT A 2n ARETES ALORS NS APPARTIENT A n CYCLES
C        -----------------------------------------------------------
C        RECHERCHE DU NOMBRE D'ARETES DE SOMMET NS
C        NS EST 1 FOIS DANS L'ARETE NA1NS
 60      NBANS = 1
         NAC1  = NA1NS
         NAC   = NAC1

C        PARCOURS DES ARETES NAC SUIVANTES POUR RETROUVER NS
 70      NAC = NSASFVP( 3, NAC )
         IF( NAC .NE. NA1NS ) THEN

            IF( NS.EQ.NSASFVP(2,NAC) .OR. NS.EQ.NSASFVP(1,NAC) ) THEN

C              NS EST RETROUVE DANS L'ARETE NAC1
               NBANS = NBANS + 1
               GOTO( 71, 72, 73, 74 ) NBANS

C              NAC1 EST L'ARETE OU NS APPARAIT POUR LA SECONDE FOIS
 71            NAC1 = NAC
               GOTO 70

C              NAC2 EST L'ARETE OU NS APPARAIT POUR LA SECONDE FOIS
 72            NAC2 = NAC
               GOTO 70

C              NAC3 EST L'ARETE OU NS APPARAIT POUR LA TROISIEME FOIS
 73            NAC3 = NAC
               GOTO 70

C              NAC4 EST L'ARETE OU NS APPARAIT POUR LA QUATRIEME FOIS
 74            NAC4 = NAC

C              IL EXISTE UN CYCLE DE L'ARETE NAC1 A NAC2 et
C              LE CYCLE INITIAL NOC EST REDUIT DE NAC4 A NAC3
C              ----------------------------------------------
C              LE CYCLE INITIAL NOC EST AMPUTE DE NAC1 A NAC2
               NSASFVP(3,NAC4) = NAC3

C              CREATION DU NOUVEAU CYCLE NAC1 A NAC2
               IF( NBCIAS .GE. MXCIAS ) THEN
                  PRINT*,'ciasfaet: AUGMENTER MXCIAS=',MXCIAS
                  GOTO 9999
               ENDIF
               NBCIAS = NBCIAS + 1
C              PREMIERE ARETE PARTANT DE NS D'ARETE NAC1
               N1CIAS( NBCIAS ) = NAC1

C              FERMETURE DU NOUVEAU CYCLE
               NSASFVP(3,NAC2) = NAC1

C              FERMETURE DU CYCLE INITIAL
               NSASFVP(3,NAC4) = NAC3

C              2 RENCONTRES DE NS EN MOINS
C              LE CYCLE DE NS REDEMMARRE AVEC L'ARETE NAC3
               NA1NS = NAC3
               GOTO 60

            ENDIF

            GOTO 70

         ENDIF

         IF( NBANS .NE. 2 ) THEN
           PRINT*,'ciasfaet: PB le SOMMET NS=',NS,' APPARTIENT A',NBANS,
     %            ' ARETES au lieu de 2 -> ABANDON'
            CALL AFCYCLE( NBCIAS, N1CIAS, NSASFVP )
            NBCIAS = 0
            IERR   = 1
            GOTO 9999
         ENDIF

C        PASSAGE AU SOMMET NS SUIVANT
         NS0 = NS

C        PASSAGE A L'ARETE SUIVANTE DE NA1NS
         NA1NS = NSASFVP( 3, NA1NS )
         IF( NA1NS .EQ. N1AC ) GOTO 50

C        LE NOUVEAU SOMMET A TESTER SOMMET DE 2n ARETES?
         NS = NSASFVP( 1, NA1NS )

         IF( NS .EQ. NS0 ) THEN
C           POUR EVITER DE TRAITER A NOUVEAU NS0
C           PASSAGE AU SECOND SOMMET:  NE DEVRAIT PAS SE PRODUIRE...
C           L'ORIENTATION DE NFETOI POSE PROBLEME ...
            NS = NSASFVP( 2, NA1NS )
            PRINT*,'ciasfaet: NE DEVRAIT PAS SE PRDODUIRE! NA1NS=',
     %              NA1NS,' NS0=',NS0,' NS=',NS
            PRINT*
         ENDIF

         GOTO 60

      ENDIF


C     AFFICHAGE DES NBCIAS CYCLES APRES DETECTION DES CYCLES INTERNES
C     ----------------------------------------------------------------
 85   CALL AFCYCLE( NBCIAS, N1CIAS, NSASFVP )


C     TRANSFORMATION DE CHAQUE CYCLE D'UN CYCLE EN UNE LISTE CHAINEE
C     SUR L'ARETE SUIVANTE ET FINIE PAR ZERO
C     --------------------------------------------------------------
      DO NOC = 1, NBCIAS

C        AFFICHAGE DU CYCLE NOC
         N1AC = N1CIAS( NOC )
         NBAC = 1
         NAC  = N1AC

C        L'ARETE SUIVANTE
 90      NAC = NSASFVP( 3, NAC )

         IF( NAC .EQ. N1AC ) THEN
C           FIN DU CYCLE NOC INITIALISE
            NSASFVP( 3, NAC ) = 0
            GOTO 100
         ENDIF

C        UNE ARETE DE PLUS POUR CE CYCLE
         NBAC = NBAC + 1
         GOTO 90

C        FIN DU CYCLE NOC
 100     IF( NBAC .LT. 3 ) THEN
            PRINT*,'ciasfaet: ERREUR CYCLE',NOC,' AVEC SEULEMENT',
     %              NBAC,' ARETES?...'
            NBCIAS = 0
            IERR   = NOC
            GOTO 9999
         ENDIF

      ENDDO

 9999 RETURN
      END
