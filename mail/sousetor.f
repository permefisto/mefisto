      SUBROUTINE SOUSETOR( XYZSOM, N1FEOC,  NFETOI,
     %                     MXSSET, NBSSET1, N1SSET, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :  DETERMINER LE NOMBRE ET LES FACES DES SOUS-ETOILES DE L'ETOILE
C ------  ISSUE DES FACES SIMPLES D'UNE ETOILE DEFINIE PAR SON POINTEUR
C         N1FEOC SUR LE TABLEAU NFETOI DE L'ETOILE
C         ET LES AJOUTER AUX SOUS-ETOILES INITIALES

C         UNE SOUS-ETOILE EST TELLE QUE TOUTE ARETE DE SES FACES
C         . APPARTIENT A 2 ET SEULEMENT 2 DE SES FACES
C         . EST PARCOURUE DANS DES SENS CONTRAIRES DANS CES 2 FACES
C         . LE VECTEUR NORMAL A CHAQUE FACE EST DIRIGEE VERS L'INTERIEUR
C           DE LA SOUS ETOILE

C          VERSION XYZSOM DES SOMMETS EN SIMPLE PRECISION

C ENTREES:
C --------
C XYZSOM : PAR POINT : X  Y  Z
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 3 LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO XYZSOM DU SOMMET 1 DE LA FACE
C          2: NUMERO XYZSOM DU SOMMET 2 DE LA FACE
C          3: NUMERO XYZSOM DU SOMMET 3 DE LA FACE
C          4: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          S1S2 x S1S3 VERS L'INTERIEUR DE L'ETOILE
C MXSSET : NOMBRE MAXIMAL DE SOUS ETOILES DECLARABLES DANS N1SSET
C NBSSET1: NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI

C SORTIES:
C --------
C N1FEOC : =0 POINTEUR SUR LA PREMIERE FACE DE L'ETOILE SI PAS DE PB
C          C-A-D QUE TOUTES LES FACES FONT PARTIE DES SOUS-ETOILES
C NBSSET1: NOMBRE FINAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET SOUS ETOILES
C IERR   : 1 SI LE TABLEAU LAPILE EST SATURE
C          2 SI N1FEOC=0 OU NBSSET1=0 EN SORTIE
C          0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& St PIERRE du PERRAY   Fevrier 2016
C2345X7..............................................................012
      PARAMETER        (MXPILE=2048)
      INTEGER           LAPILE(MXPILE)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE
      INTEGER           NFETOI(5,*), N1SSET(MXSSET), NSAR3F(2,MXPILE),
     %                  NBFASI(3), NFADSI(8,3),
     %                  NBFASM(3), NFADSM(8,3)
      REAL              XYZSOM(3,*)
      DOUBLE PRECISION  ANGLE, ANGMIN

ccc      print*,'Debut sousetor: NB SousEtoiles NBSSET0=',NBSSET0,
ccc     %       ' NB FACES de l''ETOILE NBFETO0=',NBFETO0

cccC     TRACE DES FACES DE LA SOUS-ETOILE ISSUE DE N1FEOC
cccC     SANS LES AUTRES SOUS-ETOILES
ccc      CALL TRFETO7R( XYZSOM, N1FEOC, NFETOI, 0, N1SSET,
ccc     %               0,      NSAR3F )

C     INITIALISATION DES SOUS ETOILES ISSUES DE N1FEOC
      IERR     = 0
      N1FEOC0  = N1FEOC
      NBSSET00 = NBSSET1
      NBPASS   = 0
      NBAR3F   = 0

 5    IF( N1FEOC .LE. 0 ) GOTO 9000

C     NOMBRE DE PASSAGE ICI
      NBPASS = NBPASS + 1

C     NOMBRE DE SOUS-ETOILES INITIALES DANS N1SSET-NFETOI
      NBSSET0 = NBSSET1

C     MARQUAGE SUPPRIME DES FACES DE LA SOUS-ETOILE ISSUE DE N1FEOC
C     ET NOMBRE DE SES FACES
      NBFETO = 0
      NF1    = N1FEOC
 10   IF( NF1 .GT. 0 ) THEN
         NBFETO = NBFETO + 1
         NFETOI(1,NF1) = ABS( NFETOI(1,NF1) )
         NF1 = NFETOI(5,NF1)
         GOTO 10
      ENDIF

      NBFETO0 = NBFETO

ccc      print*,'sousetor: NBPASS=',NBPASS,' N1FEOC=',N1FEOC,
ccc     %       ' NB INITIAL de Sous-Etoiles=',NBSSET0,
ccc     %       ' NB FACES de la Sous-Etoile=',NBFETO

C     INITIALISATION DE LA PILE DES FACES PAR LES ARETES
      LAPILE(1) = N1FEOC
      LHPILE = 1

 20   IF( LHPILE .GT. 0 ) THEN

C        TRAITEMENT DE LA FACE NF1 EN HAUT DE PILE
         NF1 = LAPILE( LHPILE )

C        LA FACE EST DEPILEE
         LHPILE = LHPILE - 1

C        CETTE FACE EST ELLE DEJA TRAITEE?
         IF( NFETOI(1,NF1) .LE. 0 ) GOTO 20

C        NS1-NS2 SOMMETS DE L'ARETE K-1 DE LA FACE NF1
         NS1 = NFETOI(3,NF1)
         DO 28 K=1,3

C           ARETE K-1 DE LA FACE NF1 et NS1-NS2 SES SOMMETS
            NS2 = ABS( NFETOI(K,NF1) )

C           LE TROISIEME SOMMET DE NF1
            IF( K .LT. 3 ) THEN
               KK = K+1
            ELSE
               KK = 1
            ENDIF
            NS3 = ABS( NFETOI(KK,NF1) )

C           RECHERCHE DANS L'ETOILE DES FACES D'ARETE NS1-NS2
C           PARCOURUE EN SENS INVERSE et LA PLUS PROCHE DE NF1
            NBFASI(K) = 0
            NBFASM(K) = 0

C           DEBUT DE RECHERCHE PAR LES FACES RESTANTES DE L'ETOILE
            MSSET = NBSSET0
            NF2   = N1FEOC

 22         IF( NF2 .GT. 0 ) THEN
               IF( ABS(NFETOI(1,NF2)) .GT. 0 .AND. NF2 .NE. NF1 ) THEN

C                 LA FACE NF2 N'EST PAS TRAITEE ET N'EST PAS NF1
C                 PARCOURS DES 3 ARETES DE NF2
                  NSS1 = NFETOI(3,NF2)
                  DO KK=1,3
                     NSS2 = ABS( NFETOI(KK,NF2) )
                     IF( NS1.EQ.NSS2 .AND. NS2.EQ.NSS1 ) THEN

C                       ARETE NS1-NS2 VUE DANS LE SENS OPPOSE
C                       NOMBRE DE FACES D'ARETE K PARCOURUE DANS LE SENS INVERSE
                        NBFASI(K) = NBFASI(K) + 1
                        NFADSI( NBFASI(K), K ) = NF2

                     ELSE IF( NS1.EQ.NSS1 .AND. NS2.EQ.NSS2 ) THEN

C                       ARETE NS1-NS2 VUE DANS DE LE MEME SENS
C                       NOMBRE DE FACES D'ARETE K PARCOURUE DANS LE MEME SENS
                        NBFASM(K) = NBFASM(K) + 1
                        NFADSM( NBFASM(K), K ) = NF2

                     ENDIF
                     NSS1 = NSS2
                  ENDDO
               ENDIF
               NF2 = NFETOI( 5, NF2 )
               GOTO 22
            ENDIF

            MSSET = MSSET + 1
            IF( MSSET .LE. NBSSET1 ) THEN
C              DANS LA SOUS-ETOILE MSSET NOUVELLEMENT CONSTRUITE
C              RECHERCHE DE LA FACE D'ARETE NS1-NS2 PARCOURUE EN SENS INVERSE
               NF2 = N1SSET( MSSET )
               GOTO 22
            ENDIF

ccc            print*,'sousetor:',NS1,NS2,' arete k=',k-1,' de la face',nf1,
ccc     %          ' de sommets:',(NFETOI(KK,NF1),KK=1,3),' avec NBFASI=',
ccc     %           NBFASI(K),' et NBFASM=',NBFASM(K)
ccc            print*,'sousetor: Faces en SENS INVERSE:',(NFADSI(KK,K),
ccc     %            (NFETOI(KKK,NFADSI(KK,K)),KKK=1,3),KK=1,NBFASI(K))
ccc            print*,'sousetor: Faces en MEME    SENS:',(NFADSM(KK,K),
ccc     %            (NFETOI(KKK,NFADSM(KK,K)),KKK=1,3),KK=1,NBFASM(K))

C           NOMBRE DE FACES D'ARETE K PARCOURUE DANS LE SENS INVERSE
            NB = NBFASI( K )
            IF( NB .LE. 0 ) THEN
               PRINT*,'sousetor: FACE',NF1,
     %                ' avec 0 FACE ADJACENTE PARCOURUE EN SENS INVERSE'
               GOTO 27
            ENDIF

C           NUMERO DE LA PREMIERE FACE ADJACENTE EN SENS INVERSE DE NS1-NS2
            NF = NFADSI(1,K)
            IF( NB .EQ. 1 ) THEN

               IF( NBFASI(K) .NE. NBFASM(K)+1 ) THEN
C                 AU MOINS UNE PAIRE DE FACES PARCOURUES DANS LE MEME SENS
C                 LA FACE OPPOSEE A DEJA ETE TRAITEE PARMI LES AU MOINS 4 FACES DE L'ARETE
C                 CETTE ARETE POUR CETTE FACE EST TRAITEE
C                 PASSAGE A SON ARETE SUIVANTE
                  GOTO 27
               ENDIF

C              L'UNIQUE FACE NF D'ARETE PARCOURUE EN SENS INVERSE
C              DE NS1-NS2 EST EMPILEE SI ELLE N'EST PAS DEJA TRAITEE
               IF( NFETOI(1,NF) .GT. 0 ) THEN
C                 NF EST EMPILEE CAR NON TRAITEE
                  IF( LHPILE .LT. MXPILE ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE( LHPILE ) = NF
                  ELSE
C                    TABLEAU LAPILE SATURE
                     PRINT *,'sousetor: TABLEAU LAPILE SATURE'
                     PRINT *,'AUGMENTER DANS sousetor.f MXPILE=',MXPILE
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

            ELSE

C              IL EXISTE PLUSIEURS FACES ADJACENTES A NS1-NS2 DANS LE SENS INVERSE
C              IL FAUT CHOISIR CELLE NFAMIN QUI FAIT UN ANGLE MINIMAL AVEC NF1
C              -------------------------------------------------------------------
C              UNE ARETE DE PLUS DANS AU MOINS 3 FACES
               NBAR3F = NBAR3F + 1
               NSAR3F(1,NBAR3F) = NS1
               NSAR3F(2,NBAR3F) = NS2

               ANGMIN = 16.D0
               NFAMIN = 0
               DO 25 L=1,NBFASI(K)

C                 LA FACE NF D'ARETE NS1 NS2 PARCOURUE EN SENS INVERSE
                  NF = NFADSI(L,K)

C                 RECHERCHE DU 3-EME SOMMET DE NF NON NS1 et NON NS2
                  DO M=1,3
                     NS4 = ABS( NFETOI(M,NF) )
                     IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 24
                  ENDDO

C                 CALCUL DE L'ANGLE ENTRE LES PLANS DES 2 FACES
C                 NF1 et NF D'ARETE COMMUNE NS1-NS2
 24               CALL ANG2TR3R( XYZSOM(1,NS1), XYZSOM(1,NS2),
     %                           XYZSOM(1,NS3), XYZSOM(1,NS4),
     %                           ANGLE, IERR )
                  IF( IERR .NE. 0 ) THEN
                     PRINT*,'sousetor: FACE NF1=',NF1,' St',NS1,NS2,NS3,
     %                      ' DEGENEREE'
                     GOTO 25
                  ENDIF

                  IF( ANGLE .LT. ANGMIN ) THEN
                     ANGMIN = ANGLE
                     NFAMIN = NF
                  ENDIF

 25            ENDDO

               IF( NFETOI(1,NFAMIN) .GT. 0 ) THEN
C                 NFAMIN EST EMPILEE
                  IF( LHPILE .LT. MXPILE ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE( LHPILE ) = NFAMIN
                  ELSE
C                    TABLEAU LAPILE SATURE
                     PRINT*,'sousetor: tableau LAPILE SATURE'
                     PRINT*,'AUGMENTER MXPILE=',MXPILE
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

            ENDIF

C           ARETE K SUIVANTE DE LA FACE NF1
 27         NS1 = NS2

 28      ENDDO

C        TRAITEMENT DE LA FACE NF1 VIS A VIS DES SOUS-ETOILES
C        ----------------------------------------------------
C        QUELLE EST LA FACE NF0 QUI PRECEDE NF1 DANS LE CHAINAGE
C        DES FACES DE L'ETOILE?
         NF0 = 0
         NF  = N1FEOC
 40      IF( NF .GT. 0 ) THEN
            IF( NF .EQ. NF1 ) GOTO 50
            NF0 = NF
            NF  = NFETOI( 5, NF )
            GOTO 40
         ENDIF

C        NF1 QUITTE LE CHAINAGE N1FEOC DES FACES DE L'ETOILE
 50      NF2 = NFETOI( 5, NF1 )
         IF( NF0 .GT. 0 ) THEN
            NFETOI( 5, NF0 ) = NF2
         ELSE
C           NF1 ETAIT LA PREMIERE FACE DE L'ETOILE
            N1FEOC = NF2
         ENDIF

C        LA FACE NF1 EST AJOUTEE AU DEBUT DE LA SOUS ETOILE NBSSET1
         IF( NBSSET1 .EQ. NBSSET0 ) THEN
C           NF1 DEBUTE LA SOUS-ETOILE NBSSET0+1
            NBSSET1 = NBSSET1 + 1
C           SOUS-ETOILE NBSSET1 COMPOSEE DE LA SEULE FACE NF1
            NFETOI( 5, NF1 ) = 0
         ELSE
C           LA PREMIERE FACE DE LA SOUS ETOILE NBSSET1
C           DEVIENT LA SUIVANTE DE NF1
            NFETOI( 5, NF1 ) = N1SSET( NBSSET1 )
         ENDIF

C        NF1 EST LA PREMIERE FACE DE LA SOUS ETOILE NBSSET1
         N1SSET( NBSSET1 ) = NF1

C        LA FACE NF1 EST MARQUEE COMME ETANT TRAITEE
         NFETOI(1,NF1) = -ABS( NFETOI(1,NF1) )

C        RETOUR EN HAUT DE PILE
         GOTO 20

      ENDIF


C     TOUTES LES FACES DE LA SOUS-ETOILE NBSSET1 ONT ETE TRAITEES
ccc      print *,'sousetor: N1FEOC=',N1FEOC,
ccc     %        ' NBSSET00=',NBSSET00,'  NBSSET1=',NBSSET1,
ccc     %        ' N1SSET=',(N1SSET(K),K=1,NBSSET1)
ccc      CALL TRFETO7R( XYZSOM, N1FEOC, NFETOI, NBSSET1, N1SSET,
ccc     %               NBAR3F, NSAR3F )


C     RESTE-T-IL DES FACES DANS LA SOUS-ETOILE INITIALE ISSUE DE N1FEOC?
C     ------------------------------------------------------------------
      IF( N1FEOC .GT. 0 ) THEN

C        OUI: LES FACES DE L'ETOILE SONT DEMARQUEES
         NBFETO = 0
         NF1    = N1FEOC
 70      IF( NF1 .GT. 0 ) THEN
            NBFETO = NBFETO + 1
            NFETOI(1,NF1) = ABS( NFETOI(1,NF1) )
            NF1 = NFETOI(5,NF1)
            GOTO 70
         ENDIF

C        IL RESTE NBFETO FACES NON TRAITEES DANS L'ETOILE
ccc         print*,'sousetor: sous-etoile',NBSSET1,
ccc     %          ' il reste',nbfeto,' faces a traiter'

C        LE TEMOIN DE TRAITEMENT SUR NFETOI(1,*) EST RENDU POSITIF
         DO M = NBSSET0+1, NBSSET1
            NF = N1SSET( M )
 80         IF( NF .GT. 0 ) THEN
               NFETOI(1,NF) = ABS( NFETOI(1,NF) )
               NF = NFETOI(5,NF)
               GOTO 80
            ENDIF
         ENDDO

C        RECHERCHE DE LA SOUS-ETOILE SUIVANTE
         IF( NBFETO0 .EQ. NBFETO ) THEN
C           ETAT STATIONNAIRE => BOUCLE INFINIE
         print*,'PB sousetor:',NBFETO,' NOMBRE INCHANGE de FACES N1FEOC'
            GOTO 9000
         ENDIF

C        RETOUR AU TRAITEMENT DES FACES RESTANTES DE LA SOUS-ETOILE
         GOTO 5

      ENDIF

C     ANOMALIE?
 9000 IF( N1FEOC .GT. 0 .OR. NBSSET1 .LE. NBSSET00 ) THEN
         PRINT *,'sousetor: Probleme N1FEOC0=',N1FEOC0,
     %           ' NBSSET0=',NBSSET00,
     %           ' N1FEOC1=',N1FEOC,'  NBSSET1=',NBSSET1
         IERR = 2
         RETURN
      ENDIF

C     TRACE DE TOUTES LES FACES DES SOUS-ETOILES QUI ONT ETE TRAITEES
C     ---------------------------------------------------------------
      IF( NBSSET1 .NE. NBSSET00 ) THEN
         print *,'sousetor: N1FEOC=',N1FEOC,
     %           ' NbSousEtoile00=',NBSSET00,'  NbSousEtoile1=',NBSSET1,
     %           ' N1SSET:',(N1SSET(K),K=1,NBSSET1)
         CALL TRFETO7R( XYZSOM, N1FEOC, NFETOI, NBSSET1, N1SSET,
     %                  NBAR3F, NSAR3F )
      ENDIF

 9999 RETURN
      END
