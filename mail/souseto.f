      SUBROUTINE SOUSETO( PTXYZD, NOTETR, N1FEOC00, NFETOI,
     %                    MXSSET, NBSSET1, N1SSET, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :  DETERMINER LE NOMBRE ET LES FACES DES SOUS-ETOILES DE L'ETOILE
C ------  ISSUE DES FACES SIMPLES D'UNE ETOILE DEFINIE PAR SON POINTEUR
C         N1FEOC00 SUR LE TABLEAU NFETOI DE L'ETOILE
C         ET LES AJOUTER AUX NBSSET1 SOUS-ETOILES INITIALES

C         UNE SOUS-ETOILE EST TELLE QUE TOUTE ARETE DE SES FACES
C         . APPARTIENT A 2 ET SEULEMENT 2 DE SES FACES
C         . EST PARCOURUE DANS DES SENS CONTRAIRES DANS CES 2 FACES
C         . LE VECTEUR NORMAL A CHAQUE FACE EST DIRIGEE VERS L'INTERIEUR
C           DE LA SOUS ETOILE

C          VERSION PTXYZD DES SOMMETS EN DOUBLE PRECISION

C ENTREES:
C --------
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C N1FEOC00:POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 2 LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: >0 NFETOI(2,NF)=NUMERO PTXYZD DU SOMMET 1 DE LA FACE NF
C             <0 LA FACE NF EST MARQUEE COMME ETANT TRAITEE
C                NFETOI(2,NF)=-ABS( NFETOI(2,NF) )
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2 x S1S3 VERS L'INTERIEUR DE L'ETOILE
C          5: NUMERO DANS NFETOI DE LA FACE SUIVANTE
C MXSSET : NOMBRE MAXIMAL DE SOUS ETOILES DECLARABLES DANS N1SSET

C MODIFIE :
C ---------
C NBSSET1: NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C          NOMBRE FINAL   DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI

C SORTIES:
C --------
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET1 SOUS ETOILES
C IERR   : 1 SI LE TABLEAU LAPILE EST SATURE
C          2 SI N1FEOC=0 OU NBSSET1=0 EN SORTIE
C          3 SI L'ORIENTATION VERS L'INTERIEUR DE LA NORMALE A
C            TOUTES LES FACES D'UNE SOUS-ETOILE N'EST PAS RESPECTEE
C            I.E. AU MOINS LA NORMALE D'UNE FACE EST DIRIGEE
C                 VERS L'EXTERIEUR
C          4 SI LES ARETES DES FACES SIMPLES NE SONT PAS PARCOURUUES
C            AVEC UN MEME NOMBRE DE FOIS DANS LES 2 SENS
C          0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC& St PIERRE du PERRAY Septembre 2015
C MODIFS : ALAIN PERRONNET LJLL UPMC& LJUBLJANA SLOVENIE  Decembre  2015
C MODIFS : ALAIN PERRONNET LJLL UPMC& St PIERRE du PERRAY Octobre   2017
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray         Decembre  2019
C2345X7..............................................................012
      include"./incl/trvari.inc"
      PARAMETER        (MXPILE=2048)
      INTEGER           LAPILE(MXPILE)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      INTEGER           NOTETR(8,*), NFETOI(5,*), N1SSET(MXSSET),
     %                  NSARP2F(3,MXPILE),
     %                  NOSOTR(3), NBFASI, NOFASI(8),
     %                  NBFASM, NOFASM(8)
      DOUBLE PRECISION  PTXYZD(4,*), ANGLE, ANGMIN
C     NUMERO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
C     LA NORMALE A LA FACE EST DIRIGEE VERS L'EXTERIEUR DU TETRAEDRE
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      TRACTE0 = TRACTE
      LORBITE = 1
      IERR    = 0

C     INITIALISATION DES SOUS ETOILES ISSUES DE N1FEOC00
      NBSSET00 = NBSSET1
      N1FEOC   = N1FEOC00
      NBARP2F  = 0
      NBDEFO   = 0

C     TRACE DES FACES DE LA SOUS-ETOILE ISSUE DE N1FEOC00
C     ET LES AUTRES SOUS-ETOILES
      CALL TRFETO7( PTXYZD,  N1FEOC00, NFETOI, NBSSET1, N1SSET,
     %              NBARP2F, NSARP2F )


C     DEPART DES FACES DE L'ETOILE A PARTIR DU POINTEUR N1FEOC
C     --------------------------------------------------------
 5    IF( N1FEOC .LE. 0 ) GOTO 9000

C     NOMBRE DE SOUS-ETOILES INITIALES DANS N1SSET-NFETOI
C     AVANT TRAITEMENT DES FACES ISSUES DE N1FEOC
      NBSSET0 = NBSSET1

C     MARQUAGE EVENTUEL SUPPRIME DES FACES DE LA SOUS-ETOILE
C     ISSUE DE N1FEOC ET NOMBRE DE SES FACES
      NBFETO = 0
      NF     = N1FEOC
 10   IF( NF .GT. 0 ) THEN
         NBFETO = NBFETO + 1
         NFETOI(2,NF) = ABS( NFETOI(2,NF) )
         NF = NFETOI(5,NF)
         GOTO 10
      ENDIF
      NBFETO0 = NBFETO

ccc      print*,'souseto: NBSSET0,' Sous-Etoiles INITIALES +',
ccc     %       NBFETO0,' FACES issues de N1FEOC=',N1FEOC

C     MARQUAGE EVENTUEL SUPPRIME DES FACES DES NBSSET00 SOUS-ETOILES
      DO M = 1, NBSSET00
         NF = N1SSET( M )
 15      IF( NF .GT. 0 ) THEN
            NFETOI(2,NF) = ABS( NFETOI(2,NF) )
            NF = NFETOI(5,NF)
            GOTO 15
         ENDIF
      ENDDO

C     INITIALISATION DE LA PILE DES FACES ISSUES DE N1FEOC PAR LEURS ARETES
C     =====================================================================
      LAPILE(1) = N1FEOC
      LHPILE = 1

 20   IF( LHPILE .GT. 0 ) THEN

C        TRAITEMENT DE LA FACE NF1 EN HAUT DE LAPILE C-A-D
C        RECHERCHE DE SA SOUS-ETOILE ET DE SON CHAINAGE
C        =================================================
         NF1 = LAPILE( LHPILE )

C        LA FACE EST DEPILEE
         LHPILE = LHPILE - 1

C        CETTE FACE EST ELLE DEJA TRAITEE?
         IF( NFETOI(2,NF1) .LE. 0 ) GOTO 20

ccc         PRINT*,'souseto: TRAITEMENT NF1=',NF1,' :',
ccc     %          (NFETOI(KK,NF1),KK=1,5)

C        NS1-NS2 SOMMETS DE L'ARETE K-1 DE LA FACE NF1
         NS1 = NFETOI(4,NF1)
         NS2 = NFETOI(2,NF1)
         DO 30 K=1,3

C           ARETE K-1 DE LA FACE NF1 de SOMMETS NS1-NS2
C           NS2 = NFETOI(1+K,NF1)

C           LE TROISIEME SOMMET DE NF1
            IF( K .LT. 3 ) THEN
               KK = K+1
            ELSE
               KK = 1
            ENDIF
            NS3 = NFETOI(1+KK,NF1)

C           RECENSEMENT DANS L'ETOILE DES FACES D'ARETE NS1-NS2
C           PARCOURUE EN SENS INVERSE NOFASI ou LE MEME SENS NOFASM
C           -------------------------------------------------------
            NBFASI = 0
            NBFASM = 0

C           DEBUT DE RECHERCHE PARMI LES FACES CHAINEES PAR N1FEOC
            NEWSSET = NBSSET0
            NF2     = N1FEOC

 22         IF( NF2 .GT. 0 ) THEN

C              PARCOURS DES 3 ARETES DE NF2
               NSA1 = NFETOI(4,NF2)
               DO KK=1,3
                  NSA2 = ABS( NFETOI(1+KK,NF2) )
                  IF( NS1.EQ.NSA2 .AND. NS2.EQ.NSA1 ) THEN
C                    ARETE NS1-NS2 VUE DANS LE SENS OPPOSE
C                    NOMBRE DE FACES D'ARETE K PARCOURUE EN SENS INVERSE
                     NBFASI = NBFASI + 1
                     NOFASI( NBFASI ) = NF2
                  ELSE IF( NS1.EQ.NSA1 .AND. NS2.EQ.NSA2 ) THEN
C                    ARETE NS1-NS2 VUE DANS DE LE MEME SENS
C                    NOMBRE DE FACES D'ARETE K PARCOURUE DANS LE MEME SENS
                     NBFASM = NBFASM + 1
                     NOFASM( NBFASM ) = NF2
                  ENDIF
                  NSA1 = NSA2
               ENDDO

               NF2 = NFETOI( 5, NF2 )
               GOTO 22

            ENDIF

C           RECHERCHE PARMI LES SOUS-ETOILES NOUVELLEMENT CONSTRUITES
C           A PARTIR DE N1FEOC00
            NEWSSET = NEWSSET + 1
            IF( NEWSSET .LE. NBSSET1 ) THEN
C              DANS LA SOUS-ETOILE NEWSSET NOUVELLEMENT CONSTRUITE
C              RECHERCHE DE LA FACE D'ARETE NS1-NS2 PARCOURUE EN SENS INVERSE
               NF2 = N1SSET( NEWSSET )
               GOTO 22
            ENDIF

C           TRAITEMENT SELON LE NOMBRE DE FACES D'ARETE NS1-NS2 EN SENS INVERSE
C           -------------------------------------------------------------------
C           NOMBRE DE FACES D'ARETE K PARCOURUE DANS LE SENS INVERSE
            IF( NBFASI .LE. 0 ) THEN
C              NOMBRE INCORRECT DE FACES D'ARETE NS1-NS2 PARCOURUE
C              DANS LE SENS INVERSE
               PRINT*,'souseto: FACE',NF1,' ARETE',NS1,NS2,
     %                ' avec 0 FACE ADJACENTE PARCOURUE EN SENS INVERSE'
               PRINT*,'souseto: NBSSET0=',NBSSET0,' NBSSET1=',NBSSET1,
     %                ' NEWSSET=',NEWSSET
               GOTO 9990
            ENDIF

            IF( NBFASI .NE. NBFASM ) THEN
C              PROBLEME: L'ARETE NS1-NS2 N'EST PAS PARCOURUE
C                        AUTANT DE FOIS DANS LES 2 SENS
               PRINT*,'souseto: PROBLEME: ARETE',NS1,NS2,' NON PARCOURUE
     % AUTANT DE FOIS DANS LES 2 SENS NBFASI=',NBFASI,' NBFASM=',NBFASM
               PRINT*,'souseto: NBSSET0=',NBSSET0,' NBSSET1=',NBSSET1,
     %                ' NEWSSET=',NEWSSET
               GOTO 9990
            ENDIF

C           NUMERO DE LA PREMIERE FACE ADJACENTE EN SENS INVERSE DE NS1-NS2
            NF = NOFASI(1)
            IF( NBFASI .EQ. 1 ) THEN

C              L'UNIQUE FACE NF D'ARETE PARCOURUE EN SENS INVERSE
C              DE NS1-NS2 EST EMPILEE SI ELLE N'EST PAS DEJA TRAITEE
C              -----------------------------------------------------
               IF( NFETOI(2,NF) .GT. 0 ) THEN
C                 NF EST EMPILEE CAR NON TRAITEE
                  IF( LHPILE .LT. MXPILE ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE( LHPILE ) = NF
                  ELSE
C                    TABLEAU LAPILE SATURE
                     PRINT *,'souseto: TABLEAU LAPILE SATURE'
                     PRINT *,'souseto: AUGMENTER MXPILE=',MXPILE
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

            ELSE

C              IL EXISTE PLUSIEURS FACES ADJACENTES A NS1-NS2
C              DANS LE SENS INVERSE =>
C              IL FAUT EMPILER LA FACE NFAMIN D'ANGLE MINIMAL AVEC NF1
C              -------------------------------------------------------

ccc               PRINT*
ccc               PRINT*,'souseto:',NS1,NS2,' arete de la face',NF1,
ccc     %                ' avec NBFASM=',NBFASM,' NBFASI=',NBFASI
ccc               DO KK=1,NBFASM
ccc                  M = NOFASM(KK)
ccc                PRINT*,'FACE MEME SENS   ',M,':',(NFETOI(KKK,M),KKK=1,5)
ccc               ENDDO
ccc               DO KK=1,NBFASI
ccc                  M = NOFASI(KK)
ccc                PRINT*,'FACE SENS INVERSE',M,':',(NFETOI(KKK,M),KKK=1,5)
ccc               ENDDO

C              NS1-NS2 EST UNE ARETE DANS PLUS DE 2 FACES: STOCKAGE OU NON
C              LE TABLEAU NSARP2F SERT UNIQUEMENT A VISUALISER CES ARETES
               DO L=1,NBARP2F
                  NSA1 = NSARP2F(1,L)
                  NSA2 = NSARP2F(2,L)
                  IF( NS1 .EQ. NSA1 .AND. NS2 .EQ. NSA2  .OR.
     %                NS1 .EQ. NSA2 .AND. NS2 .EQ. NSA1 ) THEN
C                    L'ARETE A DEJA ETE STOCKEE
                     GOTO 23
                  ENDIF
               ENDDO
               IF( NBARP2F .GE. MXPILE ) THEN
                  PRINT*,'souseto: AUGMENTER MXPILE car NSARP2F SATURE'
                  GOTO 23
               ENDIF
               NBARP2F = NBARP2F + 1
               NSARP2F(1,NBARP2F) = NS1
               NSARP2F(2,NBARP2F) = NS2
               NSARP2F(3,NBARP2F) = NBFASI - NBFASM
               NBDEFO = NBDEFO + ABS( NBFASI - NBFASM )

ccc               PRINT*
ccc               PRINT*,'souseto: ARETE',NBARP2F,' St:',NS1,NS2,
ccc     %               ' DANS NBFASensI=',NBFASI,' + NBFASensM=',NBFASM,
ccc     %               ' FACES'
ccc               NTEOP = NFETOI(1,NF1)
ccc               PRINT*,'souseto: NFETOI(',NF1,'): Tetraedre oppose=',
ccc     %                 NTEOP,' 3 St=',(NFETOI(M,NF1),M=2,4),
ccc     %               ' No NFETOI Triangle SUIVANT=',NFETOI(5,NF1)
ccc               IF( NTEOP .GT. 0 ) THEN
ccc                  PRINT*,'souseto: NOTETR(',NTEOP,')=',
ccc     %                   (NOTETR(M,NTEOP),M=1,8)
ccc               ENDIF

C              TRACE DES FACES DES NBSSET1 SOUS-ETOILES AVEC ARETE
C              APPARTENANT A 4 ou 6 ou 8 ... FACES
ccc               tracte = .true.
               CALL TRFETO7( PTXYZD,  N1FEOC, NFETOI, NBSSET1, N1SSET,
     %                       NBARP2F, NSARP2F )

C              LES FACES NF DE CETTE ARETE NS1-NS2 PARCOURUE EN SENS INVERSE
 23            ANGMIN = 16D0
               NFAMIN = 0
               DO 25 L=1,NBFASI

C                 LA FACE NF D'ARETE NS1 NS2 PARCOURUE EN SENS INVERSE
                  NF = NOFASI( L )

C                 RECHERCHE DU 3-EME SOMMET DE NF NON NS1 et NON NS2
                  DO M=1,3
                     NS4 = ABS( NFETOI(1+M,NF) )
                     IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 24
                  ENDDO

C                 CALCUL DE L'ANGLE DIEDRE ENTRE LES PLANS DES 2 FACES
C                 NF1 et NF D'ARETE COMMUNE NS1-NS2  I.E.
C                 ANGLE ENTRE LES PLANS DES 2 TRIANGLES S1S2S3 et S2S1S4
C                 DANS L'INTERVALLE [0, 2Pi] RADIANS
 24               CALL ANG2TR3D( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                           PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                           ANGLE, IERR )
                  IF( IERR .NE. 0 ) THEN
C                    AU MOINS UN DES 2 TRIANGLES EST DEGENERE
                     PRINT*,'souseto: FACE NF1=',NF1,' St',NS1,NS2,NS3,
     %                      ' DEGENEREE'
                     IERR = 0
                     GOTO 25
                  ENDIF

                  IF( ANGLE .LT. ANGMIN ) THEN
                     ANGMIN = ANGLE
                     NFAMIN = NF
                  ENDIF

 25            ENDDO

               IF( NFETOI(2,NFAMIN) .GT. 0 ) THEN

C                 NFAMIN EST EMPILEE CAR NON DEJA TRAITEE
                  IF( LHPILE .LT. MXPILE ) THEN
C                    LA FACE NFAMIN EST EMPILEE DANS LAPILE
                     LHPILE = LHPILE + 1
                     LAPILE( LHPILE ) = NFAMIN
                  ELSE
C                    TABLEAU LAPILE SATURE
                     PRINT*,'souseto: tableau LAPILE SATURE'
                     PRINT*,'AUGMENTER MXPILE=',MXPILE
                     IERR = 1
                     GOTO 9999
                  ENDIF

               ENDIF

            ENDIF

C           ARETE K SUIVANTE DE LA FACE NF1
            NS1 = NS2
            NS2 = NS3

 30      ENDDO


C        AJOUT DE LA FACE NF1 A LA NOUVELLE SOUS-ETOILE NBSSET1
C        ------------------------------------------------------
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
C           NF1 N'ETAIT PAS LA PREMIERE FACE DE L'ETOILE
            NFETOI( 5, NF0 ) = NF2
         ELSE
C           NF1 ETAIT LA PREMIERE FACE DE L'ETOILE
            N1FEOC = NF2
         ENDIF

C        LA FACE NF1 EST AJOUTEE AU DEBUT DE LA SOUS ETOILE NBSSET1
         IF( NBSSET1 .EQ. NBSSET0 ) THEN
C           NF1 DEBUTE LA SOUS-ETOILE NBSSET0+1
            NBSSET1 = NBSSET1 + 1
C           LA SOUS-ETOILE NBSSET1 EST COMPOSEE DE LA SEULE FACE NF1
            NFETOI( 5, NF1 ) = 0
         ELSE
C           LA PREMIERE FACE DE LA SOUS ETOILE NBSSET1
C           DEVIENT LA SUIVANTE DE NF1
            NFETOI( 5, NF1 ) = N1SSET( NBSSET1 )
         ENDIF

C        NF1 EST LA PREMIERE FACE DE LA SOUS ETOILE NBSSET1
         N1SSET( NBSSET1 ) = NF1

C        LA FACE NF1 EST MARQUEE COMME ETANT TRAITEE
         NFETOI(2,NF1) = -ABS( NFETOI(2,NF1) )

C        RETOUR EN HAUT DE PILE
         GOTO 20

      ENDIF


C     TOUTES LES FACES DE LA SOUS-ETOILE NBSSET0+1 A NBSSET1 ONT ETE TRAITEES
C     -----------------------------------------------------------------------
      IF( NBARP2F .GT. 0 .OR. NBDEFO .NE. 0 ) THEN
C        POUR VOIR LES ARETES DE PLUS DE 2 FACES OU LES ARETES N'ETANT PAS
C        PARCOURUES UN MEME NOMBRE DE FOIS DANS LES 2 SENS

ccc         tracte = .true.
ccc         print *
ccc         print *,'souseto: N1FEOC=',N1FEOC,
ccc     %           ' NBSSET00=',NBSSET00,'  NBSSET1=',NBSSET1,
ccc     %           ' N1SSET:',(N1SSET(K),K=1,NBSSET1),' NBARP2F=',NBARP2F
ccc         CALL TRFETO7( PTXYZD, N1FEOC, NFETOI,
ccc     %                 NBSSET1-NBSSET00, N1SSET(NBSSET00+1),
ccc     %                 NBARP2F, NSARP2F )

C        TRACE DES NBSSET1 SOUS-ETOILES
         CALL TRFETO7( PTXYZD,  N1FEOC, NFETOI, NBSSET1, N1SSET,
     %                 NBARP2F, NSARP2F )
      ENDIF


C     RESTE-T-IL DES FACES DANS LA SOUS-ETOILE INITIALE ISSUE DE N1FEOC?
C     ------------------------------------------------------------------
      IF( N1FEOC .GT. 0 ) THEN

C        OUI: LES FACES RESTANTES DE L'ETOILE SONT DEMARQUEES
         NBFETO = 0
         NF1    = N1FEOC
 70      IF( NF1 .GT. 0 ) THEN
            NBFETO = NBFETO + 1
            NST = NFETOI(2,NF1)
            IF( NST .LT. 0 ) THEN
C              LE MARQUAGE DE FACE TRAITEE DE NF1 EST SUPPRIME
               NFETOI(2,NF1) = -NST
            ENDIF
            NF1 = NFETOI(5,NF1)
            GOTO 70
         ENDIF

C        IL RESTE NBFETO FACES NON TRAITEES DANS L'ETOILE
ccc         print *,'souseto: sous-etoile',NBSSET1,
ccc     %  ' il reste',NBFETO,' faces a traiter a partir de N1FEOC=',N1FEOC

C        LE TEMOIN DE TRAITEMENT SUR NFETOI(2,*) EST RENDU POSITIF
         DO M = NBSSET0+1, NBSSET1
            NF = N1SSET( M )
 80         IF( NF .GT. 0 ) THEN
               NFETOI(2,NF) = ABS( NFETOI(2,NF) )
               NF = NFETOI(5,NF)
               GOTO 80
            ENDIF
         ENDDO

C        RECHERCHE DE LA SOUS-ETOILE SUIVANTE
         IF( NBFETO .EQ. NBFETO0 ) THEN
C           ETAT STATIONNAIRE => BOUCLE INFINIE
          print*,'PB souseto:',NBFETO,' NOMBRE INCHANGE de FACES N1FEOC'
            GOTO 9000
         ENDIF

C        RETOUR AU TRAITEMENT DES FACES RESTANTES DE LA SOUS-ETOILE
C        A PARTIR DE N1FEOC
C        ----------------------------------------------------------
         GOTO 5

      ENDIF


C     FIN DE CONSTRUCTION DE LA SOUS-ETOILE DES FACES ISSUES DE N1FEOC0
 9000 IF( N1FEOC .GT. 0 .OR. NBSSET1 .LE. NBSSET00 ) THEN
C        ANOMALIE?
         PRINT *,'souseto: Probleme N1FEOC00=',N1FEOC00,
     %                               ' NBSSET00=',NBSSET00,
     %           ' N1FEOC0=',N1FEOC0,' NBSSET0=',NBSSET0,
     %           ' N1FEOC1=',N1FEOC, ' NBSSET1=',NBSSET1
         IERR = 2
         GOTO 9999
      ENDIF


C     TRACE DES SOUS-ETOILES 1 A NBSSET1 ET DES ARETES DANS PLUS DE 2 FACES
C     ---------------------------------------------------------------------
      IF( NBDEFO .NE. 0  ) THEN

         PRINT*,'Les',NBSSET1,' SOUS ETOILES: %%%%%%%%%%%%%%%%%%%%%%%%%'
         PRINT *,'souseto: Sortie avec',NBARP2F,
     %           ' ARETES DANS PLUS DE 2 FACES. N1FEOC=',N1FEOC,
     %           ' NbSousEtoile0=',NBSSET0,
     %           ' NbSousEtoile1=',NBSSET1
         PRINT *,'souseto: N1SSET:',(N1SSET(K),K=1,NBSSET1)
         DO M = 1, NBSSET1
            NF = N1SSET( M )
            CALL AFETOI( NF, NFETOI )
            PRINT *
         ENDDO
         TRACTE = TRACTE0
         TRACTE = .TRUE.

      ENDIF

C     TRACE DES FACES DES NBSSET1 SOUS-ETOILES AVANT VERIFICATION
      CALL TRFETO7( PTXYZD,  N1FEOC, NFETOI, NBSSET1, N1SSET,
     %              NBARP2F, NSARP2F )
      TRACTE = TRACTE0


C     VERIFICATION DE L'ORIENTATION DES FACES DES NBSSET1 SOUS-ETOILES:
C     L'ORIENTATION DE LA FACE EST ELLE OPPOSEE A CELLE DE LA FACE
C     DU TETRAEDRE OPPOSE
C     -----------------------------------------------------------------
      DO M = 1, NBSSET1

         NBREOR = 0

C        LES FACES DE LA SOUS-ETOILE M
         NF = N1SSET( M )

 9010    IF( NF .GT. 0 ) THEN

            NST1 = ABS( NFETOI(2,NF) )
C           LE TEMOIN DE TRAITEMENT SUR NFETOI(2,*) EST RENDU POSITIF
            NFETOI(2,NF) = NST1

C           RECHERCHE DE LA FACE NF DANS LE TETRAEDRE OPPOSE
            NTEOP = NFETOI(1,NF)
            IF( NTEOP .GT. 0 ) THEN
               NOSOTR(1) = NST1
               NOSOTR(2) = NFETOI(3,NF)
               NOSOTR(3) = NFETOI(4,NF)
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTEOP), MF )
               IF( MF .GT. 0 ) THEN

C                 LA NORMALE DE LA FACE MF DU TETRAEDRE NTEOP
C                 EST DIRIGEE VERS L'EXTERIEUR DU TETRAEDRE
C                 DONC VERS L'INTERIEUR DE LA SOUS-ETOILE
                  DO K=1,3
                     NOSOTR(K) = NOTETR( NOSOFATE(K,MF), NTEOP )
                  ENDDO
C                 QUEL EST LE SOMMET 1 NST1 DE LA FACE NF NFETOI
C                 DANS LA FACE CORRECTEMENT ORIENTEE NOSOTR?
                  IF( NST1 .EQ. NOSOTR(2) ) THEN
                     NOSOTR(2) = NOSOTR(3)
                     NOSOTR(3) = NOSOTR(1)
                     NOSOTR(1) = NST1
                  ELSE IF( NST1 .EQ. NOSOTR(3) ) THEN
                     NOSOTR(3) = NOSOTR(2)
                     NOSOTR(2) = NOSOTR(1)
                     NOSOTR(1) = NST1
                  ENDIF

C                 L'ORIENTION EST BONNE SI 2 A 2 LES SOMMETS 2 ET 3
C                 SONT IDENTIQUES
                  IF( NOSOTR(2) .NE. NFETOI(3,NF) .OR.
     %                NOSOTR(3) .NE. NFETOI(4,NF) ) THEN

C                    LES ORIENTATIONS SONT OPPOSEES
C                    => IL FAUT RENVERSER L'ORIENTATION DE TOUTES
C                       LES FACES DE LA SOUS-ETOILE M
C                       SI CELA N'A PAS ETE FAIT AUPARAVANT
C                    --------------------------------------------
                     IF( NBREOR .EQ. 0 ) THEN
C                       PREMIERE RE-ORIENTATION PAR PERMUTATION DES
C                       SOMMETS 2 ET 3 DES FACES NFETOI DE LA SOUS-ETOILE M
                        NBREOR = 1
                        NF1 = N1SSET( M )
 9020                   IF( NF1 .GT. 0 ) THEN
                           K                = NFETOI( 3, NF1 )
                           NFETOI( 3, NF1 ) = NFETOI( 4, NF1 )
                           NFETOI( 4, NF1 ) = K
C                          PASSAGE A LA FACE SUIVANTE DE LA SOUS-ETOILE M
                           NF1 = NFETOI( 5, NF1 )
                           GOTO 9020
                        ENDIF
                        PRINT*,'souseto: 1-ERE RE-ORIENTATION des FACES'
C                       TRACE DES FACES DES NBSSET1 SOUS-ETOILES
ccc                        tracte = .true.
                        CALL TRFETO7( PTXYZD,  N1FEOC, NFETOI,
     %                                NBSSET1, N1SSET, NBARP2F, NSARP2F)
                     ELSE
C                       DEMANDE DE SECONDE RE-ORIENTATION:
C                       PROBLEME: LES NORMALES DE TOUTES LES FACES
C                       DE LA SOUS-ETOILE M DOIVENT ETRE DIRIGEES
C                       VERS L'INTERIEUR DE LA SOUS-ETOILE CE QUI
C                       N'EST PAS LE CAS => PROBLEME A RESOUDRE
C                       -----------------------------------------
                        PRINT*,'souseto: PROBLEME toutes les FACES de la
     % SOUS-ETOILE',M,' NE SONT PAS DIRIGEES VERS L''INTERIEUR'
                        IERR = 3
                        GOTO 9900
                     ENDIF

                  ENDIF
               ENDIF
            ENDIF

C           PASSAGE A LA FACE SUIVANTE DE LA SOUS-ETOILE M
            NF = NFETOI(5,NF)
            GOTO 9010

         ENDIF
      ENDDO

C     FIN: TRACE DES FACES DES NBSSET1 SOUS-ETOILES APRES VERIFICATION
 9900 CALL TRFETO7( PTXYZD,  N1FEOC, NFETOI, NBSSET1, N1SSET,
     %              NBARP2F, NSARP2F )
      TRACTE = TRACTE0
      GOTO 9999

C     AU MOINS UNE ARETE D'UNE FACE SIMPLE N'EST PAS PARCOURUE LE MEME
C     NOMBRE DE FOIS DANS LES 2 SENS
 9990 IERR = 4


 9999 TRACTE = TRACTE0
      RETURN
      END
