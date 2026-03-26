      SUBROUTINE TETRCYCL( QUVEVM, KTITRE, NOFOTI,
     %                     MXSOMM, NBSOMM, PTXYZD, NPSOFR,
     %                     MXFETO, N1FEVI, N1FEOC, NFETOI,
     %                     XYZ,    MXVPSI, NBVPSI,  NFVPSI0, NFVPSI,
     %                     MXCIAS, N1CIAS, MXASFVP, NSASFVP,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  EXISTE T IL DES CYCLES d'ARETES SIMPLES, BORD DES NBVPSI>0 FACES
C -----  TRIANGULABLES POUR REDUIRE LE NOMBRE DE FACES RENTRANTES
C        GENEREES PAR LES TETRAEDRES FACE VPSI + XYZ
C        ET TRIANGULATION DU CYCLE POUR FERMER LA TETRAEDRISATION DES
C        FACES NFVPSI

C ENTREES:
C --------
C QUVEVM : QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C          C-A-D AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE
C KTITRE : TITRE D'UN TRACE
C NOFOTI : =1       EXISTENCE DE LA FONCTION UTILISATEUR 'TAILLE_IDEALE'
C          =0 PAS D'EXISTENCE DE LA FONCTION UTILISATEUR 'TAILLE_IDEALE'
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS PTXYZD NPSOFR
C MXFETO : NOMBRE DE FACES DECLARABLES DANS NFETOI
C XYZ    : COORDONNEES DU POINT QUI AVEC LES FACES NFVPSI DONNE DES
C          TETRAEDRES a VOLUME>0 ET SANS INTERSECTION
C MXVPSI : NOMBRE MAXIMAL DE FACES VPSI DANS LE TABLEAU NFVPSI
C MXCIAS : NOMBRE MAXIMAL DE CYCLES DEFINIS DANS N1CIAS

C MODIFIES:
C ---------
C NBSOMM : NUMERO DU DERNIER SOMMET CREE DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE

C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C N1FEOC :  POINTEUR SUR LA PREMIERE FACE NFETOI ACTIVE DE L'ETOILE
C NFETOI : VERSION 2 LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C NFVPSI0: NUMERO NFETOI DES FACES INITIALES NFETOI + XYZ -> TETRAEDRES a V>0
C          <0 INDIQUE UNE CREATION LORS D'UNE TRIANGULATION D'UN CYCLE


C MODIFIES:
C ---------
C NBVPSI : NOMBRE DE FACES FINALES NFETOI + XYZ -> TETRAEDRES a V>0
C NFVPSI : NUMERO NFETOI DES FACES NFETOI + XYZ -> TETRAEDRES a V>0
C          <0 INDIQUE UNE CREATION LORS D'UNE TRIANGULATION D'UN CYCLE
C N1CIAS : NUMERO DANS NFVPSI DE LA PREMIERE FACE D'UN CYCLE

C SORTIES:
C --------
C NSASFVP: =1 NUMERO PTXYZD DU 1-ER  SOMMET DE L'ARETE
C          =2 NUMERO PTXYZD DU 2-EME SOMMET DE L'ARETE
C          =3 ARETE SIMPLE SUIVANTE DANS UN CYCLE
C             -1 NON ACTIVE, 0 FIN DE CYCLE
C          =4 NUMERO DE SA FACE DANS LE TABLEAU NFSVPI
C IERR   : =0 SI PAS D'ERREUR
C          =1 ETOILE AVEC MOINS DE 4 FACES
C          =2 SATURATION D'UN TABLEAU
C          =3 SOMMETS NON DANS UN TETRAEDRE
C          =4 NOMBRE INCORRECT DE FACES DANS L'ETOILE
C          =5 SOMMETS DES 4 FACES DE L'ETOILE INCORRECTS
C          =6 UNE ARETE DES FACES DE L'ETOILE APPARTIENT A MOINS DE 2
C             OU PLUS DE 2 FACES DE L'ETOILE
C          =7 ou 8 ou 9  ALGORITHME DEFAILLANT A AMELIORER...
C          =10 PLUSIEURS ARETES A PROBLEME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER       Aout 2015
C MODIFS : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER  Septembre 2017
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2020
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      DOUBLE PRECISION  PTXYZD(1:4,1:NBSOMM)
      INTEGER           NPSOFR(MXSOMM),  NFETOI(5,MXFETO)
      INTEGER           NFVPSI0(MXVPSI), NFVPSI(MXVPSI),
     %                  N1CIAS(MXCIAS),  NSASFVP(4,MXASFVP)

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  D, V, XYZ(4), XYZBAC(4),XYZBAV(3),XYZESS(3),
     %                  XYZMIMX(3), ARMIN, ARMAX, SURFTR(4)
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      TRACTE0 = TRACTE
      LORBITE = 1
ccc      tracte = .true.

C     pour forcer la construction d'une triangulation d'un cycle
      QUVEVM = 1E-5

C     PROTECTION DU TABLEAU NFVPSI DANS NFVPSI0 POUR REPRISE EVENTUELLE
      NBSOMM0 = NBSOMM
      NBVPSI0 = NBVPSI
      DO N=1,NBVPSI0
         NFVPSI0( N ) = NFVPSI( N )
      ENDDO

      NBASFVP = 0
      NBSUFV  = 0

C     TRACE EN ORANGE DES FACES+XYZ DONNANT UN VOLUME POSITIF
C     SANS INTERSECTION AVEC LES FACES SIMPLES NFETOI DE L'ETOILE
C     ET MAXIMISANT LE MINIMUM DES QUALITES DES TETRAEDRES FORMES
c     trace pour mettre au point les cycles
 10   KTITRE='tetrcycl: Recherche de CYCLES ARETES SIMPLES des faces FVP
     %SI de XYZ'
      CALL TRFETO3( KTITRE, XYZ, XYZ, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)


C     RECHERCHE DES CYCLES D'ARETES SIMPLES DES FACES NFVPSI V>0
C     ------------------------------------------------------------
C     GRANDE VALEUR POUR MISE AU POINT DE LA DETECTION DES CYCLES 
      MIARSICI = 16
      CALL CIASFAET( NBVPSI,  NFVPSI,  NFETOI,
     %               NBCIAS,  MXCIAS,  N1CIAS,
     %               MXASFVP, MIARSICI, NBASFVP, NSASFVP, IERR )
      IF( IERR .NE. 0 .OR. NBCIAS .LE. 0 ) GOTO 9999


C     ESSAIS DE TRIANGULER LES SURFACES LIMITEES PAR LES NBCIAS CYCLES
C     D'ARETES SIMPLES A FACE V>0
C     ----------------------------------------------------------------
      DO 80 NOC = 1, NBCIAS

C        NOMBRE D'ARETES DU CYCLE NOC
         NBAC = 0

C        CALCUL DU BARYCENTRE XYZBAC DES ARETES DU CYCLE NOC
         DO N=1,4
            XYZBAC( N ) = 0D0
         ENDDO
C        L'ARETE DE DEBUT DU CYCLE NOC
         NAC = N1CIAS( NOC )

C        LE SENS POUR ASSURER LA NORMALE VERS L'INTERIEUR
 25      NBAC = NBAC + 1
         NS1  = NSASFVP(1,NAC)
         DO N=1,4
            XYZBAC( N ) = XYZBAC( N ) + PTXYZD( N, NS1 )
         ENDDO

C        L'ARETE SUIVANTE DU CYCLE NOC
         NASUIV = NSASFVP(3,NAC)
         IF( NASUIV .GT. 0 ) THEN
            NAC = NASUIV
            GOTO 25
         ELSE
C           LA DERNIERE ARETE DU CYCLE NOC
            NAFIN = NAC
         ENDIF

C        XYZD DU BARYCENTRE DU CYCLE NOC
         DO N=1,4
            XYZBAC( N ) = XYZBAC( N ) / NBAC
         ENDDO
C        TAILLE_IDEALE DE L'ARETE EN XYZBAC
         CALL TAILIDEA( NOFOTI, XYZBAC, NCODEV, XYZBAC(4) )

         PRINT*,'tetrcycl:',NBAC,' ARETES a TRIANGULER du CYCLE',
     %           NOC,' parmi les',NBCIAS,' CYCLES'
         PRINT*,'tetrcycl: Barycentre XYZD=',XYZBAC


C        RECHERCHE D'UNE TRIANGULATION INTERNE DU CYCLE NOC A PARTIR
C        DE SES ARETES SANS AJOUT D'UN NOUVEAU POINT DIFFERENT DE XYZ
C        -------------------------------------------------------------
         NBFVAU = 0
         NBTRIA = 1

C        L'ARETE DE DEBUT DU CYCLE NOC
 27      NAC = N1CIAS( NOC )

C        LE SENS POUR ASSURER LA NORMALE VERS L'INTERIEUR
         NS1 = NSASFVP(1,NAC)
         NS2 = NSASFVP(2,NAC)

         DO M=1,NBAC-2

C           L'ARETE SUIVANTE NS2-NS3 DU CYCLE NOC
            NASUIV = NSASFVP(3,NAC)
            NS3    = NSASFVP(2,NASUIV)
            IF( NS3 .EQ. NS2 ) THEN
               NS3 = NSASFVP(1,NASUIV)
            ENDIF

C           LE TRIANGLE NS1 NS2 NS3 EST IL UNE FACE DE LA SOUS-ETOILE NFETOI?
            CALL TRDFETOI( NS1,NS2,NS3, N1FEOC, NFETOI, NF )
            IF( NF .GT. 0 ) THEN
C              OUI: ABANDON DE CETTE TRIANGULATION CAR
C                   CETTE FACE+XYZ A UN VOLUME<0
               PRINT*,'tetrcycl: TRIANGLE',NS1,NS2,NS3,
     %                ' FACE DEJA dans NFETOI +XYZ',XYZ
               IF( NBAC .EQ. 3 ) THEN
C                 PASSAGE AU CYCLE NOC SUIVANT
                  GOTO 80
               ELSE
C                 PASSAGE A UNE AUTRE METHODE DE TRIANGULATION
                  GOTO 30
               ENDIF
            ENDIF

C           LE VOLUME NS2 NS1 NS3 XYZ EST IL POSITIF?
            CALL QUATETD( PTXYZD(1,NS2), PTXYZD(1,NS1),
     %                    PTXYZD(1,NS3), XYZ,
     %                    ARMIN, ARMAX, SURFTR, V, Q )

            PRINT*,'tetrcycl: TRIANGLE',NS1,NS2,NS3,' +XYZ',XYZ,
     %             ' => Tetra avec V=',V,' Q=',Q

ccc            tracte = .true.
ccc            CALL TRFETO3( KTITRE(1:L), XYZ, XYZBAC, PTXYZD,
ccc     %                    N1FEOC, NFETOI, NBVPSI, NFVPSI )

            IF( Q .LE. QUVEVM ) THEN
C              TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
               IF( NBAC .EQ. 3 ) THEN
C                 PASSAGE AU CYCLE NOC SUIVANT
                  GOTO 80
               ELSE
C                 PASSAGE A UNE AUTRE METHODE DE TRIANGULATION
                  GOTO 30
               ENDIF
            ENDIF

C           LE TETRAEDRE NS2 NS1 NS3 XYZ A UN VOLUME>0 MAIS
C           INTERSECTE T IL UNE FACE DE LA SOUS-ETOILE?
            NF2 = N1FEOC
 29         IF( NF2 .GT. 0 ) THEN
ccc               IF( NF2 .NE. NF1 ) THEN

C                 INTERSECTION TRIANGLE-TETRAEDRE?
                  CALL INTRITET( PTXYZD(1,NFETOI(2,NF2)),
     %                           PTXYZD(1,NFETOI(3,NF2)),
     %                           PTXYZD(1,NFETOI(4,NF2)),
     %                           PTXYZD(1,NS2),
     %                           PTXYZD(1,NS1),
     %                           PTXYZD(1,NS3),
     %                           XYZ,   LINTER )
C                 LINTER : 0 PAS D'INTERSECTION
C                 1 SI UNE ARETE DU TRIANGLE INTERSECTE UNE FACE DU TETRAEDRE
C                 2 SI LE TRIANGLE EST INTERSECTE PAR  UNE ARETE DU TETRAEDRE
C                 3 SI UN DES SOMMETS DU TRIANGLE EST STRICTEMENT INTERNE
C                      AU TETRAEDRE
                  IF( LINTER  .NE. 0 ) THEN
C                    OUI: UN POINT D'INTERSECTION=>ABANDON DE
C                         CETTE TRIANGULATION
                     IF( NBAC .EQ. 3 ) THEN
C                       PASSAGE AU CYCLE NOC SUIVANT
                        GOTO 80
                     ELSE
C                       PASSAGE A UNE AUTRE METHODE DE TRIANGULATION
                        GOTO 30
                     ENDIF
                  ENDIF

ccc               ENDIF

C              PASSAGE A LA FACE SUIVANTE DE LA SOUS-ETOILE
               NF2 = NFETOI(5,NF2)
               GOTO 29
            ENDIF

C           POUR LE FUTUR TRIANGLE
            NS2 = NS3
            NAC = NASUIV

         ENDDO
C        TOUS LES VOLUMES NS2 NS1 NS3 XYZ SONT POSITIFS
         GOTO 32


 30      IF( NBAC .EQ. 3 ) GOTO 33

C        RECHERCHE D'UNE AUTRE TRIANGULATION DES NBAC ARETES DU CYCLE
C        EN PASSANT LE SECOND SOMMET COMME PREMIER SOMMET DU CYCLE
C        ==============================================================
         NBTRIA = NBTRIA + 1
         IF( NBTRIA .GT. NBAC-2 ) THEN
C           TOUTES LES TRIANGULATIONS ISSUES D'UN SOMMET ONT ETE TRAITEES
C           (SEULEMENT SI NBAC<=5 SINON D'AUTRES TRIANGULATIONS EXISTENT)
            GOTO 33
         ENDIF

C        PASSAGE AU SOMMET INITIAL DE LA TRIANGULATION SUIVANT
C        L'ARETE DE DEBUT DU CYCLE NOC
         NAC = N1CIAS( NOC )
C        NAFIN EST ACTUELLEMENT LA DERNIERE ARETE

C        LA SECONDE ARETE DEVIENT LA PREMIERE
         N1CIAS( NOC ) = NSASFVP( 3, NAC )

C        LA DERNIERE ARETE DEVIENT L'AVANT DERNIERE
         NSASFVP( 3, NAFIN ) = NAC

C        L'ANCIENNE PREMIERE DEVIENT LA DERNIERE
         NSASFVP( 3, NAC ) = 0
         NAFIN = NAC
         GOTO 27


C        TOUS LES VOLUMES DES TETRAEDRES NS2 NS1 NS3 XYZ SONT POSITIFS
C        CREATION EFFECTIVE DES TRIANGLES NS2 NS1 NS3 DANS NFETOI
C        POUR FERMER LA SOUS ETOILE DES NBVPSI FACES NFVPSI
C        -------------------------------------------------------------
C        L'ARETE DE DEBUT DU CYCLE
 32      NAC = N1CIAS( NOC )
C        LE SENS POUR ASSURER LA NORMALE VERS L'INTERIEUR
         NS1 = NSASFVP(1,NAC)
         NS2 = NSASFVP(2,NAC)

         DO M = 1, NBAC-2

C           L'ARETE SUIVANTE NS2-NS3 DU CYCLE NOC
            NASUIV = NSASFVP(3,NAC)
            NS3    = NSASFVP(2,NASUIV)
            IF( NS3 .EQ. NS2 ) THEN
                NS3 = NSASFVP(1,NASUIV)
            ENDIF

C           CREATION DU TRIANGLE NFETOI NS1 NS2 NS3
            NBFVAU = NBFVAU + 1

C           RECHERCHE D'UNE FACE VIDE
            NFS    = N1FEVI
            N1FEVI = NFETOI(5,N1FEVI)

C           TETRAEDRE OPPOSE INCONNU
            NFETOI(1,NFS) = -2
C           3 SOMMETS DE LA FACE NORMALE VERS L'INTERIEUR
            NFETOI(2,NFS) = NS1
            NFETOI(3,NFS) = NS2
            NFETOI(4,NFS) = NS3
C           PAS DE FACE SUIVANTE DANS LA SOUS-ETOILE CAR ETAT TRANSITOIRE
            NFETOI(5,NFS) = -3
            PRINT*,'tetrcycl: ajout TRIANGLE VPSI NFETOI(',NFS,')=',
     %             (NFETOI(N,NFS),N=1,5),' V=',V,' Q=',Q

C           AJOUT AUX FACES NFVPSI A V>0 AVEC XYZ
            NBVPSI = NBVPSI + 1
C           LE SIGNE - POUR INDIQUER UN AJOUT NON STANDARD
            NFVPSI( NBVPSI ) = -NFS

C           POUR LE FUTUR TRIANGLE
            NS2 = NS3
            NAC = NASUIV

         ENDDO

         KTITRE='tetrcycl:          FACES VPSI apres TRIANGULATION du CI
     %RCUIT'
         WRITE(KTITRE(11:15), '(I5)') NBVPSI
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD,
     %                 N1FEOC, NFETOI, NBVPSI, NFVPSI )

C        PASSAGE AU CYCLE SUIVANT
         GOTO 80


C        AUCUNE TRIANGULATION DU CYCLE NOC AVEC TOUS LES TRIANGLES+XYZ A V>0
C        et Q<CORRECT et SANS INTERSECTION AVEC UNE FACE

C        RECHERCHE D'UN POINT XYZ PROCHE DU BARYCENTRE XYZBAC DES ARETES
C        DU CYCLE NOC FORMANT AVEC LES ARETES DU CYCLE ET LE POINT XYZBAC
C        UN TETRAEDRE DE VOLUME POSITIF ET SANS INTERSECTION AVEC LES
C        AUTRES FACES NFETOI DE L'ETOILE
C        ====================================================================
C        XYZBAC+FACES NBVPSI DE VOLUME>0 ET SANS INTERSECTION AVEC LES FACES?
 33      DO K=1,NBVPSI
            NFS = ABS( NFVPSI( K ) )
            CALL QUATETD( PTXYZD(1,ABS(NFETOI(2,NFS))),
     %                    PTXYZD(1,NFETOI(3,NFS)),
     %                    PTXYZD(1,NFETOI(4,NFS)),
     %                    XYZBAC,
     %                    ARMIN, ARMAX, SURFTR, V, Q )

            PRINT*,'tetrcycl: NFETOI(',NFS,')=',(NFETOI(N,NFS),N=2,4),
     %             ' +XYZBAC',XYZBAC,
     %             ' => Tetra avec V=',V,' Q=',Q

            IF( Q .LE. QUVEVM ) THEN
C              TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
               PRINT*,'tetrcycl: XYZBAC N''EST PAS INTERNE A L''ETOILE D
     %ES NBVPSI FACES. A AMELIORER en BOUGEANT XYZBAC'
            ENDIF

C           LE TETRAEDRE NFETOI(NFS)+XYZBAC A UN VOLUME>0 MAIS
C           INTERSECTE T IL UNE FACE DE LA SOUS-ETOILE?
            NF2 = N1FEOC
 34         IF( NF2 .GT. 0 ) THEN
               IF( NF2 .NE. NFS ) THEN
C                 INTERSECTION TRIANGLE-TETRAEDRE?
                  CALL INTRITET( PTXYZD(1,NFETOI(2,NF2)),
     %                           PTXYZD(1,NFETOI(3,NF2)),
     %                           PTXYZD(1,NFETOI(4,NF2)),
     %                           PTXYZD(1,ABS(NFETOI(2,NFS))),
     %                           PTXYZD(1,NFETOI(3,NFS)),
     %                           PTXYZD(1,NFETOI(4,NFS)),
     %                           XYZBAC,  LINTER )
                  IF( LINTER  .NE. 0 ) THEN
                     PRINT*,'tetrcycl: XYZBAC+NFETOI(',NFS,') INTERSECTE
     % UNE AUTRE FACE A AMELIORER en BOUGEANT XYZBAC'
C                    ABANDON DU CYCLE NOC
                     GOTO 80
                  ENDIF
               ENDIF
C              PASSAGE A LA FACE SUIVANTE DE LA SOUS-ETOILE
               NF2 = NFETOI(5,NF2)
               GOTO 34
            ENDIF
         ENDDO

         KTITRE='tetrcycl:      FACES+XYZ TETRAEDRES + POINT XYZBAC Bary
     %centre du CYCLE'
         WRITE( KTITRE(11:15), '(I5)' ) NBVPSI
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO3( KTITRE(1:L), XYZ, XYZBAC, PTXYZD,
     %                 N1FEOC, NFETOI, NBVPSI, NFVPSI )

C        CALCUL DU BARYCENTRE DES NBVPSI FACES
         DO L=1,3
            XYZBAV(L) = 0
         ENDDO
         DO K=1,NBVPSI
            NF = ABS( NFVPSI( K ) )
            DO J=1,3
               NS = ABS( NFETOI( 1+J, NF ) )
               DO L=1,3
                  XYZBAV(L) = XYZBAV(L) + PTXYZD(L,NS)
               ENDDO
            ENDDO
         ENDDO
         DO L=1,3
            XYZBAV(L) = XYZBAV(L) / (3*NBVPSI)
         ENDDO

C        DISTANCE ENTRE XYZBAV et XYZBAC PONDEREE
         D = 0D0
         DO L=1,3
            D = D + ( XYZBAC(L) - XYZBAV(L) ) ** 2
         ENDDO
         D = SQRT( D ) / 15

         KTITRE='tetrcycl:      FACES+XYZBAV TETRAEDRES + XYZBAV->XYZ BA
     %C Barycentre des FACES V>0'
         WRITE(KTITRE(11:15), '(I5)' ) NBVPSI
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO3( KTITRE(1:L), XYZBAV, XYZBAC, PTXYZD,
     %                 N1FEOC,      NFETOI, NBVPSI, NFVPSI )

C        NOMBRE D'AUGMENTATIONS (*2) DE D
         NBAUGD = 0

C        RECHERCHE DU MEILLEUR POINT XYZESS AUTOUR DE XYZBAV
C        PARMI LES SOMMETS ET MILIEUX D'UN CUBE DE COTE 2*D
 35      QUAMIMX = -2.

         DO K3=-1,1,1
            XYZESS(3) = XYZBAV(3) + K3 * D
            DO K2=-1,1,1
               XYZESS(2) = XYZBAV(2) + K2 * D
               DO 39 K1=-1,1,1
                  XYZESS(1) = XYZBAV(1) + K1 * D

ccc            KTITRE='tetrcycl:      FACES+XYZ TETRAEDRES + POINT XYZESS  
ccc     % AUTOUR du BARYCENTRE des FACES NFVPSI'
ccc            WRITE(KTITRE(11:15), '(I5)' ) NBVPSI
ccc            CALL SANSDBL( KTITRE, L )
ccc            PRINT*, KTITRE(1:L)
ccc            CALL TRFETO3( KTITRE(1:L), XYZESS, XYZBAC, PTXYZD,
ccc     %                    N1FEOC,      NFETOI, NBVPSI, NFVPSI )

C                 TEST DES NBAC TRIANGLES ARETE+XYZESS
                  QMIN = 2.
                  NAC  = N1CIAS( NOC )
                  DO M=1,NBAC

C                    L'ARETE M DU CYCLE NOC
                     NS1 = NSASFVP(1,NAC)
                     NS2 = NSASFVP(2,NAC)

C                    LE VOLUME NS1 NS2 XYZBAC XYZESS EST IL POSITIF?
                     CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                             XYZBAC, XYZESS,
     %                             ARMIN, ARMAX, SURFTR, V, Q )

                     IF( Q .LE. QUVEVM ) THEN
C                       TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
C                       => RECHERCHE D'UN AUTRE POINT XYZESS
                        GOTO 39
                     ENDIF

C                    LE TETRAEDRE NS1 NS2 XYZESS XYZ A UN VOLUME>0
C                    INTERSECTE T IL UNE FACE DE LA SOUS-ETOILE?
                     NF2 = N1FEOC
 38                  IF( NF2 .GT. 0 ) THEN
ccc                        IF( NF2 .NE. NF1 ) THEN
C                          INTERSECTION TRIANGLE-TETRAEDRE?
                           CALL INTRITET( PTXYZD(1,NFETOI(2,NF2)),
     %                                    PTXYZD(1,NFETOI(3,NF2)),
     %                                    PTXYZD(1,NFETOI(4,NF2)),
     %                                    PTXYZD(1,NS1),
     %                                    PTXYZD(1,NS2),
     %                                    XYZBAC,
     %                                    XYZESS,  LINTER )
                           IF( LINTER  .NE. 0 ) THEN
C                             OUI: UN POINT D'INTERSECTION
C                             => ABANDON DE CE POINT XYZESS
                              GOTO 39
                           ENDIF
ccc                        ENDIF

C                       PASSAGE A LA FACE SUIVANTE DE LS SOUS-ETOILE
                        NF2 = NFETOI(5,NF2)
                        GOTO 38

                     ENDIF

                     IF( Q .LT. QMIN ) THEN
                        QMIN = Q
                     ENDIF

C                    L'ARETE SUIVANTE DE NAC
                     NAC = NSASFVP(3,NAC)

                  ENDDO

C                 RECHERCHE DU MAX DES MIN DES QUALITES DES TETRAEDRES
                  IF( QMIN .GT. QUAMIMX) THEN
                     QUAMIMX = QMIN
                     XYZMIMX(1) = XYZESS(1)
                     XYZMIMX(2) = XYZESS(2)
                     XYZMIMX(3) = XYZESS(3)
                  ENDIF

 39           ENDDO
            ENDDO
         ENDDO


         IF( QUAMIMX .EQ. -2. ) THEN

C           AUCUN POINT XYZMIMX CORRECT
C           ---------------------------
            PRINT*,'tetrcycl: PAS de POINT XYZMIMX CORRECT pour l''AUGME
     %NTATION',NBAUGD
            NBAUGD = NBAUGD + 1
            IF( NBAUGD .LE. 5 ) THEN

C              CHANGEMENT de XYZESS PAR ELARGISSEMENT DU DIAMETRE DU CUBE
               D = D * 2
               GOTO 35

            ELSE

               IF( NBSUFV .EQ. 0 ) THEN

C                 PAS DE TRIANGULATION FAVORABLE DU CYCLE
C                 ESSAI de SUPPRESSION DE L'UNE DES FACE V>0 AYANT
C                 2 ARETES CONSECUTIVES DU CYCLE
C                 ================================================
                  PRINT*,'tetrcycl: SUPPRESSION 1 FACE V>0 AYANT 2 ARETE
     %S DU CYCLE'
                  NBSUFV = 1
C                 L'ARETE DE DEBUT DU CYCLE
                  NAC = N1CIAS( NOC )

C                 LE SENS POUR ASSURER LA NORMALE VERS L'INTERIEUR
 40               NS1 = NSASFVP(1,NAC)
                  NS2 = NSASFVP(2,NAC)

C                 NUMERO NFVPSI DE LA FACE D'ARETE NAC
                  N1 = NSASFVP(4,NAC)

C                 L'ARETE SUIVANTE DU CYCLE NOC
                  NASUIV = NSASFVP(3,NAC)

C                 NUMERO NFVPSI DE LA FACE D'ARETE NASUIV
                  N2 = NSASFVP(4,NASUIV)

                  IF( N1 .EQ. N2 ) THEN
C                    LA FACE AYANT CES 2 ARETES DU CYCLE EST MARQUEE
C                    COMME SUPPRIMEE DES FACES NFVPSI
                     NFVPSI(N1) = 0
                     GOTO 45
                  ENDIF

C                 L'ARETE SUIVANTE DU CYCLE NOC
                  NASUIV = NSASFVP(3,NAC)
                  IF( NASUIV .GT. 0 ) THEN
                     NAC = NASUIV
                     GOTO 40
                  ENDIF

C                 RETRAIT DU TABLEAU NFVPSI DES FACES SUPPRIMEES
 45               N = 0
                  DO L=1,NBVPSI
                     NF = ABS( NFVPSI( L ) )
                     IF( NF .GT. 0 ) THEN
                        N = N + 1
                        NFVPSI( N ) = NFVPSI( L )
                     ENDIF
                  ENDDO

C                 RESTE T IL ASSEZ DE FACES NFVPSI POUR CONTINUER?
                  IF( NBVPSI .NE. N .AND. NBVPSI .GE. 1 ) THEN
C                    OUI: TRIANGULATION AVEC LES FACES NFVPSI RESTANTES
                     NBVPSI = N
                     GOTO 10
                  ENDIF

               ENDIF

C              LE TABLEAU NFVPSI N'EST PLUS EFFICACE
C              IMPOSSIBLE DE TROUVER UN BON POINT XYZESS
C              -----------------------------------------
               PRINT *,'tetrcycl: IMPOSSIBLE de TROUVER UN XYZESS CORREC
     %T. PAS de POINT XYZMIMX AJOUTE'
               PRINT *,'tetrcycl: RESTAURATION des',NBVPSI0,
     %                 ' TETRAEDRES V>0 INITIAUX'

C              RESTAURATION DU TABLEAU NFVPSI0 DANS NFVPSI
C              -------------------------------------------
               NBVPSI = NBVPSI0
               DO N=1,NBVPSI0
                  NFVPSI( N ) = NFVPSI0( N )
               ENDDO

            ENDIF

         ELSE

C           XYZMIMX EST LE POINT NBSOMM A AJOUTER CENTRE DES TETRAEDRES
C           CONSTRUITS A PARTIR DE CE POINT ET DES NVPSI FACES NFVPSI
C           ===========================================================
C           LE POINT CENTRE DU CYCLE
C           TAILLE_IDEALE DE L'ARETE EN XYZBAC
            CALL TAILIDEA( NOFOTI, XYZBAC, NCODEV, XYZBAC(4) )
            NBSOMM = NBSOMM + 1
            DO N=1,4
               PTXYZD( N, NBSOMM ) = XYZBAC( N )
            ENDDO
            PRINT*,'tetrcycl: Pt AJOUTE NBSOMM=',NBSOMM,
     %             ' est XYZBAC=',(PTXYZD(N,NBSOMM),N=1,4),
     %             ' CENTRE du CYCLE', NOC

C           POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
            NPSOFR( NBSOMM ) = 0

C           CONSTRUCTION DES NBAC TRIANGLES ARETE+XYZMIMX dans NFETOI
C           ---------------------------------------------------------
            NAC = N1CIAS( NOC )
            DO M = 1, NBAC

C              L'ARETE M DU CYCLE NOC
               NS1 = NSASFVP(1,NAC)
               NS2 = NSASFVP(2,NAC)

C              CREATION DU TRIANGLE NFETOI NS1 NS2 NBSOMM
C              RECHERCHE D'UNE FACE VIDE
               NFS    = N1FEVI
               N1FEVI = NFETOI(5,N1FEVI)

C              TETRAEDRE OPPOSE INCONNU
               NFETOI(1,NFS) = -1
C              3 SOMMETS DE LA FACE NORMALE VERS L'INTERIEUR
               NFETOI(2,NFS) = NS1
               NFETOI(3,NFS) = NS2
               NFETOI(4,NFS) = NBSOMM
C              PAS DE FACE SUIVANTE DANS L'ETOILE CAR ETAT TRANSITOIRE
               NFETOI(5,NFS) = -1
               PRINT*,'tetrcycl: creation TRIANGLE VPSI NFETOI(',NFS,
     %                ')=',(NFETOI(N,NFS),N=1,5),' du CYCLE',NOC

C              AJOUT AUX FACES A VOLUME>0 AVEC XYZ
               NBVPSI = NBVPSI + 1
C              LE SIGNE - POUR INDIQUER UN AJOUT NON STANDARD
               NFVPSI( NBVPSI ) = -NFS

C              L'ARETE SUIVANTE DE NAC
               NAC = NSASFVP(3,NAC)

            ENDDO

C              TRACE DES TRIANGLES AVEC NBSOMM SOMMET CENTRAL
      KTITRE='tetrcycl:         FACES VPSI apres AJOUT de XYZMIMX         
     %   '
            WRITE(KTITRE(11:15), '(I5)') NBVPSI
            WRITE(KTITRE(53:61), '(I9)') NBSOMM
            CALL SANSDBL( KTITRE, L )
            PRINT*, KTITRE(1:L)
            CALL TRFETO3( KTITRE(1:L), XYZMIMX, XYZMIMX, PTXYZD,
     %                    N1FEOC, NFETOI, NBVPSI, NFVPSI )

         ENDIF

C        FIN DU CYCLE NOC
 80   ENDDO

 9999 IF( NBVPSI0 .NE. NBVPSI .OR. NBSOMM0 .NE. NBSOMM ) THEN
         PRINT*,'tetrcycl: AJOUT de',NBVPSI-NBVPSI0,
     %     ' FACES des CYCLES  et',NBSOMM-NBSOMM0,' POINTS XYZD'
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
