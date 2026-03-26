      SUBROUTINE SIPTETO( NOAPPEL,   KTITRE,  VMOYEN, PTXYZD, NPSOFR,
     %                    N1FEOC,    NFETOI,  NBFETO,
     %                    MAXFAC,    XYZMAX,
     %                    MXVPSI,    NBVPSI,  NFVPSI, NFVPSI0,
     %                    MXCIAS,    NBCIAS,  N1CIAS,
     %                    MXASFVP,   NBASFVP, NSASFVP,
     %                    QUAMINEX,  QUAMIMX, MIARSIMI,
     %                    NOCHOIXYZ, VETXYZ,  RAVEVM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SIMULER LA TETRAEDRISATION DES FACES SIMPLES DE L'ETOILE
C -----    NFETOI ISSUES DE N1FEOC A PARTIR D'UN POINT XYZ CALCULES
C          POUR MAXIMISER LE NOMBRE DE TETRAEDRES FACE+XYZ DE VOLUMES
C          POSITIFs ET SANS INTERSECTION AVEC LES FACES

C ENTREES:
C --------
C NOAPPEL: =0 => RETOUR APRES LE TEST DE LA GRILLE AUTOUR DU BARYCENTRE
C          >0 => PAS DE TEL RETOUR
C KTITRE : TITRE DU TRACE
C VMOYEN : VOLUME MOYEN DES TETRAEDRES INITIAUX DE L'ETOILE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0  SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          5: CHAINAGE SUR LA FACE SUIVANTE OCCUPEE OU VIDE
C NBFETO : NOMBRE ACTUEL  DE FACES DE L'ETOILE ISSUE DE N1FEOC

C MAXFAC : NOMBRE MINIMUM DE FACES A VOLUME POSITIF AU DELA DUQUEL
C          ON SE CONTENTE DU MAX TROUVE POUR JOINDRE LE POINT XYZ AUX FACES
C MXVPSI : NOMBRE MAXIMAL DE FACES DU TABLEAU NFVPSI et NFVPSI0
C MXCIAS  : NOMBRE MAXIMAL DE CIRCUITS
C MXFASFVP: NOMBRE MAXIMAL DE FACES SIMPLES DES CIRCUITS

C QUAMINEX:QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE C-A-D
C          QUALITE MIN DES MAX DES QUALITES MIN REQUISE POUR RETENIR
C          UN POINT QUI AVEC UN ENSEMBLE DE FACES CONTIGUES DONNENT
C          DES TETRAEDRES DE VOLUME POSITIF SANS INTERSECTION AVEC
C          LES FACES DE L'ETOILE

C SORTIES:
C --------
C XYZMAX : POINT CHOISI POUR FORMER AVEC DES FACES DE L'ETOILE DES TETRAEDRES
C          MAXIMISANT LE NOMBRE DE TETRAEDRES FACE+XYZMAx A V>0
C NBVPSI : NOMBRE DE FACES DES TETRAEDRES DE VOLUME>0 SANS INTERSECTION
C          AVEC LES FACES
C NFVPSI : NUMERO NFETOI DES NBVPSI FACES DONNANT AVEC XYZ UN VOLUME>0
C          ET SANS INTERSECTION AVEC LES FACES DE L'ETOILE

C NBCIAS : NOMBRE DE CIRCUITS DES ARETES SIMPLES DES FACES NFVPSI DANS NFETOI
C N1CIAS : NO DANS NSASFVP DE LA PREMIERE ARETE DES NBCIAS CIRCUITS
C NBASFVP: NOMBRE DE FACES DES ARETES SIMPLES DES NBCIAS CIRCUITS
C NSASFVP: 1: NUMERO DU SOMMET 1  DE L'ARETE SIMPLE
C          2: NUMERO DU SOMMET 2  DE L'ARETE SIMPLE
C          3: NUMERO DANS NSASFVP DE L'ARETE SIMPLE SUIVANTE
C          4: NUMERO DANS NFVPSI DE LA FACE DE CETTE ARETE SIMPLE

C QUAMIMX: QUALITE MINIMALE DES NBVPSI FUTURS TETRAEDRES A VOLUME>0
C MIARSIMI:NOMBRE MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
C NOCHOIXYZ: NUMERO DU CAS DONNANT LE NOMBRE MAXIMAL DE FACES A VOLUME>0
C          SANS INTERSECTION ET OBTENUES A TRAVERS LES ARETES ADJACENTES
C VETXYZ : VOLUME DE L'ETOILE DE CENTRE XYZ DES TETRAEDRES FACE-XYZ
C RAVEVM : RAPPORT DU VOLUME DE L'ETOILE ACTUELLE AU VOLUME DE
C          NBFETO FOIS LE VOLUME MOYEN D'UN TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER   DECEMBRE 2014
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY OCTOBRE  2017
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY MAI      2018
C2345X7..............................................................012
ccc      include"./incl/darete.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      DOUBLE PRECISION  PTXYZD(4,*), XYZ(3), XYZBE(3),
     %                  XYZN(3), XYZMIL(3), VN(3), XYZMAX(3), CBTR(3),
     %                  REDUCT, D, DMIN, DMAX, COTMOY, PROSCD, COEFQP,
     %                  SURTRD, XYZMIX(3,2), XYZBAR(3), XYZBA2(3),
     %                  VETXYZ, VMOYEN, C, V, RAVEVM,
     %                  SURFTR(4), ARMIN, ARMAX

      INTEGER           NPSOFR(*), NFETOI(5,*),
     %                  NFVPSI0(MXVPSI), NFVPSI(MXVPSI),
     %                  N1CIAS(MXCIAS), NSASFVP(4,MXASFVP)
      REAL              QUAMIN, QUAMIMX
      CHARACTER*(*)     KTITRE

      TRACTE0 = TRACTE
      NBNST = 0
      NST0  = 0
      NST1  = 0

cccC     TRACE DES FACES DE L'ETOILE ET DU POINT XYZ
ccc 1    KTITRE='Entree sipteto:        FACES de NFETOI NoAPPEL=     '
ccc      WRITE(KTITRE(17:21),'(I5)') NBFETO
ccc      WRITE(KTITRE(48:48),'(I1)') NOAPPEL
ccc      CALL SANSDBL( KTITRE, L )
ccc      PRINT*, KTITRE(1:L)
ccc      CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD, N1FEOC, NFETOI,
ccc     %              0, NFVPSI )

 1    NOCHOIXYZ = 0
      NBVPSIMAX = 0
      NFVPSIMAX = 0
C     MIARSIMI: NOMBRE MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
      MIARSIMI  = 2**30
      QUAMIMX   = 0.0

C     CHOIX 1: XYZ UN SOMMET D'UNE GRILLE CUBIQUE AUTOUR DU
C              BARYCENTRE XYZBE DE TOUTES LES FACES DE L'ETOILE
C     ---------------------------------------------------------
C     CALCUL DES COORDONNEES EXTREMES ET DU POINT XYZBE DU BARYCENTRE
C     DES FACES N1FEOC A 0 DE L'ETOILE DEFINIE DANS NFETOI
      CALL XYZBF10( PTXYZD, N1FEOC, 0, NFETOI, NBFETO, XYZMIX, XYZBE )

C     CALCUL DU VOLUME DES TETRAEDRES FACES+XYZBE POUR SAVOIR
C     SI L'ETOILE EST QUASI-PLANE (RAVEVM<=1D-2)
      CALL VOTEET( PTXYZD, XYZBE, N1FEOC, NFETOI, NBFETO, V )

C     RAPPORT DU VOLUME V DE L'ETOILE ACTUELLE AU VOLUME DE
C     NBFETO VOLUMES MOYEN D'UN TETRAEDRE
      RAVEVM = V / (NBFETO*VMOYEN)

      IF( RAVEVM .LT. 1D-2 ) THEN
C        COEFFICIENT REDUCTEUR POUR UNE ETOILE QUASI-PLANE
         COEFQP = 0.6D0
      ELSE
C        COEFFICIENT REDUCTEUR POUR UNE ETOILE CORRECTEMENT VOLUMINEUSE
         COEFQP = 1D0
      ENDIF


      IF( NOAPPEL .EQ. 0 ) THEN

C        1-ER APPEL DE SIPTETO: GRILLE AUTOUR DU BARYCENTRE
C        ==================================================
C        LES PAS ENTRE SOMMETS DE LA GRILLE
ccc         NBSUB = 4
ccc         DO N=1,3
ccc            XYZN(N) = (XYZMIX(N,2)-XYZMIX(N,1)) / (2.2 * NBSUB) * COEFQP
ccc         ENDDO

         NBSUB = 3
         DO N=1,3
            XYZN(N) = (XYZMIX(N,2)-XYZMIX(N,1)) / (1.6 * NBSUB) * COEFQP
         ENDDO

C        METHODE DU MAXIMUM DE FACE+XYZ DE VOLUME>0 SANS INTERSECTION ET DE
C        MAXIMUM DU MINIMUM DES QUALITES DES TETRAEDRES NF+XYZ
C        ------------------------------------------------------------------
         MXFPVP = 0
         DO N = -NBSUB, NBSUB, 1
            XYZ(3) = XYZBE(3) + N * XYZN(3)
            DO M = -NBSUB, NBSUB, 1
               XYZ(2) = XYZBE(2) + M * XYZN(2)
               DO L = -NBSUB, NBSUB, 1
                  XYZ(1) = XYZBE(1) + L * XYZN(1)

C                 DETERMINER LE NOMBRE MAXIMAL DE TETRAEDRES XYZ+FACES
C                 SIMPLES DE VOLUME POSITIF et SANS INTERSECTION AVEC
C                 LES AUTRES FACES SIMPLES DE L'ETOILE
                  CALL FAPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI,
     %                           MXVPSI, NBVPSI, NFVPSI,
     %                           NBFETO, NBVPIN, QUAMIN )

                  IF( NBVPSI .GT. MXFPVP ) THEN
                     MXFPVP    = NBVPSI
                     QUAMIMX   = QUAMIN
                     XYZMAX(1) = XYZ(1)
                     XYZMAX(2) = XYZ(2)
                     XYZMAX(3) = XYZ(3)
                     NFVPSIMAX = NFVPSI(1)
                  ELSE IF( NBVPSI .EQ. MXFPVP   .AND.
     %                     QUAMIN .GT. QUAMIMX ) THEN
                     QUAMIMX   = QUAMIN
                     XYZMAX(1) = XYZ(1)
                     XYZMAX(2) = XYZ(2)
                     XYZMAX(3) = XYZ(3)
                     NFVPSIMAX = NFVPSI(1)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

C        LE MEILLEUR RESULTAT AVEC LE MAXIMUM DE FACES+XYZMAX A V>0
         NBVPSI    = MXFPVP
         NBVPSIMAX = MXFPVP
         MIARSIMI  = NBFETO
         NOCHOIXYZ = 0

         IF( NBVPSI .EQ. NBFETO ) THEN

C           XYZMAX TETRAEDRISE TOUTES LES FACES DE NFETOI
C           =============================================
            MIARSIMI = -1
            NOCHOIXYZ = 1
       KTITRE='sipteto:  TOUTES LES        FACES+XYZMAX -> TETRAEDRES V>
     %0 NOCHOIXYZ=   '
            WRITE(KTITRE(23:27), '(I5)' ) NBFETO
            WRITE(KTITRE(71:72), '(I2)' ) NOCHOIXYZ

         ELSE

       KTITRE='sipteto: PARMI les       FACES seules       FACES+XYZMAX 
     %-> TETRAEDRES V>0 NOCHOIXYZ=     '
            WRITE(KTITRE(20:24), '(I5)' ) NBFETO
            WRITE(KTITRE(39:43), '(I5)' ) NBVPSI
            WRITE(KTITRE(87:88), '(I2)' ) NOCHOIXYZ

         ENDIF

C        SORTIE AU PREMIER APPEL DE SIPTETO
C        RECALCULER NFVPSI POUR LE POINT XYZMAX, 
C        TABLEAU QUI A PEUT ETRE ETE ECRASE PAR LES POINTS SUIVANTS XYZMAX
         CALL FAPTVPSI( XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %                  MXVPSI, NBVPSI, NFVPSI,
     %                  NBFETO, NBVPIN, QUAMIMX )
         GOTO 9900

      ENDIF


C     METHODE DE RECHERCHE D'UN CIRCUIT DE NOMBRE MINIMUM D'ARETES
C     ------------------------------------------------------------
      NBVPSI = 0
      KTITRE='sipteto ap2:        FACES+XYZ TETRAEDRES V>0       FACES N
     %OCHOIX=  QUAMIN=              VETXYZ=                        '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(46:50), '(I5)'   ) NBFETO
      WRITE(KTITRE(66:67), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(75:87), '(G13.5)') QUAMIN
      WRITE(KTITRE(96:108),'(G13.5)') VETXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZBE, XYZBE, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)

C     PAS DE LA GRILLE DE SOMMETS A ESSAYER
      NBSUB = 3
      DO N=1,3
         XYZN(N) = (XYZMIX(N,2)-XYZMIX(N,1)) / (2.0 * NBSUB) * COEFQP
C        POUR UN TRACE CORRECT SI PAS DE VALEUR A XYZMAX
         XYZMAX(N) = XYZBE(N)
      ENDDO

C     LA COUCHE EXTERNE NBSUB N'EST PAS TESTEE
      NBSUB = 2

      DO N = -NBSUB, NBSUB, 1
         XYZ(3) = XYZBE(3) + N * XYZN(3)
         DO M = -NBSUB, NBSUB, 1
            XYZ(2) = XYZBE(2) + M * XYZN(2)
            DO L = -NBSUB, NBSUB, 1
               XYZ(1) = XYZBE(1) + L * XYZN(1)

C              SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C              RECHERCHE DES FACES CONTINUES DE L'ETOILE FORMANT AVEC XYZ
C              UN TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C              AUTRES FACES DE L'ETOILE
               CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
     %                        MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
     %                        MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
     %                        MXASFVP,  NBASFVP, NSASFVP,
     %                        VETXYZ,   QUAMIN, IERR )
C              NBVPSI: NBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C              NFVPSI: NO NFETOI DES FACES DONNANT UN VOLUME>=0
C              MIARSICI:NBRE MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
C              QUAMIN: QUALITE MINIMALE MAXIMISEE DES TETRAEDRES NF-XYZ
C              RECHERCHE DU POINT DONNANT LE MINIMUM D'ARETES SIMPLES DES CIRCUITS
               IF( NBVPSI .GE. 3 .AND. NBVPSI .GE. NBFETO/3 .AND.
     %           ( MIARSICI.LT.MIARSIMI .OR.
     %           ( MIARSICI.EQ.MIARSIMI .AND. QUAMIN.GT.QUAMIMX) )) THEN

                  NOCHOIXYZ = 1
                  NBVPSIMAX = NBVPSI
                  XYZMAX(1) = XYZ(1)
                  XYZMAX(2) = XYZ(2)
                  XYZMAX(3) = XYZ(3)
                  NFVPSIMAX = NFVPSI(1)
                  MIARSIMI  = MIARSICI
                  QUAMIMX   = QUAMIN

               ENDIF
              
            ENDDO
         ENDDO
      ENDDO

      IF( NBVPSIMAX .GT. 0 ) THEN
      KTITRE='sipteto ap3:        FACES+XYZ TETRAEDRES V>0 pour        F
     %ACES NOCHOIXYZ=   '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
      WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)
      ENDIF

      IF( NBVPSIMAX .GE. MAXFAC )THEN
         GOTO 8000
      ENDIF


cccC     CHOIX 2: XYZ=POINT SOMMET D'UNE GRILLE CUBIQUE AUTOUR DU
cccC              MILIEU D'UNE ARETE DES FACES DE L'ETOILE
cccC     --------------------------------------------------------
ccc      NBSUB = 1
ccc      NF    = N1FEOC
ccc 5    IF( NF .GT. 0 ) THEN

cccC        CALCUL DE XYZMIL MILIEU DE L'ARETE MM DE LA FACE NF
ccc         NS1 = NFETOI(4,NF)
ccc         DO MM=1,3

cccC           ARETE MM DE LA FACE NF DE SOMMETS NS1-NS2
ccc            NS2 = NFETOI(1+MM,NF)

cccC           LES 3 COORDONNEES DU MILEU DE L'ARETE MM DE NF
ccc            DO K=1,3
ccc               XYZMIL(K) = (PTXYZD(K,NS1) + PTXYZD(K,NS2)) * 0.5D0
ccc            ENDDO

ccc            D = ARETE / 5D0
ccc            DO N=-NBSUB,NBSUB,1
ccc               XYZ(3) = XYZMIL(3) + N * D
ccc               DO M=-NBSUB,NBSUB,1
ccc                  XYZ(2) = XYZMIL(2) + M * D
ccc                  DO 8 L=-NBSUB,NBSUB,1

cccC                    LE POINT A ESSAYER NE DOIT PAS ETRE LE MILIEU DE L'ARETE
ccc                     IF( L.EQ.0 .AND. M.EQ.0 .AND. N.EQ.0 ) GOTO 8
ccc                     XYZ(1) = XYZMIL(1) + L * D

cccC                    DETERMINER L'ENSEMBLE DES FACES CONTIGUES DE L'ETOILE FORMANT
cccC                    AVEC LE POINT XYZ DES TETRAEDRES DE VOLUMES>0
cccC                    SANS INTERSECTION AVEC LES FACES DE L'ETOILE TELLES QUE
cccC                    CES FACES MAXIMISENT LE MINIMUM DES QUALITES DES TETRAEDRES
ccc                     CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
ccc     %                              MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
ccc     %                              MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
ccc     %                              MXASFVP,  NBASFVP, NSASFVP,
ccc     %                              VETXYZ,   QUAMIN, IERR )
cccC                    NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
cccC                    NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0 ET
cccC                             MAXIMISENT LE MINIMUM DES QUALITES DES TETRAEDRES
cccC                     MIARSICI: NB MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
cccC                    RECHERCHE DU POINT XYZMAX DONNANT LE MAXIMUM DE FACES
cccC                    ET A EGALITE LE POINT QUI MAXIMISE LE MINIMUM DES QUALITES
ccc                     IF( NBVPSI .GT. 0 .AND. QUAMIN .GT. QUAMINEX ) THEN
ccc                        IF( NBVPSI .GT. NBVPSIMAX .OR.
ccc     %                     (NBVPSI .EQ. NBVPSIMAX .AND.
ccc     %                      QUAMIN .GT. QUAMIMX) ) THEN
ccc                           NOCHOIXYZ = 2
ccc                           NBVPSIMAX = NBVPSI
ccc                           DO K=1,3
ccc                              XYZMAX(K) = XYZ(K)
ccc                           ENDDO
ccc                           NFVPSIMAX = NFVPSI(1)
ccc                           QUAMIMX = QUAMIN
ccc                           print *,'SIpteto: nochoixyz=',nochoixyz,
ccc     %                  ' NBVPSIMAX=',NBVPSIMAX,' NFVPSIMAX=',NFVPSIMAX,
ccc     %                  ' QUAMIMX=',QUAMIMX,' XYZMAX=',XYZMAX
ccc         KTITRE='      FACES GRILLE MILIEU ARETE de V>0 et NOCHOIXYZ=2'
ccc                         CALL SANSDBL( KTITRE, L )
ccc                         CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD,
ccc     %                                 N1FEOC, NFETOI, NBVPSI, NFVPSI )
ccc                        ENDIF
ccc                     ENDIF

ccc 8                ENDDO
ccc               ENDDO
ccc            ENDDO

cccC           PASSAGE A L'ARETE SUIVANTE DE NF
ccc            NS1 = NS2
ccc         ENDDO

cccC        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
ccc         NF = NFETOI(5,NF)
ccc         GOTO 5

ccc      ENDIF

ccc      IF( NBVPSIMAX .GE. MAXFAC ) THEN
ccc         GOTO 8000
ccc      ENDIF


C     CHOIX 3: XYZ=POINT MILIEU DU SEGMENT JOIGNANT LE BARYCENTRE
C              DE 2 FACES DONT CELLE DE SURFACE MAXIMALE ET
c              SANS INTERSECTION AVEC D'AUTRES FACES
C     -----------------------------------------------------------
C     RECHERCHE DE LA FACE DE PLUS GRANDE SURFACE
      NFMAX = 0
      DMAX  = 0D0
      NF    = N1FEOC
 10   IF( NF .GT. 0 ) THEN
C        CALCUL DE LA SURFACE DE LA FACE NF
         D = SURTRD( PTXYZD(1,NFETOI(2,NF)),
     %               PTXYZD(1,NFETOI(3,NF)),
     %               PTXYZD(1,NFETOI(4,NF)) )
         IF( D .GT. DMAX ) THEN
            DMAX  = D
            NFMAX = NF
         ENDIF
         NF = NFETOI(5,NF)
         GOTO 10
      ENDIF

      IF( NFMAX .GT. 0 ) THEN

C        BARYCENTRE DE LA FACE NFMAX
         XYZBAR(1) = 0D0
         XYZBAR(2) = 0D0
         XYZBAR(3) = 0D0
         DO I=2,4
            NS = NFETOI(I,NFMAX)
            DO K=1,3
               XYZBAR(K) = XYZBAR(K) + PTXYZD(K,NS)
            ENDDO
         ENDDO
         DO K=1,3
            XYZBAR(K) = XYZBAR(K) / 3D0
         ENDDO

C        PARCOURS DES AUTRES FACES
C        CALCUL DU SEGMENT JOIGNANT LES 2 BARYCENTRES
         NF = N1FEOC
 11      IF( NF .GT. 0 ) THEN

            IF( NF .NE. NFMAX ) THEN
C              BARYCENTRE DE LA FACE NF
               XYZBA2(1) = 0D0
               XYZBA2(2) = 0D0
               XYZBA2(3) = 0D0
               DO I=2,4
                  NS = NFETOI(I,NF)
                  DO K=1,3
                     XYZBA2(K) = XYZBA2(K) + PTXYZD(K,NS)
                  ENDDO
               ENDDO
               DO K=1,3
                  XYZBA2(K) = XYZBA2(K) / 3D0
               ENDDO

C              POINT D'INTERSECTION DE XYZBAR+XYZBA2
C              AVEC UNE AUTRE FACE DE L'ETOILE
               NF1 = N1FEOC
 14            IF( NF1 .GT. 0 ) THEN
                  IF( NF1 .NE. NFMAX .AND. NF1 .NE. NF ) THEN
                     CALL INARTR( XYZBAR, XYZBA2,
     %                            PTXYZD(1,NFETOI(2,NF1)),
     %                            PTXYZD(1,NFETOI(3,NF1)),
     %                            PTXYZD(1,NFETOI(4,NF1)),
     %                            LINTER, XYZ, CBTR )
C                    LINTER=1 SI XYZBAR-XYZN INTERSECTE LA FACE NF1
C                             ENTRE CES 3 SOMMETS ET ENTRE XYZBAR-XYZBA2
                     IF( LINTER .EQ. 1 ) GOTO 16
                  ENDIF
                  NF1 = NFETOI(5,NF1)
                  GOTO 14
               ENDIF

C              XYZBAR-XYZBA2 INTERSECTE AUCUNE AUTRE FACE
C              XYZ LE POINT MILIEU DU SEGMENT XYZBAR-XYZ
               DO K=1,3
                  XYZ(K) = ( XYZBAR(K) + XYZBA2(K) ) * 0.5D0
               ENDDO

C              SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C              RECHERCHE DES FACES CONTINUES DE L'ETOILE FORMANT AVEC XYZ
C              UN TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C              AUTRES FACES DE L'ETOILE
               CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
     %                        MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
     %                        MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
     %                        MXASFVP,  NBASFVP, NSASFVP,
     %                        VETXYZ,   QUAMIN, IERR )
C              NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C              NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0
C              MIARSICI: NB MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
C              QUAMIN : QUALITE MINIMALE MAXIMISEE DES QUALITES DES TETRAEDRES NF-XYZ

C              RECHERCHE DU POINT DONNANT LE MINIMUM D'ARETES SIMPLES DES CIRCUITS
               IF( NBVPSI .GE. 3 .AND. NBVPSI .GE. NBFETO/3 .AND.
     %           ( MIARSICI.LT.MIARSIMI .OR.
     %           ( MIARSICI.EQ.MIARSIMI .AND. QUAMIN.GT.QUAMIMX) )) THEN
                  NOCHOIXYZ = 3
                  NBVPSIMAX = NBVPSI
                  XYZMAX(1) = XYZ(1)
                  XYZMAX(2) = XYZ(2)
                  XYZMAX(3) = XYZ(3)
                  NFVPSIMAX = NFVPSI(1)
                  MIARSIMI  = MIARSICI
                  QUAMIMX   = QUAMIN
               ENDIF

            ENDIF

C           PASSAGE A LA FACE SUIVANTE
 16         NF = NFETOI(5,NF)
            GOTO 11
         ENDIF

      ENDIF

      IF( NBVPSIMAX .GT. 0 ) THEN
      KTITRE='sipteto ap3:        FACES+XYZ TETRAEDRES V>0 pour        F
     %ACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
      WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(81:88), '(F8.5)' ) QUAMIN
      WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)
      ENDIF

      IF( NBVPSIMAX .GE. MAXFAC ) THEN
         GOTO 8000
      ENDIF


C     CHOIX 3': XYZ=POINT SUR LA NORMALE D'UNE FACE DE L'ETOILE
C               A MI-DISTANCE ENTRE LE BARYCENTRE DE LA FACE
C               ET LE POINT D'INTERSECTION DE LA NORMALE
C               AVEC LA PLUS PROCHE FACE OPPOSEE
C     --------------------------------------------------------
      NF = N1FEOC
 19   IF( NF .GT. 0 ) THEN

C        CALCUL DE LA LONGUEUR MOYENNE D'UN COTE, DU BARYCENTRE,
C        DE LA NORMALE DE LA FACE NF DE L'ETOILE
C        VN=S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE L'ETOILE
         CALL PTNORTR( PTXYZD, NFETOI(2,NF), COTMOY, XYZBAR, VN )

C        LA NORMALE VN POINTE T ELLE VERS LE BARYCENTRE DE L'ETOILE?
         DO K=1,3
            XYZN(K) = XYZBE(K) - XYZBAR(K)
         ENDDO
         IF( PROSCD( VN, XYZN, 3 ) .LT. 0D0 ) THEN
C           LA NORMALE A LA FACE NE VA PAS VERS LE BARYCENTRE DE L'ETOILE
C           LA FACE EST ABANDONNEE
            GOTO 40
         ENDIF

C        CALCUL DU POINT SUR LA NORMALE AU BARYCENTRE DE LA FACE NF DE L'ETOILE
         DO K=1,3
            XYZN(K) = XYZBAR(K) + 10D0 * COTMOY * VN(K)
         ENDDO

C        POINT D'INTERSECTION DE XYZBAR+NORMALE AVEC UNE AUTRE FACE
c        DE L'ETOILE
         DMIN  = 1D+111
         NFMIN = 0
         NF1  = N1FEOC
 20      IF( NF1 .GT. 0 ) THEN
            IF( NF1 .NE. NF ) THEN
               CALL INARTR( XYZBAR, XYZN,
     %                      PTXYZD(1,NFETOI(2,NF1)),
     %                      PTXYZD(1,NFETOI(3,NF1)),
     %                      PTXYZD(1,NFETOI(4,NF1)),
     %                      LINTER, XYZ, CBTR )
C              LINTER=1 SI XYZBAR-XYZN INTERSECTE LA FACE NF1
C                       ENTRE CES 3 SOMMETS ET ENTRE XYZBAR-XYZN
               IF( LINTER .EQ. 1 ) THEN
C                 DISTANCE ENTRE LE BARYCENTRE ET LE POINT D'INTERSECTION
                  D = ( XYZ(1) - XYZBAR(1) ) **2
     %              + ( XYZ(2) - XYZBAR(2) ) **2
     %              + ( XYZ(3) - XYZBAR(3) ) **2
                  IF( D .LT. COTMOY*1D-3 ) THEN
C                    XYZ=BARYCENTRE DE LA FACE NF EST SUR LA FACE NF1
C                    ABANDON DE LA FACE NF
                     GOTO 40
                  ELSE
                     IF( D .LT. DMIN ) THEN
                        DMIN  = D
                        NFMIN = NF1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            NF1 = NFETOI(5,NF1)
            GOTO 20
         ENDIF

C        TOUTES LES FACES AVEC INTERSECTION ONT ETE VUES
         IF( DMIN .NE. 1D+111 ) THEN

C           XYZ LE POINT D'INTERSECTION AVEC NFMIN EST RECALCULE
            CALL INARTR( XYZBAR, XYZN,
     %                   PTXYZD(1,NFETOI(2,NFMIN)),
     %                   PTXYZD(1,NFETOI(3,NFMIN)),
     %                   PTXYZD(1,NFETOI(4,NFMIN)),
     %                   LINTER, XYZ, CBTR )
C           LE POINT MILIEU DU SEGMENT XYZBAR-XYZ
            DO K=1,3
               XYZ(K) = ( XYZ(K) + XYZBAR(K) ) * 0.5D0
            ENDDO

C           SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C           RECHERCHE DES FACES CONTINUES DE L'ETOILE FORMANT AVEC XYZ
C           UN TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C           AUTRES FACES DE L'ETOILE
            CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
     %                     MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
     %                     MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
     %                     MXASFVP,  NBASFVP, NSASFVP,
     %                     VETXYZ,   QUAMIN, IERR )
C           NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C           NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0
C           MIARSICI: NB MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
C           QUAMIN : QUALITE MINIMALE MAXIMISEE DES QUALITES DES TETRAEDRES NF-XYZ

C           RECHERCHE DU POINT DONNANT LE MINIMUM D'ARETES SIMPLES DES CIRCUITS
            IF( NBVPSI .GE. 3 .AND. NBVPSI .GE. NBFETO/3 .AND.
     %        ( MIARSICI.LT.MIARSIMI .OR.
     %        ( MIARSICI.EQ.MIARSIMI .AND. QUAMIN.GE.QUAMIMX) ))THEN
               NOCHOIXYZ = 3
               NBVPSIMAX = NBVPSI
               XYZMAX(1) = XYZ(1)
               XYZMAX(2) = XYZ(2)
               XYZMAX(3) = XYZ(3)
               NFVPSIMAX = NFVPSI(1)
               MIARSIMI  = MIARSICI
               QUAMIMX   = QUAMIN
            ENDIF

         ENDIF

C        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
 40      NF = NFETOI(5,NF)
         GOTO 19

      ENDIF

      IF( NBVPSIMAX .GT. 0 ) THEN
      KTITRE='sipteto ap4:        FACES+XYZ TETRAEDRES V>0 pour        F
     %ACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
      WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(81:88), '(F8.5)' ) QUAMIN
      WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)
      ENDIF

      IF( NBVPSIMAX .GE. MAXFAC ) THEN
         GOTO 8000
      ENDIF

C     CHOIX 4: XYZ=POINT SUR LA NORMALE D'UNE FACE DE L'ETOILE
C              PLACEE AU MILIEU D'UN DE SES 3 COTES
C              A MI-DISTANCE ENTRE CE MILIEU
C              ET LE POINT D'INTERSECTION AVEC LA NORMALE
C              DE LA PLUS PROCHE FACE OPPOSEE
C     --------------------------------------------------------
      NF = N1FEOC
 50   IF( NF .GT. 0 ) THEN

C        CALCUL DE LA LONGUEUR MOYENNE D'UN COTE, DU BARYCENTRE,
C        DE LA NORMALE DE LA FACE NF DE L'ETOILE
C        VN=S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE L'ETOILE
         CALL PTNORTR( PTXYZD, NFETOI(2,NF), COTMOY, XYZBAR, VN )

cccC        LA NORMALE VN POINTE T ELLE VERS LE BARYCENTRE DE L'ETOILE?
ccc         DO K=1,3
ccc            XYZN(K) = XYZBE(K) - XYZBAR(K)
ccc         ENDDO
ccc         IF( PROSCD( VN, XYZN, 3 ) .LT. 0D0 ) THEN
cccC           LA NORMALE A LA FACE NE VA PAS VERS LE BARYCENTRE DE L'ETOILE
cccC           LA FACE EST ABANDONNEE
ccc            GOTO 80
ccc         ENDIF

C        CALCUL DU POINT SUR LA NORMALE AU MILIEU DE L'ARETE M
C        DE LA FACE NF DE L'ETOILE
         NS1 = NFETOI(4,NF)
         DO M=1,3
            NS2 = NFETOI(1+M,NF)
C           LES 3 COORDONNEES DU MILEU DE L'ARETE M DE NF
            DO K=1,3
               XYZMIL(K) = (PTXYZD(K,NS1) + PTXYZD(K,NS2)) * 0.5D0
               XYZN(K) = XYZMIL(K) + 2D0 * COTMOY * VN(K)
            ENDDO

C           POINT D'INTERSECTION DE MILIEU+NORMALE AVEC UNE AUTRE FACE
C           DE L'ETOILE
            DMIN = 1D+111
            NFMIN = 0
            NF1  = N1FEOC
 60         IF( NF1 .GT. 0 ) THEN
               IF( NF1 .NE. NF ) THEN
                  CALL INARTR( XYZMIL, XYZN,
     %                         PTXYZD(1,NFETOI(2,NF1)),
     %                         PTXYZD(1,NFETOI(3,NF1)),
     %                         PTXYZD(1,NFETOI(4,NF1)),
     %                         LINTER, XYZ, CBTR )
C                 LINTER=1 SI XYZMIL-XYZN INTERSECTE LA FACE NF1
C                          ENTRE CES 3 SOMMETS ET ENTRE XYZMIL-XYZN
                  IF( LINTER .EQ. 1 ) THEN
C                    DISTANCE ENTRE LE MILIEU ET LE POINT D'INTERSECTION
                     D = ( XYZ(1) - XYZMIL(1) ) **2
     %                 + ( XYZ(2) - XYZMIL(2) ) **2
     %                 + ( XYZ(3) - XYZMIL(3) ) **2
                     IF( D .LT. COTMOY*1D-3 ) THEN
C                       XYZ=MILIEU DE LA FACE NF EST SUR LA FACE NF1
C                       ABANDON DE LA FACE NF
                        GOTO 80
                     ELSE
                        IF( D .LT. DMIN ) THEN
                           DMIN  = D
                           NFMIN = NF1
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               NF1 = NFETOI(5,NF1)
               GOTO 60
            ENDIF

C           TOUTES LES FACES AVEC INTERSECTION ONT ETE VUES
            IF( DMIN .NE. 1D+111 ) THEN

C              XYZ LE POINT D'INTERSECTION AVEC NFMIN EST RECALCULE
               CALL INARTR( XYZMIL, XYZN,
     %                      PTXYZD(1,NFETOI(2,NFMIN)),
     %                      PTXYZD(1,NFETOI(3,NFMIN)),
     %                      PTXYZD(1,NFETOI(4,NFMIN)),
     %                      LINTER, XYZ, CBTR )
C              LE POINT MILIEU DU SEGMENT XYZMIL-XYZ
               DO K=1,3
                  XYZ(K) = ( XYZMIL(K) + XYZ(K) ) * 0.5D0
               ENDDO

C              SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C              RECHERCHE DES FACES CONTINUES DE L'ETOILE FORMANT AVEC XYZ
C              UN TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C              AUTRES FACES DE L'ETOILE
               CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
     %                        MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
     %                        MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
     %                        MXASFVP,  NBASFVP, NSASFVP,
     %                        VETXYZ,   QUAMIN, IERR )
C              NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C              NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0
C              MIARSICI: NB MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI
C              QUAMIN : QUALITE MINIMALE MAXIMISEE DES QUALITES DES TETRAEDRES NF-XYZ

C              RECHERCHE DU POINT DONNANT LE MINIMUM D'ARETES SIMPLES DES CIRCUITS
               IF( NBVPSI .GE. 3 .AND. NBVPSI .GE. NBFETO/3 .AND.
     %           ( MIARSICI.LT.MIARSIMI .OR.
     %           ( MIARSICI.EQ.MIARSIMI .AND. QUAMIN.GE.QUAMIMX) ))THEN
                  NOCHOIXYZ = 4
                  NBVPSIMAX = NBVPSI
                  XYZMAX(1) = XYZ(1)
                  XYZMAX(2) = XYZ(2)
                  XYZMAX(3) = XYZ(3)
                  NFVPSIMAX = NFVPSI(1)
                  MIARSIMI  = MIARSICI
                  QUAMIMX   = QUAMIN
               ENDIF

            ENDIF

C           PASSAGE A L'ARETE SUIVANTE DE NF
            NS1 = NS2
         ENDDO

C        PASSAGE A L'ARETE SUIVANTE DE NF
 80      NF = NFETOI(5,NF)
         GOTO 50

      ENDIF

      IF( NBVPSIMAX .GT. 0 ) THEN
      KTITRE='sipteto ap5:        FACES+XYZ TETRAEDRES V>0 pour        F
     %ACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
      WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(81:88), '(F8.5)' ) QUAMIN
      WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)
      ENDIF

      IF( NBVPSIMAX .GE. MAXFAC ) THEN
         GOTO 8000
      ENDIF

C     CHOIX 5: XYZ=POINT SUR LA NORMALE AU BARYCENTRE DES FACES DE L'ETOILE
C     ---------------------------------------------------------------------
C     CALCUL DU POINT SUR LA NORMALE A UNE FACE DE L'ETOILE
C     EN SON BARYCENTRE ET A UNE DISTANCE = PERIMETRE * REDUCT
C     REDUCTION = 1/2 DE LA LONGUEUR DU COTE MOYEN
      REDUCT = 0.5D0 * COEFQP
      NBVPSI = 0
      NBPAS  = 0
      C      = COTMOY * COEFQP

 100  NF = N1FEOC
 120  IF( NF .GT. 0 ) THEN

C        CALCUL DE LA LONGUEUR MOYENNE D'UN COTE, DU BARYCENTRE,
C        DE LA NORMALE DE LA FACE NF DE L'ETOILE
C        VN=S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE L'ETOILE
         CALL PTNORTR( PTXYZD, NFETOI(2,NF), COTMOY, XYZBAR, VN )

C        VN POINTE T ELLE VERS LE BARYCENTRE DE L'ETOILE?
         DO K=1,3
            XYZN(K) = XYZBE(K) - XYZBAR(K)
         ENDDO
         
         IF( PROSCD( VN, XYZN, 3 ) .LT. 0D0 ) THEN
C           LA NORMALE A LA FACE NE VA PAS VERS LE BARYCENTRE DE L'ETOILE
C           LA FACE EST ABANDONNEE
            GOTO 130
         ENDIF

C        CALCUL DU POINT SUR LA NORMALE A LA FACE NF DE L'ETOILE
         DO K=1,3
            XYZ(K) = XYZBAR(K) + REDUCT * C * VN(K)
         ENDDO

C        SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C        RECHERCHE DES FACES CONTINUES DE L'ETOILE FORMANT AVEC XYZ
C        UN TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C        AUTRES FACES DE L'ETOILE
         CALL FCPTVPSI( XYZ,    PTXYZD, N1FEOC, NFETOI, QUAMINEX,
     %                  MXVPSI, NBVPSI, NFVPSI, NFVPSI0,
     %                  MIARSICI, MXCIAS,  NBCIAS, N1CIAS,
     %                  MXASFVP,  NBASFVP, NSASFVP,
     %                  VETXYZ,   QUAMIN,  IERR )
C        NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C        NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0
C        MIARSICI: NB MINIMAL D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI

C        RECHERCHE DU POINT DONNANT LE MINIMUM D'ARETES SIMPLES DES CIRCUITS
         IF( NBVPSI .GE. 3 .AND. NBVPSI .GE. NBFETO/3 .AND.
     %     ( MIARSICI.LT.MIARSIMI .OR.
     %     ( MIARSICI.EQ.MIARSIMI .AND. QUAMIN.GE.QUAMIMX) ))THEN
            NOCHOIXYZ = 4
            NBVPSIMAX = NBVPSI
            XYZMAX(1) = XYZ(1)
            XYZMAX(2) = XYZ(2)
            XYZMAX(3) = XYZ(3)
            NFVPSIMAX = NFVPSI(1)
            MIARSIMI  = MIARSICI
            QUAMIMX   = QUAMIN
         ENDIF

C        PASSAGE A LA FACE SUIVANTE
 130     NF = NFETOI(5,NF)
         GOTO 120

      ENDIF

      IF( NBVPSIMAX .GT. 0 ) THEN
      KTITRE='sipteto ap6:        FACES+XYZ TETRAEDRES V>0 pour        F
     %ACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
      WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
      WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
      WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(81:88), '(F8.5)' ) QUAMIN
      WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      LORBITE = 1
      CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI)
      ENDIF

      IF( NBVPSIMAX .GE. MAXFAC ) THEN
         GOTO 8000
      ENDIF

      IF( NBPAS .LT. 4 ) THEN
C         ESSAI AVEC UN POINT SUR LA NORMALE PLUS PROCHE DE LA FACE
          NBPAS  = NBPAS + 1
          REDUCT = REDUCT * 0.5D0
          GOTO 100
      ENDIF


C     ==========================================
C     BILAN FINAL DE LA FACE OPTIMISEE NFVPSIMAX
C     ==========================================
 8000 IF( NBVPSIMAX .GT. 0 ) THEN

C        ICI, IL EXISTE AU MOINS UNE FACE AVEC XYZMAX
C        FORMANT UN TETRAEDRE V>0 ET SANS INTERSECTION AVEC L'ETOILE
C        LES FACES OPTIMALES POUR LE MAXIMUM DU MINIMUM DES QUALITES
C        -----------------------------------------------------------
 
C        DETERMINER L'ENSEMBLE DES NBVPSI0 FACES NFVPSI0 DE L'ETOILE FORMANT
C        AVEC LE POINT XYZ DES TETRAEDRES DE VOLUMES>0
C        SANS INTERSECTION AVEC LES FACES DE L'ETOILE SANS RECHERCHE QUALITE
         CALL F0PTVPSI( QUAMINEX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %                  MXVPSI, NBVPSI0, NFVPSI0, NBFETO, VETXYZ, QMIN )

ccc         print *,'fin sipteto: avant f1ptvpsi nochoixyz=',nochoixyz,
ccc     %           ' NBVPSIMAX=',NBVPSIMAX,' NFVPSIMAX=',NFVPSIMAX,
ccc     %           ' QMIN=',QMIN,' XYZMAX=',XYZMAX

         IF( NBVPSI0 .GT. 0 ) THEN

C           RECUPERATION DES NBVPSIMAX FACES CORRECTES
C           CAR LE CAS OPTIMAL A PU ETRE EFFACE PAR UN CAS MOINS FAVORABLE
            CALL F1PTVPSI( QUAMINEX,  XYZMAX,  PTXYZD,  N1FEOC, NFETOI,
     %                     NFVPSIMAX, NBVPSI0, NFVPSI0,
     %                     MXVPSI,    NBVPSI,  NFVPSI, VETXYZ, QUAMIMX )
            IF( NBVPSI .LE. 0 ) GOTO 8030

C           TENTATIVE D'AMELIORATION
C           LE BARYCENTRE DE CES SEULES FACES A V>0 CONTIGUES AMELIORE T IL
C           LE MAXIMUM DU MINIMUM DES QUALITES DES TETRAEDRES NFVPSIMAX?
C           ---------------------------------------------------------------
            CALL XYZBF11( PTXYZD, NFETOI, NBVPSI, NFVPSI,  XYZ )

            QMIN = 2.0
            DO NF = 1, NBVPSI

               NF1 = NFVPSI( NF )
C              QUALITE DU TETRAEDRE NF1-XYZ
               CALL QUATETD( PTXYZD(1,NFETOI(2,NF1) ),
     %                       PTXYZD(1,NFETOI(3,NF1) ),
     %                       PTXYZD(1,NFETOI(4,NF1) ),
     %                       XYZ,
     %                       ARMIN, ARMAX, SURFTR, V, QUALIT )
               IF( QUALIT .LT. QMIN ) QMIN = QUALIT

C              LA QUALITE DU TETRAEDRE NF1-XYZ DOIT ETRE >QUAMIMX
               IF( QUALIT .LT. QUAMIMX ) THEN
C                 XYZMAX DONNE LA MEILLEURE QUALITE
                  GOTO 9000
               ENDIF

            ENDDO

C           LE POINT XYZ AMELIORE LE MINIMUM DES QUALITES DES TETRAEDRES NFVPSI
            QUAMIMX   = QMIN
            XYZMAX(1) = XYZ(1)
            XYZMAX(2) = XYZ(2)
            XYZMAX(3) = XYZ(3)
            CALL F1PTVPSI( QUAMINEX,  XYZMAX,  PTXYZD,  N1FEOC, NFETOI,
     %                     NFVPSIMAX, NBVPSI0, NFVPSI0,
     %                     MXVPSI,    NBVPSI,  NFVPSI, VETXYZ, QUAMIMX )

C           RECHERCHE DES CIRCUITS D'ARETES SIMPLES DES NBVPSI FACES
C           NFVPSI DE NFETOI
            CALL CIASFAET( NBVPSI,  NFVPSI,  NFETOI,
     %                     NBCIAS,  MXCIAS,  N1CIAS,
     %                     MXASFVP, MIARSICI, NBASFVP, NSASFVP, IERR )

            PRINT*,'sipteto: QUAMIMX=',QUAMIMX,' AU BARYCENTRE DES ',
     %              NBVPSI,' FACES+XYZMAX V>0  MIARSICI=',MIARSICI
            GOTO 9100

         ENDIF

      ENDIF


C     DERNIERS ESSAIS: DEPLACEMENTS DU BARYCENTRE DES FACES DE L'ETOILE
C     -----------------------------------------------------------------
 8030 IF( NBVPSI .LE. 0 ) THEN
C        ABANDON
         GOTO 8900
      ENDIF

C     REDEPART DE ZERO
      NBPAS = 0

 8040 NBPAS     = NBPAS + 1
      NBVPSIMAX = 0
      NFVPSIMAX = 0
      QUAMIMX   = 0.0

C     COORDONNEES DU BARYCENTRE DE L'ETOILE QUASI PLANE
C     => GENERATION DE NBVPSI TETRAEDRES PLANS!...
ccc      CALL XYZBF12( NBPAS,  PTXYZD, N1FEOC, 0, NFETOI,
ccc     %              NBFETO, XYZMIX, XYZBE )

      C = COTMOY * COEFQP
      DO K=1,3

         NOCHOIXYZ = 0
         D = C / 2**K

 8050    NOCHOIXYZ = NOCHOIXYZ + 1
         GOTO( 8100, 8200, 8300, 8400, 8500, 8600, 8700 ), NOCHOIXYZ

C        DEPLACEMENT LE LONG DE L'AXE DES X>0 A PARTIR DU BARYCENTRE
 8100    XYZ(1) = XYZBE(1) + D
         XYZ(2) = XYZBE(2)
         XYZ(3) = XYZBE(3)
         GOTO 8650

C        DEPLACEMENT LE LONG DE L'AXE DES X<0 A PARTIR DU BARYCENTRE
 8200    XYZ(1) = XYZBE(1) - D
         XYZ(2) = XYZBE(2)
         XYZ(3) = XYZBE(3)
         GOTO 8650

C        DEPLACEMENT LE LONG DE L'AXE DES Y>0 A PARTIR DU BARYCENTRE
 8300    XYZ(1) = XYZBE(1)
         XYZ(2) = XYZBE(2) + D
         XYZ(3) = XYZBE(3)
         GOTO 8650

C        DEPLACEMENT LE LONG DE L'AXE DES Y<0 A PARTIR DU BARYCENTRE
 8400    XYZ(1) = XYZBE(1)
         XYZ(2) = XYZBE(2) - D
         XYZ(3) = XYZBE(3)
         GOTO 8650

C        DEPLACEMENT LE LONG DE L'AXE DES Z>0 A PARTIR DU BARYCENTRE
 8500    XYZ(1) = XYZBE(1)
         XYZ(2) = XYZBE(2)
         XYZ(3) = XYZBE(3) + D
         GOTO 8650

C        DEPLACEMENT LE LONG DE L'AXE DES Z<0 A PARTIR DU BARYCENTRE
 8600    XYZ(1) = XYZBE(1)
         XYZ(2) = XYZBE(2)
         XYZ(3) = XYZBE(3) - D
         GOTO 8650

C        SIMULATION DE LA TETRAEDRISATION DES FACES DE L'ETOILE
C        RECHERCHE DES FACES DE L'ETOILE FORMANT AVEC XYZ UN
C        TETRAEDRE DE VOLUME>0 ET SANS INTERSECTION AVEC LES
C        FACES DE L'ETOILE
 8650    CALL F0PTVPSI( QUAMINEX, XYZ,  PTXYZD, N1FEOC, NFETOI,
     %                  MXVPSI, NBVPSI, NFVPSI, NBFETO, VETXYZ, QMIN )
C        NBVPSI : NOMBRE DE TETRAEDRES DE VOLUME>=0 SANS INTERSECTION
C        NFVPSI : NUMERO NFETOI DES FACES DONNANT UN VOLUME>=0

ccc         KTITRE='      FACES+XYZ de V>0 RACCRO FINAL NOCHOIXYZ=      '
ccc         WRITE(KTITRE(46:47), '(I2)' ) NOCHOIXYZ+6
ccc         CALL SANSDBL( KTITRE, L )
ccc         CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD,
ccc     %                 N1FEOC, NFETOI, NBVPSI, NFVPSI )

         IF( NBVPSI .GT. NBFETO/3 ) THEN
            NOCHOIXYZ = NOCHOIXYZ + 6
            XYZMAX(1) = XYZ(1)
            XYZMAX(2) = XYZ(2)
            XYZMAX(3) = XYZ(3)
            GOTO 9000
         ELSE
            GOTO 8050
         ENDIF

      ENDDO

C     PROBLEME PAS DE TETRAEDRE + XYZ DE VOLUME POSITIF ET SANS INTERSECTION
C     ----------------------------------------------------------------------
 8700 IF( NBPAS .EQ. 1 ) THEN
C        CHANGEMENT DE POIDS DANS L'APPEL DE XYZBF12
C        PASSAGE AU BARYCENTRE avec POIDS des 1D0/SURFACE des FACES
         GOTO 8040
      ENDIF

C     ESSAI DE MODIFIER LES XYZ D'UN SOMMET NON CF DE L'ETOILE
      PRINT*,'sipteto : ESSAI DE MODIFIER XYZ D''UN SOMMET NON CF DE L''
     %ETOILE de',NBFETO,' FACES'

C     NOMBRE D'ESSAIS DE DEPLACEMENT D'UN SOMMET NON FRONTALIER
      NBNST = NBNST + 1
      IF( NBNST .LE . 3 ) THEN

         NF = N1FEOC

 8710    IF( NF .GT. 0 ) THEN

            DO 8720 K=1,3

               NST = NFETOI(1+K,NF)
               IF( NPSOFR(NST).LT.-4  .OR.  NPSOFR(NST).GT.0 ) GOTO 8720

C              NST N'EST PAS SUR LA FRONTIERE
               IF( NST .EQ. NST0 .OR. NST .EQ. NST1 ) GOTO 8720

C              VN=S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE L'ETOILE
               CALL PTNORTR( PTXYZD, NFETOI(2,NF), COTMOY, XYZBAR, VN )

               DO L=1,3
                  PTXYZD(L,NST) = PTXYZD(L,NST) + 5D-2 * COTMOY * VN(L)
               ENDDO

               PRINT*,'sipteto : MODIF de PTXYZD(',NST,')=',
     %                (PTXYZD(L,NST),L=1,4),' COTMOY=',COTMOY

C              RETOUR AU DEBUT AVEC CE DEPLACEMENT DU SOMMET NST
               NST0 = NST1
               NST1 = NST
               GOTO 1

 8720       ENDDO

C           PASSAGE A LA FACE SUIVANTE DE L'ETOILE
            NF = NFETOI(5,NF)
            GOTO 8710

         ENDIF

      ENDIF


C     PROBLEME PAS DE TETRAEDRE + XYZ DE VOLUME POSITIF ET SANS INTERSECTION
C     ----------------------------------------------------------------------
 8900 NOCHOIXYZ = -1
      NBVPSI    = 0
      NBVPSIMAX = 0
      XYZMAX(1) = 1D10
      XYZMAX(2) = 1D10
      XYZMAX(3) = 1D10
      PRINT*,'sipteto : AUCUN TETRAEDRE XYZMAX+FACE SANS INTERSECTION av
     %ec les',NBFETO,' FACES N''A UN VOLUME>0'


C     ==========================
C     BILAN de SORTIE de sipteto
C     ==========================
ccc 9000 PRINT *,'sipteto: NOAPPEL=',NOAPPEL,' NOCHOIXYZ=',NOCHOIXYZ,
ccc     %        ' XYZ=',XYZ,
ccc     %        ' =>',NBVPSI,' TETRAEDRES V>0 de QUALITE MIN=',QUAMIMX,
ccc     %        ' pour',NBFETO,' FACES de l''ETOILE'
ccc      PRINT *,'sipteto: les',NBVPSI,' FACES NFVPSI=',
ccc     %        (NFVPSI(K),K=1,NBVPSI)
ccc      print *,'Fin sipteto: Xmin=',XYZMIX(1,1),' Pt X=',XYZMAX(1),
ccc     %                    ' XMax=',XYZMIX(1,2)
ccc      print *,'Fin sipteto: Ymin=',XYZMIX(2,1),' Pt Y=',XYZMAX(2),
ccc     %                    ' YMax=',XYZMIX(2,2)
ccc      print *,'Fin sipteto: Zmin=',XYZMIX(3,1),' Pt Z=',XYZMAX(3),
ccc     %                    ' ZMax=',XYZMIX(3,2)


C     RECALCUL FINAL DES NBVPSI FACES DU CAS MAX OBTENU AVANT
C     =======================================================
 9000 IF( NBVPSIMAX .GT. 0 ) THEN

ccc         print *,'fin sipteto: avant f1ptvpsi nochoixyz=',nochoixyz,
ccc     %           ' NBVPSIMAX=',NBVPSIMAX,' NFVPSIMAX=',NFVPSIMAX,
ccc     %           ' QUAMIMX=',QUAMIMX,' XYZMAX=',XYZMAX
 
C        DETERMINER L'ENSEMBLE DES NBVPSI0 FACES NFVPSI0 DE L'ETOILE FORMANT
C        AVEC LE POINT XYZ DES TETRAEDRES DE VOLUMES>0
C        SANS INTERSECTION AVEC LES FACES DE L'ETOILE SANS RECHERCHE QUALITE
         CALL F0PTVPSI( QUAMINEX, XYZMAX, PTXYZD, N1FEOC, NFETOI,
     %                  MXVPSI, NBVPSI0, NFVPSI0, NBFETO, VETXYZ, QMIN)

C        RECUPERATION DES NBVPSIMAX FACES CORRECTES
C        CAR LE CAS OPTIMAL A PU ETRE EFFACE PAR UN CAS MOINS FAVORABLE
         CALL F1PTVPSI( QUAMINEX,  XYZMAX,  PTXYZD,  N1FEOC, NFETOI,
     %                  NFVPSIMAX, NBVPSI0, NFVPSI0,
     %                  MXVPSI,    NBVPSI,  NFVPSI, VETXYZ, QUAMIMX )

C        RECHERCHE DES CIRCUITS D'ARETES SIMPLES DES NBVPSI FACES
C        NFVPSI DE NFETOI
         CALL CIASFAET( NBVPSI,  NFVPSI,  NFETOI,
     %                  NBCIAS,  MXCIAS,  N1CIAS,
     %                  MXASFVP, MIARSICI, NBASFVP, NSASFVP, IERR )

      ENDIF

C     TRACE ORANGE DES FACES+XYZMAX DE VOLUME POSITIF SANS INTERSECTION
C     ET MAXIMISANT LE MINIMUM DES QUALITES DES TETRAEDRES FORMES
 9100 IF( MIARSIMI .NE. 2**30 ) THEN
      KTITRE='sipteto : MinArSiMi=          FACES+XYZMAX TETRAEDRES V>0 
     %pour        FACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
         WRITE(KTITRE(21:23), '(I3)' ) MIARSIMI
         WRITE(KTITRE(25:29), '(I5)' ) NBVPSI
      ELSE
      KTITRE='sipteto :            PAS  DE  FACES+XYZMAX TETRAEDRES V>0 
     %pour        FACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
      ENDIF

      WRITE(KTITRE(63:68), '(I5)'    ) NBFETO
      WRITE(KTITRE(87:88), '(I2)'    ) NOCHOIXYZ
      WRITE(KTITRE(96:103),'(F8.5)'  ) QUAMIMX
      WRITE(KTITRE(114:126),'(G13.5)') VETXYZ

 9900 IF( NOAPPEL .GT. 0 ) THEN
         CALL SANSDBL( KTITRE, L )
         KTITRE = KTITRE(1:L) // ' NoAppel=   '
         WRITE(KTITRE(L+10:L+11), '(I2)' ) NOAPPEL
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO3( KTITRE(1:L), XYZMAX, XYZMAX, PTXYZD,
     %                 N1FEOC, NFETOI, NBVPSI, NFVPSI )
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
