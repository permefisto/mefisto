      SUBROUTINE SISTST( NSNOFR, NSDLFR, MXSOMM, PTXYZD, NPSOFR,
     %                   MXTETR, N1TETS, NOTETR, IVOLTE, NVOLTE,
     %                   MXFETO, N1FEOC, NFETOI,
     %                   MXTENS, NBTENS, NOTENS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SIMULER DANS LA TETRAEDRISATION L'IDENTIFICATION DU SOMMET
C -----    NSNOFR NON FRONTALIER PAR LE SOMMET NSDLFR SOMMET FRONTALIER
C          TOUS LES 2 SOMMETS DE TETRAEDRES DU MAILLAGE ACTUEL
C          EN FIN IERR NON 0 INTERDIT L'IDENTIFICATION
C          L'ETOILE NFETOI EST CELLE DES TETRAEDRES DES 2 SOMMETS
C          NSNOFR ET NSDLFR

C ENTREES:
C --------
C NSNOFR : NUMERO DU SOMMET PTXYZD A IDENTIFIER AU SOMMET NSDLFR
C NSDLFR : NUMERO DU SOMMET PTXYZD
C MXSOMM : NOMBRE DE SOMMETS DECLARABLES DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C MXTENS : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTENS
C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI
C IVOLTE : =0 TABLEAU NVOLTE NON PRESENT
C          =1 TABLEAU NVOLTE     PRESENT

C MODIFIES:
C ---------
C NVOLTE : NUMERO DE VOLUME (1 a NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU (Exemple: les TETRAEDRES VIDES DE NOTETR)
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C NFETOI : AU DEBUT VERSION1  FACES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR AYANT CETTE FACE
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3: NON UTILISE ICI
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C          ENSUITE VERSION2
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST DIRIGE VERS L'INTERIEUR DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C SORTIES :
C ---------
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NBTENS : NOMBRE DE  TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C          ou LE SOMMET NSNOFR ou OPPOSE AUX FACES DES TETRAEDRES
C          DE SOMMET NSNOFR
C NOTENS : NUMERO DES TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C          ou LE SOMMET NSNOFR ou OPPOSE AUX FACES DES TETRAEDRES
C          DE SOMMET NSNOFR
C IERR   : =0 PAS D'ERREUR RENCONTREE et l'IDENTIFICATION EST ACCEPTEE
C          >0 SINON et l'IDENTIFICATION EST REFUSEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2016
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE

      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           NPSOFR(MXSOMM), NOTETR(8,MXTETR),
     %                  NOTENS(MXTENS), NFETOI(5,MXFETO),
     %                  N1TETS(MXSOMM), NVOLTE(MXTETR)

      DOUBLE PRECISION  VOLUTE0, VOLUTE1, VOLUE0, VOLUE1, VTE,COBARY(4),
     %                  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NOSOTE(4)

C     SIMULATION de L'IDENTIFICATION DU SOMMET NSNOFR PAR NSDLFR
ccc      PRINT*,'sistst: Simulation d''IDENTIFICATION du sommet NON FRONTIE
ccc     %RE', NSNOFR,' au sommet SUR LA FRONTIERE',NSDLFR

      IERR = 0

      IF( NPSOFR(NSNOFR) .GT. 0 ) THEN
         IERR = 5
         GOTO 9999
      ENDIF

C     RECHERCHE DES NBTENS TETRAEDRES DE SOMMET NSNOFR
      CALL TETR1S( NSNOFR, N1TETS, NOTETR,
     %             NBTENS, MXTENS, NOTENS, IERR )
C     NBTENS: EN ENTREE NOMBRE DE TETRAEDRES DEJA RANGES DANS NOTENS
C             EN SORTIE AJOUT DES TETRAEDRES AYANT NSNOFR COMME SOMMET
C             = 0 SI SATURATION DU TABLEAU NOTENS
C             =-1 SI UN SOMMET NSNOFR N'EST PAS UN SOMMET DE N1TETS(NS)
C             > 0 SINON
      IF( IERR .NE. 0 ) THEN
         PRINT*,'sistst: sortie tetr1s ProblemE sur les TETRAEDRES de SO
     %MMET NSNOFR=',NSNOFR,' IERR=',IERR
         TRACTE = .TRUE.
         IERR = 12
         GOTO 9999
      ENDIF

      IF( NBTENS .GT. 64 ) THEN
         PRINT*,'sistst: Sommet',NSNOFR,' dans',NBTENS,
     %' tetraedres. BEAUCOUP TROP pour une SUBSTITUTION avec le sommet',
     %NSDLFR
         IERR = 15
         GOTO 9999
       ENDIF

C     SI LE SOMMET NSNOFR DEVENANT NSDLFR SORT DE SON ETOILE
C     L'IDENTIFICATION EST REFUSEE
      DO NT = 1, NBTENS

C        TETRAEDRE DE SOMMET NSNOFR
         NTE = NOTENS( NT )

C        CALCUL DES COORDONNEES BARYCENTRIQUES DE PTXYZD(NSDLFR)
         CALL PTDSTE( PTXYZD(1,NSDLFR), PTXYZD, NOTETR(1,NTE),
     %                NONOUI, VTE, COBARY )
C        NONOUI = 0  PT N'EST PAS DANS LE TETRAEDRE NTE
C               = 1  PT   EST     DANS LE TETRAEDRE NTE
C               =-1  LE TETRAEDRE NTE EST DEGENERE

         IF( NONOUI .GT. 0 ) THEN
C           NSDLFR EST DANS UN TETRAEDRE DE SOMMET NSNOFR
C           IDENTIFICATION POSSIBLE MAIS A VERIFIER
            GOTO 5
         ENDIF

      ENDDO

C     NSDLFR EST EXTERIEUR A TOUS LES TETRAEDRES DE SOMMET NSNOFR
C     IDENTIFICATION IMPOSSIBLE
C     -----------------------------------------------------------
      IERR = 11
ccc      PRINT*,'sistst: St',NSDLFR,' EXTERIEUR A L ETOILE DU St',
ccc     %        NSNOFR,' COBARY=',COBARY
ccc      PRINT*,'sistst: Identification',NSNOFR,'->',NSDLFR,' REFUSEE'
      GOTO 9999


C     AJOUT DES TETRAEDRES DE SOMMET NSDLFR AUX TETRAEDRES NOTENS
C     -----------------------------------------------------------
 5    CALL TETR1S( NSDLFR, N1TETS,        NOTETR,
     %             NBTNS,  MXTENS-NBTENS, NOTENS(NBTENS+1), IERR )

      IF( NBTNS .GT. 64 ) THEN
         PRINT*,'sistst: Sommet',NSDLFR,' dans',NBTNS,
     %' tetraedres. BEAUCOUP TROP pour une SUBSTITUTION avec le st',
     %NSNOFR
         IERR = 16
         GOTO 9999
       ENDIF

      NBTENS = NBTENS + NBTNS

C     UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTENS
      CALL UNITABL( NOTENS, NBTENS )

      IF( IERR .NE. 0 ) THEN
         PRINT*,'sistst: sortie tetr1s ProBleme sur les TETRAEDRES de SO
     %MMET NSDLFR=',NSDLFR,' IERR=',IERR
         TRACTE = .TRUE.
         IERR = 13
         GOTO 9999
      ENDIF

      IF( IVOLTE .NE. 0 ) THEN
C        COMPARAISON DES NO DE VOLUME DES TETRAEDRES DE NSDLFR + NSNOFR
         NOVOLU = NVOLTE( NOTENS(1) )
         DO N=2,NBTENS
            NTE = NOTENS( N )
            IF( NVOLTE( NTE ) .NE. NOVOLU ) THEN
C              DEPLACEMENT MODIFIANT LES VOLUMES => INTERDIT
               IERR = 14
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

C     CONSTRUCTION NFETOI DES TRIANGLES FACES SIMPLES DES TETRAEDRES
C     DE L'ETOILE DES SOMMETS NSNOFR et NSDLFR POUR VERIFICATIONS
      CALL CRFETOI1( NBTENS, NOTENS, NOTETR,
     %               MXFETO, N1FEOC, N1FEVI, NFETOI )

C     NFETOI CONTIENT L'ETOILE DES SOMMETS NSNOFR et NSDLFR en VERSION1
C     TABLEAU NFETOI VERSION 1 => VERSION 2
      CALL NFETOI12( NOTETR, N1FEOC, NFETOI )

C     VERIFIER QUE TOUTE ARETE DE L'ETOILE APPARTIENT
C     SEULEMENT A 2n FACES DE L'ETOILE avec n=1,2
      CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )
C     NBARPB = NOMBRE D'ARETES DES FACES DE L'ETOILE N'APPARTENANT PAS A
C              2 ET SEULEMENT 2 FACES DE L'ETOILE
      IF( NBARPB .GT. 0 ) THEN
         PRINT*,'sistst: NBARPB=',NBARPB,
     %         ' FACES DE L''ETOILE MAL ORIENTEES pour cette SIMULATION'
         IERR = 9
         GOTO 9999
      ENDIF

ccc      PRINT*,'sistst: Nombre de TETRAEDRES de SOMMETS NSNOFR et NSDUC NB
ccc     %TENS=', NBTENS

C     NFETOI CONTIENT L'ETOILE DES SOMMETS NSNOFR et NSDLFR en VERSION2
C     SUPPRESSION DES TETRAEDRES D'ARETE NSNOFR-NSDLFR CAR DE VOLUME NUL
      VOLUE1 = 0D0
      VOLUE0 = 0D0
      NBTENS1= NBTENS
      DO 10 NT = 1, NBTENS

C        TETRAEDRE DE SOMMET NSNOFR ou NSDLFR
         NTE = NOTENS( NT )

         IF( NTE .LE. 0 ) THEN
            PRINT*,'sistst: PB  TETRAEDRE',NTE,' ?'
            GOTO 10
         ENDIF

         IF( NOTETR(1,NTE) .EQ. 0 ) THEN
            PRINT*,'sistst: PB du TETRAEDRE',NTE,' INACTIF',
     %      ' St:',(NOTETR(L,NTE),L=1,8)
            GOTO 10
         ENDIF

C        SON VOLUME INITIAL ET SA QUALITE INITIALE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUTE0, QUALTE0 )
         VOLUE0 = VOLUE0 + VOLUTE0

ccc         PRINT *,'sistst: ',NT,' avant NOTETR(',NTE,')=',
ccc     %           (NOTETR(M,NTE),M=1,8),' Q0=',QUALTE0,' V0=',VOLUTE0

C        BOUCLE SUR LES SOMMETS DU TETRAEDRE NTE 
         DO N0=1,4

C           LE SOMMET N0 DE NTE EST IL LE SOMMET NSNOFR?
            IF( NOTETR(N0,NTE) .EQ. NSNOFR ) THEN

C              OUI: NTE DE SOMMET NSNOFR A T IL AUSSI LE SOMMET NSDLFR?
               DO N1=1,4

                  IF( NOTETR(N1,NTE) .EQ. NSDLFR ) THEN

C                    NTE TETRAEDRE D'ARETE NSNOFR-NSDLFR =sera SUPPRIME
C                    --------------------------------------------------
ccc                     PRINT *,'sistst: ',NT,' SUppr NOTETR(',NTE,')=',
ccc     %                       (NOTETR(M,NTE),M=1,8)
ccc                     NOTENS( NT ) = -NTE
                     GOTO 10

                  ENDIF

               ENDDO

C              TETRAEDRE AVEC NSNOFR MAIS PAS NSDLFR
C              DANS NTE NSNOFR -> NSDLFR
C              -------------------------------------
               DO M=1,4
                  NOSOTE( M ) = NOTETR( M, NTE )
               ENDDO
               NOSOTE( N0 ) = NSDLFR

C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD( 1, NOSOTE(1) ),
     %                       PTXYZD( 1, NOSOTE(2) ),
     %                       PTXYZD( 1, NOSOTE(3) ),
     %                       PTXYZD( 1, NOSOTE(4) ),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE1, QUALTE1 )
               VOLUE1 = VOLUE1 + ABS( VOLUTE1 )

               IF( VOLUTE1 .LE. 0D0 ) THEN
C                 IDENTIFICATION IMPOSSIBLE
ccc                  PRINT *,'sistst: ',NT,' NOTETR(',NTE,')=',
ccc     %            NOSOTE,(NOTETR(M,NTE),M=5,8),
ccc     %           ' Q1=',QUALTE1,' V1=',VOLUTE1,
ccc     %           ' IDENTIFICATION',NSNOFR,'->',NSDLFR,' REFUSEE'
                  IERR = 11
                  GOTO 9999
               ENDIF

C              LES TETRAEDRES OPPOSES A NTE SONT AJOUTES A NOTENS
               DO 8 K = 5, 8

C                 TETRAEDRE OPPOSE A LA FACE K-4 DE NTE
                  NTEOP = NOTETR( K, NTE )

                  IF( NTEOP.GT.0 .AND. NOTETR(1,NTEOP).GT.0 ) THEN

C                    SON VOLUME ET SA QUALITE
                     CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP) ),
     %                             PTXYZD( 1, NOTETR(2,NTEOP) ),
     %                             PTXYZD( 1, NOTETR(3,NTEOP) ),
     %                             PTXYZD( 1, NOTETR(4,NTEOP) ),
     %                          ARMIN, ARMAX, SURFTR, VOLUTE0, QUALTE0 )

C                    NTEOP EST IL DEJA DANS LA LISTE?
                     DO L=1,NBTENS1
                        IF( NTEOP .EQ. NOTENS(L) ) GOTO 8
                     ENDDO
C                    NON: NTEOP EST AJOUTE A NOTENS
                     NBTENS1 = NBTENS1 + 1
                     NOTENS( NBTENS1 ) = NTEOP

                     VOLUE0 = VOLUE0 + VOLUTE0
                     VOLUE1 = VOLUE1 + ABS( VOLUTE0 )

ccc                     PRINT *,'sistst: ',NT,' ajout NOTETR(',NTEOP,')=',
ccc     %                        (NOTETR(M,NTEOP),M=1,8),
ccc     %                       ' Q0=',QUALTE0,' V0=',VOLUTE0
                  ENDIF

 8             ENDDO

ccc               PRINT *,'sistst: ',NT,' apres NOTETR(',NTE,')=',
ccc     %                  NOSOTE,(NOTETR(M,NTE),M=5,8),
ccc     %                 ' Q1=',QUALTE1,' V1=',VOLUTE1

ccc               IF( VOLUTE1 .LT. 0D0 .AND.
ccc     %             ABS(VOLUTE1) .GT. ABS(VOLUTE0)*0.025D0 ) THEN
cccC                 NSNOFR PASSE DERRIERE LA FACE OPPOSEE DE NTE
cccC                 IDENTIFICATION IMPOSSIBLE
ccc                  IERR = 10
ccc                  PRINT *,'sistst: TETRAEDRE',NTE,' V0=',VOLUTE0,
ccc     %                    ' V1=',VOLUTE1,
ccc     %                    ' VOLUMES TROP DIFFERENTS. Identification',
ccc     %                     NSNOFR,'->',NSDLFR,' REFUSEE'
ccc                  PRINT *,'sistst: NOTETR(',NTE,')=',
ccc     %                     NOSOTE,(NOTETR(M,NTE),M=5,8),
ccc     %                    ' Qual1=',QUALTE1
ccc                  GOTO 9999
ccc               ENDIF

               GOTO 10

            ENDIF

         ENDDO

C        ICI: NSNOFR N'EST PAS UN SOMMET DU TETRAEDRE NTE
C        ------------------------------------------------
C        CE TETRAEDRE NTE EST GARDE
C        CALCUL DE SON VOLUME ET DE SA QUALITE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUTE1, QUALTE1 )
         VOLUE1 = VOLUE1 + ABS( VOLUTE1 )
ccc         PRINT *,'sistst: ',NT,' egal  NOTETR(',NTE,')=',
ccc     %           (NOTETR(M,NTE),M=1,8),
ccc     %           ' Q1=',QUALTE1,' V1=',VOLUTE1

 10   ENDDO

      NBTENS = NBTENS1

C     LES VOLUMES DE L'ETOILE AVANT ET APRES SONT ILS EGAUX?
ccc      IF( ABS(VOLUE1-VOLUE0) .GT. ABS(VOLUE0)*0.025D0 ) THEN
      IF( ABS(VOLUE1-VOLUE0) .GT. ABS(VOLUE0)*0.03D0 ) THEN

C        IDENTIFICATION IMPOSSIBLE
         IERR = 11
ccc         PRINT *,'sistst: VolET0=',VOLUE0,' VolET1=',VOLUE1,
ccc     %           ' V1/V0=',VOLUE1/VOLUE0,
ccc     %           ' VOLUMES TROP DIFFERENTS. IDENTIFICATION',
ccc     %           NSNOFR,'->',NSDLFR,' REFUSEE'

      ELSE

C        IDENTIFICATION POSSIBLE
         IERR = 0

ccc         PRINT *,'sistst: VolET0=',VOLUE0,' VolET1=',VOLUE1,
ccc     %           ' V1/V0=',VOLUE1/VOLUE0,
ccc     %           ' IERR=',IERR,' IDENTIFICATION',
ccc     %           NSNOFR,'->',NSDLFR,' POSSIBLE'

      ENDIF

 9999 RETURN
      END
