      SUBROUTINE CRTEFAET( NPt,    MXSOMM, PTXYZD,  NPSOFR,
     %                     MXTETR, N1TEVI, NUDTETR, NOTETR, N1TETS,
     %                     NOVOLU, IVOLTE, NVOLTE,
     %                     MXETOI, N1FEOC, NFETOI,
     %                     NBTEET0,NBTEET1,NTETOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DES TETRAEDRES NPt-FACES SIMPLES DE
C -----    L'ETOILE DES TETRAEDRES CONTENANT LE POINT NPt
C          AUCUNE VERIFICATION SUR L'EGALITE DES VOLUMES AVANT et APRES
C          N'EST EFFECTUEE

C ENTREES:
C --------
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C MXSOMM : NOMBRE MAXIMAL DE POINTS DECLARABLES DANS PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : NUMERO DES POINTS INITIAUX
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET D'OT ET RECONNU TROP PROCHE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES STOCKABLES DANS NOTETR
C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C MXETOI : NOMBRE MAXIMAL DE FACES DU TABLEAU NFETOI et
C          NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NTETOI
C NBTEET0: NOMBRE DE TETRAEDRES DE NTETOI ETOILE DE NPt
C NOVOLU : NUMERO DE VOLUME DES NBTEET0 TETRAEDRES DE L'ETOILE
C          =0 SI IVOLTE=0

C ENTREES ET SORTIES :
C --------------------
C N1TEVI : NUMERO DU 1 PREMIER TETRAEDRE VIDE DANS LE TABLEAU NOTETR
C          LE CHAINAGE DES TETRAEDRES VIDES SE FAIT SUR NOTETR(5,.)
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE ACTIF OCCUPE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET
C NVOLTE : >0 NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE
C          -1 SI TETRAEDRE INACTIF ou SUPPRIME

C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 1
C          1: NUMERO DU TETRAEDRE DANS NOTETR
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3  NON UTILISE ICI
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C SORTIES:
C --------
C NBTEET1: NOMBRE DE TETRAEDRES DE NTETOI DE SOMMET NPt
C NTETOI : NUMERO DANS NOTETR DES NBTEET TETRAEDRES DE L'ETOILE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Auteur : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2019
C23456...............................................................012
      include"./incl/langue.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*100     KTITRE
      INTEGER           NOTETR(8,MXTETR),N1TETS(MXSOMM),NPSOFR(MXSOMM),
     %                  NTETOI(MXETOI),NFETOI(5,MXETOI),
     %                  NVOLTE(MXTETR)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM), VTE, VMIMIN, XYZBAR(4)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NSFC(3)
      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      TRACTE0 = TRACTE

C     VOLUME MINIMAL D'UN TETRAEDRE A CREER
ccc      VMIMIN = ( ( PTXYZD( 4, NPt ) ** 3 ) / 6D0 ) * 1D-4
      VMIMIN = ( ( PTXYZD( 4, NPt ) ** 3 ) / 6D0 ) * 1D-5

      NBTEET  = NBTEET0

 10   IF( NBTEET .LE. 0 ) THEN
         NBTEET1 = 0
         RETURN
      ENDIF

      IF( NBTEET .EQ. 1 .AND.
     %   ( NPSOFR(NPt) .EQ. 0 .OR. NPSOFR(NPt) .EQ. -4) ) THEN
C        DEPLACEMENT VERS LE BARYCENTRE DU POINT NPt
         DO K=1,4
            XYZBAR(K) = 0D0
         ENDDO
         NT1 = NTETOI(1)
         IF( NT1 .LE. 0 ) THEN
            PRINT*,'crtefaet: NTETOI(1)=',NT1,' A REGARDER...'
            GOTO 15
         ENDIF
         DO N=1,4
            NS = NOTETR( N, NT1 )
            DO K=1,4
               XYZBAR(K) = XYZBAR(K) + PTXYZD( K, NS )
            ENDDO
         ENDDO
C        XYZD DU BARYCENTRE DU TETRAEDRE NT1
C        DEPLACEMENT DU POINT NPt VERS LE BARYCENTRE DE NT1
         DO K=1,4
            XYZBAR( K ) = XYZBAR( K ) / 4
            PTXYZD( K, NPt ) = 0.9D0 * PTXYZD( K, NPt )
     %                       + 0.1D0 * XYZBAR( K )
         ENDDO
      ENDIF

 15   KTITRE='crtefaet: NPt=          NBTEET=         TETRAEDRES FORMENT
     % L''ETOILE'
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(33:37),'(I5)') NBTEET
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET, NTETOI )

C     CONSTRUCTION DE L'ETOILE DES FACES SIMPLES DES TETRAEDRES
C     LES FACES VUES 2 FOIS SONT ELIMINEES
C     LES FACES VUES UNE FOIS SONT AJOUTEES
C     VERSION 1 DU TABLEAU NFETOI DES FACES SIMPLES DE L'ETOILE
C     ---------------------------------------------------------
      CALL CRFETOI1( NBTEET, NTETOI, NOTETR,
     %               MXETOI, N1FEOC, N1FEVI, NFETOI )

C     TRACE DES FACES SIMPLES DE L'ETOILE du POINT NPt
      KTITRE='crtefaet: NPt=          LES FACES SIMPLES DE L''ETOILE INI
     %TIALE'
      WRITE(KTITRE(15:22),'(I8)') NPt
      CALL SANSDBL( KTITRE, NBC )
      CALL TRFETOV1(KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR, PTXYZD)

C     PASSAGE DE NFETOI VERSION 1 => VERSION 2
C     LE NO DE TETRAEDRE INTERNE ETOILE  => NO TETRAEDRE AU DELA DE LA FACE
C     LE NO FACE LOCAL DANS LE TETRAEDRE => NO 1-ER SOMMET DE LA FACE
C     LE NO INUTILISE 3                  => NO 2-ME SOMMET DE LA FACE
C     LE NO INUTILISE 4                  => NO 3-ME SOMMET DE LA FACE
C     ----------------------------------------------------------------------
      NBFETO = 0
      NF0    = 0
      NF1    = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 510  IF( NF1 .GT. 0 ) THEN

C        NFETOI VERSION 1.  UNE FACE SIMPLE DE L'ETOILE DE PLUS
C        LE NO DU TETRAEDRE NT1 INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT1   = NFETOI( 1, NF1 )
         NFNT1 = ABS( NFETOI( 2, NF1 ) )

C        PASSAGE DE NFETOI EN VERSION 2
C        LE NUMERO DU TETRAEDRE NTOP AU DELA DE LA FACE NFNT1 DE NT1
C        QUI PEUT ETRE NUL SI LA FACE EST FRONTIERE
         NTOP = NOTETR( 4+NFNT1, NT1 )

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NS1 = NOTETR( NOSOFATE(1,NFNT1), NT1 )
         NS2 = NOTETR( NOSOFATE(2,NFNT1), NT1 )
         NS3 = NOTETR( NOSOFATE(3,NFNT1), NT1 )

C        QUALITE ET VOLUME DU TETRAEDRE NF2-NPt A CREER
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         IF( NTOP .EQ. 0 .AND. VTE .LT. VMIMIN ) THEN

C           TETRAEDRE FRONTALIER DE TROP FAIBLE VOLUME POUR ETRE CREE
C           SUPPRESSION DE LA FACE NF1 DE L'ETOILE NFETOI
            IF( NF0 .EQ. 0 ) THEN
               N1FEOC = NFETOI( 5, N1FEOC )
            ELSE
               NFETOI( 5, NF0 ) = NFETOI( 5, NF1 )
            ENDIF
            GOTO 520

         ENDIF

         IF( VTE .LT. VMIMIN ) THEN
C           TETRAEDRE DE VOLUME TROP FAIBLE OU NEGATIF POUR ETRE CREE
C           LE TETRAEDRE NT1=NF1+NPt EST RETIRE DE L'ETOILE
            DO K=1,NBTEET
               IF( NT1 .EQ. NTETOI(K) ) THEN
                  DO L=K+1,NBTEET
                     NTETOI( L-1 ) = NTETOI( L )
                  ENDDO
                  NBTEET = NBTEET - 1
                  GOTO 10
               ENDIF
            ENDDO
         ENDIF

         IF( NTOP .GT. 0 ) THEN
C           RECHERCHE DU NO NFOP LOCAL A NTOP DE LA FACE NS1 NS2 NS3
            NSFC(1) = NS1
            NSFC(2) = NS2
            NSFC(3) = NS3
            CALL TRI3NO( NSFC, NSFC )
            CALL NO1F1T( NSFC, NOTETR(1,NTOP), NFOP )
            IF( NFOP .LE. 0 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*, 'crtefaet: NPt=',NPt,' 3 SOMMETS',NS1,NS2,NS3
                  PRINT*, 'NON DANS LE TETRAEDRE ',NTOP,' :',
     %                    (NOTETR(kk,NTOP),kk=1,8)
               ELSE
                  PRINT*, 'crtefaet: NPt=',NPt,' 3 VERTICES',NS1,NS2,NS3
                  PRINT*, 'NOT in TETRAHEDRON ',NTOP,' :',
     %                    (NOTETR(kk,NTOP),kk=1,8)
               ENDIF
C              LE TETRAEDRE NT1 EST RETIRE DE L'ETOILE
               DO K=1,NBTEET
                  IF( NT1 .EQ. NTETOI(K) ) THEN
                     DO L=K+1,NBTEET
                        NTETOI( L-1 ) = NTETOI( L )
                     ENDDO
                     NBTEET = NBTEET - 1
                     GOTO 10
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

C        UNE FACE DE PLUS ACCEPTEE
         NBFETO = NBFETO + 1

C        LE TETRAEDRE OPPOSE PAR LA FACE NF1
         NFETOI( 1, NF1 ) = NTOP

C        LES 3 SOMMETS DE LA FACE NF1
         NFETOI( 2, NF1 ) = NS1
         NFETOI( 3, NF1 ) = NS2
         NFETOI( 4, NF1 ) = NS3

C        PASSAGE A LA FACE SUIVANTE
         NF0 = NF1
 520     NF1 = NFETOI( 5, NF1 )
         GOTO 510

      ENDIF

      IF( NBFETO .LE. 0 ) THEN
C        AUCUN FACE SIMPLE DANS L'ETOILE
         NBTEET1 = 0
         GOTO 9999
      ENDIF


C     DESTRUCTION DES TETRAEDRES INITIAUX DE L'ETOILE DE NPt
C     ------------------------------------------------------
      DO I=1,NBTEET

         NT1 = NTETOI( I )
         IF( NT1 .GT. 0 ) THEN

C           MISE A JOUR DE N1TETS DES 4 SOMMETS AVANT NOUVEAUX TETRAEDRES
C           CELA EST NECESSAIRE POUR ELIMINER LES SOMMETS DUS A L'AJOUT
C           DES TETRAEDRES ENCOCHES
            DO K=1,4
               N1TETS( NOTETR(K,NT1) ) = 0
            ENDDO

C           NT1 DEVIENT LE PREMIER TETRAEDRE VIDE
            NOTETR(1,NT1) = 0
            NOTETR(2,NT1) = 0
            NOTETR(3,NT1) = 0
            NOTETR(4,NT1) = 0

            DO NF1 = 1, 4
C              NFOP NO LOCAL A NTOP DE LA FACE NF1 DE NT1
               CALL NOFAOP( NF1, NT1, NOTETR,  NFOP, NTOP )
               IF( NFOP .GT. 0 ) THEN
C                 NT1 EST SUPPRIME -> TETRAEDRE OPPOSE INCONNU
                  NOTETR( 4+NFOP, NTOP ) = -1
               ENDIF
               NOTETR( 4+NF1, NT1 ) = -1
            ENDDO

            NOTETR(5,NT1) = N1TEVI
            N1TEVI        = NT1

         ENDIF
      ENDDO

C     DECLARATION DES NBFETO NOUVEAUX TETRAEDRES DE SOMMET NPt-FACE SIMPLE
C     --------------------------------------------------------------------
      QUAMIN  = 2.
      QUAMOY  = 0.
      NBTEET1 = 0
      NBTEET2 = NBFETO
      NF2     = N1FEOC

C     BOUCLE SUR LES NBFETO FACES SIMPLES DE L'ETOILE
 560  IF( NF2 .GT. 0 ) THEN

C        REMPLISSAGE DU TETRAEDRE NT DE SOMMETS CEUX DE LA FACE NF2 ET NPt
         NS1 = NFETOI(2,NF2)
         NS2 = NFETOI(3,NF2)
         NS3 = NFETOI(4,NF2)

C        LE TETRAEDRE EXTERIEUR A LA FACE AU DELA DE LA FACE
C        ATTENTION: IL PEUT NE PAS EXISTER DE TEL TETRAEDRE
         NTOP = NFETOI(1,NF2)

C        QUALITE ET VOLUME DU TETRAEDRE NF2-NPt A CREER
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         IF( VTE .LE. 0D0 ) THEN
            PRINT*,'crtefaet: NPt=',NPt,' NTOP=',NTOP,
     %             ' tetraedre',ns1,ns2,ns3,npt,
     %             ' CREE MALGRE son VOLUME',VTE,' TROP FAIBLE'
ccc            tracte = .true.
ccc            GOTO 567
         ENDIF

C        QUALITE MOYENNE ET MINIMALE DES TETRAEDRES CREES
         QUAMOY = QUAMOY + QTE
         QUAMIN = MIN( QUAMIN, QTE )

C        LA FACE NF2 ET LE POINT NPt FORMENT UN NOUVEAU TETRAEDRE
         IF( N1TEVI .LE. 0 ) THEN
C           SATURATION DES TETRAEDRES NOTETR
            GOTO 9910
         ENDIF

C        MISE A JOUR DU DERNIER TETRAEDRE OCCUPE
         NUDTETR = MAX( NUDTETR, N1TEVI )
C        MISE A JOUR DU 1-ER TETRAEDRE VIDE
         NT     = N1TEVI
         NOTETR(5,N1TEVI) = ABS( NOTETR(5,N1TEVI) )
         N1TEVI = NOTETR(5,N1TEVI)

C        NS1 NS2 NS3 SOMMETS DE LA FACE DANS LE SENS DIRECT DU TETRAEDRE NT
         NOTETR(1,NT) = NS1
         NOTETR(2,NT) = NS2
         NOTETR(3,NT) = NS3
         NOTETR(4,NT) = NPt

C        LE CHAINAGE DES TETRAEDRES OPPOSES  PAR LES FACES
         IF( NTOP .GT. 0 ) THEN
C           STOCKAGE DU TETRAEDRE OPPOSE POUR RESOUDRE LES TETRAEDRES OPPOSES
            DO K=NBFETO+1,NBTEET2
               IF( NTETOI(K) .EQ. NTOP ) GOTO 565
            ENDDO
            IF( NBTEET2 .GE. MXETOI ) THEN
C              SATURATION DES TETRAEDRES DE L'ETOILE NTETOI
               GOTO 9900
            ENDIF
            NBTEET2 = NBTEET2 + 1
            NTETOI( NBTEET2 ) = NTOP
         ENDIF

 565     NOTETR(5,NT) = NTOP
         NOTETR(6,NT) = -1
         NOTETR(7,NT) = -1
         NOTETR(8,NT) = -1

         IF( NTOP .GT. 0 ) THEN
C           RECHERCHE DU NO NFOP LOCAL A NTOP DE LA FACE NS1 NS2 NS3
            CALL TRI3NO( NFETOI(2,NF2), NSFC )
            CALL NO1F1T( NSFC, NOTETR(1,NTOP), NFOP )
            IF( NFOP .GT. 0 ) THEN
               NOTETR( 4+NFOP, NTOP ) = NT
            ELSE
               PRINT*,'crtefaet: NPt=',NPt,' NT=',NT,' NTOP=',NTOP,
     %                ' CAS IMPOSSIBLE...!'
               PRINT*,'crtefaet: NOTETR(',NTOP,')=',
     %                (NOTETR(kk,NTOP),kk=1,8)
               PRINT*
            ENDIF
         ENDIF

C        MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
C        TOUS LES TETRAEDRES ONT PU DISPARAITRE
C        CAS DE L'ETOILE REDUITE AUX 2 TETRAEDRES INITIAUX
         N1TETS( NS1 ) = NT
         N1TETS( NS2 ) = NT
         N1TETS( NS3 ) = NT
         N1TETS( NPt ) = NT

C        LES TETRAEDRES CREES DANS L'ETOILE
         IF( NBTEET1 .GE. MXETOI ) THEN
C           SATURATION DES TETRAEDRES DE L'ETOILE NTETOI
            GOTO 9900
         ENDIF
         NBTEET1 = NBTEET1 + 1
         NTETOI( NBTEET1 ) = NT

C        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
ccc 567     NF2 = NFETOI( 5, NF2 )
         NF2 = NFETOI( 5, NF2 )
         GOTO 560

      ENDIF
C     QUALITE MOYENNE DES NBTEET1 TETRAEDRES CREES
      QUAMOY = QUAMOY / NBTEET1

C     SI DES TETRAEDRES N'ONT PAS ETE CREES
      IF( NBTEET1 .LT. NBFETO ) THEN
         DO K=NBFETO+1,NBTEET2
            NTETOI( K-NBFETO+NBTEET1 ) = NTETOI( K )
         ENDDO
      ENDIF

C     COMPLETION DES CHAINAGES DES TETRAEDRES CREES DANS L'ETOILE
C     -----------------------------------------------------------
      KTITRE='crtefaet: NPt=          NBTEET2=        TETRAEDRES avant C
     %OMPLETION des TETRAEDRES OPPOSES'
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(33:37),'(I5)') NBTEET2
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET2, NTETOI )

      NBTE = NBTEET2
 570  IF( NBTE .GT. 0 ) THEN

C        LE TETRAEDRE DU HAUT DE LA PILE
         NT = NTETOI( NBTE )
C        LE TETRAEDRE EST DEPILE
         NBTE = NBTE - 1
         IF( NT .LE. 0 ) GOTO 570

c        verification. A supprimer ensuite...
         if( notetr(1,NT) .EQ. 0 ) THEN
            PRINT*,'crtefaet: Anomalie NPt=',NPt,' XYZD=',
     %             (PTXYZD(mmm,NPt),mmm=1,4)
            PRINT*,'crtefaet: Anomalie NOTETR(',NT,')=',
     %             (NOTETR(mmm,NT),mmm=1,8)
         ENDIF

C        QUEL EST LE TETRAEDRE OPPOSE AUX FACES 1 2 3 4
C             ET DE SOMMETS NS1 NS2 NS3 = NSFC
         DO 590 NF=1,4

C           LE TETRAEDRE OPPOSE A LA FACE NF DU TETRAEDRE NT
 585        NTOP = NOTETR( 4+NF, NT )
            IF( NTOP .EQ. 0 ) THEN
C              FACE FRONTIERE
               GOTO 590
            ENDIF

C           LE NUMERO DES 3 SOMMETS DE LA FACE NF DU TETRAEDRE NT
            NSFC(1) = NOTETR( NOSOFATE(1,NF), NT )
            NSFC(2) = NOTETR( NOSOFATE(2,NF), NT )
            NSFC(3) = NOTETR( NOSOFATE(3,NF), NT )
C           TRI CROISSANT DE SES 3 SOMMETS
            CALL TRI3NO( NSFC, NSFC )

            IF( NTOP .LT. 0 ) THEN

C              FACE NF de NT NON TRAITEE. TETRAEDRE OPPOSE INCONNU
C              RECHERCHE de CETTE FACE PARMI LES TETRAEDRES CREES
               DO 580 I=1,NBTE

                  NT1 = NTETOI(I)
                  IF( NT1 .LE. 0 ) GOTO 580

C                 RECHERCHE DE LA FACE NSFC DANS CE TETRAEDRE NT1
                  CALL NO1F1T( NSFC, NOTETR(1,NT1), NF1 )
                  IF( NF1 .LE. 0 ) GOTO 580

C                 LA FACE NF DE NT EST LA FACE NF1 DE NT1
C                 LE NUMERO NOTETR DU TETRAEDRE A EFFACER
                  NT2 = NOTETR( 4+NF1, NT1 )

                  IF( NT2 .GT. 0 .AND. NT2 .NE. NT ) THEN
C                    ATTENTION: 1 FACE COMMUNE A 3 TETRAEDRES???
                     DO L=1,NBTEET2
                        IF( NT2 .EQ. ABS(NTETOI(L)) ) GOTO 586
                     ENDDO

C                    PROBLEME NON RESOLU!...............................
 586                 print*,'crtefaet: NPt=',NPt,' avec PROBLEME'
                     print*,'crtefaet: nbt=',nbte+1,' notetr(',nt,')=',
     %                      (notetr(kk,nt),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'crtefaet: nbt=',I,' notetr(',nt1,')=',
     %                      (notetr(kk,nt1),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'crtefaet: nbt=',L,' notetr(',nt2,')=',
     %                      (notetr(kk,nt2),kk=1,8),' AVEC FACE TRIPLE?'
                     print*
      KTITRE='crtefaet: NPt=          NBTEET2=        TETRAEDRES PENDANT
     % la COMPLETION des TETRAEDRES OPPOSES'
                    WRITE(KTITRE(15:22),'(I8)') NPt
                    WRITE(KTITRE(33:37),'(I5)') NBTEET2
                    CALL SANSDBL( KTITRE, NBC )
                    tracte = .true.
                    CALL TRACEETOILE(KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                               NBTEET2, NTETOI )
C                   PROBLEME A RESOUDRE!................................
                  ENDIF

C                 OPPOSITION DES 2 TETRAEDRES NT NT1 A TRAVERS LA FACE NSFC
                  NOTETR( 4+NF , NT  ) = NT1
                  NOTETR( 4+NF1, NT1 ) = NT

                  IF( NOTETR(5,NT1) .GE. 0 .AND.
     %                NOTETR(6,NT1) .GE. 0 .AND.
     %                NOTETR(7,NT1) .GE. 0 .AND.
     %                NOTETR(8,NT1) .GE. 0 ) THEN
C                     FACE TOTALEMENT TRAITEE =>
C                     RETRAIT DES TETRAEDRES DE RECHERCHE
                      NTETOI(I) = -NTETOI(I)
                  ENDIF

C                 PASSAGE A LA FACE NF DE NT SUIVANTE
                  GOTO 590

 580           ENDDO

C              AUCUN DES TETRAEDRES A CETTE FACE EN COMMUN
C              LA FACE NF DU TETRAEDRE EST DONC FRONTIERE
ccc            PRINT*,'crtefaet: MISE A ZERO DU TETRA OPPOSE A LA FACE',NF,
ccc     %             ' du TETRAEDRE',NT
ccc            PRINT*,'crtefaet: tetraedre',nt,':',(notetr(kk,nt),kk=1,8)

               NOTETR( 4+NF, NT ) = 0

ccc            PRINT*,'crtefaet: tetraedre',nt,':',(notetr(kk,nt),kk=1,8)
ccc            PRINT*
ccc            tracte = .true.

            ELSE

C              FACE A PRIORI TRAITEE. VERIFICATION
C              RECHERCHE DE LA FACE NSFC DANS CE TETRAEDRE NTOP
               CALL NO1F1T( NSFC, NOTETR(1,NTOP), NFOP )
               IF( NFOP .LE. 0 ) GOTO 590
               NTE = NOTETR( 4+NFOP, NTOP )
               IF( NTE .NE. NT ) THEN
C                 PROBLEME: LE TETRAEDRE OPPOSE PAR LA FACE NF DE NT
C                           N'EST PAS A JOUR
                  NOTETR( 4+NF, NT ) = -1
C                 RETOUR EN RECHERCHE DE TETRAEDRE OPPOSE
                  GOTO 585
               ENDIF

            ENDIF

 590     ENDDO

C        RETOUR EN HAUT DE PILE
         GOTO 570
      ENDIF

C     MISE A JOUR DU NUMERO POSITIF DES TETRAEDRES DE L'ETOILE de NPt
      DO I = 1,NBTEET1
         NT = ABS( NTETOI( I ) )
         NTETOI( I ) = NT
         IF( IVOLTE .NE. 0 ) THEN
            NVOLTE( NT ) = NOVOLU
         ENDIF
      ENDDO

C     VERIFIER L'ABSENCE DE TETRAEDRES OPPOSES DOUBLES ET
C     L'OPPOSITION DES TETRAEDRES DANS UNE LISTE DE TETRAEDRES
C     --------------------------------------------------------
      CALL VEOPTE( NBTEET1, NTETOI, NOTETR, PTXYZD, NBFANR )


C     TRACE DES NBTEET1 TETRAEDRES DE SOMMETS NPt
C     -------------------------------------------
      KTITRE='crtefaet: NPt=          SOMMET des NBTEET1=        TETRAED
     %RES FINAUX'
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(44:48),'(I5)') NBTEET1
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET1, NTETOI )
      GOTO 9999


C     SATURATION DES TETRAEDRES NTETOI
 9900 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'crtefaet: SATURATION DES TETRAEDRES NTETOI'
      ELSE
         PRINT*,'crtefaet: SATURATION of NTETOI TETRAHEDRA'
      ENDIF
      NBTEET1 = 0
      GOTO 9999


C     SATURATION DES TETRAEDRES NOTETR
 9910 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'crtefaet: SATURATION DES TETRAEDRES NOTETR'
      ELSE
         PRINT*,'crtefaet: SATURATION of TETRAHEDRA NOTETR'
      ENDIF
      NBTEET1 = 0


 9999 TRACTE = TRACTE0
      RETURN
      END
