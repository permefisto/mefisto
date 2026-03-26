      SUBROUTINE ETOIPT19( NPt,    QTEAME, ARETGR, NOTET0, NOTET1,
     %                     MXSOMM, PTXYZD, NPSOFR,
     %                     HEXAPAVE,NBIPAV,ECHPAV, N1SPAVE, NOPTSUIV,
     %                     MXTETR, N1TEVI, NOTETR,NUDTETR, N1TETS,
     %                     INFACO, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                     MXETOI, N1FEOC, N1FEVI, NFETOI,
     %                     NBTEET, NOTEET, VOLET0, QUAMIN, QUAMOY,
     %                     MXSTPE, NBSTPE, NOSTPE, IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DE L'ETOILE DES TETRAEDRES CONTENANT LE POINT NPtf
C -----    DANS UNE ENVELOPPE ETOILEE
C          PROGRESSION DU TETRAEDRE CONTENANT LE POINT VERS LES TETRAEDRES
C          OPPOSES AUX FACES SIMPLES DE L'ETOILE QUI GENERENT AVEC LE POINT
C          DES TETRAEDRES DE VOLUME>0
C          L'APPARTENANCE AUX VOLUMES DIFFERENTS N'EST PAS PRISE EN COMPTE
C          CAR L'AJOUT D'UN POINT SUR UNE FACE COMMUNE ENTRAINERAIT
C          LA DISCONTINUITE DU MAILLAGE

C ENTREES:
C --------
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C QTEAME : QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UN TETRAEDRE
C          EST INTERDIT ou SUPPRIME DE L'ETOILE
C ARETGR : LONGUEUR DE L'ARETE SOHAITEE POUR L'ENSEMBLE DES TETRAEDRES
C NOTET0 : NUMERO DU TETRAEDRE DE DEBUT DE RECHERCHE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DU TABLEAU PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : NUMERO DES POINTS INITIAUX
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET RECONNU TROP PROCHE
C          = -4 SI LE POINT EST SOMMET DE LA GRILLE REGULIERE ELIMINABLE
C          = -3 SI LE POINT EST SOMMET REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE OU NO DE POINT INTERNE

C HEXAPAVE: MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV  : NOMBRE D'ARETES DANS LA DIRECTION I
C ECHPAV  : ECHELLE DANS LA DIRECTION I
C N1SPAVE : NO DU 1-ER SOMMET DANS PTXYZD DU PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES

C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES STOCKABLES DANS NOTETR
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET

C INFACO : = 0 PAS DE TABLEAU LEFACO NI DE CONSERVATION DES
C              FACES FRONTIERE ( NO DE VOLUME CONNU PAR NVOLTE )
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONSERVATION DES
C              FACES DE LA FRONTIERE
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...

C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : >0 NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE
C          -1 SI TETRAEDRE INACTIF ou SUPPRIME

C MXETOI : MAXIMUM DECLARE DANS LES TABLEAUX NOTEET ET NFETOI
C MXSTPE : NOMBRE MAXIMAL DE SOMMETS DES TETRAEDRES DE L'ETOILE NPt

C ENTREES ET SORTIES:
C -------------------
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : 1: NUMERO DU TETRAEDRE DANS NOTETR
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3,4: NON UTILISES ICI
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C SORTIES:
C --------
C NOTET1 : >0 NUMERO DU TETRAEDRE CONTENANT NPt
C          =0 SI PAS DE TETRAEDRE CONTENANT NPt
C NBTEET : NOMBRE DE TETRAEDRES DE L'ETOILE DU POINT NPt
C          =0 SI AUCUN TETRAEDRE ACTUEL CONTIENT LE POINT NPt
C          =-100 LE POINT NPt REMPLACE UN SOMMET DEJA TETRAEDRISE

C NOTEET : NUMERO DANS NOTETR DES NBTEET TETRAEDRES DE L'ETOILE
C VOLET0 : LE VOLUME DE L'ETOILE FINALE DU POINT NPt
C QUAMIN : QUALITE MINIMALE DES NBTEET TETRAEDRES DE L'ETOILE FINALE
C QUAMOY : QUALITE MOYENNE  DES NBTEET TETRAEDRES DE L'ETOILE FINALE
C NBSTPE : NOMBRE DE SOMMETS PERDUS DANS LA TETRAEDRISATION DE L'ETOILE
C          AUTANT DE SOMMETS A RE-TETRAEDRISER
C NOSTPE : NUMERO DES SOMMETS DES TETRAEDRES FINAUX DE L'ETOILE
C          PUIS DES SOMMETS PERDUS
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS    du TABLEAU PTXYZD
C          3 SATURATION DES TETRAEDRES du TABLEAU NOTETR
C          4 TETRAEDRE DEGENERE
C          5 POINT NPt EXTERIEUR AUX TETRAEDRES ACTUELS
C          8 AUCUNE SPHERE CIRCONSCRITE AUX TETRAEDRES ACTUELS
C            NE CONTIENT LE POINT COURANT
C         10 AUCUN POINT DEPLACE DANS L'ETOILE NE CONVIENT
C         11 ETOILE DU POINT NPt IMPOSSIBLE A FORMER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS St PIERRE du PERRAYJanvier 2017
C23456...............................................................012
      include"./incl/langue.inc"
      PARAMETER       ( MXSMIN=1024, MXTECH=20000, MXETOI0=4096 )
      DOUBLE PRECISION  EPSDEG
      PARAMETER        (EPSDEG=1D-10)
      PARAMETER        (MXITER=16)

      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      DOUBLE PRECISION  PTXYZD(4,MXSOMM), CENTRE(4),
     %                  HEXAPAVE(3,2), ECHPAV(3),
     %                  PTPROJ(3), VTMIN, DMIN, DMAX, D, DSR,
     %                  CMAX, VTE, VOLET0, VOLET1, COBARY(4),
     %                  ARMIN, ARMAX, SURFTR(4), VOLTET,
     %                  ARETGD, ARETGD2, ARETGDM

      INTEGER           NOTETR(8,MXTETR), N1TETS(MXSOMM),NOTEET(MXETOI),
     %                  NFETOI(5,MXETOI), NPSOFR(MXSOMM),
     %                  LEFACO(1:11,0:MXFACO),
     %                  NOTECH(MXTECH), NOTERE(MXETOI0),
     %                  NOSTPE(MXSTPE), NVOLTE(MXTETR)

      INTEGER           NBIPAV(3), N1SPAVE(0:*), NOPTSUIV(1:*),
     %                  NOSOTR(3), NOSOTE(4)

      REAL              QTE, QUAMIN, QUAMOY
      CHARACTER*108     KTITRE
      LOGICAL           REPRISE

      INTEGER           NS(4)
      EQUIVALENCE      (NS(1),NS1),(NS(2),NS2),
     %                 (NS(3),NS3),(NS(4),NS4)

      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      ARETGD  = ARETGR
      ARETGDM = 0.7D0 * ARETGD ** 2
      ARETGD2 = 1.9D0 * ARETGD ** 2

CCC      ARETGD2 = 1.7D0 * ARETGD ** 2
ccc      ARETGD2 = 1.6D0 * ARETGD ** 2
ccc      ARETGD2 = 1.4D0 * ARETGD ** 2
ccc      ARETGD2 = 1.3D0 * ARETGD ** 2
CCC      ARETGD2 = 1.8D0 * ARETGD ** 2
ccc      ARETGD2 = 4.1D0 * ARETGD ** 2

      TRACTE0 = TRACTE
      REPRISE = .FALSE.
      NOTENPt = 0
      LOCTAE = 0
      NBREPR = 0
      NBTEET = 0
      NBTEET0= 0
      NBSTPE = 0

C     ======================================================
C     RECHERCHE D'UN TETRAEDRE NOTET1 CONTENANT LE POINT NPt
C     ======================================================
      CALL R1TCO1P1( NPt,     MXSOMM, PTXYZD, 1,
     %               HEXAPAVE,NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %               NOTET0,  MXTETR, NOTETR, NUDTETR, N1TETS,
     %               MXETOI,  NOTEET, MXTECH, NOTECH, 
     %               NOTET1,  COBARY, IERR)

      IF( NOTET1 .LE. 0 ) THEN
C        PAS DE TETRAEDRE CONTENANT NOTET1. ABANDON DU POINT NPt
         WRITE(IMPRIM,*) 'etoipt19 : POINT', NPt, (PTXYZD(1,NPt),L=1,3),
     %   ' EXTERIEUR AUX TETRAEDRES ACTUELS. VERIFIER ses XYZ. ABANDON'
         NBTEET0 = 0
         NBTEET  = 0
         TRACTE  = .TRUE.
         IERR    = 5
         GOTO 9999
      ENDIF

C     TRACE DU TETRAEDRE CONTENANT PTXYZD(1,NPt)
ccc      tracte = .true.
      NBTEET0 = 1
      NBTEET  = 1
      NOTEET(1) = NOTET1
      KTITRE='etoipt19 : POINT          UN TETRAEDRE LE CONTIENT'
      WRITE(KTITRE(18:25),'(I8)') NPt
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )


C     ==================================================================
C     ICI LE TETRAEDRE NOTET1 CONTIENT LE POINT NPt
C     ==================================================================
      NOTENPt = NOTET1

C     NUMERO DE VOLUME DES TETRAEDRES DE L'ETOILE DE NPt
      IF( IVOLTE .NE. 0 ) THEN
         NOVOLU = NVOLTE( NOTET1 )
      ELSE
         NOVOLU = -1
      ENDIF

      NSDMIN  = 0
      IF( NPSOFR(NPt) .GT. 0 ) THEN

C        1) SI LE POINT NPt A TETRAEDRISER EST FRONTIERE, EST IL POSSIBLE
C        DE DEPLACER LE SOMMET, NON FRONTIERE, DEJA TETRAEDRISE LE PLUS
C        PROCHE EN CE POINT NPt SANS GENERER DE TETRAEDRE DE VOLUME <=0 ?
C        ----------------------------------------------------------------
C        RECHERCHE DU SOMMET DU TETRAEDRE NOTET1 (QUI CONTIENT NPt)
C        NON FRONTALIER ET LE PLUS PROCHE DE NPt
         DMIN = 1D100
         DO I=1,4
            N = NOTETR( I, NOTET1 )
            IF( NPSOFR(N) .EQ. 0 ) THEN
C              N N'EST PAS FRONTALIER
               D = ( PTXYZD(1,NPt) - PTXYZD(1,N) ) **2
     %           + ( PTXYZD(2,NPt) - PTXYZD(2,N) ) **2
     %           + ( PTXYZD(3,NPt) - PTXYZD(3,N) ) **2
               IF( D .LT. DMIN ) THEN
                  DMIN   = D
                  NSDMIN = N
               ENDIF
            ENDIF
         ENDDO

         IF( NSDMIN .GT. 0 .AND. DMIN .LE. ARETGDM ) THEN

C           RECHERCHE DE TOUS LES TETRAEDRES DE SOMMET NSDMIN
            CALL TETR1S( NSDMIN, N1TETS, NOTETR,
     %                   NBTEET, MXETOI, NOTEET, IERR )
            NBTEETIN = NBTEET

            DO I=1,NBTEET

C              RECHERCHE DU NUMERO DE SOMMET DE NSDMIN
               NTE = NOTEET( I )
               DO K=1,4
                  IF( NOTETR(K,NTE) .EQ. NSDMIN ) THEN

C                    VOLUME DU TETRAEDRE AVEC NSDMIN DEPLACE EN NPt
                     GOTO( 51, 52, 53, 54 ), K

 51                  NOSOTE( 1 ) = NPt
                     NOSOTE( 2 ) = NOTETR(2,NTE)
                     NOSOTE( 3 ) = NOTETR(3,NTE)
                     NOSOTE( 4 ) = NOTETR(4,NTE)
                     GOTO 55

 52                  NOSOTE( 1 ) = NOTETR(1,NTE)
                     NOSOTE( 2 ) = NPt
                     NOSOTE( 3 ) = NOTETR(3,NTE)
                     NOSOTE( 4 ) = NOTETR(4,NTE)
                     GOTO 55

 53                  NOSOTE( 1 ) = NOTETR(1,NTE)
                     NOSOTE( 2 ) = NOTETR(2,NTE)
                     NOSOTE( 3 ) = NPt
                     NOSOTE( 4 ) = NOTETR(4,NTE)
                     GOTO 55

 54                  NOSOTE( 1 ) = NOTETR(1,NTE)
                     NOSOTE( 2 ) = NOTETR(2,NTE)
                     NOSOTE( 3 ) = NOTETR(3,NTE)
                     NOSOTE( 4 ) = NPt

 55                  CALL CEBOQU( 0., NOSOTE, PTXYZD, CENTRE, VTE, QTE )

                     IF( QTE .LT. QTEAME ) THEN

C                       LE TETRAEDRE SERAIT DE QUALITE MEDIOCRE
C                       NSDMIN NON CORRECT POUR REMPLACER NPt
C                       N'EST PLUS ENVISAGE

ccc                     TRACTE0 = TRACTE
ccc                     TRACTE = .TRUE.
ccc            KTITRE='etoipt19 : IMPOSSIBLE           -/>          ,         
ccc     %             TETRAEDRES'
ccc            WRITE(KTITRE(23:30),'(I8)') NSDMIN
ccc            WRITE(KTITRE(36:43),'(I8)') NPt
ccc            WRITE(KTITRE(47:54),'(I8)') NBTEET
ccc            CALL TRACEETOILE(KTITRE, PTXYZD, NOTETR, NPt, NBTEET,NOTEET)
cccccc            CALL TRFETO13( KTITRE, PTXYZD, NBTEET, NOTEET, NOTETR )
ccc            PRINT *,'etoipt19 : Sommet',NSDMIN,' IMPOSSIBLE au POINT',NPt
ccc                     TRACTE = TRACTE0

                        NSDMIN = 0
                        GOTO 60

                     ENDIF

                  ENDIF
               ENDDO

            ENDDO

C           TOUS LES TETRAEDRES RESTENT A VOLUME>0
C           LE SOMMET NSDMIN EST DEPLACE AU POINT NPt
C           DEVENANT AINSI UN SOMMET DES TETRAEDRES ACTUELS
C           -----------------------------------------------
            DO 58 I=1,NBTEET
C              RECHERCHE DU NUMERO DE SOMMET DE NSDMIN
               NTE = NOTEET( I )
               DO K=1,4
                  IF( NOTETR(K,NTE) .EQ. NSDMIN ) THEN
                      NOTETR(K,NTE) = NPt
                     GOTO 58
                  ENDIF
               ENDDO
 58         ENDDO

C           IDENTIFICATION REALISEE
            N1TETS( NPt ) = NOTEET(1)

            KTITRE='etoipt19 : REALISATION          ->          ,         
     %             TETRAEDRES'
            WRITE(KTITRE(24:31),'(I8)') NSDMIN
            WRITE(KTITRE(36:43),'(I8)') NPt
            WRITE(KTITRE(47:54),'(I8)') NBTEET
            CALL TRFETO13( KTITRE, PTXYZD, NBTEET, NOTEET, NOTETR )

ccc            PRINT*,'etoipt19 55: SOMMET',NSDMIN,' XYZ=',
ccc     %             (PTXYZD(I,NSDMIN),I=1,4),' DEPLACE au'
ccc            PRINT*,'etoipt19 55: POINT ',NPt,' XYZ=',
ccc     %             (PTXYZD(I,NPt),I=1,4)

C           TRACE FILAIRE DES ARETES DES NBTEET TETRAEDRES DE L'ETOILE
            GOTO 1100

         ENDIF

      ENDIF

C     2) LISTAGE DES TETRAEDRES INITIAUX CONTENANT LE POINT NPt PAR LA
C     RECHERCHE DU NOMBRE DE COORDONNEES BARYCENTRIQUES FAIBLES
C     POUR IDENTIFIER SI LE POINT NPt EST SUR UNE ARETE OU UNE FACE
C     ET AJOUT DES TETRAEDRES OPPOSES AUX FACES ou AUTOUR DES ARETES
C     ================================================================

C     CALCUL DES COORDONNEES BARYCENTRIQUES DU PTXYZD(1,NPt)
C     DANS LE TETRAEDRE NOTET1
 60   CALL COBATE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NOTET1),
     %             VTE, COBARY, IERR )
C     COBARY(I) EST ICI LA COORDONNEE BARYCENTRIQUE OPPOSEE A LA FACE I

ccc      PRINT*,'etoipt19 60: NPt=',NPt,' CONTENU DANS NOTETR(',NOTET1,
ccc     %       ')=', (NOTETR(I,NOTET1),I=1,8)
ccc      PRINT*,'etoipt19 60: NPt=',NPt,' NOTET1=',NOTET1,' VTE=',VTE,
ccc     %       ' COBARY=',COBARY

      IMAX = 0
      CMAX   = -1D100
      NBCBFA = 0
      DO I=1,4

CCC         IF( ABS(COBARY(I)) .LE. 0.01D0    ) THEN  24/10/2005
CCC         IF( ABS(COBARY(I)) .LE. 0.001D0   ) THEN  13/08/2014
CCC         IF(     COBARY(I)  .LE. 0.00001D0 ) THEN  17/08/2014
CCC         IF(     COBARY(I)  .LE. 0.D0      ) THEN  16/09/2014
ccc         IF( ABS(COBARY(I)) .LE. 1D-5      ) THEN  19/01/2017
ccc         IF( ABS(COBARY(I)) .LE. 1D-3      ) THEN  21/01/2017

         D = COBARY( I )
         IF( D .LE. 0.05D0 ) THEN
            NBCBFA = NBCBFA + 1
C           LE NUMERO DE LA COORDONNEE BARYCENTRIQUE QUASI NULLE
            NOTEET( NBCBFA ) = I
         ENDIF

         IF( D .GT. CMAX )THEN
            CMAX = D
            IMAX = I
         ENDIF

      ENDDO
C
C     TRAITEMENT SELON LE NOMBRE DE COORDONNEES BARYCENTRIQUES FAIBLES
 61   IF( NBCBFA .EQ. 0 ) THEN

C        PTXYZD(NPt) FRANCHEMENT INTERNE AU TETRAEDRE NOTET1
C        ---------------------------------------------------
         NBTEET0 = 1
         NBTEET  = 1
         NOTEET(1) = NOTET1

C        SI NPt EST PROCHE D'UN SOMMET (IMAX-1)
C        ESSAI D'AJOUTER LES TETRAEDRES OPPOSES AUX 3 FACES
C        NON OPPOSEE AU SOMMET DE COORDONNEE BARYCENTRIQUE MAXIMALE
C        LA FACE OPPOSEE EST IMAX
         IF( CMAX .GT. 0.6D0 ) THEN
            DO I=1,4
               IF( I .NE. IMAX ) THEN
                  NTOP = NOTETR( 4+I, NOTET1 )
                  IF( NTOP .GT. 0 ) THEN
                     IF( NOTETR(1,NTOP) .GT. 0 ) THEN
                        NBTEET = NBTEET + 1
                        NOTEET( NBTEET ) = NTOP
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         GOTO 63

      ELSE IF( NBCBFA .EQ. 1 ) THEN

C        PTXYZD(NPt) SUR UNE FACE DE TETRAEDRE OPPOSE NTOP
C        -------------------------------------------------
         NTOP = NOTETR( 4+NOTEET(1), NOTET1 )
         IF( NTOP .LE. 0 ) THEN
            NBTEET0 = 1
            NBTEET  = 1
            NOTEET(1) = NOTET1
            GOTO 63
         ENDIF

C        VOLUME ET QUALITE DU TETRAEDRE OPPOSE NTOP
         CALL QV1TET( PTXYZD, NOTETR(1,NTOP), VTE, QTE )
         IF( QTE .LT. QTEAME ) THEN

C           LE POINT NPt EST IL PROCHE D'UNE ARETE DE NOTET1?
            DO I=1,4

               IF( I .NE. NOTEET(1) ) THEN
ccc                  IF( COBARY(I) .LE. 0.0001D0 ) THEN 31/3/2016
                  IF( COBARY(I) .LE. 0.05D0 ) THEN
C                    OUI: LES TETRAEDRES AYANT CETTE ARETE SONT AJOUTES
                     IF( I .GT. NOTEET(1) ) THEN
                        NOTEET(2) = I
                     ELSE
                        NOTEET(2) = NOTEET(1)
                        NOTEET(1) = I
                     ENDIF
                     NBCBFA = 2
                     GOTO 61
                  ENDIF
               ENDIF

            ENDDO

C           LE TETRAEDRE AJOUTE NTOP EST DE MAUVAISE QUALITE
C           LUI ET SES TETRAEDRES OPPOSES SONT AJOUTES A L'ETOILE
C           POUR TENTER DE L'ELIMINER
            NBTEET0 = 2
            NBTEET  = 2
            NOTEET(1) = NOTET1
            NOTEET(2) = NTOP
            DO K=1,4
               NTE = NOTETR( 4+K, NTOP )
               IF( NTE .NE. NOTET1 .AND. NTE .GT. 0 ) THEN
                  NBTEET = NBTEET + 1
                  NOTEET( NBTEET ) = NTE
               ENDIF
            ENDDO
            GOTO 63

         ENDIF

C        LE TETRAEDRE OPPOSE A CETTE FACE EST AJOUTE A L'ETOILE
         NOTEET(1) = NOTET1
         NOTEET(2) = NTOP
         NBTEET0 = 2
         NBTEET  = 2
         GOTO 63

      ELSE IF( NBCBFA .EQ. 2 ) THEN

C        PTXYZD(NPt) SUR UNE ARETE: LES TETRAEDRES QUI PARTAGENT CETTE
C                                   ARETE SONT AJOUTES A L'ETOILE
C        -------------------------------------------------------------

C        NOTEET(3:4) LES NUMEROS DES 2 SOMMETS DE L'ARETE
         IF( NOTEET(1) .EQ. 1 ) THEN
            IF( NOTEET(2) .EQ. 2 ) THEN
               NOTEET(3) = 2
               NOTEET(4) = 3
            ELSE IF( NOTEET(2) .EQ. 3 ) THEN
               NOTEET(3) = 1
               NOTEET(4) = 3
            ELSE IF( NOTEET(2) .EQ. 4 ) THEN
               NOTEET(3) = 1
               NOTEET(4) = 2
            ENDIF
         ELSE IF( NOTEET(1) .EQ. 2 ) THEN
            IF( NOTEET(2) .EQ. 3 ) THEN
               NOTEET(3) = 3
               NOTEET(4) = 4
            ELSE IF( NOTEET(2) .EQ. 4 ) THEN
               NOTEET(3) = 2
               NOTEET(4) = 4
            ENDIF
         ELSE
            NOTEET(3) = 1
            NOTEET(4) = 4
         ENDIF

         N1 = NOTETR( NOTEET(3), NOTET1 )
         N2 = NOTETR( NOTEET(4), NOTET1 )
         CALL TETR1A( N1, N2, N1TETS, NOTETR,
     %                NBTEET, MXETOI, NOTEET, IERR )
         NBTEET0 = NBTEET

C        TRACE FILAIRE DES ARETES DES NBTEET TETRAEDRES DE L'ETOILE
         KTITRE='etoipt19 : NPt                    TETRAEDRES AUTOUR ARE
     %TE                '
         WRITE(KTITRE(16:23),'(I8)') NPt
         WRITE(KTITRE(25:28),'(I4)') NBTEET
         WRITE(KTITRE(59:66),'(I8)') N1
         WRITE(KTITRE(68:75),'(I8)') N2
         CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )

C        EXISTE T IL DES TETRAEDRES DE MAUVAISE QUALITE?
         DO L = 1, NBTEET0

            NTE = NOTEET( L )

C           VOLUME ET QUALITE DU TETRAEDRE OPPOSE
            CALL QV1TET( PTXYZD, NOTETR(1,NTE), VTE, QTE )

            IF( QTE .LT. QTEAME ) THEN

C              NTE EST DE MAUVAISE QUALITE
C              AJOUT DES TETRAEDRES OPPOSES POUR L'ELIMINER
               DO 62 K=1,4
                  NTOP = NOTETR( 4+K, NTE )
                  IF( NTOP .GT. 0 ) THEN
C                    VOLUME ET QUALITE DU TETRAEDRE OPPOSE
                     CALL QV1TET( PTXYZD, NOTETR(1,NTOP), VTE, QTE )
                     IF( QTE .GE. QTEAME ) THEN
                        DO M=1,NBTEET
                           IF( NTOP .EQ. NOTEET(M) ) THEN
                              GOTO 62
                           ENDIF
                        ENDDO
C                       NTOP EST AJOUTE
                        NBTEET = NBTEET + 1
                        NOTEET( NBTEET ) = NTOP
                     ENDIF
                  ENDIF
 62            ENDDO

            ENDIF

         ENDDO

         GOTO 63
C
      ELSE

C        PTXYZD(NPt) SUR UN SOMMET : IL EST A ABANDONNER SI POSSIBLE
C        -----------------------------------------------------------
C        NUMERO DU SOMMET DANS NOTETR A IDENTIFIER
         DO I=1,4
ccc            IF( COBARY(I) .GT. 0.949D0 ) THEN     15/4/2018
            IF( COBARY(I) .GT. 0.92D0 ) THEN
C              NUMERO LOCAL DU SOMMET DE COBARY=1
               IF( I .NE. 1 ) THEN
                  NSP = I-1
               ELSE
                  NSP = 4
               ENDIF
ccc               WRITE(IMPRIM,*)
ccc               WRITE(IMPRIM,*) 'etoipt19 : POINT ',NPt,
ccc     %         ' TRES PROCHE du SOMMET',NSP,NOTETR(NSP,NOTET1),
ccc     %                   ' DU TETRAEDRE',NOTET1
ccc               WRITE(IMPRIM,*) 'etoipt19 : COOR BARYCENTRIQUES=',COBARY
ccc               WRITE(IMPRIM,*) 'PT : XYZ',NPt,'=',(PTXYZD(L,NPt),L=1,4)
ccc               NNS = NOTETR(1,NOTET1)
ccc               WRITE(IMPRIM,*) 'ST1: XYZ',NNS,'=',PTXYZD(1,NNS),
ccc     %                    PTXYZD(2,NNS),PTXYZD(3,NNS)
ccc               NNS = NOTETR(2,NOTET1)
ccc               WRITE(IMPRIM,*) 'ST2: XYZ',NNS,'=',PTXYZD(1,NNS),
ccc     %                    PTXYZD(2,NNS),PTXYZD(3,NNS)
ccc               NNS = NOTETR(3,NOTET1)
ccc               WRITE(IMPRIM,*) 'ST3: XYZ',NNS,'=',PTXYZD(1,NNS),
ccc     %                    PTXYZD(2,NNS),PTXYZD(3,NNS)
ccc               NNS = NOTETR(4,NOTET1)
ccc               WRITE(IMPRIM,*) 'ST4: XYZ',NNS,'=',PTXYZD(1,NNS),
ccc     %                    PTXYZD(2,NNS),PTXYZD(3,NNS)

C              RECHERCHE DE TOUS LES TETRAEDRES DE SOMMET NSDMIN
               NSDMIN = NOTETR( NSP, NOTET1 )
               CALL TETR1S( NSDMIN, N1TETS, NOTETR,
     %                      NBTEET, MXETOI, NOTEET, IERR )
               IF( NBTEET .LE. 0 ) THEN
                  NBTEET = 1
                  NOTEET(1) = NOTET1
               ENDIF
               NBTEET0 = NBTEET
               GOTO 63

            ENDIF
         ENDDO

      ENDIF


C     ==================================================================
C     ICI: L'ETOILE INITIALE EST CONSTITUEE DE NBTEET0 TETRAEDRES NOTEET
C     ==================================================================
ccc 63   print*,'etoipt19 63: NOTEET(1:',NBTEET,')=',(noteet(i),i=1,nbteet)

C     SAUVEGARDE DANS NOTERE(1:NBTERE) POUR UNE REPRISE A PARTIR DE CES
C     TETRAEDRES (NOTEET(1:NBTEET0) CONTENANT PRATIQUEMENT LE POINT NPt
 63   NBTERE = NBTEET0
      DO K = 1, NBTEET0
         NOTERE( K ) = NOTEET( K )
      ENDDO

ccc      WRITE(IMPRIM,*) 'etoipt19 : POINT ',NPt,' NBTERE=',NBTERE,
ccc     %      ' TETRAEDRES INITIAUX:',(NOTERE(K),K=1,NBTERE)

C     TRACE FILAIRE DES ARETES DES NBTEET TETRAEDRES DE L'ETOILE
      KTITRE='etoipt19 : POINT                    TETRAEDRES INITIAUX DE
     % L''ETOILE'
      WRITE(KTITRE(17:24),'(I8)') NPt
      WRITE(KTITRE(27:34),'(I8)') NBTEET
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )

      IF( NBTEET .LE. 0 ) THEN

C        TRAITEMENT DES CAS A PROBLEME:
C        ------------------------------
C        LE POINT NPt PEUT IL ETRE SUPPRIME?
         NPSNPt = NPSOFR(NPt)
         IF( NPSNPt .EQ. -4 .OR. NPSNPt .EQ. -3 .OR. NPSNPt .EQ. 0 )THEN

C           OUI: LE POINT NPt EST SUPPRIME DES POINTS A TETRAEDRISER
ccc            WRITE(IMPRIM,*) 'etoipt19 : POINT ',NPt,
ccc     %                      ' SUPPRIME. NPSOFR(NPt)=',NPSNPt,
ccc     %                      '  XYZD=', (PTXYZD(I,NPt),I=1,4)
            NBTEET = -1
            NPSOFR(NPt) = -3
            N1TETS(NPt) = 0
            GOTO 9999

         ENDIF

C        NON: LE POINT NPt NE PEUT ETRE SUPPRIME
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,*) 'etoipt19 : NPt=',NPt,
     %                   ' NON SUPPRIMABLE NPSOFR(NPt)=',NPSNPt,
     %                   ' NBTEET=',NBTEET,
     %                   ' XYZ=',(PTXYZD(I,NPt),I=1,4)

         IF( NBCBFA .EQ. 2 ) THEN
            GOTO 500
         ENDIF

C        ABANDON DU POINT NPt NON SUPPRIMABLE
         NBTEET = -101
         GOTO 9999

      ENDIF


C     CONSTRUCTION DES FACES SIMPLES DES TETRAEDRES CONTENANT LE POINT NPt
C     --------------------------------------------------------------------
C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
 64   CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

      NBTEET1 = 0
      VOLET0  = 0D0
      DO I=1,NBTEET

C        NUMERO NOTETR DU TETRAEDRE DE L'ETOILE DU POINT NPt
         NTE = NOTEET( I )
         IF( NTE .GT. 0 ) THEN

C           VOLUME DE L'ETOILE DES TETRAEDRES ACTUELS
            VTE = VOLTET( PTXYZD( 1, NOTETR(1,NTE) ),
     %                    PTXYZD( 1, NOTETR(2,NTE) ),
     %                    PTXYZD( 1, NOTETR(3,NTE) ),
     %                    PTXYZD( 1, NOTETR(4,NTE) ) )

ccc            IF( QTE .LT. QTEAME ) THEN
cccC              LE TETRAEDRE NTE EST DE QUALITE MEDIOCRE MAIS CONSERVE
ccc               PRINT*,'etoipt19 : TETRAEDRE de NPt',NPt,NTE,
ccc     %         ' de VOLUME',VTE,' et de QUALITE',QTE,' est CONSERVE'
ccc            ENDIF

            VOLET0 = VOLET0 + VTE

            DO J=1,4
C              AJOUT DE LA FACE J DE NTE AUX FACES DE L'ETOILE
               CALL AJFETO( NTE, J, NOTETR, N1FEOC, N1FEVI, NFETOI, NF )
            ENDDO

            NBTEET1 = NBTEET1 + 1
            NOTEET( NBTEET1 ) = NTE

         ENDIF

      ENDDO
      NBTEET = NBTEET1

C     SUPPRESSION DE L'ETOILE DES TETRAEDRES DONT UNE FACE SIMPLE+NPt
C     CONSTITUE UN TETRAEDRE DE MAUVAISE QUALITE
C     ---------------------------------------------------------------
      NF1    = N1FEOC
      NBTESU = 0
      NBFETO = 0
      QUAMIN = 2.
      QUAMOY = 0.

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 65   IF( NF1 .GT. 0 ) THEN

C        UNE FACE DE L'ETOILE EN PLUS
         NBFETO = NBFETO + 1

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NTE = NFETOI(1,NF1)

C        PAS DE SUPPRESSION POSSIBLE DU TETRAEDRE CONTENANT NPt
         IF( NTE .EQ. NOTET1 ) GOTO 66

C        NO DE LA FACE DANS LE TETRAEDRE NTE
         I = NFETOI(2,NF1)

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NOSOTE(1) = NOTETR( NOSOFATE(1,I), NTE )
         NOSOTE(2) = NOTETR( NOSOFATE(2,I), NTE )
         NOSOTE(3) = NOTETR( NOSOFATE(3,I), NTE )
         NOSOTE(4) = NPt

C        LE TETRAEDRE A FORMER NF1-POINT NPt  POSE T IL PROBLEME ?
C        LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C        AU TETRAEDRE NT ET CARRE DE SON RAYON
         CALL CEBOQU( 0., NOSOTE, PTXYZD, CENTRE,  VTE, QTE )

         QUAMOY = QUAMOY + QTE
         IF( QTE .LT. QUAMIN ) QUAMIN = QTE

         IF( QTE .LT. QTEAME ) THEN

C           LE TETRAEDRE NF1-POINT NPt EST DE QUALITE MEDIOCRE
C           => LE TETRAEDRE NTE EST SUPPRIME DE L'ETOILE INITIALE
C              SI CE N'EST PAS LE TETRAEDRE CONTENANT NPt
            DO L = 1, NBTEET

               IF( NTE .EQ. ABS( NOTEET(L) ) ) THEN
C                 LE TETRAEDRE NTE EST SUPPRIME DE L'ETOILE
ccc                  PRINT*,'etoipt19 66: NPt=',NPt,
ccc     %                 ' RETRAIT du TETRAEDRE INITIAL NOTETR(',NTE,')=',
ccc     %                  (NOTETR(KK,NTE),KK=1,8),' VTE=',VTE
                  IF( NOTEET(L) .GT. 0 ) THEN
                     NBTESU = NBTESU + 1
                     NOTEET(L) = -NTE
                  ENDIF
                  GOTO 66
               ENDIF

            ENDDO

C           ANOMALIE SI ARRIVEE ICI
            PRINT*,'etoipt19 : Probleme TETRAEDRE',NTE,
     %             ' NON dans NOTEET'
            PRINT*,'NOTEET :',(NOTEET(L),L=1,NBTEET)

         ENDIF

C        TETRAEDRE CORRECT . PASSAGE A LA FACE SIMPLE SUIVANTE
 66      NF1 = NFETOI( 5, NF1 )
         GOTO 65

      ENDIF
      IF( NBFETO .GT. 0 ) QUAMOY = QUAMOY / NBFETO

      IF( NBTESU .GT. 0 ) THEN
C        RECONSTRUCTION DE L'ETOILE DES TETRAEDRES
         NBTEET1 = 0
         DO L=1,NBTEET
            NTE = NOTEET( L )
            IF( NTE .GT. 0 ) THEN
               NBTEET1 = NBTEET1 + 1
               NOTEET( NBTEET1 ) = NTE
            ENDIF
         ENDDO
         NBTEET = NBTEET1
         IF( NBTEET .LE. 0 ) THEN
            GOTO 500
         ELSE
            GOTO 64
         ENDIF
      ENDIF

ccc      PRINT*,'etoipt19 : AVANT EXTENSION NPt=',NPt,' NBFETO=',NBFETO,
ccc     %' NBTEET=',NBTEET,' QUAMIN=',QUAMIN,' QUAMOY=',QUAMOY
ccc      KTITRE='etoipt19 : POINT                     TETRAEDRES ETOILE INIT
ccc     %IALE avant EXTENSION'
ccc      WRITE(KTITRE(17:24),'(I8)') NPt
ccc      WRITE(KTITRE(26:32),'(I7)') NBTEET
ccc      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )


C     EXTENSION  DE L'ETOILE AUX TETRAEDRES OPPOSES PAR LES FACES DES NBTEET
C     TETRAEDRES DE L'ETOILE ET DONT LES NOUVELLES ARETES SIMPLES DE L'ETOILE
C     GENERENT AVEC LE POINT NPt DES TETRAEDRES DE VOLUME>0
C     (LE CRITERE DE LA BOULE CIRCONSCRITE N'EST PAS ICI PRIS EN COMPTE)
C     =======================================================================
 68   NF1 = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 69   IF( NF1 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NTE = NFETOI(1,NF1)
         I   = NFETOI(2,NF1)

         IF( INFACO .NE. 0 ) THEN
C           LA FACE I DE NTE EST ELLE DANS LEFACO?
            CALL NULEFT( I, NTE, NOTETR, MXFACO, LEFACO, NFLE )
            IF( NFLE .GT. 0 ) THEN

C              OUI: LA FACE EST DANS LEFACO
C              CONSERVATION DE CETTE FACE LEFACO EN NE PASSANT PAS AU DELA
               GOTO 72

cccC              AJOUT DE -NTOP A NTETOI
cccC              POUR EVITER DE LA RETROUVER LORS D'UN AUTRE PARCOURS
ccc               NBTEET = NBTEET + 1
ccc               NTETOI( NBTEET ) = - NTOP

            ENDIF
         ENDIF

         IF( IVOLTE .NE. 0 ) THEN
C           LA FACE I DE NTE EST ELLE A L'INTERFACE DE 2 VOLUMES?
C           LE TETRAEDRE OPPOSE APPARTIENT IL A UN AUTRE VOLUME
            NTEOP = NOTETR( 4+I, NTE )
            IF( NTEOP .GT. 0 ) THEN
C              COMPARAISON DES 2 VOLUMES
               IF( NVOLTE(NTE) .NE. NVOLTE(NTEOP) ) THEN
C                 CONSERVATION DE CETTE FACE
                  GOTO 72
               ENDIF
            ENDIF
         ENDIF

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NOSOTR(1) = NOTETR( NOSOFATE(1,I), NTE )
         NOSOTR(2) = NOTETR( NOSOFATE(2,I), NTE )
         NOSOTR(3) = NOTETR( NOSOFATE(3,I), NTE )

C        LE TETRAEDRE OPPOSE A CETTE FACE I DE NTE
         NTOP = NOTETR( 4+I, NTE )

C        CE TETRAEDRE NTOP PEUT IL ETRE AJOUTE A L'ETOILE DU POINT NPt?
         IF( NTOP .GT. 0 ) THEN

C           NTOP EST IL DEJA RECENSE DANS NOTEET?
            DO N=1,NBTEET
               IF( NOTEET(N) .EQ. NTOP ) THEN
                  GOTO 72
               ENDIF
            ENDDO

C           POUR VOIR SI NPt EST INTERNE A LA BOULE CIRCONSCRITE DE NTOP
            CALL PTDSBO( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NTOP),
     %                   NONOUI, VTE, DSR )

            IF( NONOUI .EQ. 0 ) THEN
C              NTOP NE VERIFIE PAS LE CRITERE DELAUNAY DE LA BOULE VIDE
               GOTO 72
            ENDIF

cccC              CE TETRAEDRE NTOP A T IL UN SOMMET TROP ELOIGNE DE NPt?
ccc               DO N=1,4
ccc                  NST = NOTETR( N, NTOP )
ccc                  D = ( PTXYZD(1,NST) - PTXYZD(1,NPt) ) **2
ccc     %              + ( PTXYZD(2,NST) - PTXYZD(2,NPt) ) **2
ccc     %              + ( PTXYZD(3,NST) - PTXYZD(3,NPt) ) **2
cccccc                  IF( D .GT. 1.6D0 * PTXYZD(4,NPt) **2 ) THEN
ccc                  IF( D .GT. ARETGD2 ) THEN
cccC                    UN SOMMET DE NTOP TROP ELOIGNE DE NPt
cccC                    LE TETRAEDRE NTOP N'EST PAS AJOUTE A L'ETOILE
ccc                     GOTO 72
ccc                  ENDIF
ccc               ENDDO

C           CE TETRAEDRE NTOP EST AJOUTE A NOTEET
C           -------------------------------------
            IF( NBTEET .GE. MXETOI ) THEN
C              LIMITATION DU NOMBRE DE TETRAEDRES DE L'ETOILE
C              JUSTE AVANT LA SATURATION DU TABLEAU NOTEET
               PRINT*,'etoipt19 : SATURATION NOTEET MXETOI=',MXETOI
               TRACTE = .TRUE.
               GOTO 73
            ENDIF
            NBTEET = NBTEET + 1
            NOTEET( NBTEET ) = NTOP

ccc            PRINT*,'etoipt19 : NPt',NPt,' opp+tetra',NTE,' St:',
ccc     %          (NOTETR(K,NTE),K=1,8),' Q=',QTE,' V=',VTE,' D+=',D0

C           AJOUT DES 4 FACES DE NTOP AUX FACES DE L'ETOILE NFETOI
            DO N=1,4
C              AJOUT DE LA FACE N DE NTOP AUX FACES DE L'ETOILE
               CALL AJFETO( NTOP, N, NOTETR, N1FEOC, N1FEVI, NFETOI,
     %                      NF )
            ENDDO

C           RETOUR AU DEBUT DU PARCOURS DES FACES SIMPLES DE L'ETOILE
            GOTO 68

         ENDIF

C        TETRAEDRE CORRECT . PASSAGE A LA FACE SIMPLE SUIVANTE
 72      NF1 = NFETOI( 5, NF1 )
         GOTO 69

      ENDIF

 73   KTITRE='etoipt19 : POINT                     TETRAEDRES apres AJOU
     %T des TETRAEDRES OPPOSES'
      WRITE(KTITRE(18:25),'(I8)') NPt
      WRITE(KTITRE(27:33),'(I7)') NBTEET
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )


C     VOLUME DE L'ETOILE DES NBTEET TETRAEDRES A ELIMINER
C     FORMATION DES FACES SIMPLES DE LEUR ETOILE
C     ===================================================
C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
 200  CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

      NBTEET1= 0
      VOLET0 = 0D0
      QUAMIN = 2.
      QUAMOY = 0.

      DO K = 1, NBTEET

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE
         NTE = NOTEET( K )
         IF( NTE .GT. 0 ) THEN

C           VOLUME ET QUALITE DU TETRAEDRE NTE
            CALL QV1TET( PTXYZD, NOTETR(1,NTE), VTE, QTE )

C           VOLUME DES NBTEET TETRAEDRES ACTUELS DE L'ETOILE
            VOLET0 = VOLET0 + VTE

            QUAMOY = QUAMOY + QTE
            IF( QTE .LT. QUAMIN ) QUAMIN = QTE

ccc            PRINT *,'etoipt19 : NPt',NPt,I,' old tetra',NTE,' St:',
ccc     %          (NOTETR(K,NTE),K=1,8),' Q=',QTE,' V=',VTE

            NBTEET1 = NBTEET1 + 1
            NOTEET( NBTEET1 ) = NTE

            DO J=1,4
C              AJOUT DE LA FACE J DE NTE AUX FACES DE L'ETOILE
               CALL AJFETO( NTE, J, NOTETR, N1FEOC, N1FEVI, NFETOI, NF )
            ENDDO

         ENDIF
      ENDDO
      NBTEET = NBTEET1
      IF( NBTEET .GT. 0 ) QUAMOY = QUAMOY / NBTEET

      IF( VOLET0 .LE. 0D0 ) THEN
C        ETOILE VIDE
ccc         WRITE(IMPRIM,*) 'etoipt19 : au POINT',NPt,' NBTEET=',NBTEET,
ccc     %   ' ETOILE de VOLUME=',VOLET0,' QUAMIN=',QUAMIN,' QUAMOY=',QUAMOY
         GOTO 500
      ENDIF

C     ICI IL RESTE AU MOINS UN TETRAEDRE CONTENANT LE POINT NPt
      VTMIN = VOLET0 / NBTEET * EPSDEG

C     BILAN : LES TETRAEDRES A FORMER FACE ETOILE-POINT NPt SERAIENT ILS
C             DEGENERES OU RENVERSES?
C             LE VOLUME DES TETRAEDRES SERAIT IL LE MEME AVANT ET APRES?
C     ==================================================================
      VOLET1 = 0D0
      NF1    = N1FEOC
      NBFETO = 0
ccc      QUAMIN = 2.
ccc      QUAMOY = 0.

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 220  IF( NF1 .GT. 0 ) THEN

C        UNE FACE DE L'ETOILE EN PLUS
         NBFETO = NBFETO + 1

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NTE = NFETOI(1,NF1)
         I   = NFETOI(2,NF1)

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NOSOTE(1) = NOTETR( NOSOFATE(1,I), NTE )
         NOSOTE(2) = NOTETR( NOSOFATE(2,I), NTE )
         NOSOTE(3) = NOTETR( NOSOFATE(3,I), NTE )
         NOSOTE(4) = NPt

C        LE TETRAEDRE A FORMER NF1-POINT NPt  POSE T IL PROBLEME ?
C        LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C        AU TETRAEDRE NT ET CARRE DE SON RAYON
         CALL CEBOQU( 0., NOSOTE, PTXYZD, CENTRE, VTE, QTE )
         IF( VTE .GT. 0D0 ) THEN
ccc         IF( QTE .GT. 0.001 ) THEN
            GOTO 240
         ENDIF

ccc         IF( VTE .GT. 0D0 .AND. QTE .GE. QTEAME ) THEN
ccc            QUAMOY = QUAMOY + QTE
ccc            IF( QTE .LT. QUAMIN ) QUAMIN = QTE
cccC           DISTANCE(NPt-CENTRE)**2 / RAYON**2
ccc            D0 = ( ( CENTRE(1) - PTXYZD(1,NPt) ) ** 2 +
ccc     %             ( CENTRE(2) - PTXYZD(2,NPt) ) ** 2 +
ccc     %             ( CENTRE(3) - PTXYZD(3,NPt) ) ** 2 ) / CENTRE(4)
ccc            IF( D0 .LE. 1.001D0 ) THEN
cccC              NPt EST INTERNE A LA BOULE DU TETRAEDRE NF1-POINT NPt
ccc               GOTO 240
ccc            ENDIF
ccc         ENDIF

C        SUPPRESSION DU TETRAEDRE NTE DE L'ETOILE
         DO K = 1, NBTEET
            IF( NOTEET(K) .EQ. NTE ) THEN
C              LE TETRAEDRE NTE EST SUPPRIME DE L'ETOILE
ccc               PRINT*,'etoipt19 : NPt=',NPt,
ccc     %                ' RETRAIT du TETRAEDRE NOTETR(',NTE,')=',
ccc     %                (NOTETR(KK,NTE),KK=1,4),' car NOSOTE=',NOSOTE,
ccc     %                ' avec VTE=',VTE,' QTE=',QTE
               NOTEET(K) = -NTE
               GOTO 200
            ENDIF
         ENDDO

C        AJOUT DU VOLUME DU TETRAEDRE AU VOLUME DE L'ETOILE
 240     VOLET1 = VOLET1 + VTE

C        TETRAEDRE CORRECT . PASSAGE A LA FACE SIMPLE SUIVANTE
         NF1 = NFETOI( 5, NF1 )
         GOTO 220

      ENDIF
      IF( NBFETO .GT. 0 ) QUAMOY = QUAMOY / NBFETO

ccc      PRINT *,'etoipt19 : NPt=',NPt,' NBFETO=',NBFETO,
ccc     %        ' VOLET0=',VOLET0,' VOLET1=',VOLET1


C     LES 2 VOLUMES DE L'ETOILE AVANT APRES DU POINT NPt SONT ILS EGAUX?
C     ------------------------------------------------------------------
      IF( ABS(VOLET1-VOLET0) .LE. ABS(VOLET0)*1D-5 ) THEN

C        OUI: CONSTRUCTION DES TETRAEDRES DE L'ETOILE DU POINT NPt
         NBTEETIN = NBTEET
         CALL TETETO( NPt,   REPRISE, QTEAME, MXSOMM, PTXYZD, NPSOFR,
ccc     %                INFACO, MXFACO, LEFACO,
     %                NOVOLU, IVOLTE, NVOLTE,
     %                MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                MXETOI, N1FEOC, NFETOI, VOLET1,
     %                MXETOI, NBTEET, NOTEET,
     %                MXSTPE, NBSTPE, NOSTPE, IERR )

         IF( IERR .NE. 0 ) THEN
            GOTO 500
         ELSE
C           TRACE DES TETRAEDRES DE L'ETOILE CORRECTE DE NPt
            GOTO 1100
         ENDIF

      ENDIF

C     =======================================================================
C     REPRISE DE L'ETOILE AVEC UNE AUTRE METHODE EN PRENANT
C     L'ETOILE MINIMALE DES NBTERE TETRAEDRES INITIAUX CONTENANT LE POINT NPt
C     LE CRITERE DE LA BOULE VIDE DE DELAUNAY N'EST PLUS APPLIQUE
C     =======================================================================
 500  REPRISE = .TRUE.
      TRACTE  = .TRUE.

ccc      PRINT *
ccc      PRINT *,'************* REPRISE NPt=',NPt,' **********************'

      IF( NPSOFR(NPt) .GT. 0 ) THEN

         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10500) NPt,(PTXYZD(I,NPt),I=1,4),NPSOFR(NPt),
     %                          NBTERE
         ELSE
            WRITE(IMPRIM,20500) NPt,(PTXYZD(I,NPt),I=1,4),NPSOFR(NPt),
     %                          NBTERE
         ENDIF
10500 FORMAT(' etoipt19 : ***** REPRISE du POINT ',I8,' XYZD=',4G15.7,
     %       ' NPSOFR=',I9,' avec',I5,' TETRAEDRES INITIAUX')
20500 FORMAT(' etoipt19 : ***** RENEWAL of POINT ',I8,' XYZD=',4G15.7,
     %       ' NPSOFR=',I9,' with',I5,' INITIAL TETRAHEDRA')

      ELSE

C        ABANDON DE LA TETRAEDRISATION DU POINT NPt NON FRONTIERE
C        CE POINT NE SERA PAS TETRAEDRISE
         NBTEET = -1
         IERR   = 11
         GOTO 9999

      ENDIF

C     REPRISE AVEC LES TETRAEDRES NOTERE INITIAUX CONTENANT LE POINT NPt
 510  NBTEET = NBTERE
      DO K = 1, NBTEET
         NOTEET( K ) = NOTERE( K )
      ENDDO

C     TRACE DES NBTEET TETRAEDRES DE REPRISE DE L'ETOILE
      tracte = .true.
      KTITRE='etoipt19 : POINT                    TETRAEDRES en ETOILE M
     %INIMALE de REPRISE'
      WRITE(KTITRE(18:25),'(I8)') NPt
      WRITE(KTITRE(28:35),'(I8)') NBTERE
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )

C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
      CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

C     FORMATION DES FACES SIMPLES DE L'ETOILE
      VOLET0 = 0D0
      DO I=NBTEET,1,-1
C        LE NUMERO DU TETRAEDRE DE FACES A TRAITER
         NTE = NOTEET( I )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

            CALL CEBOQU( 0., NOTETR(1,NTE), PTXYZD, CENTRE, VTE, QTE )

C           VOLUME DE L'ETOILE
            VOLET0 = VOLET0 + VTE

C           DISTANCE(NPt-CENTRE)**2 / RAYON**2
            D = ( ( CENTRE(1) - PTXYZD(1,NPt) ) ** 2 +
     %            ( CENTRE(2) - PTXYZD(2,NPt) ) ** 2 +
     %            ( CENTRE(3) - PTXYZD(3,NPt) ) ** 2 ) / CENTRE(4)

            PRINT*,'etoipt19 : REPRISE de NPt',NPt,' initial TETRA,',NTE
     %        ,' St:',(NOTETR(K,NTE),K=1,4),' Q=',QTE,' V=',VTE,' D1=',D

C           AJOUT DES FACES DE NTE
            DO K=1,4
               CALL AJFETO( NTE, K, NOTETR, N1FEOC, N1FEVI, NFETOI, L )
            ENDDO

         ENDIF
      ENDDO

C     LE VOLUME DES TETRAEDRES FORMES A PARTIR DE NPt ET
C     DES FACES SIMPLES EST IL TOUJOURS CORRECT ?
C     --------------------------------------------------
      NF1 = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 550  IF( NF1 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NTE = NFETOI(1,NF1)

C        PAS DE SUPPRESSION POSSIBLE DU TETRAEDRE CONTENANT NPt
         IF( NTE .EQ. NOTET1 ) GOTO 580

C        NO DE LA FACE DE NTE
         NFNT = NFETOI(2,NF1)

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NS1 = NOTETR( NOSOFATE(1,NFNT), NTE )
         NS2 = NOTETR( NOSOFATE(2,NFNT), NTE )
         NS3 = NOTETR( NOSOFATE(3,NFNT), NTE )
C
C        LE TETRAEDRE A FORMER NF1-POINT NPt POSE T IL PROBLEME ?
         CALL QUATETD( PTXYZD(1,NS1),
     %                 PTXYZD(1,NS2),
     %                 PTXYZD(1,NS3),
     %                 PTXYZD(1,NPt),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         PRINT *,'etoipt19 : REPRISE de NPt',NPt,' new+ tetra St:',
     %            NS1,NS2,NS3,NPt,' Q=',QTE,' V=',VTE

         IF( VTE .LE. 0D0 ) THEN

C           LE TETRAEDRE NF1-POINT NPt EST DE VOLUME VTE<=0
C           CAS D'UN POINT SUR UNE FACE ET TETRAEDRE OPPOSE PLAT
C           DEBORDANT D'UN COTE => CREATION DE TETRAEDRES DE VOLUMES<0
            WRITE(IMPRIM,*) 'etoipt19 : PB TETRAEDRE Face',NF1,
     %                      ' Point',NPt,' de VOLUME',VTE,' <=0'
            WRITE(IMPRIM,*)('PTXYZD(',kkk,',',NS1,')=',
     %                       PTXYZD(kkk,NS1), kkk=1,4)
            WRITE(IMPRIM,*)('PTXYZD(',kkk,',',NS2,')=',
     %                       PTXYZD(kkk,NS2), kkk=1,4)
            WRITE(IMPRIM,*)('PTXYZD(',kkk,',',NS3,')=',
     %                       PTXYZD(kkk,NS3), kkk=1,4)
            WRITE(IMPRIM,*) 'COBARY=',COBARY

C           TRACE DES FACES SIMPLES DE NFETOI VERSION 1 ET DU POINT NPt
         KTITRE='etoipt19: NPt=          LES FACES SIMPLES DE L''ETOILE'
            WRITE(KTITRE(15:22),'(I8)') NPt
            CALL SANSDBL( KTITRE, NBC )
            CALL TRFETOV1(KTITRE(1:NBC),NPt,N1FEOC,NFETOI,NOTETR,PTXYZD)

C           SUPPRESSION DU TETRAEDRE NTE DE L'ETOILE INITIALE SAUF SI UNIQUE
            IF( NBTERE .EQ. 1 ) THEN

C              TETRAEDRISATION FORCEE DE L'UNIQUE TETRAEDRE
               NOTET1 = NOTERE(1)
               CALL COBATE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NOTET1),
     %                      VTE, COBARY, IERR )
C               COBARY(I) EST ICI LA COORDONNEE BARYCENTRIQUE OPPOSEE A LA FACE I
               DO NBCBFA=1,4
                  IF( COBARY(NBCBFA) .LT. 0D0 ) GOTO 560
               ENDDO
               GOTO 570
          
C              PROJECTION DU POINT NPt SUR LA FACE PROCHE DE NOTET1
 560           CALL PRPTPLD(PTXYZD(1,NPt),
     %                      PTXYZD(1,NOTETR(NOSOFATE(1,NBCBFA),NOTET1)),
     %                      PTXYZD(1,NOTETR(NOSOFATE(2,NBCBFA),NOTET1)),
     %                      PTXYZD(1,NOTETR(NOSOFATE(3,NBCBFA),NOTET1)),
     %                      PTPROJ, IERR )
               IF( IERR .NE. 0 ) GOTO 570

C              DEPLACEMENT DU POINT NPt AU POINT PTPROJ
               NPSOFR( NPt ) = - NPSOFR( NPt )
               PTXYZD( 1, NPt ) = PTPROJ( 1 )
               PTXYZD( 2, NPt ) = PTPROJ( 2 )
               PTXYZD( 3, NPt ) = PTPROJ( 3 )

               GOTO 580

            ENDIF

C           SUPPRESSION DU TETRAEDRE NTE DE L'ETOILE INITIALE
 570        DO K=NBTERE,1,-1
               IF( NOTERE(K) .EQ. NTE ) THEN
                  DO L=K+1, NBTERE
                     NOTERE(L-1) = NOTERE(L)
                  ENDDO
                  NBTERE = NBTERE - 1
                  IF( NBTERE .LE. 0 ) THEN
                     PRINT*,'etoipt19 : ABANDON apres REPRISE MULTIPLES 
     %de NPt=',NPt,' et',NBTERE,' TETRAEDRES INITIAUX CONTENANT NPt'
C                    ABANDON DE LA TETRAEDRISATION DU POINT NPt
                     IF( NPSOFR( NPt ) .GT. 0 ) THEN
                        NBSTPE = NBSTPE + 1
                        NOSTPE( NBSTPE ) = NPt
                     ENDIF
                     NBTEET = -1
                     IERR   = 11
                     GOTO 9999
                  ENDIF
                  GOTO 510
               ENDIF
            ENDDO

         ENDIF

C        TETRAEDRE CORRECT. PASSAGE A LA FACE SIMPLE SUIVANTE
 580     NF1 = NFETOI( 5, NF1 )
         GOTO 550

      ENDIF

C     CONSTRUCTION DES TETRAEDRES DE L'ETOILE de la REPRISE du POINT NPt
C     ------------------------------------------------------------------
C     REPRISE avec TETRAEDRISATION FORCEE DANS teteto MEME SI CELA GENERE
C     DES TETRAEDRES DE VOLUME<=0 ou de QUALITE FAIBLE
      NBTEETIN = NBTEET
      CALL TETETO( NPt,   REPRISE, QTEAME, MXSOMM, PTXYZD, NPSOFR,
ccc     %             INFACO, MXFACO, LEFACO,
     %             NOVOLU, IVOLTE, NVOLTE,
     %             MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %             MXETOI, N1FEOC, NFETOI, VOLET1,
     %             MXETOI, NBTEET, NOTEET,
     %             MXSTPE, NBSTPE, NOSTPE, IERR )

ccc      PRINT*,'etoipt19 : NPt=',NPt,' VOLET0=',VOLET0,' VOLET1=',VOLET1


C     =========================================================================
C     NBTEET TETRAEDRES CONSTITUENT L'ETOILE FINALE de NPt avec ou sans REPRISE
C     =========================================================================
C     AJOUT DANS LE PAVAGE DU POINT NPt
 1100 CALL NUPAVEST( NPt,      PTXYZD,
     %               HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV )

C     ESSAI D'AMELIORER LES TETRAEDRES DE L'ETOILE PAR 2T -> 3T ou 3T -> 2T
C     ---------------------------------------------------------------------
      CALL EC2TMTNT( NPt,    PTXYZD, INFACO,  MXFACO, LEFACO,
     %               IVOLTE, NVOLTE, MXETOI,  NBTEET, NOTEET,
     %               NOTETR, N1TEVI, NUDTETR, N1TETS,
     %               NBEC2T3T, NBECMTNT )

C     BILAN FINAL: LA QUALITE DES TETRAEDRES DE L'ETOILE FINALE DU POINT NPt
C     ============
      DMAX   = 0D0
      QUAMOY = 0
      QUAMIN = 2
      NBTEVV = 0
      DO I = 1, NBTEET

         NTE = NOTEET(I)
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

            NBTEVV = NBTEVV + 1
            NOTEET( NBTEVV ) = NTE

C           POUR VOIR SI NPt INTERNE A LA BOULE CIRCONSCRITE
            CALL CEBOQU( 0., NOTETR(1,NTE), PTXYZD,  CENTRE, VTE, QTE )

C           DISTANCE(NPt-CENTRE)**2 / RAYON**2
            D=( ( CENTRE(1) - PTXYZD(1,NPt) ) ** 2 +
     %          ( CENTRE(2) - PTXYZD(2,NPt) ) ** 2 +
     %          ( CENTRE(3) - PTXYZD(3,NPt) ) ** 2 ) / CENTRE(4)

ccc            PRINT *,'etoipt19 : NPt',NPt,I,' new tetra',NTE,' St:',
ccc     %          (NOTETR(K,NTE),K=1,8),' Q=',QTE,' V=',VTE,' D=',D

            IF( D .GT. DMAX ) DMAX = D
            IF( QTE .LT. QUAMIN ) QUAMIN = QTE
            QUAMOY = QUAMOY + QTE

         ENDIF

      ENDDO
      IF( NBTEVV .GT. 0 ) QUAMOY = QUAMOY / NBTEVV

      IF( REPRISE ) THEN
       PRINT *,'etoipt19 : NPt=',NPt, NBTEETIN,' TETRAEDRES INITIAUX ->'
     %         ,NBTEVV,' FINAUX avec DMax=',DMAX,' QMIN=',QUAMIN,
     %         ' QMOY=',QUAMOY,NBEC2T3T,' Echanges 2T->3T et',
     %          NBECMTNT,' mT->2m-4T'
      ENDIF

C     TRACE FILAIRE DES ARETES DES NBTEET TETRAEDRES DE L'ETOILE FINALE
      CALL TETR1S( NPt,    N1TETS, NOTETR,
     %             NBTEET, MXETOI, NOTEET, IERR )

      IF( NBTEET .GT. 0 ) THEN
      KTITRE='etoipt19 : POINT                    TETRAEDRES pour ETOILE
     % FINALE'
         WRITE(KTITRE(18:25),'(I8)') NPt
         WRITE(KTITRE(28:35),'(I8)') NBTEET
         CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NPt, NBTEET, NOTEET )
      ENDIF

 9999 TRACTE = TRACTE0
      RETURN
      END
