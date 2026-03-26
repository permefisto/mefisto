      SUBROUTINE TETRA1ST( KNMVOLU, NPt,   NOTEIN,
     %                     MXSOMM, PTXYZD, NPSOFR,
     %                     HEXAPAVE,NBIPAV,ECHPAV, N1SPAVE,NOPTSUIV,
     %                     MXTETR, N1TEVI, NUDTETR,NOTETR, N1TETS,
     %                     INFACO, MXFACO, LEFACO, N1FASC,
     %                     IVOLTE, NVOLTE,
     %                     MXETOI, N1FEOC, N1FEVI, NFETOI,
     %                     NBTEET, NTETOI, VOLET0, VOLET1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DES FACES SIMPLES DE L'ETOILE DES NBTEET
C -----    TETRAEDRES DE BOULE CIRCONSCRITE CONTENANT LE POINT NPt
C          ET CONSTRUCTION DES NBTEET TETRAEDRES DE SOMMET NPt

C ENTREES:
C --------
C KNMVOLU: NOM DU VOLUME A TETRAEDRISER
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C NOTEIN : NUMERO DU TETRAEDRE DE DEBUT DE RECHERCHE
C MXSOMM : NOMBRE MAXIMAL DE POINTS DECLARABLES DANS PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : NUMERO DES POINTS INITIAUX
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  3 SI LE POINT EST IMPOSE PAR L'UTILISATEUR SANS
C                     OBLIGATION DE TETRAEDRISATION
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET D'OT ET RECONNU TROP PROCHE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C HEXAPAVE: MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV  : NOMBRE D'ARETES DANS LA DIRECTION I
C ECHPAV  : ECHELLE DANS LA DIRECTION I
C N1SPAVE : NO DU 1-ER SOMMET DANS PTXYZD DU PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES

C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES STOCKABLES DANS NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE ACTIF OCCUPE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET

C INFACO : = 0 PAS DE TABLEAU LEFACO NI DE CONSERVATION DES
C              FACES FRONTIERE ( NO DE VOLUME CONNU PAR NVOLTE )
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONSERVATION DES
C              FACES DE LA FRONTIERE
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3

C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3), MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS

C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : >0 NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE
C          -1 SI TETRAEDRE INACTIF ou SUPPRIME


C MXETOI : NOMBRE MAXIMAL DE FACES DU TABLEAU NFETOI et
C          NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NTETOI
C          DANS UNE ETOILE

C ENTREES ET SORTIES :
C --------------------
C N1TEVI : NUMERO DU 1 PREMIER TETRAEDRE VIDE DANS LE TABLEAU NOTETR
C          LE CHAINAGE DES TETRAEDRES VIDES SE FAIT SUR NOTETR(5,.)
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

C SORTIES :
C ---------
C NBTEET : NOMBRE DE TETRAEDRES DE L'ETOILE DU POINT NPt
C          =<0 SI PROBLEME RENCONTRE ET ETOILE NON FORMEE
C NTETOI : NUMERO DANS NOTETR DES NBTEET TETRAEDRES DE L'ETOILE
C VOLET0 : LE VOLUME DES TETRAEDRES INITIAUX DE L'ETOILE DU POINT NPt 
C VOLET1 : LE VOLUME DES TETRAEDRES FINAUX   DE L'ETOILE DU POINT NPt 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Fevrier 2019
C23456...............................................................012
ccc      PARAMETER        (QTEMED=0.001, MXTECH=1024 ) 2019/03/15
ccc      PARAMETER        (QTEMED=0.005, MXTECH=1024 ) 2020/05/22
      PARAMETER        (QTEMED=0.01, MXTECH=1024 )
      include"./incl/langue.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*160     KTITRE
      CHARACTER*24      KNMVOLU
      INTEGER           NOTETR(8,MXTETR), N1TETS(MXSOMM),NPSOFR(MXSOMM),
     %                  NTETOI(MXETOI), NFETOI(5,MXETOI),
     %                  LEFACO(1:11,0:MXFACO), N1FASC(1:MXSOMM),
     %                  NVOLTE(MXTETR)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
      INTEGER           NBIPAV(3), N1SPAVE(0:*), NOPTSUIV(1:*)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM),
     %                  VNT, VTOP, VOLET0, VOLET1, VOLTID,
     %                  D, DSR, COBARY(4)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NSFC(3), NOTECH(MXTECH)

      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      TRACTE0 = TRACTE

C     ACTUELLEMENT NPt EST UN SOMMET D'AUCUN TETRAEDRE
      N1TETS( NPt ) = 0

C     VOLUME IDEAL DU TETRAEDRE DE SOMMET NPt
C     LA TAILLE SOUHAITEE DE L'ARETE PTXYZD(4,NPt) EST SUPPOSEE INITIALISEE
      VOLTID = ( PTXYZD(4,NPt) ** 3 ) / 6
      IF( VOLTID .LE. 0D0 ) GOTO 9900

C     RECHERCHE A PARTIR DU TETRAEDRE NOTEIN D'UN TETRAEDRE NOTET1
C     CONTENANT LE POINT NPt
C     ------------------------------------------------------------
      IF( IVOLTE .EQ. 0 ) THEN
         CALL R1TCO1P0( NPt,    PTXYZD, NOTEIN, NOTETR, NUDTETR,
     %                  NOTET1, COBARY )
      ELSE
         CALL R1TCO1P1( NPt,     MXSOMM, PTXYZD, 1,
     %                  HEXAPAVE,NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                  NOTEIN,  MXTETR, NOTETR, NUDTETR, N1TETS,
     %                  MXETOI,  NTETOI, MXTECH, NOTECH,
     %                  NOTET1,  COBARY )
      ENDIF

 5    IF( NOTET1 .LE. 0 ) THEN

C        PAS DE TETRAEDRE NOTET1 RETROUVE CONTENANT NPt
C        LE POINT NPt PEUT IL ETRE SUPPRIME?
         IF( NPSOFR(NPt) .NE. -4 .AND. NPSOFR(NPt) .NE. 0 .AND.
     %       NPSOFR(NPt) .NE. 3 ) THEN

C           NON:  LE POINT NPt  NE PEUT PAS ETRE SUPPRIME
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'tetra1st: POINT',NPt,NPt,' NPSOFR=',NPSOFR(NPt),
     %                ' NE PEUT ETRE SUPPRIME'
            ELSE
               PRINT*,'tetra1st: POINT',NPt,NPt,' NPSOFR=',NPSOFR(NPt),
     %                ' CANNOT BE DELETED'
            ENDIF
            PRINT*,'tetra1st: XYZ',NPt,' =',(PTXYZD(I,NPt),I=1,4)

ccc         ELSE
cccC           OUI: LE POINT NPt EST SUPPRIME
ccc            IF( LANGAG .EQ. 0 ) THEN
ccc               PRINT*,'tetra1st: le POINT',NPt,' NPSOFR=',NPSOFR(NPt),
ccc     %                ' N''EST PAS TETRAEDRISE'
ccc            ELSE
ccc               PRINT*,'tetra1st: the POINT',NPt,' NPSOFR=',NPSOFR(NPt),
ccc     %                ' IS NOT TETRAHEDRIZED'
ccc            ENDIF

         ENDIF

         GOTO 9900

      ENDIF

C     NUMERO DU VOLUME DU TETRAEDRE NOTET1
      IF( IVOLTE .NE. 0 ) THEN
         NOVOLU = NVOLTE( NOTET1 )
      ELSE
         NOVOLU = 0
      ENDIF


C     FORMATION DE L'ETOILE DES NBTEET00 TETRAEDRES CONTENANT LE POINT NPt
C     A PARTIR DU TETRAEDRE NOTET1 DE NOTETR ISSU DE R1TCO1P
C     --------------------------------------------------------------------
      CALL ETOIL1PT( NPt,    PTXYZD, NPSOFR,
     %               NOTET1, NOTETR, N1TETS,   COBARY,
     %               NBCBA0, MXETOI, NBTEET00, NTETOI )
      IF( NBTEET00 .LE. 0 ) THEN
         NOTET1 = 0
         GOTO 5
      ENDIF

      IF( TRACTE ) THEN
      KTITRE='tetra1st: NPt=          ETOILE INITIALE de        TETRAEDR
     %ES. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(44:48),'(I5)') NBTEET00
      CALL SANSDBL( KTITRE, NBC )
ccc      PRINT*,KTITRE(1:NBC)
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET00, NTETOI )
      ENDIF


C     NPt SUPPRIMABLE EST IL TROP PROCHE DE L'UN DES SOMMETS DES NTETOI?
C     ------------------------------------------------------------------
      IF( NPSOFR(NPt) .EQ. 0 ) THEN
         DO N=1,NBTEET00
            NTE = NTETOI(N)
            DO K=1,4
               NS = NOTETR(K,NTE)
               D = SQRT( (PTXYZD(1,NS) - PTXYZD(1,NPt) ) ** 2
     %                 + (PTXYZD(2,NS) - PTXYZD(2,NPt) ) ** 2
     %                 + (PTXYZD(3,NS) - PTXYZD(3,NPt) ) ** 2 )
               IF( D .LT. PTXYZD(4,NPt)*0.333D0 ) THEN
C                 NPt TROP PROCHE DU SOMMET NS N'EST PAS TETRAEDRISE
                  GOTO 9900
               ENDIF
            ENDDO
         ENDDO
      ENDIF

C     FORMATION DES FACES SIMPLES DE L'ETOILE A PARTIR DES NBTEET00
C     TETRAEDRES CONTENANT LE POINT NPt
C     -------------------------------------------------------------
      NBFOET = 0
  8   NBFOET = NBFOET + 1

C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
      CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

      DO 10 I=1,NBTEET00

C        LE TETRAEDRE I DE L'ETOILE
         NT = NTETOI( I )
         IF( NT .LE. 0 ) GOTO 10

         DO NFNT=1,4
C           AJOUT ou SUPPRESSION de la FACE NFNT de NT de l'ETOILE NFETOI
            CALL AJFAET1( NT,     NFNT,   NOTETR,
     %                    N1FEOC, N1FEVI, NFETOI, NFS )
C           NFS>0 NUMERO DANS NFETOI DE LA FACE AJOUTEE
C              =0 SUPPRESSION DE LA FACE DE NFETOI
         ENDDO

 10   ENDDO

      KTITRE='tetra1st: NPt=          LES FACES SIMPLES DE SON ETOILE IN
     %ITIALE. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      CALL SANSDBL(  KTITRE, NBC )
      CALL TRFETOV1( KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR, PTXYZD)


C     VOLUME ET QUALITE DES TETRAEDRES INITIAUX FACE DE L'ETOILE-POINT NPt
C     --------------------------------------------------------------------
      NBTEET01 = NBTEET00
      NBFETO = 0
      QMIN01 = 2.
      NF2    = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE INITIALE
 20   IF( NF2 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NBFETO = NBFETO + 1
         NT   = NFETOI(1,NF2)
         NFNT = NFETOI(2,NF2)

         IF( NT .GT. 0 ) THEN

C           VOLUME et QUALITE DU TETRAEDRE FACE NFNT DE NT-POINT NPt
            CALL QUATETD( PTXYZD( 1, NOTETR( NOSOFATE(1,NFNT), NT ) ),
     %                    PTXYZD( 1, NOTETR( NOSOFATE(2,NFNT), NT ) ),
     %                    PTXYZD( 1, NOTETR( NOSOFATE(3,NFNT), NT ) ),
     %                    PTXYZD( 1, NPt ),
     %                    ARMIN, ARMAX, SURFTR, VNT, QNT )
            IF( QNT .LT. QMIN01 ) QMIN01 = QNT

            IF( VNT .LE. VOLTID*1D-4 .OR. QNT .LT. QTEMED ) THEN

C              FACE NFNT de NT-Point NPt est UN TETRAEDRE DEGENERE:
C              NPt A SUPPRIMER
C              ---------------------------------------------------

C              LE POINT NPt PEUT IL NE PAS ETRE TETRAEDRISE?
               IF( NPSOFR(NPt) .EQ. -4 .OR. NPSOFR(NPt) .EQ. 0 .OR.
     %             NPSOFR(NPt) .EQ. 3 ) THEN

C                 OUI: LE POINT NPt N'EST PAS TETRAEDRISE
C                 ---------------------------------------
ccc                  IF( LANGAG .EQ. 0 ) THEN
ccc                     PRINT*,'tetra1st: NPt=',NPt,
ccc     %                      ' POINT NON TETRAEDRISE'
ccc                  ELSE
ccc                     PRINT*,'tetra1st: NPt=',NPt,
ccc     %                      ' NOT TETRAHEDRIZED POINT'
ccc                  ENDIF

                  GOTO 9900

               ENDIF

C              NON: LE POINT NPt RESTE A ETRE TETRAEDRISE
C              ---- -------------------------------------
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'tetra1st: NPt=',NPt,' le TETRAEDRE SIMULE',
     %                    NOTETR( NOSOFATE(1,NFNT), NT ),
     %                    NOTETR( NOSOFATE(2,NFNT), NT ),
     %                    NOTETR( NOSOFATE(3,NFNT), NT ), NPt,
     %                  ' est DEGENERE Q=',QNT,' V=',VNT,
     %                  ' NBTEET00=',NBTEET00,' NPSOFR=',NPSOFR(NPt)
               ELSE
                  PRINT*,'tetra1st: NPt=',NPt,' the SIMULATED TETRA',
     %                    NOTETR( NOSOFATE(1,NFNT), NT ),
     %                    NOTETR( NOSOFATE(2,NFNT), NT ),
     %                    NOTETR( NOSOFATE(3,NFNT), NT ), NPt,
     %                  ' is DEGENERATED Q=',QNT,' V=',VNT,
     %                  ' NBTEET00=',NBTEET00,' NPSOFR=',NPSOFR(NPt)
               ENDIF

C              SUPPRESSION DU TETRAEDRE NT DE L'ETOILE NOTEET
C              ----------------------------------------------
               DO I=1,NBTEET01

                  IF( NTETOI( I ) .EQ. NT ) THEN
C                    SUPPRESSION DU TETRAEDRE NT DE NTETOI
                     DO II=I+1,NBTEET01
                        NTETOI( II-1 ) = NTETOI( II )
                     ENDDO
                     NBTEET01 = NBTEET01-1

                     IF( NBTEET01 .LE. 0 ) THEN
C                       NTETOI EST VIDE
C                       CREATION DES TETRAEDRES DE SOMMET NPt
C                       A PARTIR DE NOTET1
                        NBTEET00 = 1
                        NBTEET0  = 1
                        NTETOI( 1 ) = NOTET1
                        GOTO 40
                     ENDIF

C                    REPRISE DE L'ETOILE
                     NBTEET00 = NBTEET01
                     GOTO 8

                  ENDIF
               ENDDO

            ENDIF

         ENDIF

C        PASSAGE A LA FACE SIMPLE SUIVANTE
         NF2 = NFETOI( 5, NF2 )
         GOTO 20

      ENDIF

ccc      PRINT*,'tetra1st: NPt=',NPt,' QUALITE INITIALE MIN DES ',NBTEET00,
ccc     %       ' TETRAEDRES FACES-NPt QMIN01=',QMIN01

      IF( NBTEET00 .LE. 0 ) THEN
C        NTETOI EST VIDE: CREATION DES TETRAEDRES DE SOMMET NPt
C        A PARTIR DE NOTET1 POUR EVITER D'AJOUTER DES TETRAEDRES
C        DE QUALITE MEDIOCRE
C        -------------------------------------------------------
         NBTEET00 = 1
         NBTEET0  = 1
         NTETOI( 1 ) = NOTET1
         GOTO 40
      ENDIF

      NBTEET0 = NBTEET00

ccc      IF( NBCBA0 .EQ. 2 ) THEN
cccC        NPt EST SUR UNE ARETE: PAS D'EXTENSION DE L'ETOILE POUR
cccC        EVITER D'AJOUTER DES TETRAEDRES DE QUALITE MEDIOCRE
cccC        -------------------------------------------------------
ccc         GOTO 40
ccc      ENDIF


C     EXTENSION DE L'ETOILE AUX TETRAEDRES OPPOSES PAR LES FACES DES NBTEET0
C     TETRAEDRES DE L'ETOILE ET DONT LA BOULE CIRCONSCRITE CONTIENT LE POINT
C     CONSTRUCTION DU TABLEAU NTETOI DES NBTEET0 TETRAEDRES DE L'ETOILE
C     ======================================================================
 25   NF2 = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE NFETOI DE NPt
 30   IF( NF2 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT   = NFETOI(1,NF2)
         NFNT = NFETOI(2,NF2)

C        LE TETRAEDRE OPPOSE PAR LA FACE NFNT DE NT
         NTOP = NOTETR( 4+NFNT, NT )
         IF( NTOP .LE. 0 ) THEN
C           LA FACE NFNT DE NT EST SUR LA FRONTIERE. NT N'EST PAS AJOUTE
            GOTO 39
         ENDIF

         IF( INFACO .NE. 0 ) THEN
C           LA FACE EST ELLE DANS LEFACO?
            CALL NULEFT( NFNT, NT, NOTETR, MXFACO, LEFACO, NF )
C           NF>0 NUMERO LEFACO DE LA FACE NFNT DE NT RETROUVEE
C             =0 SI LA FACE N'EST PAS RETROUVEE
            IF( NF .GT. 0 ) GOTO 39
         ENDIF

         IF( IVOLTE .NE. 0 ) THEN
            IF( NVOLTE(NTOP) .NE. NOVOLU ) THEN
C              NT ET NTOP SONT DANS DES VOLUMES DIFFERENTS
C              TETRAEDRE NTOP NON AJOUTE A L'ETOILE
               GOTO 39
            ENDIF
         ENDIF

C        LE POINT NPt EST IL DANS LA BOULE DU TETRAEDRE NTOP
         CALL PTDSBO( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NTOP),
     %                NONOUI, VTOP,  DSR )
C        NONOUI=1  Pt   EST     DANS LA BOULE CIRCONSCRITE de NTOP
C               0  Pt N'EST PAS DANS LA BOULE CIRCONSCRITE
C              -1  LE TETRAEDRE EST DE VOLUME TRES FAIBLE ou NUL ou NEGATIF

         IF( NONOUI .LE. 0 ) THEN
C           TETRAEDRE NTOP NON AJOUTE A L'ETOILE de NPt
            GOTO 39
         ENDIF

C        NPt EST DANS LA BOULE CIRCONSCRITE DE NTOP
C            ET LA FACE N'EST PAS DANS LEFACO ou
C            NT ET NTOP SONT DANS LE MEME VOLUME

C        L'AJOUT DES FACES DE NTOP CREERAIT IL DES TETRAEDRES MEDIOCRES?
C        NFOP NO LOCAL A NTOP DE LA FACE NFNT DE NT
         CALL NOFAOP( NFNT, NT, NOTETR, NFOP, NTOP )

         DO 36 NF = 1, 4
            IF( NF .NE. NFOP ) THEN

C              LA FACE NF DE NTOP EST ELLE UNE FACE SIMPLE DE L'ETOILE?
               NFF2 = N1FEOC
 35            IF( NFF2 .GT. 0 ) THEN
                  IF( NFETOI(1,NFF2) .EQ. NTOP .AND.
     %                NFETOI(2,NFF2) .EQ. NF ) THEN
C                     LA FACE NF DE NTOP EST IDENTIQUE A NFF2
C                     ELLE SERA SUPPRIMEE SI NTOP EST AJOUTE
                      GOTO 36
                  ENDIF

C                 LA FACE EST DIFFERENTE. PASSAGE A LA SUIVANTE
                  NFF2 = NFETOI(5,NFF2)
                  GOTO 35
               ENDIF

C              LA FACE NF DU TETRAEDRE NTOP N'EST PAS UNE FACE SIMPLE
C              ELLE SERA AJOUTEE SI LE TETRAEDRE NTOP EST AJOUTE
C              VOLUME et QUALITE DU TETRAEDRE FACE NF-POINT NPt
               CALL QUATETD( PTXYZD( 1, NOTETR(NOSOFATE(1,NF),NTOP) ),
     %                       PTXYZD( 1, NOTETR(NOSOFATE(2,NF),NTOP) ),
     %                       PTXYZD( 1, NOTETR(NOSOFATE(3,NF),NTOP) ),
     %                       PTXYZD( 1, NPt ),
     %                       ARMIN, ARMAX, SURFTR, VNT, QNT )

               IF( VNT .LE. VOLTID*1D-4 .OR. QNT .LT. QTEMED ) THEN
C                 NTOP N'EST PAS AJOUTE A L'ETOILE
                  GOTO 39
               ENDIF

            ENDIF
 36      ENDDO

C        NPt SUPPRIMABLE TROP PROCHE DE L'UN DES SOMMETS DE NTOP?
C        --------------------------------------------------------
         IF( NPSOFR(NPt) .EQ. 0 ) THEN
            DO K=1,4
               NS = NOTETR(K,NTOP)
               D = SQRT( (PTXYZD(1,NS) - PTXYZD(1,NPt) ) ** 2
     %                 + (PTXYZD(2,NS) - PTXYZD(2,NPt) ) ** 2
     %                 + (PTXYZD(3,NS) - PTXYZD(3,NPt) ) ** 2 )
               IF( D .LT. PTXYZD(4,NPt)*0.333D0 ) THEN
C                 NPt TROP PROCHE DU SOMMET NS N'EST PAS TETRAEDRISE
                  GOTO 9900
               ENDIF
            ENDDO
         ENDIF

C        AJOUT DE NTOP AUX TETRAEDRES DE L'ETOILE DE NPt
C        -----------------------------------------------
C        NTOP EST IL DEJA UN TETRAEDRE DE L'ETOILE DE NPt?
         DO K = 1, NBTEET0
            IF( NTOP .EQ. NTETOI(K) ) GOTO 39
         ENDDO
         IF( NBTEET0 .GE. MXETOI ) THEN
            GOTO 9800
         ENDIF
         NBTEET0 = NBTEET0 + 1
         NTETOI( NBTEET0 ) = NTOP

C        LES 4 FACES DU TETRAEDRE OPPOSE NTOP PEUVENT ETRE AJOUTEES
C        OU SUPPRIMEES DE L'ETOILE DE NPt
         DO NF = 1, 4
            CALL AJFAET1( NTOP,   NF,     NOTETR,
     %                    N1FEOC, N1FEVI, NFETOI, NFS )
C           NFS>0 NUMERO DANS NFETOI DE LA FACE AJOUTEE
C              =0 SUPPRESSION DE LA FACE DE NFETOI
         ENDDO
         GOTO 25

C        PASSAGE A LA FACE SIMPLE SUIVANTE DE L'ETOILE de NPt
 39      NF2 = NFETOI( 5, NF2 )
         GOTO 30

      ENDIF

C     VERIFICATION DE LA NON EXISTENCE D'ARETES NON DOUBLE DES
C     FACES SIMPLES DE L'ETOILE NFETOI
C     --------------------------------------------------------
      CALL VEARFSET( NOTETR, N1FEOC, NFETOI, NBARNODO )

      IF( NBARNODO .GT. 0 ) THEN
         KTITRE='tetra1st: NPt=                ARETES NON DOUBLES DES FA
     %CES SIMPLES de NFETOI. Volume ' // KNMVOLU
         WRITE(KTITRE(15:22),'(I8)') NPt
         WRITE(KTITRE(24:26),'(I3)') NBARNODO
         CALL SANSDBL( KTITRE, NBC )
         CALL TRFETOV1(KTITRE(1:NBC), NPt, N1FEOC,NFETOI,NOTETR, PTXYZD)

C        ABANDON du POINT NPt
         GOTO 9900
      ENDIF

C     BILAN: QUALITE et VOLUME DES NBTEET0 TETRAEDRES DE
C     L'ETOILE ETENDUE DU POINT NPt
C     --------------------------------------------------
      QMIN0  = 2
      VOLET0 = 0D0
      DO I = 1, NBTEET0
         NT = NTETOI( I )
         IF( NT .GT. 0 ) THEN

C           VOLUME et QUALITE DU TETRAEDRE NT
            CALL QUATETD( PTXYZD( 1, NOTETR( 1, NT ) ),
     %                    PTXYZD( 1, NOTETR( 2, NT ) ),
     %                    PTXYZD( 1, NOTETR( 3, NT ) ),
     %                    PTXYZD( 1, NOTETR( 4, NT ) ),
     %                    ARMIN, ARMAX, SURFTR, VNT, QNT )

            VOLET0 = VOLET0 + VNT
            IF( QNT .LT. QMIN0 ) QMIN0 = QNT

         ENDIF
      ENDDO


C     =====================================================================
C     CONSTRUCTION DES TETRAEDRES FACES SIMPLES - POINT NPt
C     =====================================================================
C     TRACE DES FACES SIMPLES DE L'ETOILE
 40   IF( TRACTE ) THEN
      KTITRE='tetra1st: NPt=          LES FACES SIMPLES FINALES DE L''ET
     %OILE avant TETRAEDRISATION du POINT NPt. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      CALL SANSDBL(  KTITRE, NBC )
      CALL TRFETOV1( KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR, PTXYZD)
      ENDIF

C     PASSAGE DE NFETOI VERSION 1 => VERSION 2
C     LE NO DE TETRAEDRE INTERNE ETOILE  => NO TETRAEDRE OPPOSE A LA FACE
C     LE NO FACE LOCAL DANS LE TETRAEDRE => NO 1-ER SOMMET DE LA FACE
C     LE NO INUTILISE 3                  => NO 2-ME SOMMET DE LA FACE
C     LE NO INUTILISE 4                  => NO 3-ME SOMMET DE LA FACE
C     CONSTRUCTION DE LA LISTE DES TETRAEDRES OPPOSES AUX TETRAEDRES
C     ------------------------------------------------------------------
C     CALCUL DU NOMBRE DE FACES SIMPLES DE L'ETOILE DE NPt
      NBFETO  = 0
      NF      = N1FEOC
C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 45   IF( NF .GT. 0 ) THEN
C        NFETOI VERSION 1.  UNE FACE SIMPLE DE L'ETOILE DE PLUS
         NBFETO = NBFETO + 1
C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI( 5, NF )
         GOTO 45
      ENDIF

      NBTEET2 = NBTEET0 + NBFETO
      NF      = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE DE NPt
 50   IF( NF .GT. 0 ) THEN

C        NFETOI VERSION 1
C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT   = NFETOI( 1, NF )
         NFNT = ABS( NFETOI( 2, NF ) )

C        LE NUMERO DU TETRAEDRE NTOP AU DELA DE LA FACE NFNT DE NT
C        QUI PEUT ETRE NUL SI LA FACE EST FRONTIERE
C        NFOP NO LOCAL A NTOP DE LA FACE NFNT DE NT
         CALL NOFAOP( NFNT, NT, NOTETR,  NFOP, NTOP )
         IF( NFOP .GT. 0 ) THEN

C           STOCKAGE DU TETRAEDRE OPPOSE POUR RESOUDRE
C           LES TETRAEDRES OPPOSES PAR CALL  MJOPTE
C           NTOP NE PEUT ETRE UN TETRAEDRE A SUPPRIMER
            DO K = 1, NBTEET0
               IF( NTETOI(K) .EQ. NTOP ) GOTO 52
            ENDDO
C           NTOP EST PEUT ETRE DEJA STOCKE
            DO K = NBTEET0+NBFETO+1, NBTEET2
               IF( NTETOI(K) .EQ. NTOP ) GOTO 52
            ENDDO
            IF( NBTEET2 .GE. MXETOI ) THEN
C              SATURATION DES TETRAEDRES NTETOI DE L'ETOILE
               GOTO 9800
            ENDIF
            NBTEET2 = NBTEET2 + 1
            NTETOI( NBTEET2 ) = NTOP

cccC           NT VA ETRE SUPPRIME -> TETRAEDRE OPPOSE INCONNU
ccc 52         NOTETR( 4+NFOP, NTOP ) = -1
ccc            NOTETR( 4+NFNT, NT   ) = -1

         ENDIF

C        PASSAGE NFETOI VERSION 1 -> VERSION 2
C        NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NFT DU TETRAEDRE NT
 52      NFETOI( 1, NF ) = NTOP

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
C        STOCKAGE DANS NFETOI(2:4,NF)
         NFETOI( 2, NF ) = NOTETR( NOSOFATE(1,NFNT), NT )
         NFETOI( 3, NF ) = NOTETR( NOSOFATE(2,NFNT), NT )
         NFETOI( 4, NF ) = NOTETR( NOSOFATE(3,NFNT), NT )

C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI( 5, NF )
         GOTO 50

      ENDIF


C     VERIFICATION DE LA QUALITE DE LA FUTURE TETRAEDRISATION DU POINT NPt
C     --------------------------------------------------------------------
      VOLET1 = 0D0
      QMIN1  = 2.
      NF2    = N1FEOC

C     BOUCLE SUR LES FACES DE L'ETOILE
 55   IF( NF2 .GT. 0 ) THEN

C        REMPLISSAGE DU TETRAEDRE NT DE SOMMETS CEUX DE LA FACE NF2 ET NPt
         NS1 = NFETOI(2,NF2)
         NS2 = NFETOI(3,NF2)
         NS3 = NFETOI(4,NF2)

C        QUALITE ET VOLUME DU TETRAEDRE A CREER
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VNT, QNT )

         IF( VNT.LE.VOLTID*1D-4 .OR. QNT.LT.QTEMED ) THEN
C           Le TETRAEDRE NF2+NPt est TROP MEDIOCRE
C           CE CAS NE DEVRAIT PAS ARRIVER car DETECTE AUPARAVANT...
            NTOP = NFETOI(1,NF2)
            IF( NTOP .GT. 0 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'tetra1st: NPt=',NPt,' Face:',NS1,NS2,NS3,
     %                   ' NE DEVRAIT PAS EXITER car V=',VNT,' Q=',QNT
               ELSE
                  PRINT*,'tetra1st: NPt=',NPt,' Face:',NS1,NS2,NS3,
     %                   ' SHOULD NOT EXIST to V=',VNT,' Q=',QNT
               ENDIF
            ENDIF
         ENDIF

         VOLET1 = VOLET1 + ABS( VNT )
         QMIN1  = MIN( QMIN1, QNT )

C        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
         NF2 = NFETOI(5,NF2)
         GOTO 55

      ENDIF

      IF( ABS( VOLET1-VOLET0 ) .GT. 1D-5*(VOLET1+VOLET0) ) THEN
C        VOLUMES NON COMPATIBLES
         print*
         print*,'tetra1st: NPt=',NPt,' ETOILE INITIALE VOLET0=',VOLET0,
     %          ' VOLET1=',VOLET1,' QMIN0=',QMIN0,' QMIN1=',QMIN1
         print*,'tetra1st: NPt=',NPt,' POINT a SUPPRIMER NPSOFR=',
     %           NPSOFR(NPt)
         GOTO 9900
      ENDIF


C     DECLARATION DES TETRAEDRES ETOILANT LE POINT NPt
C     ------------------------------------------------
      NBTEET1= NBTEET0
      VOLET1 = 0D0
      QMIN1  = 2.
      QMOY   = 0.
      NF2    = N1FEOC

C     BOUCLE SUR LES FACES DE L'ETOILE
 60   IF( NF2 .GT. 0 ) THEN

C        LA FACE NF2 NON FRONTIERE ET LE POINT NPt FORMENT UN NOUVEAU TETRAEDRE
C        ----------------------------------------------------------------------
         IF( N1TEVI .LE. 0 ) THEN
C           SATURATION DES TETRAEDRES NOTETR
            GOTO 9700
         ENDIF

C        MISE A JOUR DU DERNIER TETRAEDRE OCCUPE
         NUDTETR = MAX( NUDTETR, N1TEVI )
C        MISE A JOUR DU 1-ER TETRAEDRE VIDE
         NT     = N1TEVI
         NOTETR(5,N1TEVI) = ABS( NOTETR(5,N1TEVI) )
         N1TEVI = NOTETR(5,N1TEVI)

C        NS1 NS2 NS3 SOMMETS DE LA FACE DANS LE SENS DIRECT DU TETRAEDRE NT
C        REMPLISSAGE DU TETRAEDRE NT DE SOMMETS CEUX DE LA FACE NF2 ET NPt
         NS1 = NFETOI(2,NF2)
         NS2 = NFETOI(3,NF2)
         NS3 = NFETOI(4,NF2)

         NOTETR(1,NT) = NS1
         NOTETR(2,NT) = NS2
         NOTETR(3,NT) = NS3
         NOTETR(4,NT) = NPt

C        RECHERCHE DES TETRAEDRES OPPOSES AUX FACES DE NT
C        NTOP LE TETRAEDRE EXTERIEUR A LA FACE 1 AU DELA DE LA FACE NF2 de NFETOI
C        ATTENTION: SI FRONTIERE IL PEUT NE PAS EXISTER UN TEL TETRAEDRE
         NTOP = NFETOI(1,NF2)
         IF( NTOP .GT. 0 ) THEN

C           RECHERCHE DU NO NFOP LOCAL A NTOP DE LA FACE NS1 NS2 NS3
            CALL TRI3NO( NFETOI(2,NF2), NSFC )
            CALL NO1F1T( NSFC, NOTETR(1,NTOP), NFOP )
            IF( NFOP .LE. 0 ) THEN

C              CE CAS NE DEVRAIT PAS ARRIVER...
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*, 'tetra1st: Probleme: 3 SOMMETS ',NS1,NS2,NS3
                  PRINT*, 'NON DANS LE TETRAEDRE ',NTOP,' :',
     %                    (NOTETR(kk,NTOP),kk=1,8)
               ELSE
                  PRINT*, 'tetra1st: Problem: 3 VERTICES ',NS1,NS2,NS3
                  PRINT*, 'NOT in TETRAHEDRON ',NTOP,' :',
     %                    (NOTETR(kk,NTOP),kk=1,8)
               ENDIF
               IERR = 12
               GOTO 9900

            ELSE

C              LE TETRAEDRE OPPOSE A LA FACE NFOP DE NTOP EST LE TETRAEDRE NT
               NOTETR( 4+NFOP, NTOP ) = NT

            ENDIF

         ENDIF

C        LE TETRAEDRE OPPOSE A LA FACE 1 DE NT EST LE TETRAEDRE NTOP
         NOTETR(5,NT) = NTOP
         NOTETR(6,NT) = -1
         NOTETR(7,NT) = -1
         NOTETR(8,NT) = -1

C        MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
C        TOUS LES TETRAEDRES ONT PU DISPARAITRE
C        CAS DE L'ETOILE REDUITE AUX 2 TETRAEDRES INITIAUX
         N1TETS(NS1) = NT
         N1TETS(NS2) = NT
         N1TETS(NS3) = NT
         N1TETS(NPt) = NT

C        QUALITE ET VOLUME DU TETRAEDRE CREE
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VNT, QNT )

         VOLET1 = VOLET1 + VNT
         QMOY   = QMOY + QNT
         QMIN1  = MIN( QMIN1, QNT )

         IF( VNT .LE. 0.0 .OR. QNT .LE. QTEMED ) THEN
C           NT TETRAEDRE DEGENERE MAIS CREE TOUT DE MEME
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'tetra1st: NPt=',NPt,' NPSOFR=',NPSOFR(NPt),
     %                ' CREATION NOTETR(',NT,')=',
     %                 (NOTETR(kk,NT),kk=1,4),
     %                ' mais DEGENERE V=',VNT,' Q=',QNT,' POURQUOI???'
            ELSE
               PRINT*,'tetra1st: NPt=',NPt,' NPSOFR=',NPSOFR(NPt),
     %                ' CREATES NOTETR(',NT,')=',
     %                 (NOTETR(kk,NT),kk=1,4),
     %                ' but DEGENERATED V=',VNT,' Q=',QNT,' WHY???'
            ENDIF
            print*
            tracte = .true.
         ENDIF

C        LE TETRAEDRE NT EST AJOUTE A L'ETOILE DE NPt
         IF( NBTEET1 .GE. MXETOI ) THEN
C           SATURATION DES TETRAEDRES NTETOI DE L'ETOILE
            GOTO 9800
         ENDIF
         NBTEET1 = NBTEET1 + 1
         NTETOI( NBTEET1 ) = NT

         IF( INFACO .NE. 0 ) THEN
C           MISE A JOUR DU NO DE TETRAEDRE DES FACES DE LEFACO et N1FASC
            DO J=1,4
C              LA FACE J DE NT EST ELLE DANS LEFACO?
               CALL NULEFT( J, NT, NOTETR, MXFACO, LEFACO,  NFLEFA )
               IF( NFLEFA .GT. 0 ) THEN
C                 OUI: LE TETRAEDRE NT CONTIENT LA FACE NFLEFA DE LEFACO
                  LEFACO( 11, NFLEFA ) = NT
C                 UN NUMERO DE FACE LEFACO AUX 3 SOMMETS
                  DO N=1,3
                     N1FASC( LEFACO(N,NFLEFA) ) = NFLEFA
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

C        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
         NF2 = NFETOI(5,NF2)
         GOTO 60

      ENDIF

C     QUALITE MOYENNE DES NBTEET1-NBTEET0 TETRAEDRES CREES
      QMOY = QMOY / (NBTEET1-NBTEET0)

      IF( QMIN1 .LE. 0 ) THEN
C        VISUALISATION DES TETRAEDRES CREES DONT 1 AU MOINS DEGENERE
         TRACTE = .TRUE.
         KTITRE='tetra1st:                 NOUVEAUX TETRAEDRES DE QUALIT
     %E MIN          . Volume: ' // KNMVOLU
         WRITE(KTITRE(11:17),'(I7)'  ) NBTEET1-NBTEET0
         WRITE(KTITRE(62:67),'(F6.3)') QMIN1
         CALL SANSDBL( KTITRE, L )
         CALL TRACEETOILE( KTITRE(1:L), PTXYZD, NOTETR, NPt,
     %                     NBTEET1-NBTEET0, NTETOI(NBTEET0+1) )
      ENDIF


C     DESTRUCTION DES TETRAEDRES INITIAUX DU TABLEAU NTETOI DE BOULE
C     CIRCONSCRITE CONTENANT LE SOMMET NPt ET ETOILANT CE POINT NPt
C     --------------------------------------------------------------
      QMIN0  = 2
      VOLET0 = 0D0
      DO I = 1, NBTEET0
         NT = NTETOI( I )
         IF( NT .GT. 0 ) THEN

C           VOLUME et QUALITE DU TETRAEDRE NT
            CALL QUATETD( PTXYZD( 1, NOTETR( 1, NT ) ),
     %                    PTXYZD( 1, NOTETR( 2, NT ) ),
     %                    PTXYZD( 1, NOTETR( 3, NT ) ),
     %                    PTXYZD( 1, NOTETR( 4, NT ) ),
     %                    ARMIN, ARMAX, SURFTR, VNT, QNT )

            VOLET0 = VOLET0 + VNT
            IF( QNT .LT. QMIN0 ) QMIN0 = QNT

C           MISE A JOUR DE N1TETS DES 4 SOMMETS AVANT NOUVEAUX TETRAEDRES
C           CELA EST NECESSAIRE POUR ELIMINER LES SOMMETS DUS A L'AJOUT
C           DES TETRAEDRES ENCOCHES
            DO K=1,4
               N1TETS( NOTETR(K,NT) ) = 0
            ENDDO

C           SUPPRESSION DU TETRAEDRE NT DE LEFACO
            CALL SUTELEFA( NT, NOTETR, INFACO, MXFACO, LEFACO )

C           DESTRUCTION DE NT DU TABLEAU NOTETR ET
C           NT DEVIENT LE PREMIER TETRAEDRE VIDE
            DO NF=1,8
               NOTETR(NF,NT) = 0
            ENDDO
            NOTETR(5,NT) = N1TEVI
            N1TEVI = NT

         ENDIF
      ENDDO

C     COMPLETION DES CHAINAGES DES TETRAEDRES OPPOSES CREES DANS L'ETOILE
C     -------------------------------------------------------------------
      IF( TRACTE ) THEN
      KTITRE='tetra1st: NPt=          NBTEET1=        TETRAEDRES avant C
     %OMPLETION des TETRAEDRES OPPOSES. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(33:37),'(I5)') NBTEET1-NBTEET0
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET1-NBTEET0, NTETOI(NBTEET0+1) )
      ENDIF

      CALL MJOPTE( NBTEET2-NBTEET0, NTETOI(NBTEET0+1), N1TETS, NOTETR,
     %             NUDTETR, N1TEVI, PTXYZD, NBFANR )
      IF( NBFANR .EQ. 0 ) THEN
C        TOUS LES TETRAEDRES OPPOSES SONT RETROUVES
         GOTO 100
      ENDIF


C-----------------------------------------------------------------------
C     A SUPPRIMER PLUS TARD SI AUCUNE EXECUTION DE CES INSTRUCTIONS...
C     IL RESTE DES FACES DE TETRAEDRE OPPOSE NON RETROUVE
C     ---------------------------------------------------
      PRINT*,'tetra1st: Probleme NPt=',NPt,' avec',NBFANR,
     %       ' FACES de TETRAEDRE OPPOSE NON RETROUVE'

C     ESSAI DE RETROUVER LES TETRAEDRES OPPOSES OUBLIES....
      NBTE = NBTEET2
 70   IF( NBTE .GT. NBTEET0 ) THEN

C        LE HAUT DE LA PILE
         NT = NTETOI( NBTE )
C        LE TETRAEDRE EST DEPILE
         NBTE = NBTE - 1
         IF( NT .LE. 0 ) GOTO 70

c        verification
         if( notetr(1,NT) .EQ. 0 ) THEN
            PRINT*,'tetra1st: Anomalie NPt=',NPt,' XYZD=',
     %             (PTXYZD(mmm,NPt),mmm=1,4)
            PRINT*,'tetra1st: Anomalie NOTETR(',NT,')=',
     %             (NOTETR(mmm,NT),mmm=1,8)
         ENDIF

C        QUEL EST LE TETRAEDRE OPPOSE AUX FACES 1 2 3 4
C             ET DE SOMMETS NS1 NS2 NS3 = NSFC
         DO 90 NF=1,4

C           LE TETRAEDRE OPPOSE A LA FACE NF DU TETRAEDRE NT
 85         NTOP = NOTETR( 4+NF, NT )
            IF( NTOP .EQ. 0 ) GOTO 90

C           LE NUMERO DES 3 SOMMETS DE LA FACE NF DU TETRAEDRE NT
            NSFC(1) = NOTETR( NOSOFATE(1,NF), NT )
            NSFC(2) = NOTETR( NOSOFATE(2,NF), NT )
            NSFC(3) = NOTETR( NOSOFATE(3,NF), NT )
C           TRI CROISSANT DE SES 3 SOMMETS
            CALL TRI3NO( NSFC, NSFC )

            IF( NTOP .LT. 0 ) THEN

C              FACE NON TRAITEE. TETRAEDRE OPPOSE INCONNU
C              RECHERCHE DE LA FACE NFSC PARMI LES TETRAEDRES CREES DE L'ETOILE
               DO 80 I=1,NBTE

                  NT1 = NTETOI(I)
                  IF( NT1 .LE. 0 ) GOTO 80
C                 RECHERCHE DE LA FACE NSFC DANS CE TETRAEDRE NT1
                  CALL NO1F1T( NSFC, NOTETR(1,NT1), NF1 )
                  IF( NF1 .LE. 0 ) GOTO 80

C                 LA FACE NF DE NT EST LA FACE NF1 DE NT1
C                 LE NUMERO NOTETR DU TETRAEDRE A EFFACER
                  NT2 = NOTETR( 4+NF1, NT1 )

                  IF( NT2 .GT. 0 .AND. NT2 .NE. NT ) THEN

C                    1 FACE COMMUNE A 3 TETRAEDRES!...
                     DO L=NBTEET0+1,NBTEET2
                        IF( NT2 .EQ. ABS(NTETOI(L)) ) GOTO 86
                     ENDDO
 86                  print*,'tetra1st: nbt=',nbte+1,' notetr(',nt,')=',
     %                      (notetr(kk,nt),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'tetra1st: nbt=',I,' notetr(',nt1,')=',
     %                      (notetr(kk,nt1),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'tetra1st: nbt=',L,'notetr(',nt2,')=',
     %                      (notetr(kk,nt2),kk=1,8),' AVEC FACE TRIPLE?'
                     print*

                     tracte = .true.
      KTITRE='tetra1st: NPt=          NBTEET2=        TETRAEDRES avant C
     %OMPLETION des TETRAEDRES OPPOSES. Volume ' // KNMVOLU
                     WRITE(KTITRE(15:22),'(I8)') NPt
                     WRITE(KTITRE(33:37),'(I5)') NBTEET2-NBTEET0
                     CALL SANSDBL( KTITRE, NBC )
                     CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR,
     %                                 NPt, NBTEET2-NBTEET0,
     %                                 NTETOI(NBTEET0+1) )

                  ENDIF

C                 NT1 et NT TETRAEDRES OPPOSES
                  NOTETR( 4+NF , NT  ) = NT1
                  NOTETR( 4+NF1, NT1 ) = NT

                  IF( NOTETR(5,NT1) .GE. 0 .AND.
     %                NOTETR(6,NT1) .GE. 0 .AND.
     %                NOTETR(7,NT1) .GE. 0 .AND.
     %                NOTETR(8,NT1) .GE. 0 ) THEN
C                     FACE TOTALEMENT TRAITEE => SUPPRESSION DE L'ETOILE
                      NTETOI(I) = -NTETOI(I)
                  ENDIF

                  GOTO 90

 80            ENDDO

            ELSE

C              FACE A PRIORI TRAITEE. VERIFICATION
C              RECHERCHE DE LA FACE NSFC DANS CE TETRAEDRE NTOP
               CALL NO1F1T( NSFC, NOTETR(1,NTOP), NFOP )
               IF( NFOP .LE. 0 ) GOTO 90
               NTE = NOTETR( 4+NFOP, NTOP )
               IF( NTE .NE. NT ) THEN
C                 PROBLEME: LE TETRAEDRE OPPOSE PAR LA FACE NF DE NT
C                           N'EST PAS A JOUR
                  NOTETR( 4+NF, NT ) = -1
C                 RETOUR EN RECHERCHE DE TETRAEDRE OPPOSE
                  GOTO 85
               ENDIF

            ENDIF

 90      ENDDO

C        RETOUR EN HAUT DE PILE
         GOTO 70
      ENDIF

C     MISE A JOUR DU NUMERO POSITIF DES TETRAEDRES DE L'ETOILE de NPt
      DO I = NBTEET0+1,NBTEET1
         NTE = ABS( NTETOI( I ) )
         NTETOI( I ) = NTE
      ENDDO

C     VERIFICATION DES TETRAEDRES OPPOSES. A SUPPRIMER PLUS TARD...
      CALL VEOPTE( NBTEET1-NBTEET0, NTETOI(NBTEET0+1), NOTETR, PTXYZD,
     %             NBFANR )
      IF( NBFANR .GT. 0 ) THEN
C        IL RESTE DES FACES DE TETRAEDRE OPPOSE NON RETROUVE
         PRINT*,'tetra1st: Probleme NPt=',NPt,' avec',NBFANR,
     %          ' FACES de TETRAEDRE OPPOSE NON RETROUVE'
      ENDIF
C-----------------------------------------------------------------------


C     ICI TOUS LES TETRAEDRES OPPOSES ONT ETE RETROUVES
C     COMPRESSION DU NUMERO NOTETR DU TABLEAU NTETOI DES TETRAEDRES CREES
C     EVENTUEL NUMERO DE VOLUME DES TETRAEDRES CREES DE SOMMET NPt
C     -------------------------------------------------------------------
 100  DO I = NBTEET0+1, NBTEET1

C        MISE EN DEBUT DE NTETOI DU NO DES TETRAEDRES CREES
         NTE                 = NTETOI( I )
         NTETOI( I-NBTEET0 ) = NTE

         IF( IVOLTE .NE. 0 ) THEN
            NVOLTE( NTE ) = NOVOLU
         ENDIF

      ENDDO
C     NOMBRE DE TETRAEDRES CREES DE NO DANS NTETOI
      NBTEET = NBTEET1 - NBTEET0


C     AMELIORER LA QUALITE DES TETRAEDRES DE L'ETOILE DU POINT NPt
C     DEFINIE PAR LES NBTEET TETRAEDRES NOTEET DU TABLEAU NOTETR
C     PAR ECHANGE 2T->3T ou mT->2m-4T des TETRAEDRES
C     ------------------------------------------------------------
      CALL EC2TMTNT( NPt,    PTXYZD, INFACO,  MXFACO, LEFACO,
     %               IVOLTE, NVOLTE, MXETOI,  NBTEET, NTETOI,
     %               NOTETR, N1TEVI, NUDTETR, N1TETS,
     %               NBEC2T3T, NBECMTNT )


C     CHAINAGE DU SOMMET NPt TETRAEDRISE DANS HEXAPAVE
C     ------------------------------------------------
      CALL NUPAVEST( NPt,      PTXYZD,
     %               HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV )

ccc      PRINT*,'tetra1st: NPt=',NPt,
ccc     %' AVANT NbTETRA=', NBTEET0, ' QMIN0=', QMIN0,' QMIN01=',QMIN01,
ccc     %' APRES NbTETRA=', NBTEET,  ' QMIN1=', QMIN1


C     TRACE DES NBTEET TETRAEDRES DE SOMMETS NPt
C     ------------------------------------------
      IF( TRACTE ) THEN
      KTITRE='tetra1st: NPt=          SOMMET des NBTEET=        TETRAEDR
     %ES FINAUX QMIN1=             . Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)'   ) NPt
      WRITE(KTITRE(43:47),'(I5)'   ) NBTEET
      WRITE(KTITRE(76:90),'(G15.6)') QMIN1
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET, NTETOI )
      ENDIF

      GOTO 9999


C     SATURATION DU TABLEAU NOTETR DES TETRAEDRES DU MAILLAGE
C     -------------------------------------------------------
 9700 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tetra1st: SATURATION DES TETRAEDRES NOTETR MXTETR=',
     %           MXTETR
      ELSE
         PRINT*,'tetra1st: SATURATION of TETRAHEDRA NOTETR MXTETR=',
     %           MXTETR
      ENDIF
      IERR = 3
      GOTO 9999


C     SATURATION DU TABLEAU NTETOI DES TETRAEDRES DE L'ETOILE
C     -------------------------------------------------------
 9800 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tetra1st: SATURATION DES TETRAEDRES NTETOI MXETOI=',
     %           MXETOI
      ELSE
         PRINT*,'tetra1st: SATURATION of STAR TETRAHEDRA NTETOI ARRAY MX
     %ETOI=',MXETOI
      ENDIF
      IERR = 3
      GOTO 9999


C     ABANDON: POINT NPt A SUPPRIMER
C     ------------------------------
 9900 NBTEET = -1
      NOTET1 = 0
C     NPt EST UN SOMMET D'AUCUN TETRAEDRE
      N1TETS( NPt ) = 0


 9999 TRACTE = TRACTE0
      RETURN
      END
