      SUBROUTINE TETRA1PT( NUVOLU, KNMVOLU, QUAMINEX, NPt, NOTEIN,
     %                     MXSOMM,  PTXYZD,NPSOFR,
     %                     HEXAPAVE,NBIPAV,ECHPAV, N1SPAVE,NOPTSUIV,
     %                     MXTETR, N1TEVI, NUDTETR,NOTETR, N1TETS,
     %                     INFACO, MXFACO, LEFACO, N1FASC,
     %                     IVOLTE, NVOLTE,
     %                     MXETOI, N1FEOC, N1FEVI, NFETOI,
     %                     NBTEET, NTETOI, VOLET0, VOLET1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DES FACES SIMPLES DE L'ETOILE DES NBTEET
C -----    TETRAEDRES DE BOULE CIRCONSCRITE CONTENANT LE POINT NPt
C          ET CONSTRUCTION DES NBTEET1 TETRAEDRES DE SOMMET NPt

C ENTREES:
C --------
C NUVOLU  : NUMERO DU VOLUME DANS LE LEXIQUE DES VOLUMES
C KNMVOLU : NOM DU VOLUME A TETRAEDRISER
C QUAMINEX: QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C           C-A-D AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C NOTEIN : NUMERO DU TETRAEDRE DE DEBUT DE RECHERCHE
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
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        Mars 1993
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray          Decembre 2018
C23456...............................................................012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))

      PARAMETER        (MXTECH=1024)
      include"./incl/langue.inc"
      PARAMETER        (RPETIT=-1E38, MXITER=16)
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
     %                  VTE, VOLTET, VOLET0, VOLET1, VTOP, D, DSR,
     %                  VTEMIN, COBARY(4)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NSFC(3), NOTECH(MXTECH)
      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      TRACTE0= TRACTE
      LORBITE= 1

ccc      if( nuvolu.eq.45 .AND. NPt.eq.480 ) then
ccc         tracte = .true.
ccc      endif

      LOCTAE = 0
      NBREPR = 0
      NBENCO = 0
      NOTET0 = NOTEIN

C     ACTUELLEMENT NPt EST UN SOMMET D'AUCUN TETRAEDRE
      N1TETS( NPt ) = 0

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
         print*,'tetra1pt: Volume ',KNMVOLU,' Numero',NUVOLU
         print*,'tetra1pt: COBARY=',COBARY,
     %          ' NPSOFR(',NPt,')=',NPSOFR(NPt),
     %          ' DANS AUCUN TETRAEDRE est il SUPPRIMABLE?'

         IF( NPSOFR(NPt) .EQ. -4 .OR. NPSOFR(NPt) .EQ. 0 ) THEN
C           OUI: LE POINT NPt EST SUPPRIME
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'tetra1pt: LE POINT',NPt,' N''EST PAS TETRAEDRISE'
            ELSE
               PRINT*,'tetra1pt: the POINT',NPt,' IS NOT TETRAHEDRIZED'
            ENDIF
         ELSE
C          NON:  LE POINT NPt  NE PEUT PAS ETRE SUPPRIME
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'tetra1pt: POINT',NPt,' NON SUPPRIMABLE'
            ELSE
               PRINT*,'tetra1pt: POINT',NPt,' CANNOT BE DELETED'
            ENDIF
           PRINT*,'tetra1pt: XYZ',NPt,' =',(PTXYZD(I,NPt),I=1,4)
         ENDIF

         GOTO 9900

      ENDIF
      NOTET0 = NOTET1

C     NUMERO DU VOLUME DU TETRAEDRE NOTET1
      IF( IVOLTE .NE. 0 ) THEN
         NOVOLU = NVOLTE( NOTET1 )
      ELSE
         NOVOLU = 0
      ENDIF


C     FORMATION DE L'ETOILE DES NBTEET TETRAEDRES CONTENANT LE POINT NPt
C     A PARTIR DU TETRAEDRE NOTET1 DE NOTETR ISSU DE R1TCO1P
C     ------------------------------------------------------------------
      CALL ETOIL1PT( NPt,    PTXYZD, NPSOFR,
     %               NOTET1, NOTETR, N1TETS,  COBARY,
     %               NBCBA0, MXETOI, NBTEET0, NTETOI )
      IF( NBTEET0 .LE. 0 ) THEN
         NOTET1 = 0
         GOTO 5
      ENDIF

      KTITRE='tetra1pt: NPt=          ETOILE INITIALE de        TETRAEDR
     %ES. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(44:48),'(I5)') NBTEET0
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET0, NTETOI )


C     NPt SUPPRIMABLE EST IL TROP PROCHE DE L'UN DES SOMMETS DES NTETOI?
C     ------------------------------------------------------------------
      IF( NPSOFR(NPt) .EQ. 0 ) THEN
         DO N=1,NBTEET0
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

C     EXTENSION  DE L'ETOILE AUX TETRAEDRES OPPOSES PAR LES FACES DES NBTEET
C     TETRAEDRES DE L'ETOILE ET DONT LA BOULE CIRCONSCRITE CONTIENT LE POINT
C     ----------------------------------------------------------------------
C     PARCOURS DES NBTEET TETRAEDRES DE L'ETOILE DU POINT NPt
      NBTEET = NBTEET0
      I = 0

 40   I = I + 1
      IF( I .LE. NBTEET ) THEN
C
C        TRAITEMENT DU TETRAEDRE I DE NTETOI
         NT = NTETOI( I )
         IF( NT .LE. 0 ) GOTO 40

         DO 50 J=1,4

C           LE TETRAEDRE OPPOSE PAR LA FACE J DE NT
            NTEOP = NOTETR( J+4, NT )

            IF( NTEOP .GT. 0 ) THEN

C              VOLUME et QUALITE DU TETRAEDRE NTEOP
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTEOP)),
     %                       PTXYZD(1,NOTETR(2,NTEOP)),
     %                       PTXYZD(1,NOTETR(3,NTEOP)),
     %                       PTXYZD(1,NOTETR(4,NTEOP)),
     %                       ARMIN, ARMAX, SURFTR, VTOP, QTOP )

               IF( VTOP .LE. 0D0 .OR. QTOP .LT. QUAMINEX ) THEN
C                 NTEOP N'EST PAS AJOUTE A L'ETOILE
                  GOTO 50
               ENDIF

C              LE TETRAEDRE OPPOSE EST IL DEJA DANS L'ETOILE?
               DO K=1,NBTEET
                  IF( NTEOP .EQ. ABS( NTETOI(K) ) ) GOTO 50
               ENDDO
C              NON: LE TETRAEDRE OPPOSE N'EST PAS DANS L'ETOILE

C              LE POINT NPt EST IL DANS LA BOULE DU TETRAEDRE NTEOP
               CALL PTDSBO( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NTEOP),
     %                      NONOUI, VTOP, DSR )
C              NONOUI=1  Pt   EST     DANS LA BOULE CIRCONSCRITE
C                     0  Pt N'EST PAS DANS LA BOULE CIRCONSCRITE
C                    -1  LE TETRAEDRE EST DE VOLUME TRES FAIBLE
C                        ou NUL ou NEGATIF

               IF( NONOUI .LE. 0 ) THEN
C                 TETRAEDRE NTEOP NON AJOUTE A L'ETOILE
                  GOTO 50
               ENDIF

C              LE POINT NPt EST INTERIEUR OU SUR LA BOULE CIRCONSCRITE
C              DU TETRAEDRE NTEOP QUI EST AJOUTE A L'ETOILE DE NPt
C              SAUF S'IL INDUIT LA PERTE D'UNE FACE FRONTIERE             
               IF( INFACO .NE. 0 ) THEN
C                 LA FACE J DE NT EST ELLE DANS LEFACO?
                  CALL NULEFT( J, NT, NOTETR, MXFACO, LEFACO, NFLE )
                  IF( NFLE .GT. 0 ) THEN
C                    LA FACE EST DANS LEFACO
C                    TETRAEDRE NTEOP NON AJOUTE A L'ETOILE
                     GOTO 50
                  ENDIF
               ENDIF

               IF( IVOLTE .NE. 0 ) THEN
                  IF( NVOLTE(NTEOP) .NE. NOVOLU ) THEN
C                    NT ET NTEOP SONT DANS DES VOLUMES DIFFERENTS
C                    TETRAEDRE NTEOP NON AJOUTE A L'ETOILE
                     GOTO 50
                  ENDIF
               ENDIF

C              NPt EST DANS LA BOULE CIRCONSCRITE DE NTEOP
C              ET LA FACE N'EST PAS DANS LEFACO ou
C              NT ET NTEOP SONT DANS LE MEME VOLUME

cccC              NPt EST IL TROP ELOIGNE DU CENTRE DE LA BOULE DE NTEOP?
ccc               CALL COBATE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NTEOP),
ccc     %                      VTOP, COBAR2, IER )
ccc               CBSOM = ABS( COBAR2(1) ) + ABS( COBAR2(2) )
ccc     %               + ABS( COBAR2(3) ) + ABS( COBAR2(4) )
ccc             test non significatif de l'eloignement
ccc               IF( CBSOM .GT. 8D0 ) THEN
cccC                 TETRAEDRE NTEOP NON AJOUTE A L'ETOILE
cccC                 CAR NTEOP EST TROP ELOIGNE DE NPt
ccc                  PRINT*,'tetra1pt: NPt=',NPt,' NTEOP=',NTEOP,
ccc                  ' COBAR2=',COBAR2,' CBSOM=',CBSOM
ccc                  GOTO 50
ccc               ENDIF

               IF( NPSOFR(NPt) .EQ. 0 ) THEN
C                 NPt SUPPRIMABLE TROP PROCHE DE L'UN DES SOMMETS DE NTEOP?
                  DO K=1,4
                     NS = NOTETR(K,NTEOP)
                     D = SQRT( (PTXYZD(1,NS) - PTXYZD(1,NPt) ) ** 2
     %                       + (PTXYZD(2,NS) - PTXYZD(2,NPt) ) ** 2
     %                       + (PTXYZD(3,NS) - PTXYZD(3,NPt) ) ** 2 )
                     IF( D .LT. PTXYZD(4,NPt)*0.333D0 ) THEN
C                       NPt TROP PROCHE DU SOMMET NS N'EST PAS TETRAEDRISE
                        GOTO 9900
                     ENDIF
                  ENDDO
               ENDIF

C              AJOUT DU TETRAEDRE NTEOP A L'ETOILE DU POINT NPt
               IF( NBTEET .GE. MXETOI ) THEN
C                 SATURATION DU TABLEAU NTETOI DES TETRAEDRES DE L'ETOILE
                  GOTO 9800
               ENDIF
               NBTEET = NBTEET + 1
               NTETOI( NBTEET ) = NTEOP
            ENDIF

 50      ENDDO
         GOTO 40

      ENDIF

C     ELIMINATION DES TETRAEDRES ISOLES HORS DE L'ETOILE
C     (SES 4 TETRAEDRES OPPOSES A SES FACES APPARTIENNENT A AUCUN
C     DES TETRAEDRES DE L'ETOILE)
C     -----------------------------------------------------------
 60   DO 70 I = NBTEET0+1, NBTEET

         NT = NTETOI( I )
         IF( NT .LE. 0 ) GOTO 70

         DO J=1,4
C           LE TETRAEDRE OPPOSE
            NTEOP = NOTETR( J+4, NT )
            IF( NTEOP .GT. 0 ) THEN
               DO K=1,NBTEET
                  IF( NTEOP .EQ. NTETOI(K) ) GOTO 70
               ENDDO
            ENDIF
         ENDDO
C        ICI: LE TETRAEDRE NT N'A PAS DE TETRAEDRE OPPOSE DE L'ETOILE

C        SUPPRESSION DE L'ETOILE DU TETRAEDRE ISOLE NT
         DO J=I+1,NBTEET
            NTETOI(J-1) = NTETOI(J)
         ENDDO
         NBTEET = NBTEET - 1

ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            PRINT*,'tetra1pt: POINT NPt=',NPt,
ccc     %             ' RETRAIT de NTETOI du TETRAEDRE ISOLE', NT
ccc         ELSE
ccc            PRINT*,'tetra1pt: POINT NPt=',NPt,
ccc     %             ' The ISOLATED TETRAHEDRON', NT,
ccc     %             ' is REMOVED of NTETOI'
ccc         ENDIF
ccc         PRINT*,'tetra1pt: TETRAEDRE(',NT,')=',(NOTETR(J,NT),J=1,8)

 70   ENDDO


C     VOLUME DE L'ETOILE DES NBTEET TETRAEDRES DE L'ETOILE A ELIMINER
C     ---------------------------------------------------------------
 75   VOLET0 = 0D0
      VTEMIN = 1D100
      DO I = 1, NBTEET
C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE
         NT = NTETOI( I )
         IF( NT .GT. 0 ) THEN
            VTE = VOLTET( PTXYZD(1,NOTETR(1,NT)),
     %                    PTXYZD(1,NOTETR(2,NT)),
     %                    PTXYZD(1,NOTETR(3,NT)),
     %                    PTXYZD(1,NOTETR(4,NT)) )
            VOLET0 = VOLET0 + VTE
            VTEMIN = MIN( VTEMIN, VTE )
         ENDIF
      ENDDO

      IF( VOLET0 .LE. 0D0 ) THEN
C        ETOILE VIDE
         PRINT*, 'tetra1pt: ETOILE du POINT',NPt,
     %           ' VOLUME NEGATIF ou NUL des NBTEET=',NBTEET,
     %           ' TETRAEDRES. VTeMin=',VTEMIN
C        VERS LA REPRISE
         GOTO 300
      ENDIF
C     ICI IL RESTE AU MOINS UN TETRAEDRE DE BOULE CONTENANT LE POINT NPt


C     FORMATION DES FACES SIMPLES DE L'ETOILE A PARTIR DES NBTEET TETRAEDRES
C     ----------------------------------------------------------------------
C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
      CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

      DO I=1,NBTEET

C        LE TETRAEDRE I DE L'ETOILE
         NTE = NTETOI( I )
         IF( NTE .GT. 0 ) THEN

            DO J=1,4
C              AJOUT ou SUPPRESSION de la FACE J de NTE de l'ETOILE
               CALL AJFAET1( NTE,    J,      NOTETR,
     %                       N1FEOC, N1FEVI, NFETOI, NF )
            ENDDO
         ENDIF

      ENDDO

      KTITRE='tetra1pt: NPt=          LES FACES SIMPLES DE SON ETOILE. V
     %olume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      CALL SANSDBL( KTITRE, NBC )
      CALL TRFETOV1( KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR, PTXYZD)


C     APRES vdtuer c-a-d SUPPRESSION DES TETRAEDRES EXTERIEURS A L'OBJET
C     ESSAI DE BOUCHER LES ENCOCHES DES FACES SIMPLES DE L'ETOILE I.E.
C     RECHERCHER LES COUPLES DE FACES SIMPLES DE L'ETOILE NFETOI
C     D'ARETE COMMUNE ET AYANT UN TETRAEDRE OPPOSE IDENTIQUE ET
C     L'AJOUTER AU TABLEAU NTETOI
C     ATTENTION: LORS DU REMPLISSAGE D'UNE ENCOCHE DES SOMMETS INTERNES NS
C     TETRAEDRISES PEUVENT DISPARAITRE ET N1TETS(NS) DOIT ETRE REMIS A ZERO
C     ---------------------------------------------------------------------
      IF( IVOLTE .EQ. 0 .OR. NBENCO .GT. 0 ) GOTO 120

      NBENCO = 1
      NBT = NBTEET

 89   NF1 = 0
      NF2 = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 90   IF( NF2 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT   = NFETOI(1,NF2)
         NFNT = NFETOI(2,NF2)

C        LE NO DU TETRAEDRE OPPOSE
         NTEOP = NOTETR( 4+NFNT, NT )

         IF( NTEOP .GT. 0 ) THEN

C           VOLUME et QUALITE DU TETRAEDRE NTEOP
            CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP) ),
     %                    PTXYZD( 1, NOTETR(2,NTEOP) ),
     %                    PTXYZD( 1, NOTETR(3,NTEOP) ),
     %                    PTXYZD( 1, NOTETR(4,NTEOP) ),
     %                    ARMIN, ARMAX, SURFTR, VTOP, QTOP )

            IF( VTOP .LE. 0D0 .OR. QTOP .LT. QUAMINEX ) GOTO 95

C           LA QUALITE DE NTEOP N'EST PAS MEDIOCRE
C           EXISTE T IL UNE AUTRE FACE SIMPLE DE TETRAEDRE OPPOSE NTEOP?
            NFF1 = 0
            NFF2 = N1FEOC
 92         IF( NFF2 .GT. 0 ) THEN

               IF( NFF2 .NE. NF2 ) THEN
C                 LE NO DU TETRAEDRE NT2 INTERNE ET NO LOCAL DE LA FACE
                  NT2   = NFETOI( 1, NFF2 )
                  NFNT2 = NFETOI( 2, NFF2 )

C                 LE NO DU TETRAEDRE OPPOSE
                  NTEOP2 = NOTETR( 4+NFNT2, NT2 )
                  IF( NTEOP2 .EQ. NTEOP ) THEN

C                    AJOUT DU TETRAEDRE NTEOP AUX TETRAEDRES
C                    DE L'ETOILE DE NPt
                     DO N=1,NBTEET
                        IF( NTETOI(N) .EQ. NTEOP ) GOTO 95
                     ENDDO
                     NBTEET = NBTEET + 1
                     NTETOI( NBTEET ) = NTEOP

                     DO J=1,4
C                       AJOUT ou SUPPRESSION de la FACE J de NTEOP de l'ETOILE
                        CALL AJFAET1( NTEOP,  J,      NOTETR,
     %                                N1FEOC, N1FEVI, NFETOI, N )
                     ENDDO

                     GOTO 89

                  ENDIF
               ENDIF

C              PASSAGE A LA FACE SIMPLE SUIVANTE
               NFF1 = NFF2
               NFF2 = NFETOI( 5, NFF2 )
               GOTO 92

            ENDIF

         ENDIF

C        PASSAGE A LA FACE SIMPLE SUIVANTE
 95      NF1 = NF2
         NF2 = NFETOI( 5, NF2 )
         GOTO 90

      ENDIF

      IF( NBT .LT. NBTEET ) THEN
ccc         tracte = .true.
ccc         PRINT*,'tetra1pt: AJOUT de',NBTEET-NBT,' TETRAEDRES POUR REMPLI
ccc     %R LES ENCOCHES'
      KTITRE='tetra1pt: NPt=          LES       TETRAEDRES APRES AJOUT D
     %ES TETRAEDRES DES ENCOCHES. Volume ' // KNMVOLU
         WRITE(KTITRE(15:22),'(I8)') NPt
         WRITE(KTITRE(29:33),'(I5)') NBTEET
         CALL SANSDBL( KTITRE, NBC )
         CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                     NBTEET, NTETOI )

C        RETOUR A LA FORMATION DE L'ETOILE DE NPt
         GOTO 75
      ENDIF

C     SUPPRESSION DES TETRAEDRES NON INITIAUX DE L'ETOILE DONT LES FACES
C     SIMPLES PRODUIRAIENT DES TETRAEDRES DE MEDIOCRE QUALITE
C     ------------------------------------------------------------------
 120  KTITRE='tetra1pt: NPt=          LES       TETRAEDRES non initiaux.
     % Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(29:33),'(I5)') NBTEET
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET, NTETOI )

      NF1 = 0
      NF2 = N1FEOC
      VOLET1 = 0D0

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 150  IF( NF2 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NTE   = NFETOI(1,NF2)
         NFNTE = NFETOI(2,NF2)

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NS1 = NOTETR( NOSOFATE(1,NFNTE), NTE )
         NS2 = NOTETR( NOSOFATE(2,NFNTE), NTE )
         NS3 = NOTETR( NOSOFATE(3,NFNTE), NTE )

C        LE TETRAEDRE A FORMER NF2-POINT NPt POSE T IL PROBLEME ?
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         IF( VTE .LE. 0D0 .OR. QTE .LT. QUAMINEX ) THEN

C           LE TETRAEDRE NF2-POINT NPt SERAIT DE MEDIOCRE QUALITE
C           NTE EST SUPPRIME SAUF SI NPt EST SUR L'OCTAEDRE GLOBAL
            NTEOP = NOTETR( 4+NFNTE, NTE )
            IF( NTEOP .EQ. 0 ) THEN

C              LE POINT NPt EST IL SUR L'OCTAEDRE ENGLOBANT ?
               CALL PTSROE( PTXYZD(1,NPt), PTXYZD, LOCTAE )

               IF( LOCTAE .NE. 0 ) THEN
C                 OUI: SUPPRESSION DE LA FACE SIMPLE QUI EST SUR
C                 L'OCTAEDRE ENGLOBANT MAIS PAS DU TETRAEDRE NTE
C                 NF1 PRECEDE NF2 DANS LE TABLEAU NFETOI
                  IF( NF1 .GT. 0 ) THEN
C                    NF2,PRECEDEE DE NF1,N'EST PAS LA PREMIERE DE L'ETOILE
                     NFETOI(5,NF1) = NFETOI(5,NF2)
C                    LA FACE NF2 DEVIENT LA PREMIERE VIDE DE L'ETOILE
                     NFETOI(5,NF2) = N1FEVI
                     N1FEVI = NF2
C                    POINTEUR POUR ARRIVER A LA FACE OCCUPEE SUIVANTE
                     NF2 = NF1
                     GOTO 190
                  ELSE
C                    LA FACE NF2=N1FEOC EST LA 1-ERE DE L'ETOILE
                     NF1    = NFETOI(5,N1FEOC)
                     NFETOI(5,N1FEOC) = N1FEVI
                     N1FEVI = N1FEOC
                     N1FEOC = NF1
C                    POINTEUR POUR REPARCOURIR L'ETOILE
                     GOTO 120
                  ENDIF
               ENDIF

            ELSE

C              IL EXISTE UN TETRAEDRE OPPOSE NTEOP PAR CETTE FACE NFNTE DE NTE
C              SUPPRESSION DU TETRAEDRE NTE DE L'ETOILE
C              HORS LES NBTEET0 TETRAEDRES INITIAUX
               DO L=NBTEET0+1,NBTEET
                  IF( NTETOI(L) .EQ. NTE ) THEN
                     DO K=L+1,NBTEET
                        NTETOI(K-1) = NTETOI(K)
                     ENDDO
                     NBTEET = NBTEET - 1

ccc                     PRINT*
ccc                     PRINT*,'tetra1pt: Pt',NPt,
ccc     %                  ' RETRAIT de son ETOILE du TETRAEDRE(',NTE,')=',
ccc     %                   (NOTETR(kkk,NTE),kkk=1,8)
ccc                PRINT*,'POUR NE PAS CREER le TETRAEDRE',NS1,NS2,NS3,NPt,
ccc     %                 ' de Volume=',VTE,' Qualite=',QTE

ccc             PRINT*,('PTXYZD(',kkk,',',NS1,')=',PTXYZD(kkk,NS1),kkk=1,4)
ccc             PRINT*,('PTXYZD(',kkk,',',NS2,')=',PTXYZD(kkk,NS2),kkk=1,4)
ccc             PRINT*,('PTXYZD(',kkk,',',NS3,')=',PTXYZD(kkk,NS3),kkk=1,4)
ccc             PRINT*,('PTXYZD(',kkk,',',NPt,')=',PTXYZD(kkk,NPt),kkk=1,4)

                     GOTO 60
                  ENDIF
               ENDDO

            ENDIF

         ENDIF

C        TETRAEDRE CORRECT. PASSAGE A LA FACE SIMPLE SUIVANTE
 190     NF1 = NF2
         NF2 = NFETOI( 5, NF2 )
         VOLET1 = VOLET1 + VTE
         GOTO 150

      ENDIF


C     ETOILE CORRECTE c-a-d VOLUMES DE L'ETOILE EGAUX AVANT et APRES?
C     ---------------------------------------------------------------
      IF( ABS(VOLET1-VOLET0) .GT. 1D-5 * ABS(VOLET0) ) THEN

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tetra1pt: POINT',NPt,
     %    ' TROP GRANDE DIFFERENCE des VOLUMES V0=',VOLET0,' V1=',VOLET1
         ELSE
            PRINT*,'tetra1pt: at POINT ',NPt,
     %    ' TOO GREAT DIFFERENCE of VOLUMES V0=',VOLET0,' V1=',VOLET1
         ENDIF
         GOTO 300

      ENDIF

C     OUI: ETOILE CORRECTE DU POINT NPt A TETRAEDRISER
      GOTO 500


C     ======================================================================
C     PROBLEME SUR L'ETOILE => REPRISE AVEC LES NBTEET0 TETRAEDRES DE NTETOI
C     ======================================================================
 300  IF( LANGAG .EQ. 0 ) THEN
         PRINT 10300, NPt,NPSOFR(NPt),(PTXYZD(I,NPt),I=1,4),NBTEET
      ELSE
         PRINT 20300, NPt,NPSOFR(NPt),(PTXYZD(I,NPt),I=1,4),NBTEET
      ENDIF
10300 FORMAT(' tetra1pt: REPRISE du POINT ',I9,' NPSOFR=',I9,
     %       ' XYZD=',4G15.7,' NBTEET=',I4)
20300 FORMAT(' tetra1pt: RENEWAL of POINT ',I9,' NPSOFR=',I9,
     %       ' XYZD=',4G15.7,' NBTEET=',I4)

C     TENTATIVE POUR TROUVER UNE ETOILE DU POINT NPt
C     SORTIE EN PRENANT L'ETOILE MINIMALE DES TETRAEDRES CONTENANT NPt
      IF( IVOLTE .EQ. 0 ) THEN
         CALL R1TCO1P0( NPt,    PTXYZD, NOTET0, NOTETR, NUDTETR,
     %                  NOTET1, COBARY )
      ELSE
         CALL R1TCO1P1( NPt,     MXSOMM, PTXYZD, 1,
     %                  HEXAPAVE,NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                  NOTET0,  MXTETR, NOTETR, NUDTETR, N1TETS,
     %                  MXETOI,  NTETOI, MXTECH, NOTECH,
     %                  NOTET1,  COBARY )
      ENDIF

      IF( NOTET1 .LE. 0 ) THEN
         GOTO 9900
      ENDIF

      CALL ETOIL1PT( NPt,    PTXYZD, NPSOFR,
     %               NOTET1, NOTETR, N1TETS, COBARY,
     %               NBCBA0, MXETOI, NBTEET, NTETOI )
      IF( NBTEET .LE. 0 ) THEN
         GOTO 9900
      ENDIF

C     REINITIALISATION A VIDE DES FACES OCCUPEES DE NFETOI
 310  CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

C     VOLUME DE L'ETOILE DES NBTEET TRIANGLES DE L'ETOILE A ELIMINER
      QMIN0  = 2.
      VOLET0 = 0D0
      DO I = 1, NBTEET

C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE
         NT = NTETOI( I )

         IF( NT .GT. 0 ) THEN

C           VOLUME et QUALITE DU TETRAEDRE NT
            CALL QUATETD( PTXYZD(1,NOTETR(1,NT)),
     %                    PTXYZD(1,NOTETR(2,NT)),
     %                    PTXYZD(1,NOTETR(3,NT)),
     %                    PTXYZD(1,NOTETR(4,NT)),
     %                    ARMIN, ARMAX, SURFTR, VTE, QTE )

            VOLET0 = VOLET0 + VTE
            IF( QTE .LT. QMIN0 ) THEN
               QMIN0 = QTE
            ENDIF

            DO J=1,4
C              AJOUT ou SUPPRESSION de la FACE J de NT de l'ETOILE
C              PAS D'ARRET AUX FACES LEFACO OU ENTRE 2 VOLUMES
               CALL AJFAET1( NT,  J, NOTETR,
     %                       N1FEOC, N1FEVI, NFETOI, NFS )
            ENDDO
         ENDIF

      ENDDO

      IF( VOLET0 .LE. 0D0 ) THEN
C        ETOILE VIDE
         PRINT*,'tetra1pt: ETOILE du POINT',NPt,
     %          ' VOLUME NEGATIF ou NUL des NBTEET=',NBTEET,
     %          ' TETRAEDRES. QMIN0=',QMIN0
C        ABANDON
         GOTO 9900
      ENDIF


C     ICI IL RESTE AU MOINS UN TETRAEDRE DE BOULE CONTENANT LE POINT NPt
      KTITRE='tetra1pt: NPt=          NBTEET=        TETRAEDRES avant VE
     %RIFICATION des TETRAEDRES A CREER. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(32:36),'(I5)') NBTEET
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET, NTETOI )

C     LE VOLUME DES TETRAEDRES FORMES A PARTIR DE NPt ET
C     DES FACES SIMPLES EST IL TOUJOURS CORRECT ?
C     --------------------------------------------------
 320  NF1 = 0
      NF2 = N1FEOC
      VOLET1 = 0D0
      QMIN1  = 2.

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 330  IF( NF2 .GT. 0 ) THEN

C        LE NO DU TETRAEDRE NT INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT   = NFETOI(1,NF2)
         NFNT = NFETOI(2,NF2)

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NS1 = NOTETR( NOSOFATE(1,NFNT), NT )
         NS2 = NOTETR( NOSOFATE(2,NFNT), NT )
         NS3 = NOTETR( NOSOFATE(3,NFNT), NT )

C        LE TETRAEDRE A FORMER NF2-POINT NPt POSE T IL PROBLEME ?
         CALL QUATETD( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                 PTXYZD( 1, NS3 ), PTXYZD( 1, NPt ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         IF( QMIN1 .LT. QTE ) THEN
            QMIN1 = QTE
         ENDIF

         IF( VTE .LE. 0D0 .OR. QTE .LT. QUAMINEX ) THEN

C           LE TETRAEDRE NF2-POINT NPt EST DEGENERE
            NTEOP = NOTETR( 4+NFNT, NT )
            IF( NTEOP .EQ. 0 ) THEN

C              LE POINT NPt EST IL SUR L'OCTAEDRE ENGLOBANT ?
               CALL PTSROE( PTXYZD(1,NPt), PTXYZD, LOCTAE )

               IF( LOCTAE .NE. 0 ) THEN
C                 OUI: SUPPRESSION DE LA FACE SIMPLE QUI EST SUR
C                 L'OCTAEDRE ENGLOBANT MAIS PAS DU TETRAEDRE NT
C                 NF1 PRECEDE NF2 DANS LE TABLEAU NFETOI
                  IF( NF1 .GT. 0 ) THEN
C                    NF2,PRECEDEE DE NF1,N'EST PAS LA PREMIERE DE L'ETOILE
                     NFETOI(5,NF1) = NFETOI(5,NF2)
C                    LA FACE NF2 DEVIENT LA PREMIERE VIDE DE L'ETOILE
                     NFETOI(5,NF2) = N1FEVI
                     N1FEVI = NF2
C                    POINTEUR POUR ARRIVER A LA FACE OCCUPEE SUIVANTE
                     GOTO 340
                  ELSE
C                    LA FACE NF2=N1FEOC EST LA 1-ERE DE L'ETOILE
                     NF1    = NFETOI(5,N1FEOC)
                     NFETOI(5,N1FEOC) = N1FEVI
                     N1FEVI = N1FEOC
                     N1FEOC = NF1
C                    POINTEUR POUR REPARCOURIR L'ETOILE
                     GOTO 320
                  ENDIF
               ENDIF

            ELSE

C              IL EXISTE UN TETRAEDRE OPPOSE PAR CETTE FACE NFNT DE NT
               PRINT*
               PRINT*,'tetra1pt: PT',NPt,' RETRAIT DE SON ETOILE DU TETR
     %AEDRE(',NT,')=',(NOTETR(kkk,NT),kkk=1,8),' NPSOFR=',NPSOFR(NPt)
               PRINT*,'POUR NE PAS CREER le TETRAEDRE',NS1,NS2,NS3,NPt,
     %                ' de Volume=',VTE,' Qualite=',QTE,' QMIN1=',QMIN1

ccc               PRINT*,('PTXYZD(',kkk,',',NS1,')=',PTXYZD(kkk,NS1),kkk=1,4)
ccc               PRINT*,('PTXYZD(',kkk,',',NS2,')=',PTXYZD(kkk,NS2),kkk=1,4)
ccc               PRINT*,('PTXYZD(',kkk,',',NS3,')=',PTXYZD(kkk,NS3),kkk=1,4)
ccc               PRINT*,('PTXYZD(',kkk,',',NPt,')=',PTXYZD(kkk,NPt),kkk=1,4)

C              TRACE DES FACES SIMPLES DE L'ETOILE
ccc               tracte = .true.
         KTITRE='tetra1pt: NPt=          PROBLEME sur les FACES SIMPLES.
     % Volume ' // KNMVOLU
               WRITE(KTITRE(15:22),'(I8)') NPt
               CALL SANSDBL( KTITRE, NBC )
               CALL TRFETOV1(KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR,
     %                       PTXYZD)

               NBREPR = NBREPR + 1
               IF( NBREPR .EQ. 1 ) THEN

C                 TENTATIVE D'ELIMINATION DE LA FACE INCORRECTE
                  K = 0
                  DO J=1,4
                     IF( ABS(COBARY(J)) .LE. 1D-3 ) K = K + 1
                  ENDDO

                  IF( K .LT. 2 ) THEN
C                    LE POINT NPt EST INTERNE OU SUR UNE FACE
                     NBTEET = 1
                     GOTO 310
                  ENDIF

               ENDIF

               IF( NPSOFR(NPt) .EQ. -4 .OR. NPSOFR(NPt) .EQ. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'tetra1pt: POINT',NPt,' SUPPRIME. Position Fron
     %tiere',NPSOFR(NPt)
                  ELSE
                  PRINT*,'tetra1pt: POINT',NPt,' DELETED. BOUNDARY LOCAT
     %ION',NPSOFR(NPt)
                  ENDIF
                  GOTO 9900
               ENDIF

               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'tetra1pt: la TETRAEDRISATION du POINT',NPt,
     %                   ' de Position Frontiere',NPSOFR(NPt),
     %                   ' CREE UN TETRAEDRE de VOLUME=',VTE,
     %                   ' et QUALITE=',QTE
               ELSE
                  PRINT*,'tetra1pt: TETRAHEDRIZATION of POINT',NPt,
     %                   ' BOUNDARY LOCATION',NPSOFR(NPt),
     %                   ' CREATES a TETRAHEDRON of VOLUME',VTE,
     %                   ' and QUALITY=',QTE
               ENDIF

            ENDIF

         ENDIF

C        TETRAEDRE CORRECT. PASSAGE A LA FACE SIMPLE SUIVANTE
 340     NF1 = NF2
         NF2 = NFETOI( 5, NF2 )
         VOLET1 = VOLET1 + VTE
         GOTO 330

      ENDIF

C     FIN DE LA REPRISE DE L'ETOILE DE NPt  ------------------------------


C     =====================================================================
C     CONSTRUCTION des TETRAEDRES de l'ETOILE du SOMMET NPt
C     =====================================================================
C     TRACE DES FACES SIMPLES DE L'ETOILE
 500  KTITRE='tetra1pt: NPt=          LES FACES SIMPLES FINALES DE L''ET
     %OILE avant TETRAEDRISATION de NPt. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      CALL SANSDBL( KTITRE, NBC )
      CALL TRFETOV1( KTITRE(1:NBC), NPt, N1FEOC, NFETOI, NOTETR, PTXYZD)

C     PASSAGE DE NFETOI VERSION 1 => VERSION 2
C     LE NO DE TETRAEDRE INTERNE ETOILE  => NO TETRAEDRE AU DELA DE LA FACE
C     LE NO FACE LOCAL DANS LE TETRAEDRE => NO 1-ER SOMMET DE LA FACE
C     LE NO INUTILISE 3                  => NO 2-ME SOMMET DE LA FACE
C     LE NO INUTILISE 4                  => NO 3-ME SOMMET DE LA FACE
C     ----------------------------------------------------------------------
      NBFETO = 0
      NF1    = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 510  IF( NF1 .GT. 0 ) THEN

C        NFETOI VERSION 1.  UNE FACE SIMPLE DE L'ETOILE DE PLUS
         NBFETO = NBFETO + 1

C        LE NO DU TETRAEDRE NT1 INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT1   = NFETOI( 1, NF1 )
         NFNT1 = ABS( NFETOI( 2, NF1 ) )

C        NFETOI VERSION 2
C        LE NUMERO DU TETRAEDRE NTEOP AU DELA DE LA FACE NFNT1 DE NT1
C        QUI PEUT ETRE NUL SI LA FACE EST FRONTIERE
         NTEOP = NOTETR( 4+NFNT1, NT1 )
         NFETOI( 1, NF1 ) = NTEOP

C        NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
         NS1 = NOTETR( NOSOFATE(1,NFNT1), NT1 )
         NS2 = NOTETR( NOSOFATE(2,NFNT1), NT1 )
         NS3 = NOTETR( NOSOFATE(3,NFNT1), NT1 )

C        STOCKAGE DANS NFETOI(2:4,NF1)
         NFETOI( 2, NF1 ) = NS1
         NFETOI( 3, NF1 ) = NS2
         NFETOI( 4, NF1 ) = NS3

C        NFOP NO LOCAL A NTEOP DE LA FACE NFNT1 DE NT1
         CALL NOFAOP( NFNT1, NT1, NOTETR,  NFOP, NTEOP )
         IF( NFOP .GT. 0 ) THEN
C           NT1 VA ETRE SUPPRIME -> TETRAEDRE OPPOSE INCONNU
            NOTETR( 4+NFOP, NTEOP ) = -1
            NOTETR( 4+NFNT1, NT1  ) = -1
         ENDIF

C        PASSAGE A LA FACE SUIVANTE
         NF1 = NFETOI( 5, NF1 )
         GOTO 510

      ENDIF

C     DESTRUCTION DES TETRAEDRES INITIAUX DU TABLEAU NTETOI DE BOULE
C     CIRCONSCRITE CONTENANT LE SOMMET NPt ET ETOILANT CE POINT NPt
C     --------------------------------------------------------------
      DO I=1,NBTEET
         NT1 = NTETOI( I )
         IF( NT1 .GT. 0 ) THEN
C           SUPPRESSION DU TETRAEDRE NT1 DE LEFACO
            CALL SUTELEFA( NT1, NOTETR, INFACO, MXFACO, LEFACO )

C           DESTRUCTION DE NT1 DU TABLEAU NOTETR ET
C           NT1 DEVIENT LE PREMIER TETRAEDRE VIDE
            DO NF=1,8
               NOTETR(NF,NT1) = 0
            ENDDO
            NOTETR(5,NT1) = N1TEVI
            N1TEVI = NT1
         ENDIF
      ENDDO

C     DECLARATION DES TETRAEDRES ETOILANT LE POINT NPt
C     ------------------------------------------------
      QUAMIN  = 2.
      QUAMOY  = 0.
      NBTEET1 = 0
      NBTEET2 = NBFETO
      NF2     = N1FEOC
C     BOUCLE SUR LES FACES DE L'ETOILE

 560  IF( NF2 .GT. 0 ) THEN

C        LA FACE NF2 ET LE POINT NPt FORMENT UN NOUVEAU TETRAEDRE
         IF( N1TEVI .LE. 0 ) THEN
C           SATURATION DES TETRAEDRES NOTETR
            GOTO 9700
         ENDIF

C        RECUPERATION ET MISE A JOUR DU 1-ER TETRAEDRE VIDE
         NT     = N1TEVI
         N1TEVI = NOTETR(5,N1TEVI)

C        MISE A JOUR DU DERNIER TETRAEDRE OCCUPE
         NUDTETR = MAX( NUDTETR, NT )

C        REMPLISSAGE DU TETRAEDRE NT DE SOMMETS CEUX DE LA FACE NF2 ET NPt
         NS1 = NFETOI(2,NF2)
         NS2 = NFETOI(3,NF2)
         NS3 = NFETOI(4,NF2)

C        NS1 NS2 NS3 SOMMETS DE LA FACE DANS LE SENS DIRECT DU TETRAEDRE NT
         NOTETR(1,NT) = NS1
         NOTETR(2,NT) = NS2
         NOTETR(3,NT) = NS3
         NOTETR(4,NT) = NPt

C        RECHERCHE DES TETRAEDRES OPPOSES AUX FACES DE NT
C        NTEOP LE TETRAEDRE EXTERIEUR A LA FACE 1 AU DELA DE LA FACE NF2 de NFETOI
C        ATTENTION: SI FRONTIERE IL PEUT NE PAS EXISTER UN TEL TETRAEDRE
         NTEOP = NFETOI(1,NF2)
         IF( NTEOP .GT. 0 ) THEN

C           STOCKAGE DU TETRAEDRE OPPOSE POUR RESOUDRE LES TETRAEDRES OPPOSES
            DO K=1,NBTEET2
               IF( NTETOI(K) .EQ. NTEOP ) GOTO 565
            ENDDO
            IF( NBTEET2 .GE. MXETOI ) THEN
C              SATURATION DES TETRAEDRES NTETOI DE L'ETOILE
               GOTO 9800
            ENDIF
            NBTEET2 = NBTEET2 + 1
            NTETOI( NBTEET2 ) = NTEOP

C           RECHERCHE DU NO NFOP LOCAL A NTEOP DE LA FACE NS1 NS2 NS3
            CALL TRI3NO( NFETOI(2,NF2), NSFC )
            CALL NO1F1T( NSFC, NOTETR(1,NTEOP), NFOP )
            IF( NFOP .LE. 0 ) THEN
C              CE CAS NE DEVRAIT PAS ARRIVER
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*, 'tetra1pt: 3 SOMMETS ',NS1,NS2,NS3
                  PRINT*, 'NON DANS LE TETRAEDRE ',NTEOP,' :',
     %                    (NOTETR(kk,NTEOP),kk=1,8)
               ELSE
                  PRINT*, 'tetra1pt: 3 VERTICES ',NS1,NS2,NS3
                  PRINT*, 'NOT in TETRAHEDRON ',NTEOP,' :',
     %                    (NOTETR(kk,NTEOP),kk=1,8)
               ENDIF
               IERR = 12
               GOTO 9900
            ELSE
               NOTETR( 4+NFOP, NTEOP ) = NT
            ENDIF

         ENDIF

 565     NOTETR(5,NT) = NTEOP
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
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )

         QUAMOY = QUAMOY + QTE
         QUAMIN = MIN( QUAMIN, QTE )

         IF( VTE .LE. 0.0 ) THEN
C           NT TETRAEDRE DEGENERE MAIS CREE TOUT DE MEME
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'tetra1pt: TETRAEDRE',NT,' DEGENERE au POINT',NPt,
     %                'QUALITE du TETRAEDRE=',QTE,' VOLUME=',VTE
            ELSE
              PRINT*,'tetra1pt: DEGENERATED TETRAHEDRON',NT,' at POINT',
     %                 NPt,' TETRAHEDRON QUALITY',QTE,' VOLUME=',VTE
            ENDIF
            PRINT*,'tetra1pt: NOTETR(',NT,')=',(NOTETR(kk,NT),kk=1,4)
         ENDIF

C        LES TETRAEDRES CREES DANS L'ETOILE
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
                  LEFACO( 11, NFLEFA ) = NTE
C                 UN NUMERO DE FACE LEFACO AUX 3 SOMMETS
                  DO N=1,3
                     N1FASC( LEFACO(N,NFLEFA) ) = NFLEFA
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

C        PASSAGE A LA FACE SUIVANTE DE L'ETOILE
         NF2 = NFETOI(5,NF2)
         GOTO 560

      ENDIF
C     QUALITE MOYENNE DES NBTEET1 TETRAEDRES CREES
      QUAMOY = QUAMOY / NBTEET1

      IF( QUAMIN .LE. 0 ) THEN
C        VISUALISATION DES NBTEET1 TETRAEDRES CREES DONT 1 AU MOINS DEGENERE
ccc         TRACTE = .TRUE.
      KTITRE='tetra1pt:                 NOUVEAUX TETRAEDRES DE  QUALITE 
     %MIN          . Volume ' // KNMVOLU
         WRITE(KTITRE(11:17),'(I7)'  ) NBTEET1
         WRITE(KTITRE(63:68),'(F6.3)') QUAMIN
         CALL SANSDBL( KTITRE, L )
         CALL TRACEETOILE( KTITRE(1:L), PTXYZD, NOTETR, NPt,
     %                     NBTEET1, NTETOI )
      ENDIF

C     COMPLETION DES CHAINAGES DES TETRAEDRES CREES DANS L'ETOILE
C     -----------------------------------------------------------
      KTITRE='tetra1pt: NPt=          NBTEET2=        TETRAEDRES avant C
     %OMPLETION des TETRAEDRES OPPOSES. Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)') NPt
      WRITE(KTITRE(33:37),'(I5)') NBTEET2
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET2, NTETOI )

      CALL MJOPTE( NBTEET2, NTETOI, N1TETS, NOTETR, NUDTETR,
     %             N1TEVI,  PTXYZD, NBFANR )

      IF( NBFANR .EQ. 0 ) THEN
C        TOUS LES TETRAEDRES OPPOSES RETROUVES
         GOTO 600
      ELSE
C        IL RESTE DES FACES DE TETRAEDRE OPPOSE NON RETROUVE
         PRINT*,'tetra1pt: PB NPt=',NPt,' avec',NBFANR,
     %          ' FACES de TETRAEDRE OPPOSE NON RETROUVE'
      ENDIF

      NBTE = NBTEET2
 570  IF( NBTE .GT. 0 ) THEN

C        LE HAUT DE LA PILE
         NT = NTETOI( NBTE )
C        LE TETRAEDRE EST DEPILE
         NBTE = NBTE - 1
         IF( NT .LE. 0 ) GOTO 570

c        verification
         if( notetr(1,NT) .EQ. 0 ) THEN
            PRINT*,'tetra1pt: Anomalie NPt=',NPt,' XYZD=',
     %             (PTXYZD(mmm,NPt),mmm=1,4)
            PRINT*,'tetra1pt: Anomalie NOTETR(',NT,')=',
     %             (NOTETR(mmm,NT),mmm=1,8)
         ENDIF

C        QUEL EST LE TETRAEDRE OPPOSE AUX FACES 1 2 3 4
C             ET DE SOMMETS NS1 NS2 NS3 = NSFC
         DO 590 NF=1,4

C           LE TETRAEDRE OPPOSE A LA FACE NF DU TETRAEDRE NT
 585        NTEOP = NOTETR( 4+NF, NT )
            IF( NTEOP .EQ. 0 ) GOTO 590

C           LE NUMERO DES 3 SOMMETS DE LA FACE NF DU TETRAEDRE NT
            NSFC(1) = NOTETR( NOSOFATE(1,NF), NT )
            NSFC(2) = NOTETR( NOSOFATE(2,NF), NT )
            NSFC(3) = NOTETR( NOSOFATE(3,NF), NT )
C           TRI CROISSANT DE SES 3 SOMMETS
            CALL TRI3NO( NSFC, NSFC )

            IF( NTEOP .LT. 0 ) THEN

C              FACE NON TRAITEE. TETRAEDRE OPPOSE INCONNU
C              RECHERCHE DE LA FACE NFSC PARMI LES TETRAEDRES CREES DE L'ETOILE
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
C                    1 FACE COMMUNE A 3 TETRAEDRES???
                     DO L=1,NBTEET2
                        IF( NT2 .EQ. ABS(NTETOI(L)) ) GOTO 586
                     ENDDO
 586                 print*,'tetra1pt: nbt=',nbte+1,' notetr(',nt,')=',
     %                      (notetr(kk,nt),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'tetra1pt: nbt=',I,' notetr(',nt1,')=',
     %                      (notetr(kk,nt1),kk=1,8),' AVEC FACE TRIPLE?'
                     print*,'tetra1pt: nbt=',L,'notetr(',nt2,')=',
     %                      (notetr(kk,nt2),kk=1,8),' AVEC FACE TRIPLE?'
                     print*
      KTITRE='tetra1pt: NPt=          NBTEET2=        TETRAEDRES avant C
     %OMPLETION des TETRAEDRES OPPOSES. Volume ' // KNMVOLU
                     WRITE(KTITRE(15:22),'(I8)') NPt
                     WRITE(KTITRE(33:37),'(I5)') NBTEET2
                     CALL SANSDBL( KTITRE, NBC )
                     CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR,
     %                                 NPt, NBTEET2, NTETOI )
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

                  GOTO 590

 580           ENDDO

            ELSE

C              FACE A PRIORI TRAITEE. VERIFICATION
C              RECHERCHE DE LA FACE NSFC DANS CE TETRAEDRE NTEOP
               CALL NO1F1T( NSFC, NOTETR(1,NTEOP), NFOP )
               IF( NFOP .LE. 0 ) GOTO 590
               NTE = NOTETR( 4+NFOP, NTEOP )
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
      NBTEET = NBTEET1
      DO I = 1,NBTEET1
         NTE = ABS( NTETOI( I ) )
         NTETOI( I ) = NTE
         IF( IVOLTE .NE. 0 ) THEN
            NVOLTE( NTE ) = NOVOLU
         ENDIF
      ENDDO


C     VERIFICATION DES TETRAEDRES OPPOSES
C     -----------------------------------
      CALL VEOPTE( NBTEET, NTETOI, NOTETR, PTXYZD, NBFANR )

      IF( NBFANR .GT. 0 ) THEN
C        IL RESTE DES FACES DE TETRAEDRE OPPOSE NON RETROUVE
         PRINT*,'tetra1pt: Probleme NPt=',NPt,' avec',NBFANR,
     %          ' FACES de TETRAEDRE OPPOSE NON RETROUVE'
      ENDIF

C     CHAINAGE DU SOMMET NPt TETRAEDRISE DANS HEXAPAVE
C     ------------------------------------------------
 600  CALL NUPAVEST( NPt,      PTXYZD,
     %               HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV )

C     TRACE DES NBTEET TETRAEDRES DE SOMMETS NPt
C     ------------------------------------------
      KTITRE='tetra1pt: NPt=          SOMMET des NBTEET=        TETRAEDR
     %ES FINAUX QUAMIN=             . Volume ' // KNMVOLU
      WRITE(KTITRE(15:22),'(I8)'   ) NPt
      WRITE(KTITRE(43:47),'(I5)'   ) NBTEET
      WRITE(KTITRE(76:90),'(G15.6)') QUAMIN
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTEET, NTETOI )
      GOTO 9999


C     SATURATION DU TABLEAU NOTETR DES TETRAEDRES DU MAILLAGE
C     -------------------------------------------------------
 9700 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tetra1pt: SATURATION DES TETRAEDRES NOTETR MXTETR=',
     %           MXTETR
      ELSE
         PRINT*,'tetra1pt: SATURATION of TETRAHEDRA NOTETR MXTETR=',
     %           MXTETR
      ENDIF
      IERR = 3
      GOTO 9999


C     SATURATION DU TABLEAU NTETOI DES TETRAEDRES DE L'ETOILE
C     -------------------------------------------------------
 9800 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tetra1pt: SATURATION DES TETRAEDRES NTETOI MXETOI=',
     %           MXETOI
      ELSE
         PRINT*,'tetra1pt: SATURATION of STAR TETRAHEDRA NTETOI ARRAY MX
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

cccC     BUT: VERIFIER SI TOUT TETRAEDRE VIDE A SON 1-ER NUMERO NOTETR(1,*) NUL
ccc      print *,'tetra1pt: NPt=',NPt,' Verif des tetraedres vides'
ccc      CALL VETEVIDE( N1TEVI, MXTETR, NOTETR, IERR )

      RETURN
      END
