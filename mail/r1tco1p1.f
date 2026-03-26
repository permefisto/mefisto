      SUBROUTINE R1TCO1P1( NPt,     MXSOMM, PTXYZD, INALLT,
     %                     HEXAPAVE,NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                     NOTEIN,  MXTETR, NOTETR, NUDTETR, N1TETS,
     %                     MXTEET,  NOTEET, MXTECH, NOTECH, 
     %                     NOTET1,  COBARY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER UN TETRAEDRE NOTET1 CONTENANT LE POINT NPt
C -----

C ENTREES:
C --------
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DU TABLEAU PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C INALLT : =1 RECHERCHE EVENTUELLE SUR TOUS LES TETRAEDRES ACTUELS
C          =0 PAS DE CETTE RECHERCHE

C HEXAPAVE: MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV  : NOMBRE D'ARETES DANS LA DIRECTION I
C ECHPAV  : ECHELLE DANS LA DIRECTION I
C N1SPAVE : NO DU 1-ER SOMMET DANS PTXYZD DU PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES

C NOTEIN : NUMERO NOTETR DU TETRAEDRE DE DEBUT DE RECHERCHE
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES STOCKABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET
C MXTEET : NOMBRE MAXIMAL DE NO DE TETRAEDRES DE L'ETOILE DU POINT NPt
C MXTECH : NOMBRE MAXIMAL DE NO DE TETRAEDRES DU CHEMIN MENANT A NOTET1

C SORTIES:
C --------
C NOTEET : NO NOTETR DES NBTEET TETRAEDRES DE L'ETOILE DU POINT NPt
C NOTECH : NO NOTETR DES NBTECH TETRAEDRES DU CHEMIN MENANT A NOTET1
C NOTET1 : >0 NUMERO NOTETR DU TETRAEDRE CONTENANT LE POINT NPt
C          =0 SI PAS DE TETRAEDRE CONTENANT LE POINT NPt
C COBARY : LES 4 COORDONNEES BARYCENTRIQUES DE POINT DANS NOTET1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Saint PIERRE du PERRAY           Decembre 2018
C23456...............................................................012
      include"./incl/langue.inc"
      PARAMETER       ( MXSMIN=1024 )
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      INTEGER           NOTETR(8,MXTETR), N1TETS(MXSOMM),
     %                  NOTEET(MXTEET),   NOTECH(MXTECH),
     %                  NOSMIN(MXSMIN),   NODSNO(MXSMIN)
      INTEGER           NBIPAV(3), N1SPAVE(0:*), NOPTSUIV(1:*)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3), DISTMIN(MXSMIN)

      DOUBLE PRECISION  PTXYZD(4,MXSOMM), COBARY(4), VTE, CBSOM, CBSOMIN

      CHARACTER*112     KTITRE
      INTEGER           NOTEAV(-6:0)
      EQUIVALENCE      (NOTEAV(-6),NOTEM6), (NOTEAV(-5),NOTEM5),
     %                 (NOTEAV(-4),NOTEM4), (NOTEAV(-3),NOTEM3),
     %                 (NOTEAV(-2),NOTEM2), (NOTEAV(-1),NOTEM1),
     %                 (NOTEAV( 0),NOTET0)

      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      TRACTE0 = TRACTE
      NESSAI  = 0
      MXCHEM  = MIN( 4096, MXTECH )
      CBSOMIN = 1D100
      NTEMIN  = 0
      NOTEDO  = 0
      NO1TED  = 0
      NOTET1  = NOTEIN

C     ======================================================
C     RECHERCHE D'UN TETRAEDRE NOTET1 CONTENANT LE POINT NPt
C     ======================================================
 10   NESSAI = NESSAI + 1
      GOTO( 100, 200, 300, 499, 499 ), NESSAI


C     LE TETRAEDRE INITIAL NOTET1 CONTIENT IL NPt?
C     --------------------------------------------
 100  IF( NOTET1 .GT. 0 .AND. NOTETR(1,NOTET1) .GT. 0 ) THEN

C        LE TETRAEDRE NOTET1 CONTIENT IL LE POINT NPt?
C        I.E. LA SOMME DES VALEURS ABSOLUES DES 4 COORDONNEES
C        BARYCENTRIQUES EST ELLE INFERIEURE A 1.0001D0?
         CALL PTDSTE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NOTET1),
     %                NONOUI, VTE, COBARY )

C        EVALUATION DE LA DISTANCE NOTET1-NPt
         CBSOM = ABS( COBARY(1) ) + ABS( COBARY(2) ) +
     %           ABS( COBARY(3) ) + ABS( COBARY(4) )

         IF( CBSOM .LT. CBSOMIN ) THEN
            CBSOMIN = CBSOM
            NTEMIN  = NOTET1
         ENDIF

         IF( NONOUI .GT. 0 ) THEN
C           OUI: LE TETRAEDRE NOTET1 CONTIENT LE POINT NPt
            GOTO 9999
         ENDIF

         IF( CBSOM .GT. 8D0 ) THEN
C           LE TETRAEDRE NOTET1 EST TROP LOIN DU POINT NPt
            GOTO 10
         ENDIF

C        ESSAI AVEC LE TETRAEDRE NOTET1
         NBTEET = 1
         NOTEET( NBTEET ) = NOTET1
         GOTO 400

      ELSE

         GOTO 10

      ENDIF


C     LES PLUS PROCHES SOMMETS TETRAEDRISES DU POINT NPt
C     A UNE DISTANCE D'UN PAVE RANGES DANS LES PAVES DE L'HEXAEDRE
C     ONT ILS UN TETRAEDRE CONTENANT NPt?
C     ------------------------------------------------------------
 200  CALL DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                PTXYZD(1,NPt), PTXYZD,
     %                MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )

C     PASSAGE AUX TETRAEDRES DU SOMMET NSDMIN LE PLUS PROCHE DE NPt
C     LES SOMMETS NODSNO SONT TRIES SELON LA DISTANCE A NPt CROISSANTE
 222  NO1TED = NO1TED + 1
      DO 226 N = NO1TED, NBSMIN

C        LE SOMMET LE PLUS PROCHE DE NPt NON LUI MEME
         NSDMIN = NOSMIN( NODSNO(N) )
         NOTET1 = N1TETS( NSDMIN )
         IF( NSDMIN .NE. NPT .AND. NOTET1.GT.0 ) THEN
            IF( NOTETR(1,NOTET1) .GT. 0 ) THEN

C              RECHERCHE DES NBTEET TETRAEDRES DE SOMMET NSDMIN
               CALL TETR1S( NSDMIN, N1TETS, NOTETR,
     %                      NBTEET, MXTEET, NOTEET, IERR )

C              pour debogger
               IF( IERR .EQ. 6 ) THEN
               PRINT*,'************************************************'
               NTEE = N1TETS(NSDMIN)
               PRINT*,'r1tco1p1: ->tetr1s  St',NSDMIN,
     %                ' NON DANS le TETRAEDRE',NTEE,' de SOMMETS',
     %                (NOTETR(J,NTEE),J=1,4)
               PRINT*,'r1tco1p1: AVANT maj N1TETS(',NSDMIN,')=',NTEE
C              MISE A JOUR DU TABLEAU N1TETS pour debog
               CALL MJN1TETS( MXTETR, NOTETR,  MXSOMM,
     %                        N1TETS, NUDTETR, NBSOMM )
               PRINT*,'r1tco1p1: APRES maj N1TETS(',NSDMIN,')=',
     %                 N1TETS(NSDMIN)
               IF( N1TETS(NSDMIN) .EQ. 0 ) THEN
                  PRINT*,'r1tco1p1: A FAIRE:  RETIRER le SOMMET',NSDMIN,
     %            ' du CHAINAGE HEXAPAVE ou POURQUOI NTEE=',NTEE,'???'
               ENDIF
               PRINT*,'************************************************'
               GOTO 226
               ENDIF

               IF( NBTEET .GT. 0 ) THEN
                  NO1TED = N
                  GOTO 400
               ENDIF

            ENDIF

         ENDIF

 226  ENDDO

      GOTO 10


C     PARCOURS EXHAUSTIF DES TETRAEDRES NUDTETR A 1 POUR TROUVER
C     UN TETRAEDRE NOTET1 CONTENANT LE POINT NPt (NON SOMMET)?
C     ----------------------------------------------------------
 300  IF( NTEMIN .GT. 0 .AND. CBSOMIN .LE. 1.001D0 ) GOTO 555

      IF( INALLT .NE. 0 ) THEN

C        RECHERCHE DANS TOUS LES TETRAEDRES ACTUELS
ccc       PRINT*,'r1tco1p1: NPt=',NPt,' NTEMIN=',NTEMIN,' CBSOMIN=',CBSOMIN
         CALL ETOIL1TET( NPt, PTXYZD(1,NPt), PTXYZD, NUDTETR, NOTETR,
     %                   NOTET1, COBARY )

         IF( NOTET1 .GT. 0 ) THEN
C           LE TETRAEDRE NOTET1 CONTIENT NPt
            NBTEET = 1
            NOTEET( 1 ) = NOTET1
            GOTO 500
         ENDIF

      ENDIF

C     PAS DE TETRAEDRE CONTENANT LE POINT NPt
      GOTO 9900


C     RECHERCHE DANS LES NBTEET TETRAEDRES NOTEET D'UN TETRAEDRE
C     NOTET1 CONTENANT LE POINT NPt ET CHEMINEMENT PAR LA FACE
C     QUI VOIT LE POINT NPt
C     ----------------------------------------------------------
 400  DO 490 NT = 1, NBTEET

C        TETRAEDRE D'UN SOMMET LE PLUS PROCHE DU POINT NPt
         NOTEDO = 0
         NBTECH = 0
         NOTET1 = NOTEET( NT )

 410     IF( NOTET1 .LE. 0 ) THEN
            GOTO 490
         ENDIF

C        EXISTE T IL UNE BOUCLE INFINIE DUE A UN CHEMIN PERIODIQUE
C        DE MOINS DE 32 TETRAEDRES?
         IF( MOD( NBTECH, 32 ) .EQ. 0 ) THEN
            DO NOTEDO = NBTECH-1, MAX(1,NBTECH-31), -1
               IF( NOTECH(NOTEDO) .EQ. NOTET1 ) THEN
C                 BOUCLE INFINIE DETECTEE
C                 NOTEDO EST LE NUMERO DANS LE CHEMIN DU TETRAEDRE
C                 ATTEINT 2 FOIS
C                 POURSUITE DE LA RECHERCHE DE NOTET1 CONTENANT NPt
                  GOTO( 10, 490, 490, 499 ), NESSAI
               ENDIF
            ENDDO
            NOTEDO = 0
         ENDIF

C        NOTET1 EST AJOUTE A LA LISTE NOTECH DES TETRAEDRES DU CHEMIN
C        ------------------------------------------------------------
         IF( NBTECH .GE. MXCHEM ) THEN

C           CHEMIN DES TETRAEDRES VERS PTXYZD(1,NPt) TROP LONG
            PRINT*,'r1tco1p1: POINT ',NPt,' XYZ=',
     %             (PTXYZD(L,NPt),L=1,3),
     %             ' CHEMIN TROP LONG de',MXCHEM,' TETRAEDRES'
            K0 = MAX(1,NBTECH-9)
            PRINT 10499,(K,NOTECH(K),K=K0,NBTECH)
10499       FORMAT(5(' NOTECH(',I7,')=',I9))

            tracte = .true.
      KTITRE='r1tco1p1: POINT          CHEMIN de          TETRAEDRES Nb 
     %ESSAIS=    '
            WRITE(KTITRE(17:24),'(I8)') NPt
            WRITE(KTITRE(36:43),'(I8)') NBTECH
            WRITE(KTITRE(66:66),'(I1)') NESSAI
            CALL SANSDBL( KTITRE, NBC )
            CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                        NBTECH, NOTECH )

            GOTO( 10, 490, 490, 499 ), NESSAI
         ENDIF

C        UN TETRAEDRE DE PLUS DANS LE CHEMIN VERS NPt
         NBTECH = NBTECH + 1
         NOTECH( NBTECH ) = NOTET1

C        LE TETRAEDRE NOTET1 CONTIENT IL LE POINT NPt?
C        I.E. LA SOMME DES VALEURS ABSOLUES DES 4 COORDONNEES
C        BARYCENTRIQUES EST ELLE INFERIEURE A 1.0001D0?
C        ----------------------------------------------------
         CALL PTDSTE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NOTET1),
     %                NONOUI, VTE, COBARY )

         IF( NONOUI .GT. 0 ) THEN
C           OUI: LE TETRAEDRE NOTET1 CONTIENT LE POINT NPt
            GOTO 9999
         ENDIF

         IF( NONOUI .LT. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'r1tco1p1: NoPOINT=',NPt,
     %              ' TETRAEDRE ',NOTET1,' DE VOLUME=',VTE,' A PROBLEME'
            ELSE
               PRINT*,'r1tco1p1: NoPOINT=',NPt,
     %          ' TETRAHEDRON ',NOTET1,' of VOLUME=',VTE,' WITH PROBLEM'
            ENDIF
            PRINT*,'r1tco1p1: NOTETR=',(NOTETR(K,NOTET1),K=1,8)
            GOTO 499
         ENDIF

C        EVALUATION DE LA DISTANCE NOTET1-NPt
         CBSOM = ABS( COBARY(1) ) + ABS( COBARY(2) ) +
     %           ABS( COBARY(3) ) + ABS( COBARY(4) )

         IF( CBSOM .LT. CBSOMIN ) THEN
C           LE TETRAEDRE LE PLUS PROCHE DE NPt
            CBSOMIN = CBSOM
            NTEMIN  = NOTET1
         ENDIF

C        NON: LE TETRAEDRE NOTET1 NE CONTIENT PAS LE PTXYZD(1,NPt)
C        RECHERCHE DE LA FACE NOFMAX DU TETRAEDRE NOTET1 QUI REGARDE
C        DIRECTEMENT LE POINT PTXYZD(1,NPt) POUR PROGRESSER VERS LUI
C        -----------------------------------------------------------

ccc         IF( METHOD .EQ. 0 ) THEN
cccC           METHODE 1: RECHERCHE DE LA FACE DE COORDONNEE
cccC                      BARYCENTRIQUE LA PLUS FAIBLE
ccc            NOFMAX = 0
ccc            CBAMAX = 1D100
ccc            DO NOFACE=1,4
ccc               IF( NOTETR( 4+NOFACE, NOTET1 ) .GT. 0 ) THEN
ccc                  D = COBARY( NOFACE )
ccc                  IF( D .LT. CBAMAX ) THEN
ccc                     CBAMAX = D
ccc                     NOFMAX = NOFACE
ccc                  ENDIF
ccc               ENDIF
ccc            ENDDO
ccc            IF( NOFMAX .GT. 0 ) THEN
cccC              PARCOURS AU DELA DE LA FACE NOFMAX
ccc               GOTO 470
ccc            ELSE
ccc               METHOD = 1
ccc               GOTO 10
ccc            ENDIF
ccc         ENDIF


C        METHODE 2: AUTRE CHEMINEMENT PLUS COUTEUX MAIS PLUS SUR
C        CALCUL DE LA FACE DU TETRAEDRE NOTET1 PAR LAQUELLE IL FAUT
C        SORTIR POUR SE RAPPROCHER DU POINT NPt
C        ----------------------------------------------------------
         CALL FAVRPT( PTXYZD(1,NPt), NOTET1, NOTETR, PTXYZD, NOFMAX )

         IF( NOFMAX .LE. 0 ) THEN
C           CHEMIN A PROBLEME VERS PTXYZD(1,NPt)
            GOTO 490
         ENDIF

C        LE TETRAEDRE OPPOSE A LA FACE NOFMAX DE NOTET1
         NTOP = NOTETR( 4+NOFMAX, NOTET1 )
         IF( NTOP .LE. 0 ) THEN
C           LA FACE NOFMAX EST FRONTIERE. PAS DE TETRAEDRE AU DELA
            GOTO 490
         ENDIF

C        LE TETRAEDRE OPPOSE SUIVANT
         NOTET1 = NOTETR( 4+NOFMAX, NOTET1 )
         GOTO 410

 490  ENDDO

C     POURSUITE DE LA RECHERCHE D'UN TETRAEDRE NOTET1 CONTENANT NPt
      GOTO( 10, 222, 499, 499 ), NESSAI


C     ICI PAS DE TETRAEDRE DETECTE CONTENANT NPt
C     TRACE DES TETRAEDRES DU CHEMIN VERS PTXYZD(1,NPt)
C     -------------------------------------------------
 499  IF( NOTEDO .EQ. 0 ) THEN
         N1CHEM = 1
         N2CHEM = NBTECH
      ELSE
         N1CHEM = NOTEDO
         N2CHEM = NBTECH
      ENDIF
 
      NBTE = N2CHEM - N1CHEM + 1
      KTITRE='r1tco1p1: POINT           TRACE du CHEMIN de     
     %  TETRAEDRES SANS EN TROUVER 1 TETRAEDRE LE CONTENANT'
      WRITE(KTITRE(18:25),'(I8)') NPt
      WRITE(KTITRE(46:53),'(I8)') NBTE
      CALL SANSDBL( KTITRE, NBC )
      CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                  NBTE, NOTECH(N1CHEM) )
      TRACTE = TRACTE0

      IF( NESSAI .LE. 2 ) THEN
C        VERS L'ESSAI SUIVANT
         GOTO 10
      ENDIF
      GOTO 9900

C     ==================================================================
C     ICI LE TETRAEDRE NOTET1 CONTIENT LE POINT NPt
C     ==================================================================
C     TRACE DES TETRAEDRES DU CHEMIN VERS PTXYZD(1,NPt)
 500  IF( TRACTE ) THEN
      KTITRE='r1tco1p1: POINT          CHEMIN de          TETRAEDRES Nb 
     %ESSAIS=    '
         WRITE(KTITRE(17:24),'(I8)') NPt
         WRITE(KTITRE(36:43),'(I8)') NBTECH
         WRITE(KTITRE(66:66),'(I1)') NESSAI
         CALL SANSDBL( KTITRE, NBC )
         CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                     NBTECH, NOTECH )
      ENDIF

C     CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT NPt DANS NOTET1
      CALL COBATE( PTXYZD(1,NPt), PTXYZD, NOTETR(1,NOTET1),
     %             VTE, COBARY, IERR )

C     EVALUATION DE LA DISTANCE NOTET1-NPt
      CBSOM = ABS( COBARY(1) ) + ABS( COBARY(2) ) +
     %        ABS( COBARY(3) ) + ABS( COBARY(4) )

      IF( CBSOM .LT. CBSOMIN ) THEN
         CBSOMIN = CBSOM
         NTEMIN  = NOTET1
      ENDIF

C     NTEMIN EST LE TETRAEDRE LE PLUS PROCHE OU CONTIENT NPt
      IF( CBSOMIN .GT. 1.001D0 ) THEN
         print*,'r1tco1p1: BIZARRE sortie avec NPt=',NPT,
     %          ' TETRAEDRE',NTEMIN,
     %          ' COBARY=',COBARY,' CBSOMIN=', CBSOMIN,' >1'
      ENDIF

C     NTEMIN EST CENSE CONTENIR LE POINT NPt
 555  NOTET1 = NTEMIN
      GOTO 9999

C     ==================================================================
C     ABANDON: PAS DE TETRAEDRE CONTENANT LE POINT NPt
C     ==================================================================
 9900 NOTET1 = 0
      
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT*,'r1tco1p1: Probleme POINT', NPt, (PTXYZD(1,NPt),L=1,4),
ccc     %  ' EXTERIEUR AU PLUS PROCHE TETRAEDRE',NTEMIN,' CBSOMIN=',CBSOMIN
ccc      ELSE
ccc         PRINT*,'r1tco1p1: Problem POINT', NPt, (PTXYZD(1,NPt),L=1,4),
ccc     %   ' OUTSIDE the NEAREST TETRAHEDRON',NTEMIN,' CBSOMIN=',CBSOMIN
ccc      ENDIF


 9999 TRACTE = TRACTE0
      RETURN
      END
