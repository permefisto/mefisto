      SUBROUTINE ETOIL1PT( NPt,    PTXYZD, NPSOFR,
     %                     NOTET1, NOTETR, N1TETS, COBARY,
     %                     NBCBA0, MXTEET, NBTEET, NOTEET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TETRAEDRE NOTET1 CONTENANT LE POINT NPt et
C -----    SELON SA POSITION INTERNE, CONSTRUCTION DES TETRAEDRES
C          PROCHES DE BOULE CIRCONSCRITE CONTENANT LE POINT

C ENTREES:
C --------
C NPt    : NUMERO DU POINT DANS PTXYZD D'ETOILE A CONSTRUIRE
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
C NOTET1 : NUMERO DU TETRAEDRE DE DEBUT DE RECHERCHE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET
C COBARY : LES 4 COORDONNEES BARYCENTRIQUES DU POINT NPt DANS
C          LE TETRAEDRE NOTEET(1) SI NBTEET>=1
C MXTEET : NOMBRE MAXIMAL DE NO DE TETRAEDRES STOCKABLES DANS NOTEET

C SORTIES:
C --------
C NBCBA0 : NOMBRE DE COORDONNEES BARYCENTRIQUES DE NPt CONSIDEREES NULLES
C NBTEET : >0  NOMBRE DE TETRAEDRES CONTENANT LE POINT NPt
C          =-1 UN TETRAEDRE DE VOLUME<0 ou
C              PAS DE NORMALE POUR RETROUVER UN POINT ou
C              POINT = SOMMET EXISTANT
C NOTEET : NUMERO DANS NOTETR DES NBTEET TETRAEDRES DE L'ETOILE DE NPt
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        Juin 1992
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2018
C23456...............................................................012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE
      INTEGER           NPSOFR(*),NOTETR(8,*), N1TETS(*),
     %                  NOTEET(MXTEET), NOCOBA(4)
      DOUBLE PRECISION  PTXYZD(1:4,1:*), COBARY(4), CBAMAX, C
ccc      CHARACTER*100     KTITRE

      TRACTE0 = TRACTE

C     LE TETRAEDRE NOTET1 CONTIENT LE POINT NPt ET
C     SES COORDONNEES BARYCENTRIQUES SONT COBARY DANS NOTET1
C     ======================================================

C     MISE A JOUR DE N1TETS POUR FAVORISER L'ACCES AU TETRAEDRE NOTET1
      DO I = 1, 4
         N1TETS( NOTETR( I, NOTET1 ) ) = NOTET1
      ENDDO

C     RECHERCHE DU NOMBRE DE COORDONNEES BARYCENTRIQUES NULLES
C     POUR CONNAITRE LA POSITION DU POINT DANS LE TETRAEDRE
      ICMAX  = 0
      CBAMAX = -1D100
      NBCBA0 = 0
      DO I=1,4

         C = ABS( COBARY(I) )
         IF( C .GT. CBAMAX ) THEN
            CBAMAX = C
            ICMAX  = I
         ENDIF

ccc         IF( C .LE. 0.01D0  ) THEN
ccc         IF( C .LE. 0.015D0 ) THEN
ccc         IF( C .LE. 0.06D0  ) THEN
ccc         IF( C .LE. 0.03D0  ) THEN

         IF( C .LE. 0.01D0 ) THEN
            NBCBA0 = NBCBA0 + 1
C           LE NUMERO DE LA COORDONNEE BARYCENTRIQUE
            NOCOBA( NBCBA0 ) = I
         ENDIF

      ENDDO
C
C     TRAITEMENT SELON LE NOMBRE DE COORDONNEES BARYCENTRIQUES NULLES
C     C-A-D LA POSITION DU POINT DANS LE TETRAEDRE
      IF( NBCBA0 .EQ. 0 ) THEN
C
C        POINT NPt FRANCHEMENT INTERNE AU TETRAEDRE NOTET1
C        -------------------------------------------------
         NBTEET = 1
         NOTEET(1) = NOTET1
C
      ELSE IF( NBCBA0 .EQ. 1 ) THEN
C
C        POINT NPt SUR UNE FACE: LE TETRAEDRE OPPOSE A CETTE FACE
C                                S'IL EXISTE EST AJOUTE A L'ETOILE
C        ---------------------------------------------------------
         NOTEET(2) = NOTETR( 4+NOCOBA(1), NOTET1 )
         IF( NOTEET(2) .GT. 0 ) THEN
            NBTEET = 2
         ELSE
            NBTEET = 1
         ENDIF
         NOTEET(1) = NOTET1
C
      ELSE IF( NBCBA0 .EQ. 2 ) THEN
C
C        POINT NPt SUR UNE ARETE: LES TETRAEDRES QUI PARTAGENT CETTE
C                                 ARETE SONT AJOUTES A L'ETOILE
C        NOCOBA(3:4) LES NUMEROS DES 2 SOMMETS DE L'ARETE
C        -----------------------------------------------------------
         IF( NOCOBA(1) .EQ. 1 ) THEN
            IF( NOCOBA(2) .EQ. 2 ) THEN
               NOCOBA(3) = 2
               NOCOBA(4) = 3
            ELSE IF( NOCOBA(2) .EQ. 3 ) THEN
               NOCOBA(3) = 1
               NOCOBA(4) = 3
            ELSE IF( NOCOBA(2) .EQ. 4 ) THEN
               NOCOBA(3) = 1
               NOCOBA(4) = 2
            ENDIF
         ELSE IF( NOCOBA(1) .EQ. 2 ) THEN
            IF( NOCOBA(2) .EQ. 3 ) THEN
               NOCOBA(3) = 3
               NOCOBA(4) = 4
            ELSE IF( NOCOBA(2) .EQ. 4 ) THEN
               NOCOBA(3) = 2
               NOCOBA(4) = 4
            ENDIF
         ELSE
            NOCOBA(3) = 1
            NOCOBA(4) = 4
         ENDIF

C        LISTAGE DES TETRAEDRES D'ARETE NS1-NS2
         NS1 = NOTETR( NOCOBA(3), NOTET1 )
         NS2 = NOTETR( NOCOBA(4), NOTET1 )
         CALL TETR1A(NS1, NS2, N1TETS, NOTETR, NBTEET,MXTEET,NOTEET,IER)
         IF( IER .NE. 0 ) THEN
            PRINT*,'etoil1pt: NPt=',NPt,' IER=',IER,' ARETE',NS1,NS2,
     %             ' DONNE',NBTEET,' TETRAEDRES'
         ENDIF
         IF( NBTEET .LE. 0 ) RETURN

      ELSE
C
C        POINT NPt PROCHE D'UN SOMMET : EST IL A ABANDONNER?
C        -----------------------------  --------------------
         IF( CBAMAX .GE. 0.97D0 ) THEN

C           NPt POINT TROP PROCHE DU SOMMET ICMAX DE NOTET1?
            IF( ICMAX .EQ. 1 ) THEN
               NSBMX = 3
            ELSE
               NSBMX = ICMAX-1
            ENDIF
            NSPR = NOTETR( NSBMX, NOTET1 )

C           DISTANCE NPt->NSPR
            C = SQRT( ( PTXYZD(1,NSPR) - PTXYZD(1,NPt) ) **2
     %              + ( PTXYZD(2,NSPR) - PTXYZD(2,NPt) ) **2
     %              + ( PTXYZD(3,NSPR) - PTXYZD(3,NPt) ) **2 )

ccc            IF( C .LT. PTXYZD(4,NSPR)/2 ) THEN

            IF( NPSOFR(NPt).EQ.-4 .OR. NPSOFR(NPt).EQ.0  .AND.
     %          C .LT. PTXYZD(4,NSPR)/5 ) THEN

C              OUI: NPt SUPPRIMABLE est TROP PROCHE du POINT NSPR
C                   => IL EST REJETE

ccc               PRINT*,
ccc               IF( LANGAG .EQ. 0 ) THEN
ccc                  PRINT*,'etoil1pt: NPt=', NPt,
ccc     %                   ' TROP PROCHE du SOMMET',NSBMX,
ccc     %                   ' du TETRAEDRE',NOTET1,':',
ccc     %                   (NOTETR(kk,NOTET1),kk=1,4)
ccc               ELSE
ccc                  PRINT*,'etoil1pt: Point', NPt,
ccc     %                   ' TOO NEAR the VERTEX',NSBMX,
ccc     %                   ' of TETRAHEDRON',NOTET1,':',
ccc     %                   (NOTETR(kk,NOTET1),kk=1,8)
ccc               ENDIF
ccc               PRINT*,'etoil1pt:',NPt,' : ',(PTXYZD(K,NPt),K=1,4)
ccc               DO I=1,4
ccc                  NS = NOTETR( I, NOTET1 )
ccc                  PRINT*,'etoil1pt:',NS,' : ',(PTXYZD(K,NS),K=1,4)
ccc               ENDDO
ccc               PRINT*,'etoil1pt: COBARY=',COBARY,' CBSom=',
ccc     %       ABS(COBARY(1))+ABS(COBARY(2))+ABS(COBARY(3))+ABS(COBARY(4))
ccc     %                ,' d(NPt,NSt)=',C,' TI(NSt)=',PTXYZD(4,NSPR)

               GOTO 9900

            ENDIF

         ENDIF

C        POINT NON TROP PROCHE D'UN SOMMET
         NBTEET = 1
         NOTEET(1) = NOTET1

      ENDIF
      GOTO 9999


C     ABANDON SUITE A UNE ERREUR OU TETRAEDRISATION NON CONVEXE
C     ---------------------------------------------------------
 9900 NBTEET = -1

 9999 TRACTE = TRACTE0
      RETURN
      END
