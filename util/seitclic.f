      SUBROUTINE SEITCLIC( NXCL, NYCL, NoTYIT, MNITEMS,
     %                     NUITDMIN, NPITEMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RECHERCHER QUEL EST L'ITEM VISIBLE LE PLUS PROCHE
C -----   DU POINT CLIQUE (NXCL,NYCL) SUR LA FENETRE AXONOMETRIQUE PIXELS

C ENTREES:
C --------
C NXCL   : ABSCISSE PIXELS DU POINT CLIQUE
C NYCL   : ORDONNEE PIXELS DU POINT CLIQUE
C NoTYIT : No DU TYPE D'ITEMS
C          =6 SOMMETS
C          =7 ARETES
C          =8 FACES
C MNITEMS: ADRESSE MCN DU TABLEAU DES ITEMS VISIBLES
C          (SOMMETS ou ARETES ou FACES)

C SORTIES:
C --------
C NUITDMIN: >0 NUMERO DE L'ITEM LE PLUS PROCHE DU POINT CLIQUE (NXCL,NYCL)
C           =0 PAS D'ITEM PROCHE
C NPITEMIN: >0 NUMERO DU SOMMET ou ARETE ou FACE DE L'ITEM LE PLUS PROCHE
C           =0 SI PAS D'ITEM PROCHE (NPITEMIN=MCN(MNITEMS+MOITEM*NUITDMIN+3))
C              SI ITEMS SOMMETS le NO EST L'INDICE DANS LE TABLEAU XYZSOM
C              SI ITEMS ARETES  le NO EST L'INDICE DANS LE TABLEAU LARETE
C              SI ITEMS FACES   le NO EST L'INDICE DANS LE TABLEAU NOSOEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C MODIFS : Alain PERRONNET  Saint PIERRE du PERRAY          Janvier 2021
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/darete.inc"
      include"./incl/pp.inc"
      REAL          RMCN(MOTMCN)
      COMMON         MCN(MOTMCN)
      EQUIVALENCE  (RMCN(1),MCN(1))
      REAL          XYZIT(3), XYZAIT(3)

C     DISTANCE PTV-OEIL POUR CETTE VISEE
C     ----------------------------------
      DPTVOEIL = SQRT( ( AXOEIL(1) - AXOPTV(1) ) **2
     %               + ( AXOEIL(2) - AXOPTV(2) ) **2
     %               + ( AXOEIL(3) - AXOPTV(3) ) **2 )

C     CARACTERISTIQUES DU TABLEAU DES ITEMS VISIBLES
C     ----------------------------------------------
C     NOMBRE DE MOTS PAR ITEM (=4 ou PLUS)
      MOITEM = MCN( MNITEMS )

C     NOMBRE ACTUEL D'ITEMS
      NBITEM = MCN( MNITEMS+2 )

C     RECHERCHE DE L'ITEM VISIBLE LE PLUS PROCHE DU CLIC NXCL NYCL
      NUITDMIN = 0
      NPITEMIN = 0
      DISTPXMIN = 1E28
      DISTMOEIL = 1E28

C     PARCOURS DES ITEMS VISIBLES
C     ---------------------------
C     ADRESSE DE DEBUT DES ITEMS PRESENTS SUR L'ECRAN
      MNIT = MNITEMS
      DO 50 NIT = 1, NBITEM

C        L'ADRESSE MCN DU DEBUT DE L'ITEM NIT
         MNIT = MNIT + MOITEM

C        ABSCISSE,ORDONNEE,COTE XYZ DU BARYCENTRE DE L'ITEM NIT
         XYZIT(1) = RMCN( MNIT   )
         XYZIT(2) = RMCN( MNIT+1 )
         XYZIT(3) = RMCN( MNIT+2 )

C        LE NUMERO dans le type d'ITEM: SOMMET ( VALEUR DANS XYZSOM )
C                                       ARETE  ( VALEUR DANS LARETE )
C                                       FACE   ( VALEUR DANS NOSOEF )
         NPITEM = MCN( MNIT+3 )

C        LES COORDONNEES AXONOMETRIQUES ET PIXELS de l'ITEM NIT
         IF( NDIMLI .EQ. 2 ) THEN
C           ABSCISSE ORDONNEE PIXELS DE L'ITEM
            NXIT = NUPXEX( XYZIT(1) )
            NYIT = NUPXEY( XYZIT(2) )
         ELSE
C           COORDONNEES AXONOMETRIQUES de l'ITEM NIT
            CALL XYZAXO( XYZIT, XYZAIT )
C           ABSCISSE ORDONNEE PIXELS DE L'ITEM
            NXIT = NUPXEX( XYZAIT(1) )
            NYIT = NUPXEY( XYZAIT(2) )
         ENDIF

C        DISTANCE**2 PIXELS ENTRE (NXCL,NYCL) et (NXIT,NYIT) de l'ITEM NIT
         DISTITPX = (NXCL-NXIT) ** 2 + (NYCL-NYIT) ** 2

         IF( DISTITPX .GT. 1600 ) THEN
C           CLIC NON DANS LE CERCLE DE RAYON 40 PIXELS
C           POINT CLIQUE TROP ELOIGNE DE L'ITEM NIT DANS LE PLAN des PIXELS
            GOTO 50
         ENDIF

C        L'ITEM NIT EST DANS LE CERCLE DE RAYON 200 PIXELS DE CENTRE NXCL NYCL
         IF( DISTITPX .LT. DISTPXMIN ) THEN

C           L'ITEM NIT EST LE PLUS PROCHE EN PIXELS DU POINT CLIQUE
            IF( NDIMLI .EQ. 3 ) THEN

C              DISTANCE A L'OEIL DU BARYCENTRE DE L'ITEM NIT
               DISTAROE = DPTVOEIL - XYZAIT(3)

C              L'ITEM NIT EST IL TRES ELOIGNE DE L'OEIL?
C              ET SURTOUT LOIN DERRIERE LE PRECEDENT POINT IT PROCHE DU CLIC?
               IF( NUITDMIN .GT. 0 ) THEN
C                 IL EXISTE UN ITEM NUITDMIN PRECEDENT POINT PROCHE DU CLIC
                  IF( DISTAROE .GE. DISTMOEIL ) THEN
C                    L'ITEM NIT EST TROP LOIN DE L'OEIL
                     GOTO 50
                  ENDIF
               ENDIF

            ENDIF

C           NIT EST LE NOUVEL ITEM LE PLUS PROCHE DU POINT CLIQUE (NXCL,NYCL)
            NUITDMIN = NIT
            NPITEMIN = NPITEM
            DISTPXMIN = DISTITPX
            DISTMOEIL = DISTAROE

         ENDIF

 50   ENDDO

      IF( NUITDMIN .LT. 0 ) THEN
C        JAMAIS TRACE
         MNITMI = MNITEMS + MOITEM * NUITDMIN
         PRINT*,'sestclic: l''ITEM PROCHE du CLIC', NXCL, NYCL,
     %          ' de TYPE',NoTYIT,' est',NUITDMIN,' No SAF=',NPITEMIN,
     %          ' XYZ=',(RMCN(MNITMI+K),K=0,2)
C        TRACE '*' AU BARYCENTRE DE L'ITEM
         CALL SYMBOLE3D( NoTYIT-5, RMCN(MNITMI), '*' )
      ENDIF

      RETURN
      END
