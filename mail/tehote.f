      SUBROUTINE TEHOTE( NOFOTI, NUTYSU, RAP2P3,
     %                   NBARPI, MXSOMM, NBSOMM, PXYD,
     %                   ORIGIN, D2D3,   XMIN,   YMIN,  XYMXMI,
     %                   COMXMI, ARETMX,
     %                   LETREE, MXQUEU, LAQUEU,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    HOMOGENEISATION DE L'ARBRE DES TE A UN SAUT DE TAILLE AU PLUS
C -----    PRISE EN COMPTE DES DISTANCES SOUHAITEES AUTOUR DES SOMMETS INITIAUX
C
C ENTREES:
C --------
C NOFOTI : NUMERO DANS LE LX DES FONCTIONS DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C          >0 SI ELLE EXISTE ET 0 SINON
C NUTYSU : NUMERO D'OPTION DE GENERATION DU MAILLAGE DE CETTE SURFACE
C          8 APPEL DIRECT
C          1 APPEL PAR LE SP SUEX01-QUNSAL (QUADRANGLE ALGEBRIQUE)
C          6 APPEL PAR LE SP SUEX06-TRNSAL (TRIANGLE   ALGEBRIQUE)
C          3 APPEL PAR LE SP SUEX03 (SURFACE B-SPLINE D'INTERPOLATION)
C RAP2P3 : RAPPORT DU PERIMETRE DE L'ENVELOPPE DU DOMAINE PLAN SUR
C                  LE PERIMETRE DE L'ENVELOPPE DU DOMAINE DE R**3
C          APPROXIMATION POUR CALCULER RAPIDEMENT LA TAILLE DE L'ARETE
C          DANS LE PLAN CONNAISSANT SA TAILLE DANS R**3
C NBARPI : NOMBRE DE SOMMETS DE LA FRONTIERE + NOMBRE DE POINTS INTERNES
C          IMPOSES PAR L'UTILISATEUR
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION  ET TE
C MXQUEU : NOMBRE D'ENTIERS UTILISABLES DANS LAQUEU
C ORIGIN : POINT ORIGINE DANS LE PLAN DE MAILLAGE
C D2D3   : MATRICE DE PASSAGE DES COORDONNEES 3D EN COORDONNEES 2D
C XMIN   : ABSCISSE MINIMALE APRES TRANSFORMATION 3D->2D
C YMIN   : ORDONNEE MINIMALE APRES TRANSFORMATION 3D->2D
C XYMXMI : ECART MAXIMAL ENTRE XMAX-XMIN ET YMAX-YMIN
C COMXMI : MINIMUM ET MAXIMUM DES COORDONNEES DE L'OBJET
C ARETMX : LONGUEUR MAXIMALE DES ARETES DES TRIANGLES EQUILATERAUX
C PERMTR : PERIMETRE DE LA LIGNE ENVELOPPE DANS LE PLAN
C          AVANT MISE A L'ECHELLE A 2**20
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE DE SOMMETS APRES IDENTIFICATION
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C LETREE : ARBRE-4 DES TRIANGLES EQUILATERAUX (TE) FOND DE LA TRIANGULATION
C          LETREE(0,0) : NO DU 1-ER TE VIDE DANS LETREE
C          LETREE(1,0) : MAXIMUM DU 1-ER INDICE DE LETREE (ICI 8)
C          LETREE(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LETREE (ICI MXTREE)
C          LETREE(0:8,1) : RACINE DE L'ARBRE  (TRIANGLE SANS SUR TRIANGLE)
C          SI LETREE(0,.)>0 ALORS
C             LETREE(0:3,J) : NO (>0) LETREE DES 4 SOUS-TRIANGLES DU TRIANGLE J
C          SINON
C             LETREE(0:3,J) :-NO PXYD DES 1 A 4 POINTS INTERNES AU TRIANGLE J
C                             0  SI PAS DE POINT
C                             ( J EST ALORS UNE FEUILLE DE L'ARBRE )
C          LETREE(4,J) : NO LETREE DU SUR-TRIANGLE DU TRIANGLE J
C          LETREE(5,J) : 0 1 2 3 NO DU SOUS-TRIANGLE J POUR SON SUR-TRIANGLE
C          LETREE(6:8,J) : NO PXYD DES 3 SOMMETS DU TRIANGLE J
C
C AUXILIAIRE :
C ------------
C LAQUEU : MXQUEU ENTIERS SERVANT DE QUEUE POUR LE PARCOURS DE LETREE
C
C SORTIES:
C --------
C IERR   :  0 SI PAS D'ERREUR
C          51 SI SATURATION LETREE DANS TE4STE
C          52 SI SATURATION PXYD   DANS TE4STE
C          >0 SI AUTRE ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC      AVRIL 1997
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/ampli.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION

      DOUBLE PRECISION  PXYD(3,MXSOMM), D2D3(3,3), D, D2, ARETM2
      REAL              ORIGIN(3), XYZ(3), COMXMI(3,2)
      INTEGER           LETREE(0:8,0:*)
C
      INTEGER           LAQUEU(1:MXQUEU),LEQUEU
C                       LEQUEU : ENTREE DANS LA QUEUE
C                       LHQUEU : LONGUEUR DE LA QUEUE
C                       GESTION CIRCULAIRE
C
      INTEGER           NUSTE(3)
      EQUIVALENCE      (NUSTE(1),NS1),(NUSTE(2),NS2),(NUSTE(3),NS3)
      DOUBLE PRECISION  RAP2P3
C
      IERR = 0

      TRATRI = .TRUE.
      CALL EFFACEMEMPX
C
C     LE CARRE DE LA LONGUEUR DE L'ARETE DE TRIANGLES EQUILATERAUX
C     SOUHAITEE POUR LE FOND DE LA TRIANGULATION
      ARETM2 = ( ARETMX * AMPLI ) ** 2
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE' DES ARETES
C     AUTOUR DU POINT.  ICI LA CARTE EST SUPPOSEE ISOTROPE
C     ==========================================================
C     ATTENTION: SI LA FONCTION TAILLE_IDEALE EXISTE
C                ALORS PXYD(3,*) EST LA TAILLE_IDEALE DANS L'ESPACE INITIAL
C                SINON PXYD(3,*) EST LA DISTANCE CALCULEE DANS LE PLAN PAR
C                PROPAGATION A PARTIR DES TAILLES DES ARETES DE LA FRONTIERE
C
      IF( NOFOTI .GT. 0 ) THEN
C
C        LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C        ---------------------------------------
C        INITIALISATION DE LA DISTANCE SOUHAITEE AUTOUR DES POINTS 1 A NBSOMM
         DO I=1,NBSOMM
C           CALCUL DE PXYZD(1:2,I) DANS LE REPERE INITIAL => XYZ(1:3)
            CALL TETAID( NUTYSU, NOFOTI, RAP2P3, PXYD(1,I), PXYD(2,I),
     %                   ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                   XYZ,    PXYD(3,I), IERR )
            IF( IERR .NE. 0 ) GOTO 9990
         ENDDO
C
      ELSE
C
C        LA FONCTION TAILLE_IDEALE(X,Y,Z) N'EXISTE PAS
C        ---------------------------------------------
C        PRISE EN COMPTE DES DISTANCES SOUHAITEES DANS LE PLAN
C        AUTOUR DES POINTS FRONTALIERS ET DES POINTS INTERNES IMPOSES
C        TOUTES LES AUTRES DISTANCES SOUHAITEES ONT ETE MIS A ARETMX
C        LORS DE L'EXECUTION DU SP TEQINI
         DO I=1,NBARPI

C           LE SOMMET I N'EST PAS UN SOMMET DE LETREE => SOMMET FRONTALIER
C           RECHERCHE DU SOUS-TRIANGLE MINIMAL FEUILLE CONTENANT LE POINT I
            NTE = 1
C
 2          NTE = NOTRPT( PXYD(1,I), PXYD, NTE, LETREE )
C           LA DISTANCE AU SOMMET LE PLUS ELOIGNE EST ELLE INFERIEURE
C           A LA DISTANCE SOUHAITEE?
            NS1 = LETREE(6,NTE)
            NS2 = LETREE(7,NTE)
            NS3 = LETREE(8,NTE)
            D2  = MAX( ( PXYD(1,I)-PXYD(1,NS1) )**2 +
     %                 ( PXYD(2,I)-PXYD(2,NS1) )**2
     %               , ( PXYD(1,I)-PXYD(1,NS2) )**2 +
     %                 ( PXYD(2,I)-PXYD(2,NS2) )**2
     %               , ( PXYD(1,I)-PXYD(1,NS3) )**2 +
     %                 ( PXYD(2,I)-PXYD(2,NS3) )**2 )

            IF( D2 .GT. PXYD(3,I)**2 ) THEN
C              LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C              EN 4 SOUS-TRIANGLES
               CALL TE4STE( NBSOMM, MXSOMM, PXYD, NTE, LETREE,
     &                      IERR )
               IF( IERR .NE. 0 ) THEN
                  GOTO 9999
               ELSE
                  GOTO 2
               ENDIF
            ENDIF

         ENDDO
      ENDIF

C     LE SOUS-TRIANGLE CENTRAL DE LA RACINE EST DECOUPE SYSTEMATIQUEMENT
C     ==================================================================
      NTE = 2
      IF( LETREE(0,2) .LE. 0 ) THEN
C        LE SOUS-TRIANGLE CENTRAL DE LA RACINE N'EST PAS SUBDIVISE
C        IL EST DONC DECOUPE EN 4 SOUSTRIANGLES
         NBSOM0 = NBSOMM
         CALL TE4STE( NBSOMM, MXSOMM, PXYD, NTE, LETREE,
     %                IERR )
         IF( IERR .NE. 0 ) THEN
            GOTO 9999
         ENDIF
         DO I=NBSOM0+1,NBSOMM
C           MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS DE TE
            CALL TETAID( NUTYSU, NOFOTI, RAP2P3, PXYD(1,I), PXYD(2,I),
     %                   ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                   XYZ,    PXYD(3,I), IERR )
            IF( IERR .NE. 0 ) GOTO 9990
         ENDDO
      ENDIF

C     TOUT TE CONTENU DANS LE RECTANGLE ENGLOBANT DOIT AVOIR UN
C     COTE < ARETMX ET ETRE DE MEME TAILLE QUE LES TE VOISINS
C     S'IL CONTIENT UN POINT; SINON UN SEUL SAUT DE TAILLE EST PERMIS
C     ===============================================================
C     LE RECTANGLE ENGLOBANT POUR SELECTIONNER LES TE "INTERNES"
C     LE NUMERO DES 3 SOMMETS DU TE ENGLOBANT RACINE DE L'ARBRE DES TE
      NS1 = LETREE(6,1)
      NS2 = LETREE(7,1)
      NS3 = LETREE(8,1)
      A   = ARETMX * 0.01

C     ABSCISSE DU MILIEU DE L'ARETE GAUCHE DU TE 1
      S      = REAL( ( PXYD(1,NS1) + PXYD(1,NS3) ) / 2 )
      XRMIN  = MIN( S, COMXMI(1,1) - ARETMX ) - A

C     ABSCISSE DU MILIEU DE L'ARETE DROITE DU TE 1
      S      = REAL( ( PXYD(1,NS2) + PXYD(1,NS3) ) / 2 )
      XRMAX  = MAX( S, COMXMI(1,2) + ARETMX ) + A

      YRMIN  = COMXMI(2,1) - ARETMX
C     ORDONNEE DE LA DROITE PASSANT PAR LES MILIEUS DES 2 ARETES
C     DROITE GAUCHE DU TE 1
      S      = REAL( ( PXYD(2,NS1) + PXYD(2,NS3) ) / 2 )
      YRMAX  = MAX( S, COMXMI(2,2) + ARETMX ) + A

C     CAS PARTICULIER DE 3 OU 4 OU PEU D'ARETES FRONTALIERES
      IF( NBARPI .LE. 8 ) THEN
C        TOUT LE TRIANGLE ENGLOBANT (RACINE) EST A PRENDRE EN COMPTE
         XRMIN = REAL( PXYD(1,NS1) - A )
         XRMAX = REAL( PXYD(1,NS2) + A )
         YRMIN = REAL( PXYD(2,NS1) - A )
         YRMAX = REAL( PXYD(2,NS3) + A )
      ENDIF
C
      NBS0   = NBSOMM
      NBITER = -1
C
C     INITIALISATION DE LA QUEUE DES TRIANGLES EQUILATERAUX
  5   NBITER = NBITER + 1
      LEQUEU = 1
      LHQUEU = 0
C     LA RACINE DE LETREE INITIALISE LA QUEUE
      LAQUEU(1) = 1
C
C     TANT QUE LA LONGUEUR DE LA QUEUE EST >=0 TRAITER LE DEBUT DE QUEUE
C     ==================================================================
 10   IF( LHQUEU .GE. 0 ) THEN
C
C        LE TRIANGLE TE A TRAITER
         I   = LEQUEU - LHQUEU
         IF( I .LE. 0 ) I = MXQUEU + I
         NTE = LAQUEU( I )
C        LA LONGUEUR DE LA QUEUE EST REDUITE
         LHQUEU = LHQUEU - 1
C
C        NTE EST IL UN SOUS-TRIANGLE FEUILLE MINIMAL ?
 15      IF( LETREE(0,NTE) .GT. 0 ) THEN
C
C           NON LES 4 SOUS-TRIANGLES SONT MIS DANS LA QUEUE
            IF( LHQUEU + 4 .GE. MXQUEU ) THEN
               WRITE(IMPRIM,*) 'TEHOTE: SATURATION DE LA QUEUE'
               IERR = 7
               GOTO 9999
            ENDIF
            DO 20 I=3,0,-1
C              AJOUT DU SOUS-TRIANGLE I
               LHQUEU = LHQUEU + 1
               LEQUEU = LEQUEU + 1
               IF( LEQUEU .GT. MXQUEU ) LEQUEU = LEQUEU - MXQUEU
               LAQUEU( LEQUEU ) = LETREE( I, NTE )
 20         CONTINUE
            GOTO 10
C
         ENDIF
C
C        ICI NTE EST UN TRIANGLE MINIMAL NON SUBDIVISE
C        ---------------------------------------------
C        LE TE EST IL DANS LE CADRE ENGLOBANT DE L'OBJET ?
         NS1 = LETREE(6,NTE)
         NS2 = LETREE(7,NTE)
         NS3 = LETREE(8,NTE)
         IF( PXYD(1,NS1) .GT. PXYD(1,NS2) ) THEN
            DMIN = REAL( PXYD(1,NS2) )
            DMAX = REAL( PXYD(1,NS1) )
         ELSE
            DMIN = REAL( PXYD(1,NS1) )
            DMAX = REAL( PXYD(1,NS2) )
         ENDIF
         IF( (XRMIN .LE. DMIN .AND. DMIN .LE. XRMAX) .OR.
     %       (XRMIN .LE. DMAX .AND. DMAX .LE. XRMAX) ) THEN
            IF( PXYD(2,NS1) .GT. PXYD(2,NS3) ) THEN
               DMIN = REAL( PXYD(2,NS3) )
               DMAX = REAL( PXYD(2,NS1) )
            ELSE
               DMIN = REAL( PXYD(2,NS1) )
               DMAX = REAL( PXYD(2,NS3) )
            ENDIF
            IF( (YRMIN .LE. DMIN .AND. DMIN .LE. YRMAX) .OR.
     %          (YRMIN .LE. DMAX .AND. DMAX .LE. YRMAX) ) THEN

C              NTE EST UN TE FEUILLE ET INTERNE AU RECTANGLE ENGLOBANT
C              =======================================================
C              LE CARRE DE LA LONGUEUR DE L'ARETE DU TE DE NUMERO NTE
               D2 = (PXYD(1,NS1)-PXYD(1,NS2)) ** 2
     %            + (PXYD(2,NS1)-PXYD(2,NS2)) ** 2

               IF( NOFOTI .EQ. 0 ) THEN

C                 IL N'EXISTE PAS DE FONCTION 'TAILLE_IDEALE'
C                 -------------------------------------------
C                 SI LA TAILLE EFFECTIVE DE L'ARETE DU TE EST SUPERIEURE A ARETM
C                 ALORS LE TE EST DECOUPE
                  IF( D2 .GT. ARETM2 ) THEN
C                    LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C                    EN 4 SOUS-TRIANGLES
                     CALL TE4STE( NBSOMM, MXSOMM, PXYD,
     %                            NTE,    LETREE, IERR )
                     IF( IERR .NE. 0 ) THEN
                        GOTO 9999
                     ELSE
                        GOTO 15
                     ENDIF
                  ENDIF

               ELSE

C                 IL EXISTE ICI UNE FONCTION 'TAILLE_IDEALE'
C                 ------------------------------------------
C                 SI LA TAILLE EFFECTIVE DE L'ARETE DU TE EST SUPERIEURE AU MINI
C                 DES 3 TAILLES_IDEALES AUX SOMMETS  ALORS LE TE EST DECOUPE
                  DO 28 I=1,3
                     D = ( PXYD(3,NUSTE(I)) * AMPLI ) ** 2
                     IF( D2 .GT. D ) THEN
C                       LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C                       EN 4 SOUS-TRIANGLES
                        NBSOM0 = NBSOMM
                        CALL TE4STE( NBSOMM, MXSOMM, PXYD,
     &                               NTE,    LETREE, IERR )
                        IF( IERR .NE. 0 ) THEN
                           GOTO 9999
                        ENDIF
                        DO 27 J=NBSOM0+1,NBSOMM
C                          MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS DE
                           CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                                  PXYD(1,J), PXYD(2,J),
     %                                  ORIGIN, D2D3, XMIN,YMIN, XYMXMI,
     %                                  XYZ,    PXYD(3,J), IERR )
                           IF( IERR .NE. 0 ) GOTO 9990
 27                     CONTINUE
                        GOTO 15
                     ENDIF
 28               CONTINUE
               ENDIF

C              RECHERCHE DU NOMBRE DE NIVEAUX ENTRE NTE ET LES TE VOISINS PAR SE
C              SI LA DIFFERENCE DE SUBDIVISIONS EXCEDE 1 ALORS LE PLUS GRAND DES
C              =================================================================
               DO 30 I=1,3

C                 NOTEVA TRIANGLE VOISIN DE NTE PAR L'ARETE I
                  CALL N1TRVA( NTE, I, LETREE, NOTEVA, NIVEAU )
                  IF( NOTEVA .LE. 0 ) GOTO 30
C                 IL EXISTE UN TE VOISIN
                  IF( NIVEAU .GT. 0 ) GOTO 30
C                 NTE A UN TE VOISIN PLUS PETIT OU EGAL
                  IF( LETREE(0,NOTEVA) .LE. 0 ) GOTO 30
C                 NTE A UN TE VOISIN NOTEVA SUBDIVISE AU MOINS UNE FOIS

                  IF( NBITER .GT. 0 ) THEN
C                    LES 2 SOUS TRIANGLES VOISINS SONT-ILS SUBDIVISES?
                     NS2 = LETREE(I,NOTEVA)
                     IF( LETREE(0,NS2) .LE. 0 ) THEN
C                       NS2 N'EST PAS SUBDIVISE
                        NS2 = LETREE(NOSUI3(I),NOTEVA)
                        IF( LETREE(0,NS2) .LE. 0 ) THEN
C                          LES 2 SOUS-TRIANGLES NE SONT PAS SUBDIVISES
                           GOTO 30
                        ENDIF
                     ENDIF
                  ENDIF

C                 SAUT>1 => LE TRIANGLE NTE DOIT ETRE SUBDIVISE EN 4 SOUS-TRIANG
C                 --------------------------------------------------------------
                  NBSOM0 = NBSOMM
                  CALL TE4STE( NBSOMM, MXSOMM, PXYD, NTE, LETREE,
     &                         IERR )
                  IF( IERR .NE. 0 ) THEN
                     GOTO 9999
                  ENDIF
                  IF( NOFOTI .GT. 0 ) THEN
                     DO 32 J=NBSOM0+1,NBSOMM
C                       MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS DE TE
                        CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                               PXYD(1,J), PXYD(2,J),
     %                               ORIGIN, D2D3, XMIN,YMIN, XYMXMI,
     %                               XYZ,    PXYD(3,J), IERR )
                        IF( IERR .NE. 0 ) GOTO 9990
 32                  CONTINUE
                  ENDIF
                  GOTO 15

 30            CONTINUE
            ENDIF
         ENDIF
         GOTO 10
      ENDIF
      IF( NBS0 .LT. NBSOMM ) THEN
         NBS0 = NBSOMM
         GOTO 5
      ENDIF

      GOTO 9999

C     PB DANS LE CALCUL DE LA FONCTION TAILLE_IDEALE
 9990 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'PROBLEME dans le CALCUL de TAILLE_IDEALE(x,y,z)'
      ELSE
         KERR(1) = 'COMPUTATION of EDGE_LENGTH(x,y,z) with PROBLEM'
      ENDIF
      CALL LEREUR

 9999 TRATRI = .FALSE.
      RETURN
      END
