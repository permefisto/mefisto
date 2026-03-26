      SUBROUTINE DVHOTE( NUTYSU, NBSTEC, NOSTEC, NOSTE1,
     %                   MXSOMM, NBSOMM, PXYD,
     %                   ORIGIN, D2D3,   XMIN,   YMIN,  XYMXMI,
     %                   COMXMI, ARETMX,
     %                   LETREE, NOLESO, MXQUEU, LAQUEU,
     %                   NBS   , NOPXYD, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    HOMOGENEISATION DE L'ARBRE DES TEE
C -----    CONSTRUCTION DU TABLEAU DES NUMEROS DES SOMMETS A TRIANGULER
C          NOLESO DES SOMMETS DES TE TROP PROCHES CHANGE DE SIGNE
C
C ENTREES:
C --------
C NUTYSU : NUMERO D'OPTION DE GENERATION DU MAILLAGE DE CETTE SURFACE
C          9 APPEL DIRECT
C          1 APPEL PAR LE SP SUEX01-QUNSAL (QUADRANGLE ALGEBRIQUE)
C          6 APPEL PAR LE SP SUEX06-TRNSAL (TRIANGLE   ALGEBRIQUE)
C NBSTEC : NOMBRE DE SOMMETS DE L'ENVELOPPE CONVEXE
C NOSTEC : NOSTEC(1,NA) NUMERO DU SOMMET DE L'ARETE NA DE L'ENVELOPPE
C          NOSTEC(2,NA) NUMERO DE L'ARETE SUIVANTE DANS L'ENVELOPPE
C          NA=1 EST LA PREMIERE ARETE DE L'ENVELOPPE PARCOURUE SELON
C          LE SENS DES AIGUILLES D'UNE MONTRE
C          LE CHAINAGE EST CIRCULAIRE
C NOSTE1 : NUMERO PXYD DU PREMIER SOMMET DU TE ENGLOBANT (RACINE)
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION  ET TE
C MXQUEU : NOMBRE D'ENTIERS UTILISABLES DANS LAQUEU
C ORIGIN : POINT ORIGINE DANS LE PLAN DE MAILLAGE
C D2D3   : MATRICE DE PASSAGE DES COORDONNEES 3D EN COORDONNEES 2D
C XMIN   : ABSCISSE MINIMALE APRES TRANSFORMATION 3D->2D
C YMIN   : ORDONNEE MINIMALE APRES TRANSFORMATION 3D->2D
C XYMXMI : ECART MAXIMAL ENTRE XMAX-XMIN ET YMAX-YMIN
C COMXMI : MINIMUM ET MAXIMUM DES COORDONNEES DE L'OBJET
C ARETMX : LONGUEUR MAXIMALE DES ARETES DES TRIANGLES EQUILATERAUX
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE DE SOMMETS APRES IDENTIFICATION
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C LETREE : ARBRE-4 DES TRIANGLES EQUILATERAUX (TE) FOND DE LA TRIANGULATION
C          LETREE(0,0) : NO DU 1-ER TE VIDE DANS LETREE
C          LETREE(0,1) : MAXIMUM DU 1-ER INDICE DE LETREE (ICI 8)
C          LETREE(0,2) : MAXIMUM DECLARE DU 2-EME INDICE DE LETREE (ICI MXTREE)
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
C NOLESO : TABLEAU DU NUMERO DANS LETREE D'UN SOMMET DE TE
C          >0 +NO LETREE DU TRIANGLE DU SOMMET ET SOMMET DE LA TRIANGULATION
C          <0 -NO LETREE DU TRIANGLE DU SOMMET ELIMINE   DE LA TRIANGULATION
C          =0 SI SOMMET N'APPARTENANT PAS AUX SOMMETS DE LETREE
C
C AUXILIAIRE :
C ------------
C LAQUEU : MXQUEU ENTIERS SERVANT DE QUEUE POUR LE PARCOURS DE LETREE
C
C SORTIES:
C --------
C NBS    : NOMBRE DE SOMMETS A TRIANGULER
C NOPXYD : NUMERO DANS PXYD DES NBS SOMMETS A TRIANGULER
C IERR   : 0 SI PAS D'ERREUR , >0 SI SATURATION DE TABLEAU OU AUTRE PB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JUILLET 1994
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
CCC      PARAMETER        (PROXPI=0.9)  ESSAYE MAIS MOINS BON QUE 0.8
      PARAMETER        (PROXPI=0.8)
C
C     SI SUEX03 A APPELE SUEX99 (NUTYSU=3) ALORS LES VARIABLES DU
C     COMMON QUI SUIT DOIVENT ETRE INITIALISEES
      COMMON / S09S03 / LDEXSB, LRX, MNRX,
     %                  LDEYSB, LRY, MNRY, MNSPQB
C
      include"./incl/pp.inc"
      COMMON            RMCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      DOUBLE PRECISION  PXYD(3,MXSOMM), D2D3(3,3)
      REAL              ORIGIN(3), XYZ(3), XY(2), XY0(2),
     %                  XYZSTE(3,3), HST(3), COMXMI(3,2)
      INTEGER           LETREE(0:8,0:*),
     %                  NOLESO(MXSOMM),
     %                  NOPXYD(0:*),
     %                  NOSTEC(1:2,*)
C
      INTEGER           LAQUEU(1:MXQUEU),LEQUEU
C                       LEQUEU : ENTREE DANS LA QUEUE
C                       LHQUEU : LONGUEUR DE LA QUEUE
C                       GESTION CIRCULAIRE
C
      INTEGER           NUSTE(3)
      EQUIVALENCE      (NUSTE(1),NS),(NUSTE(2),NSTE),(NUSTE(3),NS3)
      DOUBLE PRECISION  DPARAF(3), DTAILL
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE' DES ARETES
C     AUTOUR DU POINT.  ICI LA CARTE EST SUPPOSEE ISOTROPE
C     ==========================================================
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) THEN
         XYZ(1) = ORIGIN(1)
         XYZ(2) = ORIGIN(2)
         XYZ(3) = ORIGIN(3)
C        PASSAGE AUX COORDONNEES (X,Y) DE L'ORIGIN
C        DANS LE PLAN APRES MISE A L"ECHELLE
         CALL DV3D2D( XYZ, ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                XY0 )
      ENDIF
C
C     PRISE EN COMPTE DES DISTANCES SOUHAITEES AUTOUR DES POINTS FRONTALIERS
C     ======================================================================
      DO 4 I=1,NOSTE1-1
C
         IF( NOLESO(I) .EQ. 0 ) THEN
C
C           LE SOMMET I N'EST PAS UN SOMMET DE LETREE => SOMMET FRONTALIER
C           RECHERCHE DU SOUS-TRIANGLE MINIMAL CONTENANT LE POINT I
            NTE = 1
 2          NTE  = NOTRPT( PXYD(1,I), PXYD, NTE, LETREE )
C           LA LONGUEUR DE SON COTE EST ELLE INFERIEURE A DISTANCE SOUHAITEE?
            NS   = LETREE(6,NTE)
            NSTE = LETREE(7,NTE)
            D2   = REAL( ( PXYD(1,NS)-PXYD(1,NSTE) )**2 +
     %                   ( PXYD(2,NS)-PXYD(2,NSTE) )**2 )
C
            IF( D2 .GT. 6.75*(PXYD(3,I)**2) ) THEN
C
C              TOUT POINT INTERNE AU TE D'ARETE DE LONGUEUR L EST AU PIRE
C              A LA DISTANCE L/RACINE(3) DU PLUS PROCHE DE SES 3 SOMMETS
C              ICI NOUS VOULONS QUE TOUT POINT INTERNE AU TE SOIT AU PLUS
C              A 1.5 FOIS LA DISTANCE SOUHAITEE AUTOUR DE CE POINT =>
C              L/RACINE(3) < 1.5 L => L**2 < 27/4 DS**2
C
C              LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C              EN 4 SOUS-TRIANGLES
               CALL DV4STE( NBSOMM, MXSOMM, PXYD, NTE, LETREE, NOLESO,
     &                      IERR )
               IF( IERR .NE. 0 ) RETURN
               GOTO 2
C
            ENDIF
         ENDIF
 4    CONTINUE
C
C     LONGUEUR DE L'ARETE DE TRIANGLES EQUILATERAUX SOUHAITEE
C     POUR LE FOND DE LA TRIANGULATION
      ARETE2 = ARETMX * ARETMX
C
C     TOUT TE CONTENU DANS LE RECTANGLE ENGLOBANT DOIT AVOIR UN
C     COTE < ARETMX ET ETRE DE MEME TAILLE QUE LES TE VOISINS
C     S'IL CONTIENT UN POINT; SINON UN SEUL SAUT DE TAILLE EST PERMIS
C     MARQUAGE DES SOMMETS DE TE (NOLESO(NS)<0) TROP PROCHES D'UN POINT
C     =================================================================
C     LE RECTANGLE ENGLOBANT
      XRMIN = COMXMI(1,1) - ARETMX
      XRMAX = COMXMI(1,2) + ARETMX
      YRMIN = COMXMI(2,1) - ARETMX
      YRMAX = COMXMI(2,2) + ARETMX
C
      NBS0   = NBSOMM
      NBITER = -1
C
C     INITIALISATION DE LA QUEUE
  5   NBITER = NBITER + 1
      LEQUEU = 1
      LHQUEU = 0
C     LA RACINE DE LETREE INITIALISE LA QUEUE
      LAQUEU(1) = 1
C
C     TANT QUE LA LONGUEUR DE LA QUEUE EST >=0 TRAITER LE DEBUT DE QUEUE
 10   IF( LHQUEU .GE. 0 ) THEN
C
C        LE TRIANGLE TE A TRAITER
         I   = LEQUEU - LHQUEU
         IF( I .LE. 0 ) I = MXQUEU + I
         NTE = LAQUEU( I )
C        LA LONGUEUR EST REDUITE
         LHQUEU = LHQUEU - 1
C
C        NTE EST IL UN SOUS-TRIANGLE MINIMAL ?
 15      IF( LETREE(0,NTE) .GT. 0 ) THEN
C           NON LES 4 SOUS-TRIANGLES SONT MIS DANS LA QUEUE
            IF( LHQUEU + 4 .GE. MXQUEU ) THEN
               WRITE(IMPRIM,*) 'DVHOTE: SATURATION DE LA QUEUE'
               IERR = 1
               RETURN
            ENDIF
            DO 20 I=3,0,-1
C              AJOUT DU SOUS-TRIANGLE I
               LHQUEU = LHQUEU + 1
               LEQUEU = LEQUEU + 1
               IF( LEQUEU .GT. MXQUEU ) LEQUEU = LEQUEU - MXQUEU
               LAQUEU( LEQUEU ) = LETREE( I, NTE )
 20         CONTINUE
            GOTO 10
         ENDIF
C
C        ICI NTE EST UN TRIANGLE MINIMAL NON SUBDIVISE
C        ---------------------------------------------
C        LE SOUS-TRIANGLE CENTRAL DE LA RACINE EST DECOUPE SYSTEMATIQUEMENT
         IF( NTE .EQ. 2 ) THEN
            CALL DV4STE( NBSOMM, MXSOMM, PXYD, NTE, LETREE, NOLESO,
     &                   IERR )
            IF( IERR .NE. 0 ) RETURN
            GOTO 15
         ENDIF
C
C        LE TE EST IL DANS LE CADRE ENGLOBANT DE L'OBJET ?
         NS   = LETREE(6,NTE)
         NSTE = LETREE(7,NTE)
         NS3  = LETREE(8,NTE)
         IF( PXYD(1,NS) .GT. PXYD(1,NSTE) ) THEN
            DMIN = REAL( PXYD(1,NSTE) )
            DMAX = REAL( PXYD(1,NS) )
         ELSE
            DMIN = REAL( PXYD(1,NS) )
            DMAX = REAL( PXYD(1,NSTE) )
         ENDIF
         IF( (XRMIN .LE. DMIN .AND. DMIN .LE. XRMAX) .OR.
     %       (XRMIN .LE. DMAX .AND. DMAX .LE. XRMAX) ) THEN
            IF( PXYD(2,NS) .GT. PXYD(2,NS3) ) THEN
               DMIN = REAL( PXYD(2,NS3) )
               DMAX = REAL( PXYD(2,NS) )
            ELSE
               DMIN = REAL( PXYD(2,NS) )
               DMAX = REAL( PXYD(2,NS3) )
            ENDIF
            IF( (YRMIN .LE. DMIN .AND. DMIN .LE. YRMAX) .OR.
     %          (YRMIN .LE. DMAX .AND. DMAX .LE. YRMAX) ) THEN
C
C              TE MINIMAL ET INTERNE AU RECTANGLE ENGLOBANT
C              --------------------------------------------
C              LE CARRE DE LA LONGUEUR DE L'ARETE DU TE MINIMAL
               D2 = REAL( (PXYD(1,NS)-PXYD(1,NSTE))**2 +
     %                    (PXYD(2,NS)-PXYD(2,NSTE))**2 )
C              LE CARRE DE LA LONGUEUR DES ARETES SOUHAITEES ICI
               ARETS2 = REAL( MIN( DBLE(ARETE2), PXYD(3,NS)**2,
     %                        PXYD(3,NSTE)**2, PXYD(3,NS3)**2 ) )
               IF( D2 .GT. ARETS2 ) THEN
C                 LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C                 EN 4 SOUS-TRIANGLES
                  CALL DV4STE(NBSOMM,MXSOMM, PXYD, NTE, LETREE, NOLESO,
     &                        IERR )
                  IF( IERR .NE. 0 ) RETURN
                  GOTO 15
               ENDIF
C
C              PRISE EN COMPTE DE LA FONCTION "TAILLE_IDEALE" SI ELLE EXISTE
C              -------------------------------------------------------------
               IF( NOFOTI .GT. 0 ) THEN
C
C                 IL EXISTE ICI UNE FONCTION 'TAILLE_IDEALE'
C                 AJUSTEMENT DE LA TAILLE DE L'ARETE AVEC CELUI DE TAILLE_IDEALE
C                 AUTOUR DU SOMMET 1 DU TE MINIMAL INTERNE (NS=LETREE(6,NTE))
C
C                 CALCUL DE PXYZD(1:2,NS) DANS LE REPERE INITIAL => XYZ(1:3)
                  XY(1) = REAL((PXYD(1,NS)+PXYD(1,NSTE)+PXYD(1,NS3))/ 3)
                  XY(2) = REAL((PXYD(2,NS)+PXYD(2,NSTE)+PXYD(2,NS3))/ 3)
                  CALL DV2D3D( XY, ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                         XYZ )
C
C                 CHOIX DU DECOUPAGE OU NON DU TE SELON
C                 NUTYSU : NUMERO D'OPTION DE GENERATION DU MAILLAGE DE CETTE SU
C                          9 APPEL DIRECT DE SUEX99
C                          1 APPEL PAR LE SP SUEX01-QUNSAL (QUADRANGLE ALGEBRIQU
C                          6 APPEL PAR LE SP SUEX06-TRNSAL (TRIANGLE   ALGEBRIQU
C                          3 APPEL PAR LE SP SUEX03 (SURFACE B-SPLINE D'INTERPOL
C
                  IF( NUTYSU .EQ. 9 ) THEN
C
C                    APPEL DIRECT DE SUEX99
C                    ----------------------
C                    LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
                     DPARAF(1) = XYZ(1)
                     DPARAF(2) = XYZ(2)
                     DPARAF(3) = XYZ(3)
                     CALL FONVAL( NOFOTI, 3, DPARAF,  NCODEV, DTAILL )
                     IF( NCODEV .GT. 0 ) THEN
C
C                       TAILLE EST INITIALISEE -> VECTEUR ORIGIN + TAILLE
                        DTAILL = ABS( DTAILL )
                        XYZ(1) = REAL( ORIGIN(1) + DTAILL )
                        XYZ(2) = ORIGIN(2)
                        XYZ(3) = ORIGIN(3)
C                       PASSAGE AUX COORDONNEES (X,Y) DANS LE PLAN APRES
C                       MISE A L'ECHELLE
                        CALL DV3D2D( XYZ, ORIGIN, D2D3,
     %                               XMIN, YMIN, XYMXMI,  XY )
C                       TAILLE_IDEALE DANS LE PLAN APRES MISE A L'ECHELLE
                        TAILL2 = (XY(1)-XY0(1))**2 + (XY(2)-XY0(2))**2
                        TAILLE = SQRT( TAILL2 )
C
C                       MISE A JOUR DE LA TAILLE SOUHAITEE AUTOUR DES 3 SOMMETS
                        PXYD(3,NS  )= MIN( PXYD(3,NS  ), DBLE(TAILLE) )
                        PXYD(3,NSTE)= MIN( PXYD(3,NSTE), DBLE(TAILLE) )
                        PXYD(3,NS3 )= MIN( PXYD(3,NS3 ), DBLE(TAILLE) )
                        DO 25 I=0,3
C                          LE NUMERO EVENTUEL DU POINT I APPARTENANT A CE TE
                           NPT = ABS( LETREE(I,NTE) )
                           IF( NPT .GT. 0 ) THEN
                              PXYD(3,NPT)=MIN(PXYD(3,NPT),DBLE(TAILLE))
                           ENDIF
 25                     CONTINUE
C
                        IF( D2 .GT. 1.69*TAILL2 ) THEN
C                          LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C                          EN 4 SOUS-TRIANGLES
                           CALL DV4STE( NBSOMM, MXSOMM, PXYD,
     &                                  NTE, LETREE, NOLESO,
     &                                  IERR )
                           IF( IERR .NE. 0 ) RETURN
                           GOTO 15
                        ENDIF
                     ENDIF
C
                  ELSE IF( NUTYSU .EQ. 1 .OR. NUTYSU .EQ. 6 .OR.
     %                     NUTYSU .EQ. 3 ) THEN
C
C                    APPEL DE SUEX99 A PARTIR DE SUEX01 OU SUEX03 OU SUEX06
C                    ------------------------------------------------------
C                    CALCUL DES COORDONNEES DES 3 SOMMETS DU TE
C                    ET DES 3 TAILLE_IDEALE SOUHAITEE EN CES SOMMETS
C                    SI LA TAILLE D'UNE ARETE EST SUPERIEURE A 1.5 FOIS
C                    LE MINIMUM DES 3 TAILLES SOUHAITEES ALORS LE TE EST DECOUPE
C
                     DO 26 I=1,3
C
C                       RETOUR AUX COORDONNEES (R,S) DANS LE POLYGONE PLAN
                        XY(1) = REAL( PXYD( 1, NUSTE(I) ) )
                        XY(2) = REAL( PXYD( 2, NUSTE(I) ) )
                        CALL DV2D3D( XY, ORIGIN, D2D3,
     %                               XMIN, YMIN, XYMXMI, XYZ )
                        R = XYZ(1)
                        S = XYZ(2)
C
                        IF( NUTYSU .EQ. 1 ) THEN
C
C                          QUADRANGLE ALGEBRIQUE TRANSFINI
C                          ...............................
                           CALL XYZTRF( R, S, XYZSTE(1,I), IERR )
                           IF( IERR .NE. 0 ) RETURN
C
                        ELSE IF( NUTYSU .EQ. 3 ) THEN
C
C                          SURFACE B-SPLINE D'INTERPOLATION
C                          ................................
C                          LES VARIABLES DE L'APPEL APPARTENANT AU
C                          COMMON /S09S03/ DOIVENT ETRE INITIALISEES
                           CALL XYZBST( R, S,
     %                                  LDEXSB, LRX, RMCN(MNRX),
     %                                  LDEYSB, LRY, RMCN(MNRY),
     %                                  RMCN(MNSPQB),
     %                                  XYZSTE(1,I) )
C
                        ELSE IF( NUTYSU .EQ. 6 ) THEN
C
C                          TRIANGLE ALGEBRIQUE
C                          ...................
                           CALL XYZTTF( R, S, XYZSTE(1,I) )
C
                        ENDIF
C
C                       LA TAILLE_IDEALE EN CE SOMMET
C                       LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
                        DPARAF(1) = XYZSTE(1,I)
                        DPARAF(2) = XYZSTE(2,I)
                        DPARAF(3) = XYZSTE(3,I)
                        CALL FONVAL( NOFOTI, 3, DPARAF, NCODEV, DTAILL )
                        IF( NCODEV .LE. 0 ) GOTO 29
C                       LE CARRE DE LA TAILLE SOUHAITEE AUTOUR DU SOMMET I DU TR
                        HST(I) = REAL( DTAILL ** 2 )
C
 26                  CONTINUE
C
C                    LA TAILLE EFFECTIVE DES 3 ARETES DANS LE QUADRANGLE TRANSFI
                     DO 27 I=1,3
                        IF( I .LT. 3 ) THEN
                           I1 = I + 1
                        ELSE
                           I1 = 1
                        ENDIF
C                       LE CARRE DE LA LONGUEUR DE L'ARETE DROITE DU TRIANGLE
C                       SUR LE QUADRANGLE TRANSFINI
                        D2 = (XYZSTE(1,I1)-XYZSTE(1,I))**2 +
     %                       (XYZSTE(2,I1)-XYZSTE(2,I))**2 +
     %                       (XYZSTE(3,I1)-XYZSTE(3,I))**2
C
                        IF( D2 .GT. 1.69*HST( I) .OR.
     %                      D2 .GT. 1.69*HST(I1) ) THEN
C                          LE TRIANGLE NTE TROP GRAND DOIT ETRE SUBDIVISE
C                          EN 4 SOUS-TRIANGLES
                           CALL DV4STE( NBSOMM, MXSOMM, PXYD,
     &                                  NTE, LETREE, NOLESO, IERR )
                           IF( IERR .NE. 0 ) RETURN
                           GOTO 15
                        ENDIF
 27                  CONTINUE
C
                  ELSE
C
C                    NUTYSU INCORRECT
                     WRITE(IMPRIM,*) 'ERREUR DVHOTE: NUTYSU=',NUTYSU,
     %                               ' INCONNU'
                     IERR = 10
                     RETURN
                  ENDIF
               ENDIF
C
C              RECHERCHE DU NOMBRE DE NIVEAUX ENTRE NTE ET LES TE VOISINS
C              PAR SES ARETES
 29            DO 30 I=1,3
C
C                 NOTEVA TRIANGLE VOISIN DE NTE PAR L'ARETE I
                  CALL N1TRVA( NTE, I, LETREE, NOTEVA, NIVEAU )
                  IF( NOTEVA .LE. 0 ) GOTO 30
C                 IL EXISTE UN TE VOISIN
                  IF( NIVEAU .GT. 0 ) GOTO 30
C                 NTE A UN TE VOISIN PLUS PETIT OU EGAL
                  IF( LETREE(0,NOTEVA) .LE. 0 ) GOTO 30
C                 NTE A UN TE VOISIN NOTEVA SUBDIVISE AU MOINS UNE FOIS
C
                  IF( NBITER .GT. 0 ) THEN
C
C                    LES 2 SOUS TRIANGLES VOISINS SONT-ILS SUBDIVISES?
                     NSTE = LETREE(I,NOTEVA)
                     IF( LETREE(0,NSTE) .LE. 0 ) THEN
C                       NSTE N'EST PAS SUBDIVISE
                        NSTE = LETREE(NOSUI3(I),NOTEVA)
                        IF( LETREE(0,NSTE) .LE. 0 ) THEN
C                          LES 2 SOUS-TRIANGLES NE SONT PAS SUBDIVISES
                           GOTO 30
                        ENDIF
                     ENDIF
                  ENDIF
C
C                 LE TRIANGLE NTE DOIT ETRE SUBDIVISE EN 4 SOUS-TRIANGLES
                  CALL DV4STE(NBSOMM,MXSOMM, PXYD, NTE, LETREE, NOLESO,
     &                        IERR )
                  IF( IERR .NE. 0 ) RETURN
                  GOTO 15
C
 30            CONTINUE
            ENDIF
         ENDIF
         GOTO 10
      ENDIF
      IF( NBS0 .LT. NBSOMM ) THEN
         NBS0 = NBSOMM
         GOTO 5
      ENDIF
C
C     ====================================================
C     RENUMEROTATION FINALE DES SOMMETS NON ELIMINES DU TE
C     PAR ORDRE DECROISSANT DE LEUR TAILLE
C     ====================================================
      DO 34 I=0,MXSOMM
         NOPXYD(I) = 0
 34   CONTINUE
C
C     PARCOURS DE L'ARBRE PAR NIVEAUX
      NBS    = 0
      LEQUEU = 1
      LHQUEU = 0
C     LA RACINE DE LETREE INITIALISE LA QUEUE
      LAQUEU(1) = 1
C
C     TANT QUE LA LONGUEUR DE LA QUEUE EST >=0 TRAITER LE DEBUT DE QUEUE
 40   IF( LHQUEU .GE. 0 ) THEN
C
C        LE TRIANGLE TE A TRAITER
         I   = LEQUEU - LHQUEU
         IF( I .LE. 0 ) I = MXQUEU + I
         NTE = LAQUEU( I )
C        LA LONGUEUR EST REDUITE
         LHQUEU = LHQUEU - 1
C
C        NTE EST IL UN SOUS-TRIANGLE MINIMAL ?
         IF( LETREE(0,NTE) .GT. 0 ) THEN
C           NON LES 4 SOUS-TRIANGLES SONT MIS DANS LA QUEUE
            IF( LHQUEU + 4 .GE. MXQUEU ) THEN
               WRITE(IMPRIM,*) 'DVHOTE: SATURATION DE LA QUEUE'
               IERR = 1
               RETURN
            ENDIF
            DO 42 I=3,0,-1
C              AJOUT DU SOUS-TRIANGLE I
               LHQUEU = LHQUEU + 1
               LEQUEU = LEQUEU + 1
               IF( LEQUEU .GT. MXQUEU ) LEQUEU = LEQUEU - MXQUEU
               LAQUEU( LEQUEU ) = LETREE( I, NTE )
 42         CONTINUE
         ENDIF
C
C        PARCOURS DES 3 SOMMETS DU TE
         DO 44 I=6,8
            NS = LETREE(I,NTE)
            IF( NOLESO(NS) .GT. 0 ) THEN
C              NS EST UN SOMMET DE TE NON ELIMINE
               IF( NOPXYD(NS) .EQ. 0 ) THEN
C                 SOMMET DE TE NON ELIMINE NON ENCORE INSCRIT
                  NBS = NBS + 1
                  NOPXYD(NS) = NBS
               ENDIF
            ENDIF
 44      CONTINUE
         GOTO 40
      ENDIF
C
C     SUPPRESSION DES SOMMETS TROP PROCHES DES POINTS UTILISATEUR
C     NUMEROTATION FINALE DES POINTS INTERNES ET SOMMETS DES ARETES
C     DES LIGNES FERMEES  (LES SOMMETS DE L'ENVELOPPE CONVEXE SONT TRIANGULES)
C     APRES LES SOMMETS DE TE RETENUS POR ETRE TRAITES DANS
C     UN FOND DE MAILLAGE PREPARE
C     ========================================================================
      DO 50 I=NOSTE1+3,NBSOMM
         IF( NOLESO(I) .EQ. 0 ) THEN
C           AJOUT DU POINT
            NBS = NBS + 1
            NOPXYD(I) = NBS
C           SUPPRESSION DES SOMMETS TROP PRES
            CALL PPVTR2( PXYD(1,I), PXYD(3,I)*PROXPI,
     %                   PXYD, LETREE, NOLESO, NS, NT )
         ENDIF
 50   CONTINUE
C
C     CONSTRUCTION DU TABLEAU NOPXYD
C     ==============================
C     RENUMEROTATION DES SOMMETS DES TE NON SUPPRIMES
      NS = 0
      DO 60 I=NOSTE1+3,NBSOMM
         IF( NOLESO(I) .GT. 0 ) THEN
            NS = NS + 1
            NOPXYD(I) = NS
         ELSE
            NOPXYD(I) = 0
         ENDIF
 60   CONTINUE
C
C     LES NOSTE1-2 SOMMETS DE LA PREMIERE LIGNE SONT AJOUTES ENSUITE
      DO 70 I=1,NOSTE1-2
         NOPXYD(I) = I + NS
 70   CONTINUE
C
C     SUPPRESSION DES SOMMETS DE L'ENVELOPPE CONVEXE (MAILLES DIRECTEMENT)
      NA = 1
      DO 80 I=1,NBSTEC
         NOPXYD( NOSTEC(1,NA) ) = 0
         NA = NOSTEC(2,NA)
 80   CONTINUE
C
C     LE BARYCENTRE ET LES 3 SOMMETS DU TE RACINE
      DO 82 I=NOSTE1-1,NOSTE1+2
         NOPXYD(I) = 0
 82   CONTINUE
C
C     NUMEROTATION DES SOMMETS INITIAUX HORS ENVELOPPE CONVEXE
      DO 90 I=1,NOSTE1-2
         IF( NOPXYD(I) .GT. 0 ) THEN
            NS = NS + 1
            NOPXYD(I) = NS
         ENDIF
 90   CONTINUE
C
C     INVERSION DU TABLEAU NOPXYD
C     ===========================
C     NOPXYD( NOST ) = ORDRE   =>  NOPXYD( ORDRE ) = NOST
      DO 100 I=1,NBSOMM
         LAQUEU(I) = NOPXYD(I)
         NOPXYD(I) = 0
 100  CONTINUE
      DO 110 I=1,NBSOMM
         NOPXYD( LAQUEU(I) ) = I
 110  CONTINUE
      NOPXYD(0) = 0
C
C     MINIMISATION A ARETMX DES DISTANCES SOUHAITEES TROP GRANDES DES SOMMETS DE
C     ==========================================================================
      DO 190 I=NOSTE1+3, NBSOMM
         IF( NOLESO(I) .GT. 0 ) THEN
C           SOMMET DE TE NON ELIMINE
            PXYD(3,I) = MIN( PXYD(3,I), DBLE(ARETMX) )
         ENDIF
 190  CONTINUE
C
C     LE NOMBRE DE SOMMETS A TRIANGULER
      NBS = NS
      END
