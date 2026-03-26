      SUBROUTINE TEAMQS( NOFOTI, NUTYSU, RAP2P3,
     %                   NOARST, MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR,
     %                   MXTRCF, NOTRCF, NOSTBO,
     %                   N1ARCF, NOARCF, LARMIN,
     %                   ORIGIN, D2D3,   XMIN,   YMIN,   XYMXMI,
     %                   NBDPFR, NBDPII, NBSOMM,
     %                   PXYD, NSLIGN,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    UNE ITERATION DE BARYCENTRAGE DES POINTS INTERNES
C -----    MODIFICATION DE LA TOPOLOGIE POUR AVOIR 4 OU 5 OU 6 TRIANGLES
C          POUR CHAQUE SOMMET DE LA TRIANGULATION
C          MISE EN TRIANGULATION DELAUNAY
C
C ENTREES:
C --------
C NOFOTI : NUMERO DANS LE LX DES FONCTIONS DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C          >0 SI ELLE EXISTE ET 0 SINON
C NUTYSU : NUMERO D'OPTION DE GENERATION DU MAILLAGE DE CETTE SURFACE
C          9 APPEL DIRECT A PARTIR DU SP SUEX09
C          1 APPEL PAR LE SP SUEX01-QUNSAL (QUADRANGLE ALGEBRIQUE)
C          3 APPEL PAR LE SP SUEX03        (B-SPLINE POLYNOMIALE )
C          6 APPEL PAR LE SP SUEX06-TRNSAL (TRIANGLE   ALGEBRIQUE)
C RAP2P3 : RAPPORT DU PERIMETRE DE L'ENVELOPPE DU DOMAINE PLAN SUR
C                  LE PERIMETRE DE L'ENVELOPPE DU DOMAINE DE R**3
C          APPROXIMATION POUR CALCULER RAPIDEMENT LA TAILLE DE L'ARETE
C          DANS LE PLAN CONNAISSANT SA TAILLE DANS R**3
C NOARST : NOARST(I) NUMERO D'UNE ARETE DE SOMMET I
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C MXSOAR : NOMBRE MAXIMAL D'ARETES FRONTALIERES DECLARABLES
C N1SOAR : NUMERO DE LA PREMIERE ARETE VIDE DANS LE TABLEAU NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARTR : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOARTR
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C MXTRCF : NOMBRE MAXIMAL DE TRIANGLES EMPILABLES
C NBDPFR : NUMERO DU DERNIER SOMMET FRONTALIER
C NBDPII : NUMERO DU DERNIER SOMMET INTERNE IMPOSE
C          =NBDPFR SI PAS DE POINT INTERNE IMPOSE
C NSLIGN : >0 => NS NUMERO DU POINT DANS LE LEXIQUE POINT SI INTERNE IMPOSE
C          OU => 1 000 000 * N + NS1
C              OU N   EST LE NUMERO (1 A NBLFTR) DE LA LIGNE DE CE POINT
C                 NS1 EST LE NUMERO DU POINT DANS SA LIGNE
C          = 0 SI LE POINT EST INTERNE NON IMPOSE PAR L'UTILISATEUR
C          =-1 SI LE SOMMET EST EXTERNE AU DOMAINE
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS DE LA TRIANGULATION
C          (CERTAINS SOMMETS INTERNES ONT ETE DESACTIVES OU AJOUTES)
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C
C AUXILIAIRES:
C ------------
C NOTRCF : TABLEAU ( MXTRCF ) AUXILIAIRE D'ENTIERS
C          NUMERO DANS NOARTR DES TRIANGLES DE SOMMET NS
C NOSTBO : TABLEAU ( MXTRCF ) AUXILIAIRE D'ENTIERS
C          NUMERO DANS PXYD DES SOMMETS DES ARETES SIMPLES DE LA BOULE
C N1ARCF : TABLEAU (0:MXTRCF) AUXILIAIRE D'ENTIERS
C NOARCF : TABLEAU (3,MXTRCF) AUXILIAIRE D'ENTIERS
C LARMIN : TABLEAU ( MXTRCF ) AUXILIAIRE D'ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE    UPMC PARIS  MAI     1997
C MODIFS : ALAIN PERRONNET LABORATOIRE JL LIONS UPMC PARIS  OCTOBRE 2006
C....................................................................012
      PARAMETER        (LCHAIN=6)
      include"./incl/ampli.inc"
      include"./incl/langue.inc"
      INCLUDE"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  RAP2P3, PXYD(3,*)
      DOUBLE PRECISION  PONDER, PONDE1, XBAR, YBAR, X, Y, D, DMIN, DMAX
      DOUBLE PRECISION  SURTD2, XXX, YYY, ZZZ, S0, S1, S
      DOUBLE PRECISION  D2D3(3,3)
      REAL              ORIGIN(3), XYZ(3)
      INTEGER           NOARTR(MOARTR,MXARTR),
     %                  NOSOAR(MOSOAR,MXSOAR),
     %                  NOARST(*),
     %                  NOTRCF(MXTRCF),
     %                  NSLIGN(*),
     %                  NOSTBO(*),
     %                  N1ARCF(0:MXTRCF),
     %                  NOARCF(3,MXTRCF),
     %                  LARMIN(MXTRCF)
      INTEGER           NOSOTR(3)
C
C     LE NOMBRE D'ITERATIONS POUR AMELIORER LA QUALITE
      NBITAQ = 6
      IER    = 0
      NOAR0  = 0
C
C     LE CHAINAGE DES ARETES A RENDRE DELAUNAY EST MARQUE A -1
      DO 1 NA=1,MXSOAR
         NOSOAR(LCHAIN,NA) = -1
 1    CONTINUE
C
C     INITIALISATION DU PARCOURS
      NBS1 = NBSOMM
      NBS2 = NBDPII + 1
C     => PAS DE TRAITEMENT SUR LES POINTS DES LIGNES DE LA FRONTIERE
      NBS3 = -1
C
      DO 5000 ITER=1,NBITAQ
C
C        LE NOMBRE DE SOMMETS SUPPRIMES
         NBSTSU = 0
C
C        LES COMPTEURS DE PASSAGE SUR LES DIFFERENTS CAS
         NBST4 = 0
         NBST5 = 0
CCC         NBST8 = 0
C
C        COEFFICIENT DE PONDERATION CROISSANT AVEC LES ITERATIONS
         PONDER = 0.1D0 + ITER * 0.9D0 / NBITAQ
CCC 27 MARS    2008 PONDER = 0.1D0 + ITER * 0.9D0 / NBITAQ
CCC 10 OCTOBRE 2006 PONDER = MIN( 1D0, 0.1D0 + ITER * 0.9D0 / NBITAQ )
CCC  9 MARS    2006 PONDER = MIN( 1D0, ( 50 + (50*ITER)/NBITAQ ) * 0.01D0 )
         PONDE1 = 1D0 - PONDER
C
C        L'ORDRE DU PARCOURS DANS LE SENS CROISSANT OU DECROISSANT
         NT   = NBS1
         NBS1 = NBS2
         NBS2 = NT
C        ALTERNANCE DU PARCOURS
         NBS3 = -NBS3
C
         DO 1000 NS = NBS1, NBS2, NBS3
C
C           LE SOMMET EST IL INTERNE AU DOMAINE?
            IF( NSLIGN(NS) .NE. 0 ) GOTO 1000
C
C           TRAITEMENT D'UN SOMMET INTERNE NON IMPOSE PAR L'UTILISATEUR
C           ===========================================================
C           EXISTE-T-IL UNE ARETE DE SOMMET NS ?
            NOAR = NOARST( NS )
            IF( NOAR .LE. 0 ) GOTO 1000
C
C           LE 1-ER TRIANGLE DE L'ARETE NOAR
            NT = NOSOAR( 4, NOAR )
            IF( NT .LE. 0 ) GOTO 1000
C
C           RECHERCHE DES TRIANGLES DE SOMMET NS
C           ILS DOIVENT FORMER UN CONTOUR FERME DE TYPE ETOILE
            CALL TRP1ST( NS,     NOARST, MOSOAR, NOSOAR,
     %                   MOARTR, MXARTR, NOARTR,
     %                   MXTRCF, NBTRCF, NOTRCF )
            IF( NBTRCF .LE. 2 ) GOTO 1000
C
C           BOUCLE SUR LES TRIANGLES QUI FORMENT UNE BOULE AUTOUR DU SOMMET NS
            NBSTBO = 0
C           CHAINAGE DES ARETES SIMPLES DE LA BOULE A RENDRE DELAUNAY
C
C           CHAINAGE DES ARETES SIMPLES DE LA BOULE A RENDRE DELAUNAY
C           REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA BOULE NS
            NOAR1 = NOAR0
 16         IF( NOAR1 .GT. 0 ) THEN
C              PROTECTION DU NO DE L'ARETE SUIVANTE
               NA = NOSOAR(LCHAIN,NOAR1)
C              L'ARETE INTERNE EST REMISE A -1
               NOSOAR(LCHAIN,NOAR1) = -1
C              L'ARETE SUIVANTE
               NOAR1 = NA
               GOTO 16
            ENDIF
            NOAR0 = 0
C
            S0 = 0D0
            DO 40 I=1,NBTRCF
C
C              LE NUMERO DE L'ARETE DU TRIANGLE NT NE CONTENANT PAS LE SOMMET NS
               NT = NOTRCF(I)
C              LE NUMERO DES 3 SOMMETS DANS LE SENS DIRECT
               CALL NUSOTR( NT, MOSOAR, NOSOAR,
     %                      MOARTR, NOARTR, NOSOTR )
               S = SURTD2( PXYD(1,NOSOTR(1)),
     %                     PXYD(1,NOSOTR(2)),
     %                     PXYD(1,NOSOTR(3)) )
               IF( S .LE. 0 ) THEN
                   PRINT *,'TEAMQS: TRIANGLE',NT,' D AIRE NEGATIVE'
                   PRINT *,'TEAMQS: ABANDON DU TRAITEMENT DU SOMMET',NS
                   GOTO 1000
               ENDIF
               S0 = S0 + ABS(S)
C
               DO 20 NA=1,3
C                 LE NUMERO DE L'ARETE NA DANS LE TABLEAU NOSOAR
                  NOAR = ABS( NOARTR(NA,NT) )
                  IF( NOSOAR(1,NOAR) .NE. NS   .AND.
     %                NOSOAR(2,NOAR) .NE. NS ) GOTO 25
 20            CONTINUE
C
C              CONSTRUCTION DE LA LISTE DES SOMMETS DES ARETES SIMPLES
C              DE LA BOULE DES TRIANGLES DE SOMMET NS
C              -------------------------------------------------------
 25            DO 35 NA=1,2
                  NS1 = NOSOAR(NA,NOAR)
                  DO 30 J=NBSTBO,1,-1
                     IF( NS1 .EQ. NOSTBO(J) ) GOTO 35
 30               CONTINUE
C                 NS1 EST UN NOUVEAU SOMMET A AJOUTER
                  NBSTBO = NBSTBO + 1
                  NOSTBO(NBSTBO) = NS1
 35            CONTINUE
C
C              NOAR EST UNE ARETE POTENTIELLE A RENDRE DELAUNAY
               IF( NOSOAR(3,NOAR) .EQ. 0 ) THEN
C                 ARETE NON FRONTALIERE
                  NOSOAR(LCHAIN,NOAR) = NOAR0
                  NOAR0 = NOAR
               ENDIF
C
 40         CONTINUE
C
C           CALCUL DES 2 COORDONNEES DU BARYCENTRE DE LA BOULE DU SOMMET NS
C           CALCUL DE L'ARETE DE TAILLE MAXIMALE ET MINIMALE ISSUE DE NS
C           ---------------------------------------------------------------
            XBAR = 0D0
            YBAR = 0D0
            DMIN = 1D28
            DMAX = 0D0
            DO 50 I=1,NBSTBO
               X    = PXYD(1,NOSTBO(I))
               Y    = PXYD(2,NOSTBO(I))
               XBAR = XBAR + X
               YBAR = YBAR + Y
               D    = (X-PXYD(1,NS)) ** 2 + (Y-PXYD(2,NS)) ** 2
               IF( D .GT. DMAX ) THEN
                  DMAX = D
                  IMAX = I
               ENDIF
               IF( D .LT. DMIN ) THEN
                  DMIN = D
                  IMIN = I
               ENDIF
 50         CONTINUE
C
C           PAS DE MODIFICATION DE LA TOPOLOGIE LORS DE LA DERNIERE ITERATION
C           =================================================================
            IF( ITER .GE. NBITAQ ) GOTO 200
C
ccc            ccc le 11/4/2008
cccC           SI LE POINT A SUPPRIMER A SON PLUS PROCHE VOISIN SOUS LA DISTANCE
cccC           SOUHAITEE EN NS ALORS IL N'EST PAS DETRUIT  (SEPTEMBRE 2006)
cccC           =================================================================
ccc            D = MIN( PXYD(3,NS), PXYD(3,NOSTBO(IMIN)) )
ccc            IF( DMIN .LT. (AMPLI*D)**2 ) GOTO 200
C
C           SI CE POINT NS EST DANS UNE BOULE DE L'UN DES POINTS INTERNES
C           IMPOSES (RAYON=3 DIST SOUHAITEE) ALORS IL N'EST PAS DETRUIT
            CALL PTBPTI( NS, NBDPFR, NBDPII, PXYD, NONOUI )
            IF( NONOUI .NE. 0 ) GOTO 200
C
C           SI LA BOULE DE NS CONTIENT 3 ou 4 TRIANGLES LE SOMMET NS EST DETRUIT
C           ====================================================================
            IF( NBTRCF .EQ. 3 .OR. NBTRCF .EQ. 4 ) THEN
C
C              REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA BOULE NS
               NOAR = NOAR0
 60            IF( NOAR .GT. 0 ) THEN
C                 PROTECTION DU NO DE L'ARETE SUIVANTE
                  NA = NOSOAR(LCHAIN,NOAR)
C                 L'ARETE INTERNE EST REMISE A -1
                  NOSOAR(LCHAIN,NOAR) = -1
C                 L'ARETE SUIVANTE
                  NOAR = NA
                  GOTO 60
               ENDIF
               CALL TE1STM( NS,     NBDPII, PXYD,   NOARST,
     %                      MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                      MOARTR, MXARTR, N1ARTR, NOARTR,
     %                      MXTRCF, N1ARCF, NOARCF,
     %                      LARMIN, NOTRCF, NOSTBO,
     %                      IERR )
               IF( IERR .EQ. -543 ) THEN
                  IERR = 0
                  GOTO 1000
               ELSE IF( IERR .LT. 0 ) THEN
C                 LE SOMMET NS EST EXTERNE DONC NON SUPPRIME
C                 OU BIEN LE SOMMET NS EST LE CENTRE D'UN CF DONT TOUTES
C                 LES ARETES SIMPLES SONT FRONTALIERES
C                 DANS LES 2 CAS LE SOMMET NS N'EST PAS SUPPRIME
                  IERR = 0
                  GOTO 200
               ELSE IF( IERR .EQ. 0 ) THEN
                  NBST4  = NBST4 + 1
                  NBSTSU = NBSTSU + 1
               ELSE
C                 ERREUR IRRECUPERABLE
               PRINT *,'TEAMQS: ERREUR1 IRRECUPERABLE EN SORTIE TE1STM'
                  GOTO 9999
               ENDIF
               GOTO 1000
C
            ENDIF
C
C           SI LA BOULE DE NS CONTIENT 5 TRIANGLES ET A UN SOMMET VOISIN
C           SOMMET DE 5 TRIANGLES ALORS L'ARETE JOIGNANT CES 2 SOMMETS
C           EST TRANSFORMEE EN UN SEUL SOMMET DE 6 TRIANGLES
C           ============================================================
            IF( NBTRCF .EQ. 5 ) THEN
C
               DO 80 I=1,5
C                 LE NUMERO DU SOMMET DE L'ARETE I ET DIFFERENT DE NS
                  NS1 = NOSTBO(I)
                  IF( NS1 .LE. NBDPII ) GOTO 80
C
C                 LA LISTE DES TRIANGLES DE SOMMET NS1
                  CALL TRP1ST( NS1, NOARST,
     %                         MOSOAR, NOSOAR, MOARTR, MXARTR, NOARTR,
     %                         MXTRCF-5, NBTRC1, NOTRCF(6) )
                  IF( NBTRC1 .EQ. 5 ) THEN
C
C                    L'ARETE DE SOMMETS NS-NS1 DEVIENT UN POINT
C                    PAR SUPPRESSION DU SOMMET NS
C
C                    REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA BOUL
                     NOAR = NOAR0
 70                  IF( NOAR .GT. 0 ) THEN
C                       PROTECTION DU NO DE L'ARETE SUIVANTE
                        NA = NOSOAR(LCHAIN,NOAR)
C                       L'ARETE INTERNE EST REMISE A -1
                        NOSOAR(LCHAIN,NOAR) = -1
C                       L'ARETE SUIVANTE
                        NOAR = NA
                        GOTO 70
                     ENDIF
C
C                    LE POINT NS1 DEVIENT LE MILIEU DE L'ARETE NS-NS1
                     X = PXYD(1,NS1)
                     Y = PXYD(2,NS1)
                     D = PXYD(3,NS1)
                     DO 75 J=1,3
                        PXYD(J,NS1) = (PXYD(J,NS) + PXYD(J,NS1)) * 0.5D0
 75                  CONTINUE
C
                     IF( NOFOTI .GT. 0 ) THEN
C                       LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                       CALCUL DE PXYZD(3,NS1) DANS LE REPERE INITIAL => XYZ(1:3
                        CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                               PXYD(1,NS1), PXYD(2,NS1),
     %                               ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                               XYZ,    PXYD(3,NS1), IER )
                     ENDIF
C
C                    SUPPRESSION DU POINT NS ET MISE EN DELAUNAY
                     CALL TE1STM( NS,     NBDPII, PXYD,   NOARST,
     %                            MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                            MOARTR, MXARTR, N1ARTR, NOARTR,
     %                            MXTRCF, N1ARCF, NOARCF,
     %                            LARMIN, NOTRCF, NOSTBO,
     %                            IERR )
                     IF( IERR .LT. 0 ) THEN
C                       RESTAURATION DU SOMMET NS1 A SON ANCIENNE PLACE
                        PXYD(1,NS1) = X
                        PXYD(2,NS1) = Y
                        PXYD(3,NS1) = D
                        IERR = 0
                        GOTO 1000
                     ELSE IF( IERR .EQ. 0 ) THEN
                        NBST5  = NBST5 + 1
                        NBSTSU = NBSTSU + 1
                        GOTO 1000
                     ELSE
C                       ERREUR IRRECUPERABLE
               PRINT *,'TEAMQS: ERREUR2 IRRECUPERABLE EN SORTIE TE1STM'
                        GOTO 9999
                     ENDIF
                  ENDIF
 80            CONTINUE
            ENDIF
CCCC
CCCC           SI LA BOULE DE NS CONTIENT AU MOINS 8 TRIANGLES
CCCC           ALORS UN TRIANGLE INTERNE EST AJOUTE + 3 TRIANGLES (1 PAR ARETE)
CCCC           ================================================================
CCC            IF( NBTRCF .GE. 8 ) THEN
CCCC
CCCC              MODIFICATION DES COORDONNEES DU SOMMET NS
CCCC              IL DEVIENT LE BARYCENTRE DU TRIANGLE NOTRCF(1)
CCC               CALL NUSOTR( NOTRCF(1), MOSOAR, NOSOAR,
CCC     %                      MOARTR, NOARTR, NOSOTR )
CCC               DO 110 I=1,3
CCC                  PXYD(I,NS) = ( PXYD(I,NOSOTR(1))
CCC     %                         + PXYD(I,NOSOTR(2))
CCC     %                         + PXYD(I,NOSOTR(3)) ) / 3D0
CCC 110           CONTINUE
CCCC
CCC               IF( NOFOTI .GT. 0 ) THEN
CCCC                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
CCCC                 CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(1:3
CCC                  CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
CCC     %                         PXYD(1,NS), PXYD(2,NS),
CCC     %                         ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
CCC     %                         XYZ,    PXYD(3,NS), IER )
CCC               ENDIF
CCCC
CCCC              AJOUT DES 2 AUTRES SOMMETS COMME BARYCENTRES DES TRIANGLES
CCCC              NOTRCF(1+NBTRCF/3) ET NOTRCF(1+2*NBTRCF/3)
CCC               NBT1 = ( NBTRCF + 1 ) / 3
CCC               DO 140 N=1,2
CCCC
CCCC                 LE TRIANGLE TRAITE
CCC                  NT = NOTRCF(1 + N * NBT1 )
CCCC
CCCC                 LE NUMERO PXYD DE SES 3 SOMMETS
CCC                  CALL NUSOTR( NT, MOSOAR, NOSOAR,
CCC     %                         MOARTR, NOARTR, NOSOTR )
CCCC
CCCC                 AJOUT DU NOUVEAU BARYCENTRE
CCC                  IF( NBSOMM .GE. MXSOMM ) THEN
CCC                     WRITE(IMPRIM,*) 'SATURATION DU TABLEAU PXYD'
CCCC                    ABANDON DE L'AMELIORATION
CCC                     GOTO 1100
CCC                  ENDIF
CCC                  NBSOMM = NBSOMM + 1
CCC                  DO 120 I=1,3
CCC                     PXYD(I,NBSOMM) = ( PXYD(I,NOSOTR(1))
CCC     %                                + PXYD(I,NOSOTR(2))
CCC     %                                + PXYD(I,NOSOTR(3)) ) / 3D0
CCC 120              CONTINUE
CCCC
CCC                  IF( NOFOTI .GT. 0 ) THEN
CCCC                    LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
CCCC                    CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(
CCC                     CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
CCC     %                            PXYD(1,NBSOMM), PXYD(2,NBSOMM),
CCC     %                            ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
CCC     %                            XYZ,    PXYD(3,NBSOMM), IER )
CCC                  ENDIF
CCCC
CCCC                 SOMMET INTERNE A LA TRIANGULATION
CCC                  NSLIGN(NBSOMM) = 0
CCCC
CCCC                 LES 3 ARETES DU TRIANGLE NT SONT A RENDRE DELAUNAY
CCC                  DO 130 I=1,3
CCC                     NOAR = ABS( NOARTR(I,NT) )
CCC                     IF( NOSOAR(3,NOAR) .EQ. 0 ) THEN
CCCC                       ARETE NON FRONTALIERE
CCC                        IF( NOSOAR(LCHAIN,NOAR) .LT. 0 ) THEN
CCCC                          ARETE NON ENCORE CHAINEE
CCC                           NOSOAR(LCHAIN,NOAR) = NOAR0
CCC                           NOAR0 = NOAR
CCC                        ENDIF
CCC                     ENDIF
CCC 130              CONTINUE
CCCC
CCCC                 TRIANGULATION DU TRIANGLE DE BARYCENTRE NBSOMM
CCCC                 PROTECTION A NE PAS MODIFIER SINON ERREUR!
CCC                  CALL TR3STR( NBSOMM, NT,
CCC     %                         MOSOAR, MXSOAR, N1SOAR, NOSOAR,
CCC     %                         MOARTR, MXARTR, N1ARTR, NOARTR,
CCC     %                         NOARST,
CCC     %                         NOSOTR, IERR )
CCC                  IF( IERR .NE. 0 ) THEN
CCC               PRINT *,'TEAMQS: ERREUR IRRECUPERABLE EN SORTIE TR3STR'
CCC                     GOTO 9999
CCC                  ENDIF
CCC 140           CONTINUE
CCCC
CCC               NBST8  = NBST8 + 1
CCCC
CCCC              LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
CCC               GOTO 300
CCCC
CCC            ENDIF
C
C           NBTRCF EST SUPERIEUR A 5 => BARYCENTRAGE SIMPLE
C           ===============================================
C           LES 2 COORDONNEES DU BARYCENTRE DES SOMMETS DES ARETES
C           SIMPLES DE LA BOULE DU SOMMET NS
 200        XBAR = XBAR / NBSTBO
            YBAR = YBAR / NBSTBO
C
C           PONDERATION POUR EVITER LES DEGENERESCENSES AVEC PROTECTION
            XXX = PXYD(1,NS)
            YYY = PXYD(2,NS)
            ZZZ = PXYD(2,NS)
C
            PXYD(1,NS) = PONDE1 * PXYD(1,NS) + PONDER * XBAR
            PXYD(2,NS) = PONDE1 * PXYD(2,NS) + PONDER * YBAR
C
C DEBUT AJOUT 10 OCTOBRE 2006
          IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C              CALCUL DE PXYZD(3,NS) DANS LE REPERE INITIAL => XYZ(1:3
               CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                      PXYD(1,NS), PXYD(2,NS),
     %                      ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                      XYZ,    PXYD(3,NS), IER )
            ENDIF
C
C           CALCUL DES SURFACES APRES DEPLACEMENT DE NS
            NS2 = 0
            NS3 = 0
            S1  = 0D0
            DO 210 I=1,NBTRCF
C              NUMERO DE L'ARETE DU TRIANGLE NT NE CONTENANT PAS LE SOMMET NS
               NT = NOTRCF(I)
               DO 204 NA=1,3
C                 LE NUMERO DE L'ARETE NA DANS LE TABLEAU NOSOAR
                  NOAR = ABS( NOARTR(NA,NT) )
                  IF( NOSOAR(1,NOAR) .NE. NS   .AND.
     %                NOSOAR(2,NOAR) .NE. NS ) THEN
                     NS2 = NOSOAR(1,NOAR)
                     NS3 = NOSOAR(2,NOAR)
                     GOTO 206
                  ENDIF
 204           CONTINUE
C              AIRE SIGNEE DU TRIANGLE
 206           S1 = S1 + ABS(SURTD2(PXYD(1,NS),PXYD(1,NS2),PXYD(1,NS3)))
 210        CONTINUE
C
            IF( ABS(S0-S1) .GT. 1D-10*ABS(S0) ) THEN
C
C              RETOUR A LA POSITION INITIALE
C              CAR LE POINT EST PASSE AU DELA D'UNE ARETE DE SON ETOILE
               WRITE(IMPRIM,*) 'TEAMQS: BARYCENTRAGE SOMMET',NS,
     %              '   XB=',PXYD(1,NS),' YB=',PXYD(2,NS),
     %              ' => ETOILE D AIRE INCORRECTE'
               WRITE(IMPRIM,*) 'TEAMQS: ST REMIS EN X =',XXX,
     %              ' Y =',YYY,' ARETE=',ZZZ
               PXYD(1,NS) = XXX
               PXYD(2,NS) = YYY
               PXYD(3,NS) = ZZZ
C              LA PONDERATION EST REDUITE  10 OCTOBRE 2006
               PONDER = MAX( 0.1D0, PONDER*0.5D0 )
               PONDE1 = 1D0 - PONDER
CCC               PRINT *,'PONDERATION DU BARYCENTRE PONDER=',PONDER
               GOTO 1000
            ENDIF
C
C FIN AJOUT 10 OCTOBRE 2006
C
C           LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
            CALL TEDELA( PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOAR0,
     %                   MOARTR, MXARTR, NOARTR, MODIFS )
C
 1000    CONTINUE
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,11000) ITER, NBST4, NBST5
         ELSE
            WRITE(IMPRIM,21000) ITER, NBST4, NBST5
         ENDIF
C                                            , NBST8
11000    FORMAT(' TEAMQS: ITERATION',I2,':',I7,' SOMMETS de 4Tr',
     %           I7,' SOMMETS 5Tr+5Tr')
21000    FORMAT(' TEAMQS: ITERATION',I2,':',I7,' VERTICES of 4Tr',
     %           I7,' VERTICES of 5Tr+5Tr')
CCC     %        I7,' SOMMETS >7T' )
C
cccC        TRACE DE LA TRIANGULATION ACTUELLE ET CALCUL DE LA QUALITE
ccc         CALL TETRMA( 1+MOD(ITER,14), NCNOIR,
ccc     %                COMXMI, PXYD,
ccc     %                MOSOAR, MXSOAR, NOSOAR,
ccc     %                MOARTR, MXARTR, NOARTR,
ccc     %                NBTRIA, QUAMOY, QUAMIN )
C
C        MISE A JOUR POUR NE PAS OUBLIER LES NOUVEAUX SOMMETS
         IF( NBS1 .GT. NBS2 ) THEN
            NBS1 = NBSOMM
            NBS2 = NBDPII + 1
         ELSE
            NBS1 = NBDPII + 1
            NBS2 = NBSOMM
         ENDIF
C
 5000 CONTINUE

 9999 RETURN
      END
