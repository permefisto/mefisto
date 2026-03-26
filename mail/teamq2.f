      SUBROUTINE TEAMQ2( NOFOTI, NUTYSU,
     %                   NOARST, MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR,
     %                   MXTRCF, NOTRCF, NOSTBO,
     %                   N1ARCF, NOARCF, LARMIN,
     %                   ORIGIN, D2D3,   XMIN,   YMIN,  XYMXMI,
     %                   NBDPII, NBSOMM,
     %                   nbarpi, PXYD,   NSLIGN,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: LES ARETES CROISSENT ET LA TOPOLOGIE AUGMENTE AVEC DES POINTS EN -
C ---- BOULE 3, 4, 5+5 TRIANGLES => SUPPRESSION D'UN POINT
C      DMOYENNE<0.5 TAILLE_IDEALE => SUPPRESSION DU POINT
C      PUIS MISE EN DELAUNAY DE LA TRIANGULATION
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
C NBDPII : NUMERO DU DERNIER SOMMET FRONTALIER OU INTERNE IMPOSE
C NSLIGN : TABLEAU DU NUMERO DE SOMMET DANS SA LIGNE POUR CHAQUE
C          SOMMET FRONTALIER
C          NUMERO DU POINT DANS LE LEXIQUE POINT SI INTERNE IMPOSE
C          0 SI LE POINT EST INTERNE NON IMPOSE PAR L'UTILISATEUR
C         -1 SI LE SOMMET EST EXTERNE AU DOMAINE
C COMXMI : MIN ET MAX DES COORDONNEEES DES SOMMETS DU MAILLAGE
C NBARPI : NUMERO DU DERNIER POINT INTERNE IMPOSE PAR L'UTILISATEUR
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
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1997
C....................................................................012
      PARAMETER        (LCHAIN=6)
      include"./incl/ampli.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PXYD(3,*)
      DOUBLE PRECISION  PONDER, PONDE1, XBAR, YBAR, X, Y, D, DMOY
      DOUBLE PRECISION  D2D3(3,3)
      REAL              ORIGIN(3), XYZ(3)
      INTEGER           NOARTR(MOARTR,*),
     %                  NOSOAR(MOSOAR,*),
     %                  NOARST(*),
     %                  NOTRCF(MXTRCF),
     %                  NSLIGN(*),
     %                  NOSTBO(*),
     %                  N1ARCF(0:MXTRCF),
     %                  NOARCF(3,MXTRCF),
     %                  LARMIN(MXTRCF)
C
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION
CCC      TRATRI = .TRUE.
      TRATRI = .FALSE.
C
C     TRACE DE LA TRIANGULATION ACTUELLE ET CALCUL DE LA QUALITE
      CALL TETRMA( NCBLAN, NCNOIR,
     %             COMXMI, PXYD,
     %             MOSOAR, MXSOAR, NOSOAR,
     %             MOARTR, MXARTR, NOARTR,
     %             NBTRIA, QUAMOY, QUAMIN )
C
C     LE NOMBRE D'ITERATIONS POUR AMELIORER LA QUALITE
      NBITAQ = 5
      IER    = 0
C
C     INITIALISATION DU PARCOURS
      NBS1 = NBSOMM
      NBS2 = NBDPII + 1
      NBS3 = -1
C
      DO 5000 ITER=1,NBITAQ
C
C        LE NOMBRE DE SOMMETS SUPPRIMES
         NBSTSU = 0
C
CCCC        LES COMPTEURS DE PASSAGE SUR LES DIFFERENTS CAS
CCC         NBST4 = 0
CCC         NBST5 = 0
CCC         NBST8 = 0
C
C        COEFFICIENT DE PONDERATION CROISSANT AVEC LES ITERATIONS
         PONDER = MIN( 1D0, ( 50 + (50*ITER)/NBITAQ ) * 0.01D0 )
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
            CALL TRP1ST( NS, NOARST, MOSOAR, NOSOAR,
     %                   MOARTR, MXARTR, NOARTR,
     %                   MXTRCF, NBTRCF, NOTRCF )
            IF( NBTRCF .LE. 0 ) GOTO 1000
C
C           MISE A JOUR DE LA DISTANCE SOUHAITEE
            IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C              CALCUL DE PXYZD(3,NS) DANS LE REPERE INITIAL => XYZ(1:3)
               CALL DVTAIL( NUTYSU, NOFOTI,
     %                      PXYD(1,NS), PXYD(2,NS),
     %                      ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                      XYZ,    D, IER )
               IF( IER .EQ. 0 ) THEN
C                 LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
                  R = REAL( XMIN + D )
                  PXYD(3,NS) = MISECH( R, XMIN, XYMXMI )
               ENDIF
            ENDIF
C
C           BOUCLE SUR LES TRIANGLES QUI FORMENT UNE BOULE AUTOUR DU SOMMET NS
            NBSTBO = 0
C           CHAINAGE DES ARETES SIMPLES DE LA BOULE A RENDRE DELAUNAY
            NOAR0  = 0
            DO 40 I=1,NBTRCF
C
C              LE NUMERO DE L'ARETE DU TRIANGLE NT NE CONTENANT PAS LE SOMMET NS
               NT = NOTRCF(I)
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
C           CALCUL DE LA LONGUEUR MOYENNE DES ARETES ISSUES DU SOMMET NS
C           ---------------------------------------------------------------
            XBAR = 0D0
            YBAR = 0D0
            DMOY = 0D0
            DO 50 I=1,NBSTBO
               X    = PXYD(1,NOSTBO(I))
               Y    = PXYD(2,NOSTBO(I))
               XBAR = XBAR + X
               YBAR = YBAR + Y
               DMOY = DMOY + SQRT( (X-PXYD(1,NS))**2+(Y-PXYD(2,NS))**2 )
 50         CONTINUE
            DMOY = DMOY / NBSTBO
C
C           REGULARISATION PAR BARYCENTRAGE SEULE POUR LES 2 DERNIERES ITERATION
            IF( ITER .GT. NBITAQ-2 ) GOTO 200
CCCC
CCCC           SI LA BOULE DE NS CONTIENT 3 OU 4 TRIANGLES LE SOMMET NS EST DETR
CCCC           =================================================================
CCC            IF( NBTRCF .LE. 4 ) THEN
CCCC
CCCC              REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA BOULE N
CCC               NOAR = NOAR0
CCC 60            IF( NOAR .GT. 0 ) THEN
CCCC                 PROTECTION DU NO DE L'ARETE SUIVANTE
CCC                  NA = NOSOAR(LCHAIN,NOAR)
CCCC                 L'ARETE INTERNE EST REMISE A -1
CCC                  NOSOAR(LCHAIN,NOAR) = -1
CCCC                 L'ARETE SUIVANTE
CCC                  NOAR = NA
CCC                  GOTO 60
CCC               ENDIF
CCC               CALL TE1STM( NS,     nbarpi, PXYD,   NOARST,
CCC     %                      MOSOAR, MXSOAR, N1SOAR, NOSOAR,
CCC     %                      MOARTR, MXARTR, N1ARTR, NOARTR,
CCC     %                      MXTRCF, N1ARCF, NOARCF,
CCC     %                      LARMIN, NOTRCF, NOSTBO,
CCC     %                      IERR )
CCC               IF( IERR .GT. 0 ) THEN
CCC                  GOTO 9999
CCC               ELSE IF( IERR .LT. 0 ) THEN
CCC                  IERR = 0
CCC                  GOTO 200
CCC               ENDIF
CCC               NBST4  = NBST4 + 1
CCC               NBSTSU = NBSTSU + 1
CCC               GOTO 1000
CCCC
CCC            ENDIF
CCCC
CCCC           SI LA BOULE DE NS CONTIENT 5 TRIANGLES ET A UN SOMMET VOISIN
CCCC           SOMMET DE 5 TRIANGLES ALORS L'ARETE JOIGNANT CES 2 SOMMETS
CCCC           EST TRANSFORMEE EN UN SEUL SOMMET DE 6 TRIANGLES
CCCC           ============================================================
CCC            IF( NBTRCF .EQ. 5 ) THEN
CCCC
CCC               DO 80 I=1,5
CCCC                 LE NUMERO DU SOMMET DE L'ARETE I ET DIFFERENT DE NS
CCC                  NS1 = NOSTBO(I)
CCCC                 LA LISTE DES TRIANGLES DE SOMMET NS1
CCC                  CALL TRP1ST( NS1, NOARST,
CCC     %                         MOSOAR, NOSOAR, MOARTR, MXARTR, NOARTR,
CCC     %                         MXTRCF-5, NBTRC1, NOTRCF(6) )
CCC                  IF( NBTRC1 .EQ. 5 ) THEN
CCCC
CCCC                    L'ARETE DE SOMMETS NS-NS1 DEVIENT UN POINT
CCCC                    PAR SUPPRESSION DU SOMMET NS
CCCC
CCCC                    REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA B
CCC                     NOAR = NOAR0
CCC 70                  IF( NOAR .GT. 0 ) THEN
CCCC                       PROTECTION DU NO DE L'ARETE SUIVANTE
CCC                        NA = NOSOAR(LCHAIN,NOAR)
CCCC                       L'ARETE INTERNE EST REMISE A -1
CCC                        NOSOAR(LCHAIN,NOAR) = -1
CCCC                       L'ARETE SUIVANTE
CCC                        NOAR = NA
CCC                        GOTO 70
CCC                     ENDIF
CCCC
CCCC                    LE POINT NS1 DEVIENT LE MILIEU DE L'ARETE NS-NS1
CCC                     DO 75 J=1,3
CCC                        PXYD(J,NS1) = (PXYD(J,NS) + PXYD(J,NS1)) * 0.5D0
CCC 75                  CONTINUE
CCCC
CCC                     IF( NOFOTI .GT. 0 ) THEN
CCCC                       LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
CCCC                       CALCUL DE PXYZD(3,NS1) DANS LE REPERE INITIAL => XYZ(
CCC                        CALL DVTAIL( NUTYSU, NOFOTI,
CCC     %                               PXYD(1,NS1), PXYD(2,NS1),
CCC     %                               ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
CCC     %                               XYZ,    D, IER )
CCC                        IF( IER .EQ. 0 ) THEN
CCCC                          LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
CCC                           R = XMIN + D
CCC                           PXYD(3,NS1) = MISECH( R, XMIN, XYMXMI )
CCC                        ENDIF
CCC                     ENDIF
CCCC
CCCC                    SUPPRESSION DU POINT NS ET MISE EN DELAUNAY
CCC                     CALL TE1STM( NS,     nbarpi, PXYD,   NOARST,
CCC     %                            MOSOAR, MXSOAR, N1SOAR, NOSOAR,
CCC     %                            MOARTR, MXARTR, N1ARTR, NOARTR,
CCC     %                            MXTRCF, N1ARCF, NOARCF,
CCC     %                            LARMIN, NOTRCF, NOSTBO,
CCC     %                            IERR )
CCC                     IF( IERR .GT. 0 ) THEN
CCC                         GOTO 9999
CCC                     ELSE IF( IERR .LT. 0 ) THEN
CCC                        IERR = 0
CCC                        GOTO 200
CCC                     ENDIF
CCC                     NBSTSU = NBSTSU + 1
CCC                     NBST5  = NBST5 + 1
CCC                     GOTO 1000
CCC                  ENDIF
CCC 80            CONTINUE
CCC            ENDIF
C
C           SI LA TAILLE DE L'ARETE MINIMALE EST <0.5 TAILLE SOUHAITEE
C           ALORS SUPPRESSION DU SOMMET NS
C           ===========================================================
            IF( DMOY .LT. AMPLI2*PXYD(3,NS) ) THEN
C              REMISE A -1 DU CHAINAGE DES ARETES PERIPHERIQUES DE LA BOULE NS
               NOAR = NOAR0
 90            IF( NOAR .GT. 0 ) THEN
C                 PROTECTION DU NO DE L'ARETE SUIVANTE
                  NA = NOSOAR(LCHAIN,NOAR)
C                 L'ARETE INTERNE EST REMISE A -1
                  NOSOAR(LCHAIN,NOAR) = -1
C                 L'ARETE SUIVANTE
                  NOAR = NA
                  GOTO 90
               ENDIF
C
C              SI CE POINT EST DANS UNE BOULE DE L'UN DES POINTS INTERNES
C              IMPOSES (RAYON=3 DIST SOUHAITEE) ALORS IL N'EST PAS DETRUIT
               CALL PTBPTI( NS, NBDPFR, NBDPII, PXYD, NONOUI )
               IF( NONOUI .NE. 0 ) GOTO 200
C
               CALL TE1STM( NS,     nbarpi, PXYD,   NOARST,
     %                      MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                      MOARTR, MXARTR, N1ARTR, NOARTR,
     %                      MXTRCF, N1ARCF, NOARCF,
     %                      LARMIN, NOTRCF, NOSTBO,
     %                      IERR )
               IF( IERR .GT. 0 ) GOTO 9999
               IF( IERR .LT. 0 ) THEN
                  IERR = 0
                  GOTO 200
               ENDIF
               NBSTSU = NBSTSU + 1
               GOTO 1000
C
            ENDIF
C
C           NBTRCF EST COMPRIS ENTRE 5 ET 7 => BARYCENTRAGE SIMPLE
C           ======================================================
C           LES 2 COORDONNEES DU BARYCENTRE DES SOMMETS DES ARETES
C           SIMPLES DE LA BOULE DU SOMMET NS
 200        XBAR = XBAR / NBSTBO
            YBAR = YBAR / NBSTBO
C
C           PONDERATION POUR EVITER LES DEGENERESCENSES
            PXYD(1,NS) = PONDE1 * PXYD(1,NS) + PONDER * XBAR
            PXYD(2,NS) = PONDE1 * PXYD(2,NS) + PONDER * YBAR
C
            IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C              CALCUL DE PXYZD(3,NS) DANS LE REPERE INITIAL => XYZ(1:3)
               CALL DVTAIL( NUTYSU, NOFOTI,
     %                      PXYD(1,NS), PXYD(2,NS),
     %                      ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                      XYZ,    D, IER )
               IF( IER .EQ. 0 ) THEN
C                 LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
                  R = REAL( XMIN + D )
                  PXYD(3,NS) = MISECH( R, XMIN, XYMXMI )
               ENDIF
            ENDIF
C
C           LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
            CALL TEDELA( PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOAR0,
     %                   MOARTR, MXARTR, NOARTR, MODIFS )
C
 1000    CONTINUE
C
         WRITE(IMPRIM,11000) NBSTSU
11000 FORMAT( I6,' SOMMETS SUPPRIMES' )
C
CCC         WRITE(IMPRIM,11000) NBST4, NBST5, NBSTSU
CCC11000 FORMAT( I6,' SOMMETS CENTRES DE 4 TRIANGLES'/
CCC     %        I6,' SOMMETS CENTRES DE 5 TRIANGLES + ST DE 5 TRIANGLES'/
CCC     %        I6,' SOMMETS SUPPRIMES' )
C
C        TRACE DE LA TRIANGULATION ACTUELLE ET CALCUL DE LA QUALITE
         CALL TETRMA( 1+MOD(ITER,14), NCNOIR,
     %                COMXMI, PXYD,
     %                MOSOAR, MXSOAR, NOSOAR,
     %                MOARTR, MXARTR, NOARTR,
     %                NBTRIA, QUAMOY, QUAMIN )
C
C        MISE A JOUR POUR NE PAS OUBLIER LES NOUVEAUX SOMMETS
         IF( NBS1 .GT. NBS2 ) THEN
            NBS1 = NBSOMM
         ELSE
            NBS2 = NBSOMM
         ENDIF
C
 5000 CONTINUE
C
 9999 RETURN
      END
