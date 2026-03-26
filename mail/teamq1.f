      SUBROUTINE TEAMQ1( NOFOTI, NUTYSU,
     %                   NOARST, MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR,
     %                   MXTRCF, NOTRCF, NOSTBO,
     %                   ORIGIN, D2D3,   XMIN,   YMIN,  XYMXMI,
     %                   COMXMI, NBARPI, NBSOMM, MXSOMM, PXYD, NSLIGN,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LES ARETES DECROISSENT ET LA TOPOLOGIE DIMINUE DES POINTS EN +
C -----    BOULE>7TRIANGLES + DELAUNAY + BARYCENTRAGE DES 3 SOMMETS
C          DMAX>TAILLE_IDEALE => BARYCENTRE DU PLUS GRAND TRIANGLE
C                             + DELAUNAY + BARYCENTRE DU PLUS GRAND TRIANGLE
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
C NBARPI : NUMERO DU DERNIER SOMMET FRONTALIER OU INTERNE IMPOSE
C NSLIGN : TABLEAU DU NUMERO DE SOMMET DANS SA LIGNE POUR CHAQUE
C          SOMMET FRONTALIER
C          NUMERO DU POINT DANS LE LEXIQUE POINT SI INTERNE IMPOSE
C          0 SI LE POINT EST INTERNE NON IMPOSE PAR L'UTILISATEUR
C         -1 SI LE SOMMET EST EXTERNE AU DOMAINE
C COMXMI : MIN ET MAX DES COORDONNEEES DES SOMMETS DU MAILLAGE
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
C LARMIN : TABLEAU ( MXTRCF ) AUXILIAIRE D'ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1997
C....................................................................012
      PARAMETER        (LCHAIN=6)
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PXYD(3,*)
      DOUBLE PRECISION  PONDER, PONDE1, XBAR, YBAR, X, Y, D, DMAX
      DOUBLE PRECISION  D2D3(3,3), SURTD2
      REAL              ORIGIN(3), XYZ(3)
      INTEGER           NOARTR(MOARTR,*),
     %                  NOSOAR(MOSOAR,*),
     %                  NOARST(*),
     %                  NOTRCF(MXTRCF),
     %                  NSLIGN(*),
     %                  NOSTBO(*)
      INTEGER           NOSOTR(3,2)
      REAL              COMXMI(3,2)
C
C     LE NOMBRE D'ITERATIONS POUR AMELIORER LA QUALITE
      NBITAQ = 5
      IER    = 0
C
C     INITIALISATION DU PARCOURS
      NBS1 = NBSOMM
      NBS2 = NBARPI + 1
      NBS3 = -1
C
      DO 5000 ITER=1,NBITAQ
C
C        LE NOMBRE DE BARYCENTRES AJOUTES
         NBBAAJ = 0
C
C        LES COMPTEURS DE PASSAGE SUR LES DIFFERENTS CAS
         NBST8 = 0
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
C           CALCUL DE L'ARETE DE TAILLE MAXIMALE ET MINIMALE ISSUE DE NS
C           ---------------------------------------------------------------
            XBAR = 0D0
            YBAR = 0D0
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
 50         CONTINUE
C
C           REGULARISATION PAR BARYCENTRAGE SEULE POUR LES 2 DERNIERES ITERATION
            IF( ITER .GT. NBITAQ-2 ) GOTO 200
C
C           SI LA BOULE DE NS CONTIENT AU MOINS 8 TRIANGLES
C           ALORS UN TRIANGLE INTERNE EST AJOUTE + 3 TRIANGLES (1 PAR ARETE)
C           ================================================================
            IF( NBTRCF .GE. 8 ) THEN
C
C              MODIFICATION DES COORDONNEES DU SOMMET NS
C              IL DEVIENT LE BARYCENTRE DU TRIANGLE NOTRCF(1)
               CALL NUSOTR( NOTRCF(1), MOSOAR, NOSOAR,
     %                      MOARTR, NOARTR, NOSOTR )
               DO 110 I=1,3
                  PXYD(I,NS) = ( PXYD(I,NOSOTR(1,1))
     %                         + PXYD(I,NOSOTR(2,1))
     %                         + PXYD(I,NOSOTR(3,1)) ) / 3D0
 110           CONTINUE
C
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                 CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(1:3)
                  CALL DVTAIL( NUTYSU, NOFOTI,
     %                         PXYD(1,NS), PXYD(2,NS),
     %                         ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                         XYZ,    D, IER )
                  IF( IER .EQ. 0 ) THEN
C                    LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
                     R = REAL( XMIN + D )
                     PXYD(3,NS) = MISECH( R, XMIN, XYMXMI )
                  ENDIF
               ENDIF
C
C              AJOUT DES 2 AUTRES SOMMETS COMME BARYCENTRES DES TRIANGLES
C              NOTRCF(1+NBTRCF/3) ET NOTRCF(1+2*NBTRCF/3)
               NBT1 = ( NBTRCF + 1 ) / 3
               DO 140 N=1,2
C
C                 LE TRIANGLE TRAITE
                  NT = NOTRCF(1 + N * NBT1 )
C
C                 LE NUMERO PXYD DE SES 3 SOMMETS
                  CALL NUSOTR( NT, MOSOAR, NOSOAR,
     %                         MOARTR, NOARTR, NOSOTR )
C
C                 AJOUT DU NOUVEAU BARYCENTRE
                  IF( NBSOMM .GE. MXSOMM ) THEN
                     WRITE(IMPRIM,*)'TEAMQ1: SATURATION DU TABLEAU PXYD'
C                    ABANDON DE L'AMELIORATION DU SOMMET NS
                     GOTO 1000
                  ENDIF
                  NBSOMM = NBSOMM + 1
                  DO 120 I=1,3
                     PXYD(I,NBSOMM) = ( PXYD(I,NOSOTR(1,1))
     %                                + PXYD(I,NOSOTR(2,1))
     %                                + PXYD(I,NOSOTR(3,1)) ) / 3D0
 120              CONTINUE
C
                  IF( NOFOTI .GT. 0 ) THEN
C                    LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                    CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(1:3
                     CALL DVTAIL( NUTYSU, NOFOTI,
     %                            PXYD(1,NBSOMM), PXYD(2,NBSOMM),
     %                            ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                            XYZ,    D, IER )
                     IF( IER .EQ. 0 ) THEN
C                       LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
                        R = REAL( XMIN + D )
                        PXYD(3,NBSOMM) = MISECH( R, XMIN, XYMXMI )
                     ENDIF
                  ENDIF
C
C                 SOMMET INTERNE A LA TRIANGULATION
                  NSLIGN(NBSOMM) = 0
C
C                 LES 3 ARETES DU TRIANGLE NT SONT A RENDRE DELAUNAY
                  DO 130 I=1,3
                     NOAR = ABS( NOARTR(I,NT) )
                     IF( NOSOAR(3,NOAR) .EQ. 0 ) THEN
C                       ARETE NON FRONTALIERE
                        IF( NOSOAR(LCHAIN,NOAR) .LT. 0 ) THEN
C                          ARETE NON ENCORE CHAINEE
                           NOSOAR(LCHAIN,NOAR) = NOAR0
                           NOAR0 = NOAR
                        ENDIF
                     ENDIF
 130              CONTINUE
C
C                 TRIANGULATION DU TRIANGLE DE BARYCENTRE NBSOMM
C                 PROTECTION A NE PAS MODIFIER SINON ERREUR!
                  CALL TR3STR( NBSOMM, NT,
     %                         MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                         MOARTR, MXARTR, N1ARTR, NOARTR,
     %                         NOARST,
     %                         NOSOTR, IERR )
                  IF( IERR .NE. 0 ) GOTO 9999
 140           CONTINUE
C
               NBST8  = NBST8 + 1
C
C              LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
               GOTO 900
C
            ENDIF
C
C           SI L'ARETE MAXIMALE EST SUPERIEURS A LA TAILLE SOUHAITEE
C           ALORS LE BARYCENTRE DU TRIANGLE DE PLUS GRANDE SURFACE EST AJOUTE
C           =================================================================
            IF( DMAX .GT. PXYD(3,NS)**2 ) THEN
C
               IMAX = 0
               DMAX = 0D0
               DO 150 I=1,NBTRCF
C                 RECHERCHE DU PLUS GRAND TRIANGLE EN SURFACE
                  CALL NUSOTR( NOTRCF(I), MOSOAR, NOSOAR,
     %                         MOARTR, NOARTR, NOSOTR )
                  D  = SURTD2( PXYD(1,NOSOTR(1,1)),
     %                         PXYD(1,NOSOTR(2,1)),
     %                         PXYD(1,NOSOTR(3,1)) )
                  IF( D .GT. DMAX ) THEN
                     DMAX = D
                     IMAX = I
                  ENDIF
 150           CONTINUE
C
C              AJOUT DU BARYCENTRE DU TRIANGLE NOTRCF(IMAX)
               NT = NOTRCF( IMAX )
               CALL NUSOTR( NT, MOSOAR, NOSOAR,
     %                      MOARTR, NOARTR, NOSOTR )
               IF( NBSOMM .GE. MXSOMM ) THEN
                  WRITE(IMPRIM,*)'TEAMQ1: SATURATION DU TABLEAU PXYD'
C                 ABANDON DE L'AMELIORATION DU SOMMET NS
                  GOTO 1000
               ENDIF
               NBSOMM = NBSOMM + 1
               DO 160 I=1,3
                  PXYD(I,NBSOMM) = ( PXYD(I,NOSOTR(1,1))
     %                             + PXYD(I,NOSOTR(2,1))
     %                             + PXYD(I,NOSOTR(3,1)) ) / 3D0
 160           CONTINUE
C
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                 CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(1:3)
                  CALL DVTAIL( NUTYSU, NOFOTI,
     %                         PXYD(1,NBSOMM), PXYD(2,NBSOMM),
     %                         ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                         XYZ,    D, IER )
                  IF( IER .EQ. 0 ) THEN
C                    LA TAILLE IDEALE EST REMISE A L'ECHELLE 2**20
                     R = REAL( XMIN + D )
                     PXYD(3,NBSOMM) = MISECH( R, XMIN, XYMXMI )
                  ENDIF
               ENDIF
C
C              SOMMET INTERNE A LA TRIANGULATION
               NSLIGN(NBSOMM) = 0
C
C              LES 3 ARETES DU TRIANGLE NT SONT A RENDRE DELAUNAY
               DO 170 I=1,3
                  NOAR = ABS( NOARTR(I,NT) )
                  IF( NOSOAR(3,NOAR) .EQ. 0 ) THEN
C                    ARETE NON FRONTALIERE
                     IF( NOSOAR(LCHAIN,NOAR) .LT. 0 ) THEN
C                       ARETE NON ENCORE CHAINEE
                        NOSOAR(LCHAIN,NOAR) = NOAR0
                        NOAR0 = NOAR
                     ENDIF
                  ENDIF
 170           CONTINUE
C
C              TRIANGULATION DU TRIANGLE DE BARYCENTRE NBSOMM
C              PROTECTION A NE PAS MODIFIER SINON ERREUR!
               CALL TR3STR( NBSOMM, NT,
     %                      MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                      MOARTR, MXARTR, N1ARTR, NOARTR,
     %                      NOARST,
     %                      NOSOTR, IERR )
               IF( IERR .NE. 0 ) GOTO 9999
C
C              UN BARYCENTRE AJOUTE DE PLUS
               NBBAAJ = NBBAAJ + 1
C
C              LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
               GOTO 900
C
            ENDIF
C
C           NBTRCF EST < 8 => BARYCENTRAGE SIMPLE
C           =====================================
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
 900        CALL TEDELA( PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOAR0,
     %                   MOARTR, MXARTR, NOARTR, MODIFS )
C
 1000    CONTINUE
C
         WRITE(IMPRIM,11000) NBST8,NBBAAJ
11000 FORMAT( I6,' TRIANGLES AJOUTES DANS >7 TRIANGLES'/
     %        I6,' BARYCENTRES AJOUTES')
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
