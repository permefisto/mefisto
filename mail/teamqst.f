      SUBROUTINE TEAMQST( NOFOTI, NUTYSU, RAP2P3,
     %                    NOARST, MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                    MOARTR, MXARTR, N1ARTR, NOARTR,
     %                    MXTRCF, NOTRCF, NOSTBO,
     %                    N1ARCF, NOARCF, LARMIN,
     %                    ORIGIN, D2D3,   XMIN,   YMIN,   XYMXMI,
     %                    NBDPFR, NBDPII, NBSOMM, MXSOMM,
     %                    PXYD,   NSLIGN, QUMIST, NOANCI,
     %                    IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    SUR LES SOMMETS DE PLUS MAUVAISE QUALITE DES TRIANGLES
C ----    SI LA TAILLE DE L'ARETE MOYENNE EST >AMPLI*TAILLE SOUHAITEE
C         ALORS AJOUT D'UN SOMMET BARYCENTRE DU PLUS GRAND TRIANGLE
C               DE SOMMET NS
C         SI LA TAILLE DE L'ARETE MOYENNE EST <AMPLI/2*TAILLE SOUHAITEE
C         ALORS SUPPRESSION DU SOMMET NS
C         SINON LE SOMMET NS DEVIENT LE BARYCENTRE PONDERE DE SES VOISINS
C
C         REMARQUE: AMPLI EST DEFINI DANS $MEFISTO/mail/tehote.f
C         ET DOIT AVOIR LA MEME VALEUR POUR EVITER TROP DE MODIFICATIONS
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
C NSLIGN : TABLEAU DU NUMERO DE SOMMET DANS SA LIGNE POUR CHAQUE
C          SOMMET FRONTALIER
C          NUMERO DU POINT DANS LE LEXIQUE POINT SI INTERNE IMPOSE
C          0 SI LE POINT EST INTERNE NON IMPOSE PAR L'UTILISATEUR
C         -1 SI LE SOMMET EST EXTERNE AU DOMAINE
C ORIGIN : POINT ORIGINE DANS LE PLAN DE MAILLAGE
C D2D3   : MATRICE DE PASSAGE DES COORDONNEES 3D EN COORDONNEES 2D
C XMIN   : ABSCISSE MINIMALE APRES TRANSFORMATION 3D->2D
C YMIN   : ORDONNEE MINIMALE APRES TRANSFORMATION 3D->2D
C XYMXMI : ECART MAXIMAL ENTRE XMAX-XMIN ET YMAX-YMIN
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
C QUMIST : TABLEAU (MXSOMM) QUALITE MINIMALE DES TRANGLES DES SOMMETS
C NOANCI : TABLEAU (MXSOMM) NUMERO ANCIEN APRES TRI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY     MARS 2008
C....................................................................012
      include"./incl/ampli.inc"
      PARAMETER        (LCHAIN=6)
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  RAP2P3, PXYD(3,*)
      DOUBLE PRECISION  SURTD2, XXX, YYY, ZZZ
      DOUBLE PRECISION  PONDER, PONDE1, XBAR, YBAR, X, Y, D, DMOY,
     %                  DMAX, DMIN, DNS
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
     %                  LARMIN(MXTRCF),
     %                  NOANCI(MXSOMM)
      REAL              QUMIST(MXSOMM)
      INTEGER           NOSOTR(3)
C
ccc      LOGICAL           TRATRI
ccc      COMMON / DV2DCO / TRATRI
cccC     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION
ccc      TRATRI = .TRUE.
cccC
cccC     TRACE DE LA TRIANGULATION ACTUELLE ET CALCUL DE LA QUALITE
ccc      CALL TETRMA( NCBLAN, NCNOIR,
ccc     %             COMXMI, PXYD,
ccc     %             MOSOAR, MXSOAR, NOSOAR,
ccc     %             MOARTR, MXARTR, NOARTR,
ccc     %             NBTRIA, QUAMOY, QUAMIN )
C
C     QUALITE MINIMALE DES TRIANGLES DE CHAQUE SOMMET
      DO 5 NS=1,NBSOMM
         QUMIST( NS ) = 2
         NOANCI( NS ) = NS
 5    CONTINUE
C
      DO 8 NT=1,MXARTR
         IF( NOARTR(1,NT) .EQ. 0 ) GOTO 8
C
C        LE NUMERO DES 3 SOMMETS DU TRIANGLE NT
         CALL NUSOTR( NT, MOSOAR, NOSOAR, MOARTR, NOARTR,
     %                NOSOTR )
C
C        LA QUALITE DU TRIANGLE NOSOTR
         CALL QUTR2D( PXYD(1,NOSOTR(1)), PXYD(1,NOSOTR(2)),
     %                PXYD(1,NOSOTR(3)), QUALIT )
C
         DO 6 NS=1,3
            N = NOSOTR(NS)
            QUMIST( N ) = MIN( QUMIST( N ), QUALIT )
 6       CONTINUE
 8    CONTINUE
C
C     TRI CROISSANT DES SOMMETS SELON LEUR QUALITE
      CALL TRITRP( NBSOMM, QUMIST, NOANCI )
C
C     LE NOMBRE D'ITERATIONS POUR AMELIORER LA QUALITE DES MAUVAIS SOMMETS
      NBITAQ = 6
      IER    = 0
      QUAMAX = 0.5
      PASQUA = (1.0-QUAMAX) / NBITAQ
C
      DO 5000 ITER=1,NBITAQ
C
C        LE NOMBRE DE SOMMETS SUPPRIMES
         NBSTSU = 0
         NBBAAJ = 0
C
C        COEFFICIENT DE PONDERATION CROISSANT AVEC LES ITERATIONS
         PONDER = 0.1D0 + ITER * 0.6D0 / NBITAQ
CCC 21 NOVEMBRE 2006 PONDER = 0.1D0 + ITER * 0.5D0 / NBITAQ
CCC 10 OCTOBRE  2006 PONDER = MIN( 1D0, 0.1D0 + ITER * 0.9D0 / NBITAQ )
CCC  9 MARS     2006 PONDER = MIN( 1D0, ( 50 + (50*ITER)/NBITAQ ) * 0.01D0 )
         PONDE1 = 1D0 - PONDER
C
C        QUALITE MAX A ENVISAGER POUR CETTE ITERATION
         QUAMAX = QUAMAX + PASQUA
C
         DO 1000 NST = 1, NBSOMM
C
C           LA QUALITE EST ELLE SUFFISANTE?
            IF( QUMIST(NST) .GT. QUAMAX ) GOTO 1010
C
C           LE NUMERO DU SOMMET DE QUALITE MINIMALE TRIEE
            NS = NOANCI( NST )
C
C           LE SOMMET EST IL INTERNE AU DOMAINE?
            IF( NSLIGN(NS) .NE. 0 ) GOTO 1000
C
ccc            print *,'Sommet visite',NS,' qualite=',QUMIST(NST)
C
C           LE SOMMET EST IL INTERNE A UNE BOULE DE POINT INTERNE ET
C           DE RAYON 5 FOIS LA DISTANCE SOUHAITEE?
            CALL PTBPTI( NS,  NBDPFR, NBDPII, PXYD, NONOUI )
            IF( NONOUI .NE. 0 ) GOTO 1000
C
C           EXISTE-T-IL UNE ARETE DE SOMMET NS ?
            NOAR = NOARST( NS )
            IF( NOAR .LE. 0 ) GOTO 1000
            IF( NOSOAR(1,NOAR) .LE. 0 ) GOTO 1000
C
C           LES TRIANGLES DE L'ARETE NOAR
            IF( NOSOAR( 4, NOAR ) .LE. 0 .OR.
     %          NOSOAR( 5, NOAR ) .LE. 0 ) GOTO 1000
C
C           RECHERCHE DES TRIANGLES DE SOMMET NS
C           ILS DOIVENT FORMER UN CONTOUR FERME DE TYPE ETOILE
            CALL TRP1ST( NS,     NOARST, MOSOAR, NOSOAR,
     %                   MOARTR, MXARTR, NOARTR,
     %                   MXTRCF, NBTRCF, NOTRCF )
            IF( NBTRCF .LE. 0 ) GOTO 1000
C
C           MISE A JOUR DE LA DISTANCE SOUHAITEE
            IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C              CALCUL DE PXYZD(3,NS) DANS LE REPERE INITIAL => XYZ(1:3)
               CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                      PXYD(1,NS), PXYD(2,NS),
     %                      ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                      XYZ,    PXYD(3,NS), IER )
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
            DMAX = 0D0
            DMIN = 1D124
            DNS  = 1D133
            DO 50 I=1,NBSTBO
               NSTB = NOSTBO(I)
               X    = PXYD(1,NSTB)
               Y    = PXYD(2,NSTB)
               XBAR = XBAR + X
               YBAR = YBAR + Y
               D    = SQRT( (X-PXYD(1,NS))**2 + (Y-PXYD(2,NS))**2 )
               DMOY = DMOY + D
               DMAX = MAX( DMAX, D )
               DMIN = MIN( DMIN, D )
               DNS  = MIN( DNS, PXYD(3,NSTB) )
 50         CONTINUE
            XBAR = XBAR / NBSTBO
            YBAR = YBAR / NBSTBO
            DMOY = DMOY / NBSTBO
C
cccC           PAS DE MODIFICATION DE LA TOPOLOGIE LORS DE LA DERNIERE ITERATION
cccC           =================================================================
ccc            IF( ITER .EQ. NBITAQ ) GOTO 200
C
C           SI LA TAILLE DE L'ARETE MAXIMALE EST >AMPLI*MIN des TAILLES SOUHAITE
C           ALORS AJOUT D'UN SOMMET BARYCENTRE DU PLUS GRAND TRIANGLE DE SOMMET
C           ====================================================================
CCC19/10/2006            IF( DMOY .GT. AMPLI*PXYD(3,NS) ) THEN
            IF( DMAX .GT. AMPLI*DNS  ) THEN
C
               IMAX = 0
               DMOY = 0D0
               DO 150 I=1,NBTRCF
C                 RECHERCHE DU PLUS GRAND TRIANGLE EN SURFACE
                  CALL NUSOTR( NOTRCF(I), MOSOAR, NOSOAR,
     %                         MOARTR, NOARTR, NOSOTR )
                  D  = SURTD2( PXYD(1,NOSOTR(1)),
     %                         PXYD(1,NOSOTR(2)),
     %                         PXYD(1,NOSOTR(3)) )
                  IF( D .GT. DMOY ) THEN
                     DMOY = D
                     IMAX = I
                  ENDIF
 150           CONTINUE
C
C              AJOUT DU BARYCENTRE DU TRIANGLE NOTRCF(IMAX)
               NT = NOTRCF( IMAX )
               CALL NUSOTR( NT, MOSOAR, NOSOAR,
     %                      MOARTR, NOARTR, NOSOTR )
               IF( NBSOMM .GE. MXSOMM ) THEN
                  WRITE(IMPRIM,*) 'TEAMQST: SATURATION DU TABLEAU PXYD'
C                 ABANDON DE L'AMELIORATION DU SOMMET NS
                  IERR = 53
                  GOTO 9999
               ENDIF
               NBSOMM = NBSOMM + 1
               DO 160 I=1,3
                  PXYD(I,NBSOMM) = ( PXYD(I,NOSOTR(1))
     %                             + PXYD(I,NOSOTR(2))
     %                             + PXYD(I,NOSOTR(3)) ) / 3D0
 160           CONTINUE
C
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                 CALCUL DE PXYZD(3,NBSOMM) DANS LE REPERE INITIAL => XYZ(1:3)
                  CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                         PXYD(1,NBSOMM), PXYD(2,NBSOMM),
     %                         ORIGIN, D2D3,   XMIN, YMIN, XYMXMI,
     %                         XYZ,    PXYD(3,NBSOMM), IER )
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
C           SI LA TAILLE DE L'ARETE MAXIMALE EST <AMPLI/2*TAILLE SOUHAITEE
C           ALORS SUPPRESSION DU SOMMET NS
C           ==============================================================
CCC19/10/2006            IF( DMOY .LT. AMPLI2*PXYD(3,NS) ) THEN
            IF( DMAX .LT. AMPLI2*DNS ) THEN
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
               CALL TE1STM( NS,     NBDPII, PXYD,   NOARST,
     %                      MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                      MOARTR, MXARTR, N1ARTR, NOARTR,
     %                      MXTRCF, N1ARCF, NOARCF,
     %                      LARMIN, NOTRCF, NOSTBO,
     %                      IERR )
               if( ierr .eq. -543 ) then
                  ierr = 0
                  goto 1000
               else if( ierr .lt.    0 ) then
C                 LE SOMMET NS EST EXTERNE DONC NON SUPPRIME
C                 OU BIEN LE SOMMET NS EST LE CENTRE D'UN CF DONT TOUTES
C                 LES ARETES SIMPLES SONT FRONTALIERES
C                 DANS LES 2 CAS LE SOMMET NS N'EST PAS SUPPRIME
                  IERR = 0
                  GOTO 200
               ELSE IF( IERR .GT. 0 ) THEN
C                 ERREUR IRRECUPERABLE
                  GOTO 9999
               ENDIF
               NBSTSU = NBSTSU + 1
               GOTO 1000
C
            ENDIF
C
C           LES 2 COORDONNEES DU BARYCENTRE DES SOMMETS DES ARETES
C           SIMPLES DE LA BOULE DU SOMMET NS
C           ======================================================
C DEBUT AJOUT 10 OCTOBRE 2006
C           PONDERATION POUR EVITER LES DEGENERESCENSES AVEC PROTECTION
C           SI UN TRIANGLE DE SOMMET NS A UNE AIRE NEGATIVE APRES BARYCENTRAGE
C           ALORS LE SOMMET NS N'EST PAS BOUGE
C
C           PROTECTION DES XY DU POINT INITIAL
 200        XXX = PXYD(1,NS)
            YYY = PXYD(2,NS)
            ZZZ = PXYD(3,NS)
C
C           PONDERATION POUR EVITER LES DEGENERESCENSES
            PXYD(1,NS) = PONDE1 * PXYD(1,NS) + PONDER * XBAR
            PXYD(2,NS) = PONDE1 * PXYD(2,NS) + PONDER * YBAR
C
            IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C              CALCUL DE PXYZD(3,NS) DANS LE REPERE INITIAL => XYZ(1:3)
               CALL TETAID( NUTYSU, NOFOTI, RAP2P3,
     %                      PXYD(1,NS), PXYD(2,NS),
     %                      ORIGIN, D2D3, XMIN, YMIN, XYMXMI,
     %                      XYZ,    PXYD(3,NS), IER )
            ENDIF
C
CCC         WRITE(IMPRIM,*)'TEAMQST: NS=',NS,' ANCIEN =',XXX,YYY
CCC         WRITE(IMPRIM,*)'TEAMQST: NS=',NS,' NOUVEAU=',PXYD(1,NS),PXYD(2,NS)
C
            NS2 = 0
            NS3 = 0
            DO 240 I=1,NBTRCF
C              LE NUMERO DE L'ARETE DU TRIANGLE NT NE CONTENANT PAS LE SOMMET NS
               NT = NOTRCF(I)
               DO 220 NA=1,3
C                 LE NUMERO DE L'ARETE NA DANS LE TABLEAU NOSOAR
                  NOAR = ABS( NOARTR(NA,NT) )
                  IF( NOSOAR(1,NOAR) .NE. NS   .AND.
     %                NOSOAR(2,NOAR) .NE. NS ) THEN
                     IF( NOARTR(NA,NT) .GE. 0 ) THEN
                        NS2 = NOSOAR(1,NOAR)
                        NS3 = NOSOAR(2,NOAR)
                     ELSE
                        NS3 = NOSOAR(1,NOAR)
                        NS2 = NOSOAR(2,NOAR)
                     ENDIF
                     GOTO 225
                  ENDIF
 220           CONTINUE
C              AIRE SIGNEE DU TRIANGLE NT
 225           D = SURTD2( PXYD(1,NS), PXYD(1,NS2), PXYD(1,NS3) )
               IF( D .LE. 0D0 ) THEN
ccc                  WRITE(IMPRIM,*) 'TEAMQST: BARYCENTRAGE POINT NS=',NS,
ccc     %            '   XB=',PXYD(1,NS),' YB=',PXYD(2,NS),
ccc     %            ' => NOUVEAU TRIANGLE AVEC AIRE<0'
ccc                  WRITE(IMPRIM,*) 'TEAMQST: => PT REMIS EN X =',XXX,
ccc     %            ' Y =',YYY,' LONGUEUR ARETE=',ZZZ
                  PXYD(1,NS) = XXX
                  PXYD(2,NS) = YYY
                  PXYD(3,NS) = ZZZ
C                 LA PONDERATION EST REDUITE  10 OCTOBRE 2006
                  PONDER = MAX( 0.1D0, PONDER*0.5D0 )
                  PONDE1 = 1D0 - PONDER
ccc               WRITE(IMPRIM,*)'PONDERATION DU BARYCENTRE PONDER=',PONDER
                  GOTO 1000
               ENDIF
 240        CONTINUE
C
C FIN AJOUT 10 OCTOBRE 2006
C
C           LES ARETES CHAINEES DE LA BOULE SONT RENDUES DELAUNAY
 900        CALL TEDELA( PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOAR0,
     %                   MOARTR, MXARTR, NOARTR, MODIFS )
C
 1000    CONTINUE
C
 1010    IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,11010) ITER, NBSTSU, NBBAAJ, NST
         ELSE
            WRITE(IMPRIM,21010) ITER, NBSTSU, NBBAAJ, NST
         ENDIF
11010 FORMAT(' TEAMQST: ITERATION',I2,':',I7,' SOMMETS SUPPRIMES ',
     %        I7,' BARYCENTRES AJOUTES',I7,' SOMMETS VISITES' )
21010 FORMAT(' TEAMQST: ITERATION',I2,':',I7,' DELETED VERTICES',
     %        I7,' ADDED BARYCENTRES',I7,' VISITED VERTICES' )
C
cccC        TRACE DE LA TRIANGULATION ACTUELLE ET CALCUL DE LA QUALITE
ccc         CALL TETRMA( 1+MOD(ITER,14), NCNOIR,
ccc     %                COMXMI, PXYD,
ccc     %                MOSOAR, MXSOAR, NOSOAR,
ccc     %                MOARTR, MXARTR, NOARTR,
ccc     %                NBTRIA, QUAMOY, QUAMIN )
C
 5000 CONTINUE
C
ccc 9999 TRATRI = .FALSE.

 9999 RETURN
      END
