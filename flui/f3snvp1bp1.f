      SUBROUTINE F3SNVP1BP1( NEF00,  NBEF,     NUSOTE,
     %                       MOFACE, MXFACE,   LFACES,
     %                       DT,     Coeftr,   XYZSOM, NDDLNO,
     %                       NTDLHB, VXYZPNtn, VITMXtn,
     %                       BE,     IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU VECTEUR TRANSPORT NAVIER STOKES DU TETRAEDRE BREZZI-FORTIN
C ----- INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 1 + BULLE P4 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TETRAEDRE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C       INTEGRATION EN 15 POINTS DU TETRAEDRE P5-EXACTE
C       Coeftr Integrale tP1B ui(tn,X(tn;tn+1,x)) dx

C ENTREES:	
C --------
C NEF00  : NUMERO DE L'EF A TRAITER
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NUSOTE : NUMERO DES 4 SOMMETS DES TETRAEDRES

C MOFACE : NOMBRE DE MOTS DE CHAQUE FACE DU TABLEAU LFACES
C MXFACE : NOMBRE MAXIMAL DEFACES DU TABLEAU LFACES
C LFACES : TABLEAU DES FACES DU MAILLAGE  (cf hachag.f)
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C                       1 CUBE EST UN TETRA ou PENTA ou HEXAEDRE
C          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                       0 SI PAS DE 1-ER  CUBE
C          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                       ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
C          LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE
C                       CETTE INFORMATION EXISTE SEULEMENT POUR UN OBJET
C                       PAS POUR UNE SURFACE OU UN VOLUME!

C DT     : LE PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C CoefTr : COEFFICIENT DE LA VITESSE TRANSPORTEE
C XYZSOM : 3 COORDONNEES DES SOMMETS DES TETRAEDRES
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU
C NTDLHB : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET HORS BARYCENTRES SI BF
C          (NDIM+1) NBNOHB (SANS BARYCENTRES)

C VXYZPNtn: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLHB) au TEMPS tn
C          SUIVI DES 3 VITESSES AU BARYCENTRE DES TETRAEDRES
C VITMXtn: NORME DE LA VITESSE MAXIMALE AU TEMPS tn

C SORTIES:
C --------
C BE     : VECTEUR ELEMENTAIRE DE 3x5 COEFFICIENTS + 4 PRESSIONS=0 ICI
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C          =4 NOMBRE DE TETRAEDRES PARCOURUS TROP GRAND
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Mars 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Mai 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C MODIFS : ALAIN PERRONNET St Pierre du Perray              Octobre 2020
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)

      REAL               DT, XYZSOM(3,*)
      INTEGER            MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                   NEF00, NBEF, NUSOTE(NBEF,4), NTDLHB,
     %                   NDDLNO(0:*), IERR
      DOUBLE PRECISION   VXYZPNtn(1:*), VITMXtn
      DOUBLE PRECISION   BE(19)

      INTEGER            I, K, L, ITERT, IK, NDL, NS,
     %                   NEF, NEF0, NBSSPATE, LINTER, NOFAMX

      DOUBLE PRECISION   DBLE, ABS, DETM33
      DOUBLE PRECISION   Coeftr, RDPOIDS, D, PRCBL1(5), PRCBL0(5),
     %                   VITMXNEF, V(3), VA(3), Xe(3,4), DELTAe1,DELTAe,
     %                   CB1BL, CB2BL, CB3BL, CB4BL, BULLE, X1, Y1, Z1,
     %                   DDT, CBTR(3), PTITR(3)
ccc     DOUBLE PRECISION VNORM

      DOUBLE PRECISION   XEBL0, YEBL0, ZEBL0, XEBL1, YEBL1, ZEBL1,
     %                   XYZ0(3), XYZ1(3)
      EQUIVALENCE      (XYZ0(1),XEBL0), (XYZ0(2),YEBL0), (XYZ0(3),ZEBL0)
      EQUIVALENCE      (XYZ1(1),XEBL1), (XYZ1(2),YEBL1), (XYZ1(3),ZEBL1)

C     LE NUMERO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER            NOSOFATE(3,4)
      DATA               NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      INTRINSIC          NINT, SQRT

C     POIDS ET POINTS D'INTEGRATION SUR LE TETRAEDRE UNITE P5 EXACTE
C     DOUBLE PRECISION   POTEP5( NBPTI=15 ), COTEP5( 3, NBPTI )
      include"./incl/tep5poco.inc"

C     INITIALISATION DU VECTEUR ELEMENTAIRE
      IERR = 0
      DO I=1,19
         BE(I) = 0D0
      ENDDO

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP1BP13D( DT, XYZSOM, NEF00, NBEF, NUSOTE, NDDLNO,
     %                 VITMXtn, VXYZPNtn,   VITMXNEF, NBSSPATE, DDT )

C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF
      DO I=1,4
         NS = NUSOTE(NEF00,I)
         DO K=1,3
            Xe(K,I) = DBLE( XYZSOM( K, NS ) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF INITIAL
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      Z1 = Xe(3,1)
      DELTAe1 = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                  Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                  Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )

      IF( DELTAe1 .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'f3snvp1bp1: ATTENTION EF',NEF00,
     %                     ' de VOLUME*6=',DELTAe1,' NON PRIS EN COMPTE'
         ELSE
            WRITE(IMPRIM,*) 'f3snvp1bp1: ATTENTION FE',NEF00,
     %                      ' of VOLUME*6=',DELTAe1,' is NOT COMPUTED'
         ENDIF
         GOTO 9999
      ENDIF

C     CALCUL DU VECTEUR tPchapeau(bl) ( Pchapeau(X(tn;tn+1,Fe(bl)) Vitesse(k) )
C     =========================================================================
C     Boucle sur les points d'integration dans le tetraedre NEF00
      DO 100 L=1,NBPTI

C        INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
         IERR = 0
         NEF0 = 0
         NEF  = NEF00

C        RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF
         DO I=1,4
            NS = NUSOTE(NEF,I)
            DO K=1,3
               Xe(K,I) = DBLE( XYZSOM( K, NS ) )
            ENDDO
         ENDDO

C        CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
         X1 = Xe(1,1)
         Y1 = Xe(2,1)
         Z1 = Xe(3,1)
         DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                    Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                    Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )
         IF( DELTAe .LE. 0D0 ) THEN
C           JACOBIEN de l'EF NEF NUL ou NEGATIF
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'f3snvp1bp1: attention EF',NEF,
     %      ' de VOLUME*6=',DELTAe,' Be NON CALCULE pour L=',L
            ELSE
               WRITE(IMPRIM,*) 'f3snvp1bp1: attention FE',NEF,
     %      ' of VOLUME*6=',DELTAe,' Be is NOT COMPUTED for L=',L
            ENDIF
            GOTO 100
         ENDIF

C        COORDONNEES BARYCENTRIQUES DU POINT bl D'INTEGRATION DANS EF REFERENCE
C        LES 3 DERNIERES COORDONNEES BARYCENTRIQUES SONT AUSSI X Y Z
C        DANS L'ELEMENT DE REFERENCE
         CB2BL = COTEP5( 1, L )
         CB3BL = COTEP5( 2, L )
         CB4BL = COTEP5( 3, L )
         CB1BL = 1.D0 - CB2BL - CB3BL - CB4BL

C        COORDONNEES (XEBL1,YEBL1,ZEBL1)=Fe(bl) SUR LE TETRAEDRE NEF a tn+1
C       (XEBL1,YEBL1,ZEBL1)=Fe(CB2BL,CB3BL,CB4BL) SUR LE TETRAEDRE NEF
         XEBL1=CB1BL*Xe(1,1)+CB2BL*Xe(1,2)+CB3BL*Xe(1,3)+CB4BL*Xe(1,4)
         YEBL1=CB1BL*Xe(2,1)+CB2BL*Xe(2,2)+CB3BL*Xe(2,3)+CB4BL*Xe(2,4)
         ZEBL1=CB1BL*Xe(3,1)+CB2BL*Xe(3,2)+CB3BL*Xe(3,3)+CB4BL*Xe(3,4)

C        CALCUL DES Pchapeau(bl) POLYNOMES LAGRANGE DE DEGRE 1 + BULLE P4
         BULLE = 64D0 * CB1BL * CB2BL * CB3BL * CB4BL
         PRCBL1(1) = CB1BL - BULLE
         PRCBL1(2) = CB2BL - BULLE
         PRCBL1(3) = CB3BL - BULLE
         PRCBL1(4) = CB4BL - BULLE
         PRCBL1(5) = 4D0 * BULLE

C        CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL,CB4BL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         V(3) = 0D0
         DO I=1,4
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUSOTE(NEF,I) - 1 )
            V(1) = V(1) + PRCBL1(I) * VXYZPNtn( NDL+1 )
            V(2) = V(2) + PRCBL1(I) * VXYZPNtn( NDL+2 )
            V(3) = V(3) + PRCBL1(I) * VXYZPNtn( NDL+3 )
         ENDDO
C        LE NO DU DERNIER DL VITESSE DU BARYCENTRE
         NDL  = NTDLHB + 3 * NEF
         V(1) = V(1) + PRCBL1(5) * VXYZPNtn( NDL-2 )
         V(2) = V(2) + PRCBL1(5) * VXYZPNtn( NDL-1 )
         V(3) = V(3) + PRCBL1(5) * VXYZPNtn( NDL   )

         IF( NBSSPATE .LE. 0 ) GOTO 90
C        SI VITESSE QUASI NULLE PAS DE REMONTEE DE LA CARACTERISTIQUE

C        -----------------------------------------------------------------
C        CALCUL DU POINT X(tn;tn+1,Fe(bl))=(XEBL0,YEBL0,ZEBL0)  c-a-d
C        du POINT a L'INSTANT tn QUI SERA AU POINT D'INTEGRATION Fe(Bl)
C        a L'INSTANT tn+1 par integration retrograde de la caracteristique
C        -----------------------------------------------------------------
C        COORDONNEES du POINT INITIAL sur l'EF NEF pour les ITERATIONS en temps
         XEBL0 = XEBL1
         YEBL0 = YEBL1
         ZEBL0 = ZEBL1

C        NBSSPATE SOUS PAS DE TEMPS DDT ENTRE tn+1-DT=tn et tn+1
C        pour calculer X(tn;tn+1,(XEBL1,YEBL1,ZEBL1))

         DO 10 ITERT = 1, NBSSPATE
C
C           SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
            VA(1) = V(1)
            VA(2) = V(2)
            VA(3) = V(3)

C           CALCUL DES Pchapeau(XCB234BL)   POLYNOMES de BREZZI-FORTIN P1+BULLE
            BULLE = 64D0 * CB1BL * CB2BL * CB3BL * CB4BL
            PRCBL0(1) = CB1BL - BULLE
            PRCBL0(2) = CB2BL - BULLE
            PRCBL0(3) = CB3BL - BULLE
            PRCBL0(4) = CB4BL - BULLE
            PRCBL0(5) = 4D0 * BULLE
C
C           CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL)) dans l'EF NEF
            V(1) = 0D0
            V(2) = 0D0
            V(3) = 0D0
            DO I=1,4
C              NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C              AU POINT X(tn;tn+1,Fe(bl))
               NDL  = NDDLNO( NUSOTE(NEF,I) - 1 )
               V(1) = V(1) + PRCBL0(I) * VXYZPNtn( NDL+1 )
               V(2) = V(2) + PRCBL0(I) * VXYZPNtn( NDL+2 )
               V(3) = V(3) + PRCBL0(I) * VXYZPNtn( NDL+3 )
            ENDDO
C           LE NO DU DERNIER DL VITESSE DU BARYCENTRE
            NDL  = NTDLHB + 3 * NEF
            V(1) = V(1) + PRCBL0(5) * VXYZPNtn( NDL-2 )
            V(2) = V(2) + PRCBL0(5) * VXYZPNtn( NDL-1 )
            V(3) = V(3) + PRCBL0(5) * VXYZPNtn( NDL   )


C           LE POINT (XEBL0,YEBL0,ZEBL0) SERA TRANSPORTE en (XEBL1,YEBL1,ZEBL1)
C           SUR L'EF COURANT NEF00 de l'instant tn a tn+1 par pas ddt
C           -------------------------------------------------------------------
            XEBL0 = XEBL0 - DDT * V(1)
            YEBL0 = YEBL0 - DDT * V(2)
            ZEBL0 = ZEBL0 - DDT * V(3)

C           SAUVEGARDE DU NO DU DERNIER EF PARCOURU
            NEF0 = NEF

C           DANS QUEL ELEMENT FINI NEF du MAILLAGE EST CE POINT XYZEBL0?
            CALL XYZDSTE( NBEF,   4,      NUSOTE,
     %                    MOFACE, MXFACE, LFACES,
     %                    XYZSOM, NEF0,   XEBL0, YEBL0, ZEBL0,
     %                    NEF,    NOFAMX, IERR )
C NEF    : SI IERR=0 NUMERO DU TETRAEDRE CONTENANT LE POINT XEBL0, YEBL0, ZEBL0
C          SI IERR>0 NUMERO DU DERNIER TETRAEDRE PARCOURU
C NOFAMX : NUMERO (1 a 4) DE LA FACE EN CAS DE POINT EXTERIEUR AU MAILLAGE
C IERR   : 0 PAS D'ERREUR LE POINT XEBL0, YEBL0, ZEBL0 EST DANS LE TETRAEDRE NEF
C          1 POINT EXTERIEUR AU MAILLAGE
C            PAS de TETRAEDRE derriere LA FACE NOFAMX de
C            NEF LE DERNIER TETRAEDRE DU PARCOURS POUR ATTEINDRE XYZEBL0
C          2 POINT PROCHE D'UNE FACE NON RETROUVEE DANS LE MAILLAGE
C          3 RENCONTRE D'UNE FACE DANS 3 TETRAEDRES CE QUI EST INCORRECT
C          4 NOMBRE DE TETRAEDRES PARCOURUS TROP GRAND pour DETECTION de BOUCLE

            IF( IERR .GE. 2 ) THEN

C              VITESSE AU NOEUD NONOEU AVANT SORTIE DES TETRAEDRES
               PRINT*,'f3snvp1bp1: POINT:', XEBL0, YEBL0, ZEBL0

               IF( IERR .EQ. 4 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT 10000, V, VA
                  ELSE
                     PRINT 20000, V, VA
                  ENDIF
10000 FORMAT('f3snvp1bp1:CAUSE sortie des TETRAEDRES: VITESSE CALCULEE='
     %             ,3G15.6,' REMPLACEE par VITESSE PRECEDENTE=', 3G15.6)
20000 FORMAT('f3snvp1bp1: Before TETRAHEDRA EXIT: COMPUTED VELOCITY='
     %                ,3G15.6,' REPLACED by PREVIOUS VELOCITY=', 3G15.6)
               ENDIF

C              RETOUR A LA VITESSE ANTERIEURE
               V(1) = VA(1)
               V(2) = VA(2)
               V(3) = VA(3)
               IERR = 0
               GOTO 80

            ENDIF

            IF( NEF0 .NE. NEF ) THEN

C              RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF
               DO I=1,4
                  NS = NUSOTE( NEF, I )
                  DO K=1,3
                     Xe(K,I) = DBLE( XYZSOM( K, NS ) )
                  ENDDO
               ENDDO

C              CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
               X1 = Xe(1,1)
               Y1 = Xe(2,1)
               Z1 = Xe(3,1)
               DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                          Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                          Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )
               IF( DELTAe .LE. 0D0 ) THEN
C                 JACOBIEN de l'EF NEF NUL ou NEGATIF
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,*) 'f3snvp1bp1: Probleme EF',NEF,
     %      ' de VOLUME*6=',DELTAe,' Be NON CALCULE pour L=',L
                  ELSE
                     WRITE(IMPRIM,*) 'f3snvp1bp1: Problem FE',NEF,
     %      ' of VOLUME*6=',DELTAe,' Be is NOT COMPUTED for L=',L
                  ENDIF
                  GOTO 100
               ENDIF

            ENDIF

            IF( IERR .EQ. 1 ) THEN

C              LE POINT XYZ0 EST EXTERIEUR AUX TETRAEDRES
C              IL N'Y A PAS DE TETRAEDRE derriere LA FACE NOFAMX de
C              NEF LE DERNIER TETRAEDRE DU PARCOURS DANS xyzdste
C              ----------------------------------------------------
C              CALCUL DU POINT D'INTERSECTION DU SEGMENT
C              XEBL1,YEBL1,ZEBL1-XEBL0,YEBL0,ZEBL0 ET DE LA FACE
C              FRONTIERE NOFAMX DU TETRAEDRE NEF
               CALL INARTR( XYZ1, XYZ0,
     %                      Xe( 1, NOSOFATE(1,NOFAMX) ),
     %                      Xe( 1, NOSOFATE(2,NOFAMX) ),
     %                      Xe( 1, NOSOFATE(3,NOFAMX) ),
     %                      LINTER, PTITR, CBTR )
C              LINTER : -3 SI XYZ1-XYZ0 PAS DE CALCUL DE PTINT
C                       -2 SI XYZ1-XYZ0 EST DANS  LE PLAN DU TRIANGLE
C                       -1 SI XYZ1-XYZ0 PARALLELE AU PLAN DU TRIANGLE SANS ETRE DEDANS
C                        0 SI XYZ1-XYZ0 N'INTERSECTE PAS LE TRIANGLE ENTRE SES 3 St
C                        1 SI XYZ1-XYZ0   INTERSECTE     LE TRIANGLE ENTRE SES 3 St
C                          OU SUR UNE DE SES 3 ARETES
C              PTITR  : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=0
C              CBTR   : LES 3 COORDONNEES BARYCENTRIQUES DE PT DANS LE TRIANGLE

               IF( LINTER .EQ. 1 ) THEN
C
C                 IL EXISTE UN TEL POINT D'INTERSECTION SUR LA FACE NOFAMX de NEF
C                 ---------------------------------------------------------------
                  XEBL0 = PTITR(1)
                  YEBL0 = PTITR(2)
                  ZEBL0 = PTITR(3)

ccc            WRITE(IMPRIM,*)'f3snvp1bp1: Point frontiere intersecte',
ccc     %                      XEBL0, YEBL0, ZEBL0,' Vitesse',V,' EF',NEF
ccc            WRITE(IMPRIM,*)
C
               ELSE
C
C                 IL N'EXISTE PAS DE POINT PROJETE SUR UNE FACE FRONTIERE
C                 => POINT EXTERIEUR AU DOMAINE
C                 -------------------------------------------------------
C                 VITESSE AU POINT AVANT SORTIE DES TETRAEDRES
                  V(1) = VA(1)
                  V(2) = VA(2)
                  V(3) = VA(3)
                  GOTO 80

               ENDIF
C
            ENDIF

C           CALCUL DES 4 COORDONNEES BARYCENTRIQUES DU
C           POINT (XEBL0,YEBL0,ZEBL0) DANS LE TETRAEDRE NEF QUI LE CONTIENT
C           CE SONT ABSCISSE,ORDONNEE,COTE DU POINT DANS L'EF DE REFERENCE

C           Lambda 1
            CB1BL = DETM33( Xe(1,2)-XEBL0, Xe(1,3)-XEBL0, Xe(1,4)-XEBL0,
     %                      Xe(2,2)-YEBL0, Xe(2,3)-YEBL0, Xe(2,4)-YEBL0,
     %                      Xe(3,2)-ZEBL0, Xe(3,3)-ZEBL0, Xe(3,4)-ZEBL0)
     %            / DELTAe
C           Lambda 2
            CB2BL = DETM33( Xe(1,3)-XEBL0, Xe(1,1)-XEBL0, Xe(1,4)-XEBL0,
     %                      Xe(2,3)-YEBL0, Xe(2,1)-YEBL0, Xe(2,4)-YEBL0,
     %                      Xe(3,3)-ZEBL0, Xe(3,1)-ZEBL0, Xe(3,4)-ZEBL0)
     %            / DELTAe
C           Lambda 3
            CB3BL = DETM33( Xe(1,4)-XEBL0, Xe(1,1)-XEBL0, Xe(1,2)-XEBL0,
     %                      Xe(2,4)-YEBL0, Xe(2,1)-YEBL0, Xe(2,2)-YEBL0,
     %                      Xe(3,4)-ZEBL0, Xe(3,1)-ZEBL0, Xe(3,2)-ZEBL0)
     %            / DELTAe
C           Lambda 4
            CB4BL = DETM33( Xe(1,1)-XEBL0, Xe(1,3)-XEBL0, Xe(1,2)-XEBL0,
     %                      Xe(2,1)-YEBL0, Xe(2,3)-YEBL0, Xe(2,2)-YEBL0,
     %                      Xe(3,1)-ZEBL0, Xe(3,3)-ZEBL0, Xe(3,2)-ZEBL0)
     %            / DELTAe

C           Verification (XEBL0,YEBL0,ZEBL0) INTERNE au TETRAEDRE NEF?
            D = ABS(CB1BL) + ABS(CB2BL) + ABS(CB3BL) + ABS(CB4BL)
            IF( D .GT. 1.01D0 ) THEN

ccc               print*,'f3snvp1bp1: Pb NEF=',NEF,' L=',L,' CB1234=',D,
ccc     %                ' CB1BL=',CB1BL,' CB2BL=',CB2BL,' CB3BL=',CB3BL,
ccc     %                ' CB4BL=',CB4BL,'->  XYZBL0 NON INTERNE a NEF'
ccc               print*,'f3snvp1bp1: XYZBL0=',XYZ0,
ccc     %                ' est PROJETE dans NEF.  LINTER=',LINTER

C              PROJECTION DU POINT XYZBL0 DANS LE TETRAEDRE NEF
               IF( CB1BL .LT. 0D0 ) CB1BL=0D0
               IF( CB2BL .LT. 0D0 ) CB2BL=0D0
               IF( CB3BL .LT. 0D0 ) CB3BL=0D0
               IF( CB4BL .LT. 0D0 ) CB4BL=0D0
               D = CB1BL + CB2BL + CB3BL + CB4BL
               CB1BL = CB1BL / D
               CB2BL = CB2BL / D
               CB3BL = CB3BL / D
               CB4BL = CB4BL / D

C         (XEBL0,YEBL0,ZEBL0)=Fe(CB2BL,CB3BL,CB4BL) SUR LE TETRAEDRE NEF
           XEBL0=CB1BL*Xe(1,1)+CB2BL*Xe(1,2)+CB3BL*Xe(1,3)+CB4BL*Xe(1,4)
           YEBL0=CB1BL*Xe(2,1)+CB2BL*Xe(2,2)+CB3BL*Xe(2,3)+CB4BL*Xe(2,4)
           ZEBL0=CB1BL*Xe(3,1)+CB2BL*Xe(3,2)+CB3BL*Xe(3,3)+CB4BL*Xe(3,4)

ccc            print*,'f3snvp1bp1: XYZBL0=',XYZ0,' DANS NEF'

               IF( LINTER .EQ. 1 ) GOTO 30
            ENDIF

 10      ENDDO

ccc      print*,'XYZtn+1=',XEBL1,YEBL1,ZEBL1,' XYZtn=',XEBL0,YEBL0,ZEBL0,
ccc     %       '  V=',V,' EF=',NEF,' CB=',CB2BL,CB3BL,CB4BL


C        NBSSPATE SOUS PAS DE TEMPS DDT INTEGRES et X(tn;tn+1,Fe(bl)) CALCULE
C        LES 4 COORDONNEES BARYCENTRIQUES CB1234BL ONT ETE CALCULEES DANS NEF
C        ====================================================================
C        CALCUL DES Pchapeau(CB1234BL) P1+BULLE
 30      BULLE = 64D0 * CB1BL * CB2BL * CB3BL * CB4BL
         PRCBL0(1) = CB1BL - BULLE
         PRCBL0(2) = CB2BL - BULLE
         PRCBL0(3) = CB3BL - BULLE
         PRCBL0(4) = CB4BL - BULLE
         PRCBL0(5) = 4D0 * BULLE

C        CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL,CB4BL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         V(3) = 0D0
         DO I=1,4
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUSOTE(NEF,I) - 1 )
            V(1) = V(1) +  PRCBL0(I) * VXYZPNtn( NDL+1 )
            V(2) = V(2) +  PRCBL0(I) * VXYZPNtn( NDL+2 )
            V(3) = V(3) +  PRCBL0(I) * VXYZPNtn( NDL+3 )
         ENDDO
C        LE NO DU DERNIER DL VITESSE DU BARYCENTRE
         NDL  = NTDLHB + 3 * NEF
         V(1) = V(1) +  PRCBL0(5) * VXYZPNtn( NDL-2 )
         V(2) = V(2) +  PRCBL0(5) * VXYZPNtn( NDL-1 )
         V(3) = V(3) +  PRCBL0(5) * VXYZPNtn( NDL   )

 80      CONTINUE

cccC        PAS D'ECRETAGE: PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
cccC        -------------------------------------------------------------
ccc 80      VNORM = SQRT( V(1)**2 + V(2)**2 + V(3)**2 )
ccc         D     = VNORM / VITMXtn
ccc         IF( D .GT. 1.5D0 ) THEN
ccc            print*,'f3snvp1bp1: NEF=',NEF,
ccc     %             ' CBT=',ABS(CB1BL)+ABS(CB2BL)+ABS(CB3BL)+ABS(CB4BL),
ccc     %             ' CBBL=',CB1BL,CB2BL,CB3BL,CB4BL
ccc            print*,'f3snvp1bp1: PRCBL0T=',
ccc     %              PRCBL0(1)+PRCBL0(2)+PRCBL0(3)+PRCBL0(4)+PRCBL0(5),
ccc     %             ' PRCBL0=',PRCBL0
ccc            print*,'f3snvp1bp1: Attention ECRETAGE de VNORM=',VNORM,
ccc     %             ' en 1.5 VITMXtn=',1.5D0*VITMXtn,' V=',V
ccc            D = 1.5D0 * VITMXtn / VNORM
ccc            V(1) = V(1) * D
ccc            V(2) = V(2) * D
ccc            V(3) = V(3) * D
ccc         ENDIF

C        CALCUL DU VECTEUR BE DE L'EF NEF00
C        ----------------------------------
C        POIDS * JACOBIEN AU POINT L D'INTEGRATION * Coeftr
 90      RDPOIDS = POTEP5(L) * DELTAe1 * Coeftr
C
C        V(1:3) EST ICI LA VITESSE AU POINT X(tn;tn+1,Fe(bl))
C        CONSTRUCTION DU BLOC ASSOCIE A Vk DU SECOND MEMBRE ELEMENTAIRE
         IK = 0
         DO K=1,3
C           V(k) COMPOSANTE k de la VITESSE AU POINT X(tn;tn+1,Fe(bl))
            D = RDPOIDS * V(K)
            DO I=1,5
C              COMPOSANTE IK du VECTEUR
C              Coeftr Poids Jacobien tP(bl) Vk( X(tn;tn+1,Fe(bl)) )
               BE(IK+I) = BE(IK+I) + PRCBL1(I) * D
            ENDDO
            IK = IK + 5
         ENDDO

C        FIN DE L'INTEGRATION DU POINT L
 100  ENDDO

 9999 RETURN
      END
