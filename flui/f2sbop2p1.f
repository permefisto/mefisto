      SUBROUTINE F2SBOP2P1( NEF00,   NTYEF,  NBEF, NUNOTR,
     %                      MOARET,  MXARET, MNLARE,
     %                      DT,      CoefTr, XYZNOE, NDDLNO,
     %                      TempouVit, TEMPER, VXYZPN, VITMAX,
     %                      BE,      IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU VECTEUR TRANSPORT NAVIER STOKES DU TRIANGLE TAYLOR HOOD
C ----- INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 2 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TRIANGLE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C       INTEGRATION EN 7 POINTS DES MEDIANES DU TRIANGLE P5-EXACTE de
C       Si TempouVit=0 calcul de
C          CoefTr Integrale tP2 TEMPER(tn,X(tn;tn+1,x) dx  (TEMPERATURE)
C       Sinon Si TempouVit=1 calcul de
C          CoefTr Integrale tP2 ui(tn,X(tn;tn+1,x) dx   pour i=1,2
C       Sinon Si TempouVit=2 calcul de   ( NON PROGRAMME )
C          CoefTr Integrale tP2 ( TEMPER(tn,x)
C                                - dt u1(tn,x) d TEMPER(tn,x)/dx1
C                                - dt u2(tn,x) d TEMPER(tn,x)/dx2 ) dx
C
C ENTREES:	
C --------
C NEF00  : NUMERO DE L'EF A TRAITER
C NTYEF  : NUMERO DANS LES NO DE TYPE D'EF (SI UN SEUL => =1)
C NBEF   : NOMBRE DE TRIANGLES DU MAILLAGE
C NUNOTR : NUMERO DES 6 NOEUDS DE CHAQUE TRIANGLE P2
C MOARET : NOMBRE DE MOTS DE CHAQUE ARETE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DU TABLEAU LARETE
C MNLARE : ADRESSE MCN DU TABLEAU LARETE cf hachag.f
C DT     : LE PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C CoefTr : DENSITE SURFACIQUE DE MASSE DE L'ELEMENT FINI
C XYZNOE : 3 COORDONNEES DES 6 NOEUDS DES TRIANGLES
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU

C TempouVit: =0 calcul de CoefTr Integrale tP2 TEMPER(tn,X(tn;tn+1,x) dx
C            =1 CoefTr Integrale tP2 Som ui(tn,x) d TEMPER(tn)/dxi dx

C TEMPER : VECTEUR GLOBAL des DL TEMPERATURE au TEMPS tn
C VXYZPN : VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLVP) PAR NOEUDS
C          NOEUD PAR NOEUD au TEMPS tn
C VITMAX : NORME DE LA VITESSE MAXIMALE AU TEMPS tn
C
C SORTIES:
C --------
C BE     : VECTEUR ELEMENTAIRE DE 6 ou 2x6 + 3 COEFFICIENTS selon TempouVit
C          LES 3 DERNIERS COEFFICIENTS BE(13:15) SONT NULS
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Aout 2009
C AUTEUR : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN          Decembre 2009
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C MODIFS : ALAIN PERRONNET Veulettes sur mer              Septembre 2022
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)

      INTEGER            NEF00, NTYEF, NBEF, MOARET, MXARET, MNLARE,IERR

      REAL               DT, XYZNOE(3,*)
      INTEGER            NUNOTR(NBEF,6)
      INTEGER            NDDLNO(0:*), TempouVit, MOTBE

      DOUBLE PRECISION   TEMPER(1:*), VXYZPN(1:*), VITMAX, TEMPERA
      DOUBLE PRECISION   BE(*)
      DOUBLE PRECISION   Xe(2,3), DELTAe1, DELTAe
      DOUBLE PRECISION   CoefTr, RDPOIDS, D, PRCBL1(6), PRCBL0(6),
     %                   V(2), VA(2), VIT0MX, XRCBL, YRCBL
      DOUBLE PRECISION   XEBL0, YEBL0, XEBL1, YEBL1, X1, Y1
      DOUBLE PRECISION   PTI(2,3), DMIN, DDT
      INTEGER            LINTER(3)
      INTEGER            I, K, L, ITERT, IK, NDL, NS,
     %                   I0, NEF, NEF0, NBSSPATE
      INTRINSIC          NINT, SQRT, DBLE

      include"./incl/trp5poco.inc"
C     POIDS ET POINTS D'INTEGRATION SUR LE TRIANGLE UNITE P5 EXACTE
C     DOUBLE PRECISION   POTRP5( NBPTIP5 ), COTRP5( 2, NBPTIP5 )

C     INITIALISATION A ZERO DU VECTEUR ELEMENTAIRE BE
      IF( TempouVit .EQ. 0 ) THEN
         MOTBE = 6
      ELSE
         MOTBE = 12
      ENDIF
      DO I=1,MOTBE
         BE(I) = 0D0
      ENDDO
      IERR = 0

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP2P12D( DT,     XYZNOE, NEF00,  NBEF,  NUNOTR, NDDLNO,
     %                VITMAX, VXYZPN, VIT0MX, NBSSPATE, DDT )

C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
      DO K = 1, 3
         NS = NUNOTR(NEF00,K)
         DO I = 1, 2
            Xe( I, K ) = DBLE( XYZNOE( I, NS ) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF00
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      DELTAe1 = ABS( (Xe(1,2) - X1) * (Xe(2,3) - Y1)
     %             - (Xe(1,3) - X1) * (Xe(2,2) - Y1) )
      IF( DELTAe1 .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'f2sbop2p1: ATTENTION EF',NEF00,
     %                    ' de SURFACE*2=',DELTAe1,' NON PRIS EN COMPTE'
         ELSE
            WRITE(IMPRIM,*) 'f2sbop2p1: ATTENTION FE',NEF00,
     %                       ' of SURFACE*2=',DELTAe1,' is NOT COMPUTED'
         ENDIF
         GOTO 9999
      ENDIF

C     CALCUL DU VECTEUR tPchapeau(bl) ( Pchapeau(X(tn;tn+1,Fe(bl)) Vitesse(k) )
C     =========================================================================
C     Boucle sur les points d'integration dans le triangle NEF00
      DO 100 L=1,NBPTIP5

C        INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
         IERR = 0
         NEF0 = 0
         NEF  = NEF00

C        RECUPERATION DES COORDONNEES Xe(2,3) DES SOMMETS DU TRIANGLE NEF
         DO I=1,3
            NS = NUNOTR(NEF,I)
            DO K=1,2
               Xe(K,I) = DBLE( XYZNOE( K, NS ) )
            ENDDO
         ENDDO

C        CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
         X1 = Xe(1,1)
         Y1 = Xe(2,1)
         DELTAe = (Xe(1,2) - X1) * (Xe(2,3) - Y1)
     %          - (Xe(1,3) - X1) * (Xe(2,2) - Y1)
         IF( DELTAe .LE. 0D0 ) GOTO 100

C        COORDONNEES DU POINT bl D'INTEGRATION DANS EF REFERENCE
         XRCBL = COTRP5( 1, L )
         YRCBL = COTRP5( 2, L )

C        COORDONNEES de (XEBL1,YEBL1)=Fe(bl) SUR LE TRIANGLE COURANT a tn+1
         D     = 1.D0 - XRCBL - YRCBL
         XEBL1 = D * Xe(1,1) + XRCBL * Xe(1,2) + YRCBL * Xe(1,3)
         YEBL1 = D * Xe(2,1) + XRCBL * Xe(2,2) + YRCBL * Xe(2,3)

C        CALCUL DES Pchapeau(bl)   POLYNOMES LAGRANGE DE DEGRE 2
         PRCBL1(1) = D    * ( 2.D0 * D - 1.D0 )
         PRCBL1(2) = 2.D0 * XRCBL * XRCBL - XRCBL
         PRCBL1(3) = 2.D0 * YRCBL * YRCBL - YRCBL
         D         = 4.D0 * D
         PRCBL1(4) = D    * XRCBL
         PRCBL1(5) = 4.D0 * XRCBL * YRCBL
         PRCBL1(6) = D    * YRCBL

C        CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         DO I=1,6
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUNOTR(NEF,I) - 1 )
            V(1) = V(1) +  PRCBL1(I) * VXYZPN( NDL+1 )
            V(2) = V(2) +  PRCBL1(I) * VXYZPN( NDL+2 )
         ENDDO

         IF( NBSSPATE .LE. 0 ) THEN
C           SI VITESSE QUASI NULLE PAS DE REMONTEE DE LA CARACTERISTIQUE
C           Les  Pchapeau(XYRCBL) des POLYNOMES LAGRANGE de DEGRE 2
C           SONT DEJA CALCULES DANS PRCBL1
            DO I=1,6
               PRCBL0(I) = PRCBL1(I)
            ENDDO
            GOTO 80
         ENDIF

C        ---------------------------------------------------------------
C        CALCUL DU POINT X(tn;tn+1,Fe(bl)) c-a-d du POINT a L'INSTANT tn
C        QUI SERA AU POINT D'INTEGRATION Fe(Bl) a L'INSTANT tn+1
C        (XRCBL,YRCBL) = Fe-1(XEBL1,YEBL1)
C        ---------------------------------------------------------------
C        COORDONNEES du POINT INITIAL sur l'EF e pour les ITERATIONS n
         XEBL0 = XEBL1
         YEBL0 = YEBL1
C
C        NBSSPATE PAS DE TEMPS DDT ENTRE t-DT et t pour calculer
C        X(tn;tn+1,(XEBL1,YEBL1,ZEBL1))
         DO 10 ITERT = 1, NBSSPATE

C           SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
            VA(1) = V(1)
            VA(2) = V(2)

C           CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
            D         = 1.D0 - XRCBL - YRCBL
            PRCBL0(1) = D    * ( 2.D0 * D - 1.D0 )
            PRCBL0(2) = 2.D0 * XRCBL * XRCBL - XRCBL
            PRCBL0(3) = 2.D0 * YRCBL * YRCBL - YRCBL
            D         = 4.D0 * D
            PRCBL0(4) = D    * XRCBL
            PRCBL0(5) = 4.D0 * XRCBL * YRCBL
            PRCBL0(6) = D    * YRCBL
C
C           CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF
            V(1) = 0D0
            V(2) = 0D0
            DO I=1,6
C              NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C              AU POINT X(tn;tn+1,Fe(bl))
               NDL  = NDDLNO( NUNOTR(NEF,I) - 1 )
               V(1) = V(1) +  PRCBL0(I) * VXYZPN( NDL+1 )
               V(2) = V(2) +  PRCBL0(I) * VXYZPN( NDL+2 )
            ENDDO
C
C           LE POINT (XEBL0,YEBL0) SERA TRANSPORTE en (XEBL1,YEBL1)
C           SUR L'EF COURANT NEF00 de l'instant tn a tn+1 par pas ddt
            XEBL0 = XEBL0 - DDT * V(1)
            YEBL0 = YEBL0 - DDT * V(2)
C
C           SAUVEGARDE DU DERNIER EF PARCOURU
            NEF0 = NEF
C
C           DANS QUEL ELEMENT FINI NEF du MAILLAGE EST CE POINT (XEBL0,YEBL0)?
            CALL XYDSTR( NTYEF,  NBEF,   6, NUNOTR,
     %                   MOARET, MXARET, MNLARE,
     %                   XYZNOE, NEF0, XEBL0, YEBL0, NEF, IERR )
C           IERR=0 PAS D'ERREUR LE POINT XEBL0, YEBL0 EST DANS LE TRIANGLE NEF
C                1 POINT PROCHE D'UNE ARETE NON RETROUVEE DANS LE MAILLAGE
C                2 ARETE APPARTENANT A 3 TRIANGLES
C                3 POINT SORTANT. PAS de TRIANGLE derriere LA PLUS PROCHE ARETE
C                4 NOMBRE DE TRIANGLES PARCOURUS TROP GRAND (BOUCLE?)
            IF( IERR .NE. 0 ) THEN
               WRITE(IMPRIM,*) 'f2sbop2p1: sortie xydstr: ierr=',ierr,
     %                         ' No EF=',NEF00,' Pt integration=',L
               WRITE(IMPRIM,*) 'f2sbop2p1: AUCUN TRIANGLE contenant ',
     %                          XEBL0, YEBL0
               WRITE(IMPRIM,*) 'f2sbop2p1: le CAS est TRAITE ENSUITE'
            ENDIF

            IF( IERR .GE. 2 ) THEN
C              VITESSE AU NOEUD NONOEU AVANT SORTIE DES TRIANGLES
               IF( IERR .EQ. 4 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10000) V, VA
                  ELSE
                     WRITE(IMPRIM,20000) V, VA
                  ENDIF
10000 FORMAT('f2sbop2p1: Avant sortie des TRIANGLES: VITESSE CALCULEE='
     %                   ,2G15.6,' VITESSE IMPOSEE =', 2G15.6)
20000 FORMAT('f2sbop2p1: Before TRIANGLES EXIT: COMPUTED VELOCITY='
     %                   ,2G15.6,' IMPOSED VELOCITY=', 2G15.6)
               ENDIF
               V(1) = VA(1)
               V(2) = VA(2)
               IERR = 0
               GOTO 80
            ENDIF
C
            IF( NEF .NE. NEF0 ) THEN
C              RECUPERATION DES COORDONNEES Xe(2,3) DES SOMMETS DU TRIANGLE NEF
               DO I=1,3
                  NS = NUNOTR(NEF,I)
                  DO K=1,2
                     Xe(K,I) = DBLE( XYZNOE( K, NS ) )
                  ENDDO
               ENDDO
C
C              CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
               X1 = Xe(1,1)
               Y1 = Xe(2,1)
               DELTAe = (Xe(1,2) - X1) * (Xe(2,3) - Y1)
     %                - (Xe(1,3) - X1) * (Xe(2,2) - Y1)
               IF( DELTAe .LE. 0D0 ) GOTO 100
            ENDIF
C
            IF( IERR .EQ. 1 ) THEN
C
C              LE POINT EST EXTERIEUR AU DOMAINE ET DU TRIANGLE NEF
C              CALCUL DU POINT D'INTERSECTION DU SEGMENT
C              XEBL1,YEBL1-XEBL0,YEBL0 ET DE L'ARETE
C              SUR LA FRONTIERE DU TRIANGLE NEF
               I0 = 3
               DO I=1,3
                  CALL INTARSE( XEBL1, YEBL1, XEBL0, YEBL0,
     %                          Xe(1,I0), Xe(2,I0), Xe(1,I), Xe(2,I),
     %                          LINTER(I), PTI(1,I), PTI(2,I) )
C                 LINTER : -1 SI S3-S4 PARALLELE A S1-S2
C                           0 SI S3-S4 N'INTERSECTE PAS S1-S2 ENTRE CES 2 SOMMETS
C                           1 SI S3-S4   INTERSECTE     S1-S2 ENTRE CES 2 SOMMETS
C                 PTI(I) :  2 COORDONNEES DU POINT D'INTERSECTION DU SEGMENT X3Y3-X4Y4
C                           SUR L'ARETE X1Y1-X2Y2
                  I0 = I
               ENDDO
C
C              UNE OU PLUSIEURS ARETES D'INTERSECTION
C              RECHERCHE DU POINT D'INTERSECTION LE PLUS PROCHE DE XYZ0
               I0   = 0
               DMIN = 1D100
               DO I=1,3
                  IF( LINTER(I) .EQ. 1 ) THEN
                     D = (XEBL0-PTI(1,I))**2 + (YEBL0-PTI(2,I))**2
                     IF( D .LT. DMIN ) THEN
                        I0   = I
                        DMIN = D
                     ENDIF
                  ENDIF
               ENDDO
C              LE POINT D'INTERSECTION SUR L'ARETE LA PLUS PROCHE DE XEBL0,YEBL0
               IF( I0 .NE. 0 ) THEN
                  IERR = 0
                  XEBL0 = PTI(1,I0)
                  YEBL0 = PTI(2,I0)
ccc                  WRITE(IMPRIM,*) 'f2sbop2p1: Point frontiere Calcule',
ccc     %                             XEBL0, YEBL0,' EF',NEF
               ELSE

C                 IL N'EXISTE PAS DE POINT PROJETE SUR UNE FACE FRONTIERE
C                 => POINT EXTERIEUR AU DOMAINE
C                 -------------------------------------------------------
C                 VITESSE AU POINT AVANT SORTIE DES TRIANGLES
                  V(1) = VA(1)
                  V(2) = VA(2)
                  IERR = 0
                  GOTO 80
               ENDIF

            ENDIF

C           CALCUL DES 2 FINALES COORDONNEES BARYCENTRIQUES DU
C           POINT (XEBL0,YEBL0) DANS LE TRIANGLE NEF QUI LE CONTIENT
C           CE SONT L'ABSCISSE ET ORDONNEE DU POINT DANS L'EF DE REFERENCE
            XRCBL = ( (Xe(1,3)-XEBL0) * (Y1     -YEBL0)
     %              - (X1     -XEBL0) * (Xe(2,3)-YEBL0) ) / DELTAe
            YRCBL = ( (X1     -XEBL0) * (Xe(2,2)-YEBL0)
     %              - (Xe(1,2)-XEBL0) * (Y1     -YEBL0) ) / DELTAe

 10      ENDDO

ccc      print*,'XYtn+1=',XEBL1,YEBL1,' XYtn=',XEBL0,YEBL0,
ccc     %'  V=',V,' EF=',NEF,' CB=',XRCBL,YRCBL

C        Calcul de X(tn;tn+1,Fe(bl)) fait:
C        Les Pchapeau(X(tn;tn+1,Fe(bl))) POLYNOMES LAGRANGE DE DEGRE 2
C        DANS LE TRIANGLE NEF contenant (XEBL0,YEBL0) sont calcules
C        -------------------------------------------------------------
C        CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
         D         = 1.D0 - XRCBL - YRCBL
         PRCBL0(1) = D    * ( 2.D0 * D - 1.D0 )
         PRCBL0(2) = 2.D0 * XRCBL * XRCBL - XRCBL
         PRCBL0(3) = 2.D0 * YRCBL * YRCBL - YRCBL
         D         = 4.D0 * D
         PRCBL0(4) = D    * XRCBL
         PRCBL0(5) = 4.D0 * XRCBL * YRCBL
         PRCBL0(6) = D    * YRCBL

C        CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         DO I=1,6
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUNOTR(NEF,I) - 1 )
            V(1) = V(1) +  PRCBL0(I) * VXYZPN( NDL+1 )
            V(2) = V(2) +  PRCBL0(I) * VXYZPN( NDL+2 )
         ENDDO

cccC        PAS de PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
ccc 80      D = SQRT( V(1)**2 + V(2)**2 ) / VITMAX
ccc         IF( D .GT. 1.001D0 ) THEN
ccc            V(1) = V(1) / D
ccc            V(2) = V(2) / D
ccc         print*,'f2sbop2p1: d=',d,' V=',V
ccc         ENDIF

C        CALCUL DU VECTEUR BE DE L'EF NEF00
C        ----------------------------------
C        POIDS * JACOBIEN AU POINT L D'INTEGRATION * CoefTr
 80      RDPOIDS = POTRP5(L) * DELTAe1 * CoefTr

         IF( TempouVit .EQ. 0 ) THEN

C           Calcul  CoefTr RPOIDS Jacob tP2(bl) TEMPER(tn,X(tn;tn+1,bl) (TEMPERATURE)
C           ICI LA TEMPERATURE AU POINT X(tn;tn+1,Fe(bl))
            TEMPERA = 0D0
            DO I=1,6
C              NO DU DL I DANS L'EF NEF DE LA TEMPERATURE AU TEMPS tn
C              AU POINT X(tn;tn+1,Fe(bl))
               NDL = NUNOTR(NEF,I)
               TEMPERA = TEMPERA + PRCBL0(I) * TEMPER(NDL)
            ENDDO

            D = RDPOIDS * TEMPERA
            DO I=1,6
C              CoefTr Poids Jacobien tP(bl) TEMPER(X(tn;tn+1,Fe(bl))
C              AU POINT X(tn;tn+1,Fe(bl))
               BE(I) = BE(I) + D * PRCBL1(I)
            ENDDO

         ELSE IF( TempouVit .EQ. 1 ) THEN

C           Calcul  CoefTr RPOIDS Jacob tP2(bl) ui(tn,X(tn;tn+1,bl)   pour i=1,2
C           V(1:2) EST ICI LA VITESSE AU POINT X(tn;tn+1,Fe(bl))
C           CONSTRUCTION DU BLOC ASSOCIE A Vk DU SECOND MEMBRE ELEMENTAIRE
            IK = 0
            DO K=1,2
C              V(k) COMPOSANTE k de la VITESSE
               D = RDPOIDS * V(K)
               DO I=1,6
C                 COMPOSANTE IK du VECTEUR
C                 CoefTr Poids Jacobien tP(bl) Vk(X(tn;tn+1,Fe(bl)) )
                  BE(IK+I) = BE(IK+I) + D * PRCBL1(I)
               ENDDO
               IK = IK + 6
            ENDDO


C        ELSE IF( TempouVit .EQ. 2 ) THEN
C           CoefTr Integrale tP2 ( TEMPER(tn,x)
C                                - dt u1(tn,x) d TEMPER(tn,x)/dx1
C                                - dt u2(tn,x) d TEMPER(tn,x)/dx2 ) dx
C           N'EST PAS CALCULE ICI.   A FAIRE SI NECESSAIRE


         ENDIF

C        FIN DE L'INTEGRATION DU POINT L
 100  ENDDO

 9999 RETURN
      END
