      SUBROUTINE CONVECTH2( NONOEU,  NEF00,  NBEF,   NUNOTR,
     %                      MOARET,  MXARET, MNLARE, DT, XYZNOE, NDDLNO,
     %                      VXYZPN0, VITMAX0,
     %                      V,       IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE LA VITESSE ui(tn,X(tn;tn+1,NONOEU)) CONVECTEE EN NONOEU
C ----- PAR REMONTEE DE LA CARACTERISTIQUE en ce NOEUD sur UN TEMPS DT
C       PROBLEME de NAVIER STOKES sur des TRIANGLES TAYLOR HOOD
C       INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 2 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TRIANGLE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C
C ENTREES:	
C --------
C NONOEU : NUMERO DU NOEUD DE VITESSE CONVECTEE A CALCULER EN REMONTANT
C          LA CARACTERISTIQUE D'UN TEMPS DT
C NEF00  : NUMERO D'UN EF CONTENANT LE NOEUD NONOEU
C NBEF   : NOMBRE DE TRIANGLES DU MAILLAGE
C NUNOTR : NUMERO DES 6 NOEUDS DE CHAQUE TRIANGLE
C MOARET : NOMBRE DE MOTS DE CHAQUE ARETE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DU TABLEAU LARETE
C MNLARE : ADRESSE MCN DU TABLEAU LARETE cf hachag.f
C DT     : LE PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C XYZNOE : 3 COORDONNEES DES NOEUDS DES TRIANGLES
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU
C VXYZPN0: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLVP) au TEMPS tn
C VITMAX0: NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT
C
C SORTIES:
C --------
C V      : VITESSE CONVECTEE de DT en NONOEU  ui(tn,X(tn;tn+1,NONOEU))
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      INTEGER            NEF00, NBEF, MOARET, MXARET, MNLARE, IERR
C
      REAL               DT, XYZNOE(3,*)
      INTEGER            NONOEU, NUNOTR(NBEF,6)
      INTEGER            NDDLNO(0:*)
      DOUBLE PRECISION   VXYZPN0(1:*), VITMAX0, V(2), VA(2)
C
      DOUBLE PRECISION   Xe(2,3), X21, Y21, X31, Y31, DELTAe
      DOUBLE PRECISION   D, PRCBL0(6), XRCBL, YRCBL
      DOUBLE PRECISION   XEBL0, YEBL0, XEBL1, YEBL1, X1, Y1
      DOUBLE PRECISION   VIT0MX, PTI(2,3), DMIN, DDT
      INTEGER            LINTER(3)
      INTRINSIC          NINT, SQRT, DBLE
C
      INTEGER            I, K, ITERT, NDL, NS,
     %                   I0, NEF0, NEF1, NBSSPATE
C
      IERR = 0
C     NO DU DL VITESSE AU NOEUD NONOEU
      NDL  = NDDLNO( NONOEU - 1 )
C     VITESSE AU NOEUD NONOEU NON CONVECTEE
      V(1) = VXYZPN0(NDL+1)
      V(2) = VXYZPN0(NDL+2)
C
C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
      DO K = 1, 3
         NS = NUNOTR(NEF00,K)
         DO I = 1, 2
            Xe( I, K ) = XYZNOE( I, NS )
         ENDDO
      ENDDO
C
C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      X1  = Xe(1,1)
      X21 = Xe(1,2) - X1
      X31 = Xe(1,3) - X1

      Y1  = Xe(2,1)
      Y21 = Xe(2,2) - Y1
      Y31 = Xe(2,3) - Y1
C
      DELTAe = ABS( X21*Y31 - X31*Y21 )
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'CONVECTH2: ATTENTION EF',NEF00,
     %                    ' de SURFACE*2=',DELTAe,' NON PRIS EN COMPTE'
     %                   ,' VITESSE INCHANGEE AU NOEUD ',NONOEU        
         ELSE
            WRITE(IMPRIM,*) 'CONVECTH2: ATTENTION FE',NEF00,
     %                       ' of SURFACE*2=',DELTAe,' is NOT COMPUTED'
     %                   ,' UNCHANGED VELOCITY AT NODE ',NONOEU
         ENDIF
         WRITE(IMPRIM,*) 'convecth2: V(',NONOEU,')=',V
         GOTO 9900
      ENDIF

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP2P12D( DT,      XYZNOE,  NEF00,  NBEF,  NUNOTR, NDDLNO,
     %                VITMAX0, VXYZPN0, VIT0MX, NBSSPATE, DDT )
      IF( NBSSPATE .LE. 0 ) GOTO 9900
C
C     QUEL EST LE NO LOCAL A NEF00 DU NOEUD NONOEU?
C     --------------------------------------------
      DO K=1,6
         IF( NONOEU .EQ. NUNOTR(NEF00,K) ) GOTO 10
      ENDDO
C
C     COORDONNEES BARYCENTRIQUES 2 et 3         DU NOEUD K
C     COORDONNEES DANS LE TRIANGLE DE REFERENCE DU NOEUD K
 10   GOTO(11,12,13,14,15,16),K
 11   XRCBL = 0D0
      YRCBL = 0D0
      GOTO 20
 12   XRCBL = 1D0
      YRCBL = 0D0
      GOTO 20
 13   XRCBL = 0D0
      YRCBL = 1D0
      GOTO 20
 14   XRCBL = 0.5D0
      YRCBL = 0D0
      GOTO 20
 15   XRCBL = 0.5D0
      YRCBL = 0.5D0
      GOTO 20
 16   XRCBL = 0D0
      YRCBL = 0.5D0
C
C     POINT INITIAL = XYZNOE(NONOEU)
 20   XEBL1 = XYZNOE(1,NONOEU)
      YEBL1 = XYZNOE(2,NONOEU)
C
C     INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
      NEF0 = 0
      NEF1 = NEF00
C
C     CALCUL DU POINT X(tn;tn+1,NONOEU) c-a-d du POINT a L'INSTANT tn
C     QUI SERA AU NOEUD NONOEU=XYZNOE(1,NONOEU)=(XEBL1,YEBL1) a L'INSTANT tn+1
C     (XRCBL,YRCBL) = Fe-1(XEBL1,YEBL1)
C     ------------------------------------------------------------------------
C     COORDONNEES du POINT INITIAL sur l'EF NEF1 pour les ITERATIONS itert
      XEBL0 = XEBL1
      YEBL0 = YEBL1
C
C     NBSSPATE PAS DE TEMPS ENTRE t-DT et t pour calculer
C     X(tn;tn+1,(XEBL1,YEBL1,ZEBL1))
      DO 30 ITERT = 1, NBSSPATE
C
C        SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
         VA(1) = V(1)
         VA(2) = V(2)
C
C        CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
         D         = 1.D0 - XRCBL - YRCBL
         PRCBL0(1) = D    * ( 2.D0 * D - 1.D0 )
         PRCBL0(2) = 2.D0 * XRCBL * XRCBL - XRCBL
         PRCBL0(3) = 2.D0 * YRCBL * YRCBL - YRCBL
         D         = 4.D0 * D
         PRCBL0(4) = D    * XRCBL
         PRCBL0(5) = 4.D0 * XRCBL * YRCBL
         PRCBL0(6) = D    * YRCBL

ccc         SP =PRCBL0(1)+PRCBL0(2)+PRCBL0(3)+PRCBL0(4)+PRCBL0(5)+PRCBL0(6)
ccc         IF( SP .GT. 1.00001D0 ) then
ccc            print*,'convecth2: somme P1 a P6=',sp
ccc         ENDIF
C
C        CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
         V(1) = 0D0
         V(2) = 0D0
         DO I=1,6
C           NO DU DL I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,NONOEU)
            NDL  = NDDLNO( NUNOTR(NEF1,I) - 1 )
            V(1) = V(1) +  PRCBL0(I) * VXYZPN0( NDL+1 )
            V(2) = V(2) +  PRCBL0(I) * VXYZPN0( NDL+2 )
         ENDDO
C
C        LE POINT (XEBL0,YEBL0) SERA TRANSPORTE en (XEBL1,YEBL1)
C        SUR L'EF COURANT NEF1 de l'instant tn a tn+1 par pas ddt
C        --------------------------------------------------------
         XEBL0 = XEBL0 - DDT * V(1)
         YEBL0 = YEBL0 - DDT * V(2)
C
C        LISTE DES 2 DERNIERS EF PARCOURUS NEF0 NEF1
         NEF0 = NEF1
C
C        DANS QUEL ELEMENT FINI NEF1 du MAILLAGE EST CE POINT (XEBL0,YEBL0)?
         CALL XYDSTR( 1,      NBEF,   6, NUNOTR,
     %                MOARET, MXARET, MNLARE,
     %                XYZNOE, NEF0,   XEBL0, YEBL0, NEF1, IERR )
C        IERR=0 PAS D'ERREUR LE POINT XEBL0, YEBL0 EST DANS LE TRIANGLE NEF1
C             1 POINT PROCHE D'UNE ARETE NON RETROUVEE DANS LE MAILLAGE
C             2 ARETE APPARTENANT A 3 TRIANGLES
C             3 POINT SORTANT. PAS de TRIANGLE derriere LA PLUS PROCHE ARETE
C             4 NOMBRE DE TRIANGLES PARCOURUS TROP GRAND (BOUCLE?)
ccc         IF( IERR .NE. 0 ) THEN
ccc            WRITE(IMPRIM,*) 'convecth2: sortie xydstr: ierr=',ierr
ccc            WRITE(IMPRIM,*) 'convecth2: AUCUN TRIANGLE contenant ',
ccc     %                       XEBL0, YEBL0
ccc            WRITE(IMPRIM,*) 'convecth2: CAS TRAITE ENSUITE'
ccc         ENDIF

         IF( IERR .GE. 2 ) THEN
C           VITESSE AU NOEUD NONOEU AVANT SORTIE DES TRIANGLES
            V(1) = VA(1)
            V(2) = VA(2)
            IERR = 0
            GOTO 9900
         ENDIF

         IF( NEF0 .NE. NEF1 ) THEN
C
C           RECUPERATION DES COORDONNEES Xe(2,3) DES SOMMETS DU TRIANGLE NEF1
            DO I=1,3
               NS = NUNOTR(NEF1,I)
               DO K=1,2
                  Xe(K,I) = DBLE( XYZNOE( K, NS ) )
               ENDDO
            ENDDO
C
C           CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF1
            X1 = Xe(1,1)
            Y1 = Xe(2,1)
            DELTAe = (Xe(1,2) - X1) * (Xe(2,3) - Y1)
     %             - (Xe(1,3) - X1) * (Xe(2,2) - Y1)
            IF( DELTAe .LE. 0D0 ) GOTO 9900
         ENDIF
C
         IF( IERR .EQ. 1 ) THEN
C
C           LE POINT EST EXTERIEUR AU DOMAINE ET DU TRIANGLE NEF1
C           CALCUL DU POINT D'INTERSECTION DU SEGMENT
C           XEBL1,YEBL1 - XEBL0,YEBL0  ET DE L'ARETE
C           SUR LA FRONTIERE DU TRIANGLE NEF1
            I0 = 3
            DO I=1,3
               CALL INTARSE( XEBL1, YEBL1, XEBL0, YEBL0,
     %                       Xe(1,I0), Xe(2,I0), Xe(1,I), Xe(2,I),
     %                       LINTER(I), PTI(1,I), PTI(2,I) )
C              LINTER :-1 SI S3-S4 PARALLELE A S1-S2
C                       0 SI S3-S4 N'INTERSECTE PAS S1-S2 ENTRE CES 2 SOMMETS
C                       1 SI S3-S4   INTERSECTE     S1-S2 ENTRE CES 2 SOMMETS
C              PTI(I) : 2 COORDONNEES DU PT D'INTERSECTION DU SEGMENT X3Y3-X4Y4
C              SUR L'ARETE X1Y1-X2Y2
               I0 = I
            ENDDO
C
C           UNE OU PLUSIEURS ARETES D'INTERSECTION
C           RECHERCHE DU POINT D'INTERSECTION LE PLUS PROCHE DE XYZ0
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
C           LE POINT D'INTERSECTION SUR L'ARETE LA PLUS PROCHE DE XEBL0,YEBL0
            IF( I0 .NE. 0 ) THEN
               XEBL0 = PTI(1,I0)
               YEBL0 = PTI(2,I0)
ccc            WRITE(IMPRIM,*) 'convecth2: Point frontiere Calcule',
ccc  %                          XEBL0, YEBL0,' EF',NEF1
            ELSE
C
C              IL N'EXISTE PAS DE POINT PROJETE SUR UNE ARETE FRONTIERE.
C              POINT EXTERIEUR AU DOMAINE
ccc               print*,'convecth2.f: Pas d''intersection avec ef NEF1=',NEF1
ccc               print*,'convecth2.f: NOEUD NON CONVECTE'
cccC              NO DU DL VITESSE AU NOEUD NONOEU
ccc               NDL  = NDDLNO( NONOEU - 1 )
cccC              VITESSE AU NOEUD NONOEU NON CONVECTEE
ccc               V(1) = VXYZPN0(NDL+1)
ccc               V(2) = VXYZPN0(NDL+2)
C              VITESSE AU NOEUD NONOEU AVANT SORTIE DES TRIANGLES
               V(1) = VA(1)
               V(2) = VA(2)
               IERR = 0
               GOTO 9900
            ENDIF
C
         ENDIF
C
C        CALCUL DES 2 DERNIERES COORDONNEES BARYCENTRIQUES DU
C        POINT (XEBL0,YEBL0) DANS LE TRIANGLE NEF1 QUI LE CONTIENT
C        CE SONT L'ABSCISSE ET ORDONNEE DU POINT DANS L'EF DE REFERENCE
         XRCBL = ( (Xe(1,3)-XEBL0) * (Y1     -YEBL0)
     %           - (X1     -XEBL0) * (Xe(2,3)-YEBL0) ) / DELTAe
         YRCBL = ( (X1     -XEBL0) * (Xe(2,2)-YEBL0)
     %           - (Xe(1,2)-XEBL0) * (Y1     -YEBL0) ) / DELTAe
C
 30   CONTINUE
ccc
ccc      print*,'XYtn+1=',XEBL1,YEBL1,' XYtn=',XEBL0,YEBL0,
ccc     %'  V=',V,' EF=',NEF1,' CB=',XRCBL,YRCBL
C
C     Calcul de X(tn;tn+1,NONOEU) fait:
C     Les Pchapeau(X(tn;tn+1,NONOEU)) POLYNOMES LAGRANGE DE DEGRE 2
C     DANS LE TRIANGLE NEF1 contenant (XEBL0,YEBL0) sont calcules
C     =============================================================
C     CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
      D         = 1.D0 - XRCBL - YRCBL
      PRCBL0(1) = D    * ( 2.D0 * D - 1.D0 )
      PRCBL0(2) = 2.D0 * XRCBL * XRCBL - XRCBL
      PRCBL0(3) = 2.D0 * YRCBL * YRCBL - YRCBL
      D         = 4.D0 * D
      PRCBL0(4) = D    * XRCBL
      PRCBL0(5) = 4.D0 * XRCBL * YRCBL
      PRCBL0(6) = D    * YRCBL
C
C     CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
      V(1) = 0D0
      V(2) = 0D0
      DO I=1,6
C        NO DU NOEUD I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C        AU POINT X(tn;tn+1,NONOEU)
         NDL  = NDDLNO( NUNOTR(NEF1,I) - 1 )
         V(1) = V(1) +  PRCBL0(I) * VXYZPN0( NDL+1 )
         V(2) = V(2) +  PRCBL0(I) * VXYZPN0( NDL+2 )
      ENDDO
C
C     PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
 9900 D = SQRT( V(1)**2 + V(2)**2 ) / VITMAX0
      IF( D .GT. 1.001D0 ) THEN
         V(1) = V(1) / D
         V(2) = V(2) / D
ccc      print*,'convecth2: nonoeu=',nonoeu,' d=',d,' V=',V
      ENDIF

      RETURN
      END
