      SUBROUTINE F2EX1P1BP1( DT,     NOEUD,   NBSOM,  XYZNOE,
     %                       NEF1,   NBEF,    NUSOTR,
     %                       MOARET, MXARET,  MNLARE,
     %                       NBNOVI, VITESSE, VITMAX,
     %                       V,      IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU POINT X(tn;tn+1,NOEUD) c-a-d du POINT a L'INSTANT tn
C ----  QUI SERA AU NOEUD a L'INSTANT tn+1 TRANSPORTE PAR LA VITESSE(tn)
C       DU PROBLEME DE NAVIER STOKES POUR UN FLUIDE INCOMPRESSIBLE
C       INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       TRIANGLE BREZZI-FORTIN
C       POLYNOME LAGRANGE DE DEGRE 1 + BULLE P3 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TRIANGLE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C
C ENTREES:	
C --------
C DT     : PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C NOEUD  : NUMERO DU NOEUD VITESSE DANS LE MAILLAGE DES VITESSES
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C XYZNOE : 3 COORDONNEES DES 3 SOMMETS DES TRIANGLES
C
C NEF1   : NUMERO D'UN TRIANGLE CONTENANT CE NOEUD (SOMMET OU BARYCENTRE)
C NBEF   : NOMBRE DE TRIANGLES DU MAILLAGE
C NUSOTR : NUMERO DES 3 SOMMETS DES TRIANGLES
C MOARET : NOMBRE DE MOTS DE CHAQUE ARETE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DU TABLEAU LARETE
C MNLARE : ADRESSE MCN DU TABLEAU LARETE cf hachag.f
C
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DES DL DE LA VITESSE
C VITESSE: VECTEUR VITESSE(1:NBNOVI, 1:2) (COMPOSANTES suivant X et Y)
C VITMAX : NORME DE LA VITESSE MAXIMALE
C
C SORTIES:
C --------
C V      : 2 COMPOSANTES DE LA VITESSE EN X(tn;tn+1,NOEUD)
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      REAL               DT,    XYZNOE(3,*)
      INTEGER            NOEUD, NBSOM, NEF1, NBEF, NUSOTR(NBEF,3),
     %                   MOARET, MXARET,  MNLARE, NBNOVI
      DOUBLE PRECISION   VITESSE(NBNOVI,2), VITMAX
C
      INTEGER            IERR, I0, I, K, ITERT, NS, NS1, NS2,
     %                   NEF, NEF00, NEF01, NBSSPATE, NTYEF
      DOUBLE PRECISION   DBLE, MAX
      DOUBLE PRECISION   D, PREFND1(4), PREFND0(4)
      DOUBLE PRECISION   VIT0MX, S, ARETMX
      DOUBLE PRECISION   V(2), Xe(2,4), DELTAe1, DELTAe,
     %                   CB1, BULLE, X1, Y1, XREFND, YREFND
      DOUBLE PRECISION   XEFND0, YEFND0, XEFND1, YEFND1
      DOUBLE PRECISION   DDT, PTI(2,3), DMIN
      INTEGER            LINTER(3)
      INTRINSIC          NINT, SQRT
C
      INTEGER            NOSOAR(2,3)
      DATA               NOSOAR / 1,2,  2,3,  3,1 /
C
      IERR = 0
C
C     NTYEF : NUMERO DU TYPE DE L'EF DANS LES TYPES DU MAILLAGE
C             DE 1 A NBTYEF. SI UN SEUL TYPE D'EF => NTYEF=1
      NTYEF = 1
C
C     RECUPERATION DES COORDONNEES Xe(2,4) DES NOEUDS DU TRIANGLE NEF1
C     LES 3 SOMMETS D'ABORD ET LE BARYCENTRE ENSUITE
      DO I=1,3
         NS = NUSOTR(NEF1,I)
         DO K=1,2
            Xe(K,I) = DBLE( XYZNOE( K, NS ) )
         ENDDO
      ENDDO
      NS = NBSOM + NEF1
      DO K=1,2
         Xe(K,4) = DBLE( XYZNOE( K, NS ) )
      ENDDO
C
C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      DELTAe1 = ( Xe(1,2)-X1 ) * ( Xe(2,3)-Y1 )
     %        - ( Xe(2,2)-Y1 ) * ( Xe(1,3)-X1 )
C
C     CALCUL DE LA LONGUEUR MAXIMALE DES 3 ARETES DU TRIANGLE NEF1
      ARETMX = 0D0
      DO I=1,3
         NS1 = NOSOAR(1,I)
         NS2 = NOSOAR(2,I)
         S = ( Xe(1,NS2) - Xe(1,NS1) ) ** 2
     %     + ( Xe(2,NS2) - Xe(2,NS1) ) ** 2
         IF( S .GT. ARETMX ) ARETMX = S
      ENDDO
      ARETMX = SQRT( ARETMX )
C
C     VITESSE MAXIMALE AU TEMPS t - DT AUX 3 SOMMETS DU TRIANGLE NEF1
      VIT0MX = 0D0
      DO I=1,3
C        NO DU SOMMET DU TRIANGLE
         NS = NUSOTR(NEF1,I)
         S  = VITESSE(NS,1)**2 + VITESSE(NS,2)**2
         IF( S .GT. VIT0MX ) VIT0MX = S
      ENDDO
      VIT0MX = SQRT( VIT0MX )
C
C     LE SOUS PAS DE TEMPS POUR L'INTEGRATION RETROGRADE EN TEMPS
C     -----------------------------------------------------------
C     LE NOMBRE DE SOUS PAS DE TEMPS DANS LE PAS DE TEMPS DT
      I = NINT( 8 * VIT0MX / VITMAX )
C     4 SOUS PAS DE TEMPS POUR PARCOURIR L'ARETE MAXIMALE A
C     LA VITESSE VIT0MX DURANT LE TEMPS DT
      K = NINT( 4 * VIT0MX*DT / ARETMX )
      NBSSPATE = MAX( 2, I, K )
      NBSSPATE = MIN( 16, NBSSPATE )
C
C     SOUS PAS DE TEMPS ENTRE T-DT et T
      DDT = DT / NBSSPATE
C
C     INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
C     ------------------------------------------------------------
      NEF00  = 0
      NEF01  = 0
      NEF    = NEF1
      DELTAe = DELTAe1
C
C     LES 2 COORDONNEES DU NOEUD SUR LE TRIANGLE NEF1 a tn+1
      XEFND1 = DBLE( XYZNOE( 1, NOEUD ) )
      YEFND1 = DBLE( XYZNOE( 2, NOEUD ) )
C
C     CALCUL DES 2 DERNIERES COORDONNEES BARYCENTRIQUES DU
C     POINT (XEFND1,YEFND1) DANS LE TRIANGLE NEF1 QUI LE CONTIENT
C     CE SONT ABSCISSE,ORDONNEE DU POINT DANS L'EF DE REFERENCE
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
C     Lambda 2
      XREFND = ( (Xe(1,3)-XEFND1) * (Y1     -YEFND1)
     %         - (X1     -XEFND1) * (Xe(2,3)-YEFND1) ) / DELTAe1
C     Lambda 3
      YREFND = ( (X1     -XEFND1) * (Xe(2,2)-YEFND1)
     %         - (Xe(1,2)-XEFND1) * (Y1     -YEFND1) ) / DELTAe1
C
C     CALCUL DES P1Bchapeau(NOEUD)  POLYNOMES LAGRANGE DE DEGRE 1 + BULLE P3
      CB1   = 1.D0 - XREFND - YREFND
      BULLE = 9.D0 * CB1 * XREFND * YREFND
      PREFND1(1) = CB1    - BULLE
      PREFND1(2) = XREFND - BULLE
      PREFND1(3) = YREFND - BULLE
      PREFND1(4) = 3D0 * BULLE
C
C     CALCUL DU POINT X(tn;tn+1,NOEUD) c-a-d du POINT a L'INSTANT tn
C     QUI SERA AU NOEUD a L'INSTANT tn+1 TRANSPORTE PAR LA VITESSE(tn)
C     (XREFND,YREFND) = Fe-1(XEFND1,YEFND1)
C     ---------------------------------------------------------------
C     COORDONNEES du POINT INITIAL sur l'EF e pour les ITERATIONS n
      XEFND0 = XEFND1
      YEFND0 = YEFND1
C
C     NBSSPATE SOUS PAS DE TEMPS ENTRE t-DT et t pour calculer
C     X(tn;tn+1,(XEFND1,YEFND1))
      DO 10 ITERT = 1, NBSSPATE
C
C        CALCUL DES P1Bchapeau(XYREFND) POLYNOMES LAGRANGE DE DEGRE 1+BULLE P3
         CB1   = 1.D0 - XREFND - YREFND
         BULLE = 9.D0 * CB1 * XREFND * YREFND
         PREFND0(1) = CB1    - BULLE
         PREFND0(2) = XREFND - BULLE
         PREFND0(3) = YREFND - BULLE
         PREFND0(4) = 3.D0 * BULLE
C
C        CALCUL DE LA VITESSE V(tn,(XREFND,YREFND)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         DO I=1,3
C           NO DU SOMMET I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,SOMMET)
            NS = NUSOTR(NEF,I)
            V(1) = V(1) + PREFND0(I) * VITESSE( NS, 1 )
            V(2) = V(2) + PREFND0(I) * VITESSE( NS, 2 )
         ENDDO
C        LE NO DU BARYCENTRE DU TRIANGLE NEF
         NS = NBSOM + NEF
         V(1) = V(1) + PREFND0(4) * VITESSE( NS, 1 )
         V(2) = V(2) + PREFND0(4) * VITESSE( NS, 2 )
C
C        LE POINT (XEFND0,YEFND0) SERA TRANSPORTE en (XEFND1,YEFND1)
C        SUR L'EF COURANT NEF1 de l'instant tn a tn+1 par pas ddt
         XEFND0 = XEFND0 - DDT * V(1)
         YEFND0 = YEFND0 - DDT * V(2)
C
C        LISTE DES 3 DERNIERS EF PARCOURUS
         NEF00 = NEF01
         NEF01 = NEF
C
C        DANS QUEL ELEMENT FINI NEF du MAILLAGE EST CE POINT XYEFND0?
         CALL XYDSTR( NTYEF,  NBEF,   3, NUSOTR,
     %                MOARET, MXARET, MNLARE,
     %                XYZNOE, NEF01,  XEFND0, YEFND0, NEF, IERR )
C        IERR=0 PAS D'ERREUR LE POINT  XEFND0, YEFND0 EST DANS LE TRIANGLE NEF
C             1 POINT PROCHE D'UNE ARETE NON RETROUVEE DANS LE MAILLAGE
C             2 POINT SORTANT PAS de TRIANGLE derriere LA PLUS PROCHE ARETE
         IF( IERR .NE. 0 ) THEN
            WRITE(IMPRIM,*) 'f2snvp1bp1: sortie xydstr: ierr=',ierr
            WRITE(IMPRIM,*) 'f2snvp1bp1: AUCUN TRIANGLE contenant ',
     %                       XEFND0, YEFND0
         ENDIF
C
C        RECUPERATION DES 2 COORDONNEES Xe(2,3) DES 3 SOMMETS DU TRIANGLE NEF
         IF( NEF .NE. NEF01 ) THEN
            DO I=1,3
               NS = NUSOTR(NEF,I)
               DO K=1,2
                  Xe(K,I) = DBLE( XYZNOE( K, NS ) )
               ENDDO
            ENDDO
C
C           CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
            X1 = Xe(1,1)
            Y1 = Xe(2,1)
            DELTAe = ( Xe(1,2)-X1 ) * ( Xe(2,3)-Y1 )
     %             - ( Xe(2,2)-Y1 ) * ( Xe(1,3)-X1 )
         ENDIF
C
         IF( IERR .EQ. 3 ) THEN
C
C           LE POINT EST EXTERIEUR AU DOMAINE ET DU TRIANGLE NEF
C           CALCUL DU POINT D'INTERSECTION DU SEGMENT
C           XEFND1,YEFND1-XEFND0,YEFND0 ET DE L'ARETE
C           SUR LA FRONTIERE DU TRIANGLE NEF
            I0 = 3
            DO I=1,3
               CALL INTARSE( XEFND1, YEFND1, XEFND0, YEFND0,
     %                       Xe(1,I0), Xe(2,I0), Xe(1,I), Xe(2,I),
     %                       LINTER(I), PTI(1,I), PTI(2,I) )
C LINTER : -1 SI S3-S4 PARALLELE A S1-S2
C           0 SI S3-S4 N'INTERSECTE PAS S1-S2 ENTRE CES 2 SOMMETS
C           1 SI S3-S4   INTERSECTE     S1-S2 ENTRE CES 2 SOMMETS
C X,Y    :  2 COORDONNEES DU POINT D'INTERSECTION DU SEGMENT X3Y3-X4Y4
C           SUR L'ARETE X1Y1-X2Y2
                  I0 = I
            ENDDO
C
C           UNE OU PLUSIEURS ARETES D'INTERSECTION
C           RECHERCHE DU POINT D'INTERSECTION LE PLUS PROCHE DE XYEFND0
            I0   = 0
            DMIN = 1D100
            DO I=1,3
               IF( LINTER(I) .EQ. 1 ) THEN
                  D = (XEFND0-PTI(1,I))**2 + (YEFND0-PTI(2,I))**2
                  IF( D .LT. DMIN ) THEN
                     I0   = I
                     DMIN = D
                  ENDIF
               ENDIF
            ENDDO
C           LE POINT D'INTERSECTION SUR L'ARETE LA PLUS PROCHE DE XEFND0,YEFND0
            IF( I0 .NE. 0 ) THEN
               XEFND0 = PTI(1,I0)
               YEFND0 = PTI(2,I0)
               WRITE(IMPRIM,*) 'f2snvp1bp1: Point frontiere',
     %                          XEFND0, YEFND0,' EF',NEF
            ELSE
C
C              IL N'EXISTE PAS DE POINT PROJETE SUR UNE ARETE FRONTIERE.
C              POURQUOI????
               print *,'Pas d''intersection avec ef NEF=',NEF
               print *,'f2ex1p1bp1.f trouver un autre algo...'
               read *,I
            ENDIF
C
         ENDIF
C
C        CALCUL DES 3 DERNIERES COORDONNEES BARYCENTRIQUES DU
C        POINT (XEFND0,YEFND0) DANS LE TRIANGLE NEF QUI LE CONTIENT
C        CE SONT ABSCISSE,ORDONNEE,COTE DU POINT DANS L'EF DE REFERENCE
C        Lambda 2
         XREFND = ( (Xe(1,3)-XEFND0) * (Y1     -YEFND0)
     %            - (X1     -XEFND0) * (Xe(2,3)-YEFND0) ) / DELTAe
C        Lambda 3
         YREFND = ( (X1     -XEFND0) * (Xe(2,2)-YEFND0)
     %            - (Xe(1,2)-XEFND0) * (Y1     -YEFND0) ) / DELTAe
C
         IF( IERR .EQ. 3 ) GOTO 20
C        SI LE POINT EST EXTERIEUR AU MAILLAGE,
C        LE POINT ACTUEL EST SUR UNE ARETE DE LA FRONTIERE
C
 10   CONTINUE
ccc
ccc   print*,'XYtn+1=',XEFND1,YEFND1,' XYtn=',XEFND0,YEFND0,
ccc  %'  V=',V,' EF=',NEF,' CB=',XREFND,YREFND
C
C     calcul de X(tn;tn+1,NOEUD) effectue:
C     Les Pchapeau(X(tn;tn+1,NOEUD)) POLYNOMES BREZZI-FORTIN
C     DANS LE TRIANGLE NEF contenant (XEFND0,YEFND0) ont ete calcules
C     ======================================================================
C     CALCUL DES Pchapeau(XYREFND) P1+BULLE
 20   CB1   = 1.D0 - XREFND - YREFND
      BULLE = 9.D0 * CB1 * XREFND * YREFND
      PREFND0(1) = CB1    - BULLE
      PREFND0(2) = XREFND - BULLE
      PREFND0(3) = YREFND - BULLE
      PREFND0(4) = 3.D0 * BULLE
C
C     CALCUL DE LA VITESSE V(tn,(XREFND,YREFND)) dans l'EF NEF
      V(1) = 0D0
      V(2) = 0D0
      DO I=1,3
C        NO DU SOMMET I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C        AU POINT X(tn;tn+1,SOMMET)
         NS   = NUSOTR(NEF,I)
         V(1) = V(1) + PREFND0(I) * VITESSE( NS, 1 )
         V(2) = V(2) + PREFND0(I) * VITESSE( NS, 2 )
      ENDDO
C     LE NO DU BARYCENTRE DU TRIANGLE NEF
      NS = NBSOM + NEF
      V(1) = V(1) + PREFND0(4) * VITESSE( NS, 1 )
      V(2) = V(2) + PREFND0(4) * VITESSE( NS, 2 )
C
C     V(1:2) EST ICI LA VITESSE AU POINT X(tn;tn+1,NOEUD)
      RETURN
      END
