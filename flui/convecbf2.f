      SUBROUTINE CONVECBF2( NOSOMT,  NEF00,   NBEF,   NUSOTR,
     %                      MOARET,  MXARET,  MNLARE,
     %                      DT,      XYZSOM,  NDDLNO,
     %                      NTDLHB,  VXYZPN0, VITMAX0,
     %                      V,       IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU VECTEUR VITESSE CONVECTEE ui(tn,X(tn;tn+1,SOMMETS)
C ----- POUR LE PROBLEME DE NAVIER STOKES AVEC TRIANGLES BREZZI-FORTIN
C       INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 1 + BULLE P3 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TRIANGLE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C
C ENTREES:	
C --------
C NOSOMT : NUMERO DU SOMMET DE VITESSE CONVECTEE A CALCULER EN REMONTANT
C          LA CARACTERISTIQUE D'UN TEMPS DT
C NEF00  : NUMERO DE L'EF A TRAITER
C NBEF   : NOMBRE DE TRIANGLES DU MAILLAGE
C NUSOTR : NUMERO DES 3 SOMMETS DES TRIANGLES
C MOARET : NOMBRE DE MOTS DE CHAQUE ARETE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DU TABLEAU LARETE
C MNLARE : ADRESSE MCN DU TABLEAU LARETE cf hachag.f
C DT     : LE PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C Rho    : DENSITE VOLUMIQUE DE MASSE DE L'ELEMENT FINI
C XYZSOM : 3 COORDONNEES DES SOMMETS DES TRIANGLES (Z=0)
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU
C NTDLHB : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET HORS BARYCENTRES
C VXYZPN0: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLHB) au TEMPS tn
C          SUIVI DES 3 VITESSES AU BARYCENTRE DES TRIANGLES
C VITMAX0: NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT
C
C SORTIES:
C --------
C V      : VITESSE CONVECTEE de DT en NOSOMT  ui(tn,X(tn;tn+1,NOSOMT))
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      REAL               DT, XYZSOM(3,*)
      INTEGER            NDDLNO(0:*), IERR,
     %                   NOSOMT, NEF00, NBEF, NUSOTR(NBEF,3),
     %                   MOARET, MXARET, MNLARE, NTDLHB
      DOUBLE PRECISION   VXYZPN0(1:*), VITMAX0, V(2), VA(2)
C
      INTEGER            I, K, ITERT, IMIN, NDL, NS,
     %                   NEF0, NEF1, NBSSPATE
      DOUBLE PRECISION   DDT, D, PRCBL0(4), VIT0MX, PTI(2,3), DMIN,
     %                   Xe(2,3), DELTAe, CB1, BULLE, X1, Y1,
     %                   XRCBL, YRCBL, XEBL0, YEBL0, XEBL1, YEBL1
      INTEGER            LINTER(3)
      INTRINSIC          NINT, SQRT, DBLE, ABS

      IERR  = 0

C     NO DU DL VITESSE AU SOMMET NOSOMT
      NDL  = NDDLNO( NOSOMT - 1 )
C     VITESSE AU NOEUD NOSOMT NON CONVECTEE
      V(1) = VXYZPN0(NDL+1)
      V(2) = VXYZPN0(NDL+2)

C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
C     RECUPERATION DES COORDONNEES Xe(3,2) DES SOMMETS DU TRIANGLE NEF00
      DO I=1,3
         NS = NUSOTR(NEF00,I)
         DO K=1,2
            Xe(K,I) = DBLE( XYZSOM( K, NS ) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      DELTAe = (Xe(1,2) - X1) * (Xe(2,3) - Y1)
     %       - (Xe(1,3) - X1) * (Xe(2,2) - Y1)
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'convecbf2: ATTENTION EF',NEF00,
     %                 ' de SURFACE*2=',DELTAe,' NON PRIS EN COMPTE'
         ELSE
            WRITE(IMPRIM,*) 'convecbf2: ATTENTION FE',NEF00,
     %                      ' of SURFACE*2=',DELTAe,' is NOT COMPUTED'
         ENDIF
         WRITE(IMPRIM,*) 'convecbf2: V(',NOSOMT,')=',V
         GOTO 9900
      ENDIF

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP1BP12D( DT,      XYZSOM,  NEF00,  NBEF,  NUSOTR,
     %                 NDDLNO,  NTDLHB,
     %                 VITMAX0, VXYZPN0, VIT0MX, NBSSPATE, DDT )
      IF( NBSSPATE .LE. 0 ) GOTO 9900

C     QUEL EST LE NO LOCAL A NEF00 DU SOMMET NOSOMT?
C     ----------------------------------------------
      DO K=1,3
         IF( NOSOMT .EQ. NUSOTR(NEF00,K) ) GOTO 10
      ENDDO

C     COORDONNEES BARYCENTRIQUES 2 et 3          DU SOMMET K
C     COORDONNEES DANS LE TRIANGLE DE REFERENCE DU SOMMET K
 10   GOTO(11,12,13),K
 11   XRCBL = 0D0
      YRCBL = 0D0
      GOTO 21
 12   XRCBL = 1D0
      YRCBL = 0D0
      GOTO 21
 13   XRCBL = 0D0
      YRCBL = 1D0

C     POINT INITIAL = XYZSOM(NOSOMT)
 21   XEBL1 = XYZSOM(1,NOSOMT)
      YEBL1 = XYZSOM(2,NOSOMT)

C     CALCUL DU POINT X(tn;tn+1,NOSOMT) c-a-d du POINT a L'INSTANT tn
C     QUI SERA AU SOMMET NOSOMT a L'INSTANT tn+1
C     ---------------------------------------------------------------
C     INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
      NEF0 = 0
      NEF1 = NEF00

C     NBSSPATE PAS DE TEMPS ENTRE t-DT et t pour calculer
C     X(tn;tn+1,(XEBL1,YEBL1))
      DO 30 ITERT = 1, NBSSPATE
C
C        SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
         VA(1) = V(1)
         VA(2) = V(2)

C        SAUVEGARDE DU POINT AVANT NOUVELLE REMONTEE DE LA CONVECTION
         XEBL0 = XEBL1
         YEBL0 = YEBL1
C
C        CALCUL DES Pchapeau(XYRCBL)   POLYNOMES P1 + BULLE
         CB1   = 1.D0 - XRCBL - YRCBL
         BULLE = 9D0 * CB1 * XRCBL * YRCBL
         PRCBL0(1) = CB1   - BULLE
         PRCBL0(2) = XRCBL - BULLE
         PRCBL0(3) = YRCBL - BULLE
         PRCBL0(4) = 3D0 * BULLE

C        CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
         V(1) = 0D0
         V(2) = 0D0
         DO I=1,3
C           NO DU DL I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,NOSOMT)
            NDL  = NDDLNO( NUSOTR(NEF1,I) - 1 )
            V(1) = V(1) + PRCBL0(I) * VXYZPN0( NDL+1 )
            V(2) = V(2) + PRCBL0(I) * VXYZPN0( NDL+2 )
         ENDDO
C        LE NO DU DERNIER DL VITESSE DU BARYCENTRE
         NDL  = NTDLHB + 2 * NEF1
         V(1) = V(1) + PRCBL0(4) * VXYZPN0( NDL-1 )
         V(2) = V(2) + PRCBL0(4) * VXYZPN0( NDL   )

C        LE POINT (XEBL0,YEBL0) A ETE TRANSPORTE de (XEBL1,YEBL1)
C        AU SOUS-PAS DDT
         XEBL1 = XEBL0 - DDT * V(1)
         YEBL1 = YEBL0 - DDT * V(2)

C        SAUVEGARDE DU DERNIER EF PARCOURU
         NEF0 = NEF1

C        DANS QUEL ELEMENT FINI NEF du MAILLAGE EST CE POINT (XEBL1,YEBL1)?
         CALL XYDSTR( 1, NBEF, 3, NUSOTR, MOARET, MXARET, MNLARE,
     %                XYZSOM, NEF0, XEBL1, YEBL1, NEF1, IERR )
C        IERR=0 PAS D'ERREUR LE POINT XEBL1, YEBL1 EST DANS LE TRIANGLE NEF1
C             1 POINT PROCHE D'UNE ARETE NON RETROUVEE DANS LE MAILLAGE
C             2 ARETE APPARTENANT A 3 TRIANGLES
C             3 POINT SORTANT. PAS de TRIANGLE derriere LA PLUS PROCHE ARETE
C             4 NOMBRE DE TRIANGLES PARCOURUS TROP GRAND (BOUCLE?)
ccc      IF( IERR .NE. 0 ) THEN
ccc         WRITE(IMPRIM,*) 'convecbf2: sortie xydstr: ierr=',ierr
ccc         WRITE(IMPRIM,*) 'convecbf2: AUCUN TRIANGLE contenant ',
ccc                          XEBL1, YEBL1
ccc         WRITE(IMPRIM,*) 'convecbf2: CAS TRAITE ENSUITE'
ccc      ENDIF

         IF( IERR .GE. 2 ) THEN
C           VITESSE AU NOEUD NONOEU AVANT SORTIE DES TRIANGLES
            V(1) = VA(1)
            V(2) = VA(2)
            IERR = 0
            GOTO 9900
         ENDIF

         IF( NEF0 .NE. NEF1 ) THEN
C           RECUPERATION DES COORDONNEES Xe(2,3) DES SOMMETS DU TRIANGLE NEF1
            DO I=1,3
               NS = NUSOTR(NEF1,I)
               DO K=1,2
                  Xe(K,I) = DBLE( XYZSOM( K, NS ) )
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

         IF( IERR .EQ. 1 ) THEN
C
C           LE POINT EST EXTERIEUR AU DOMAINE ET DU TRIANGLE NEF1
C           CALCUL DU POINT D'INTERSECTION DU SEGMENT
C           XEBL1,YEBL1-XEBL0,YEBL0 ET DE L'ARETE
C           SUR LA FRONTIERE DU TRIANGLE NEF1
            IMIN = 3
            DO I=1,3
               CALL INTARSE( XEBL1, YEBL1, XEBL0, YEBL0,
     %                       Xe(1,IMIN), Xe(2,IMIN), Xe(1,I), Xe(2,I),
     %                       LINTER(I), PTI(1,I), PTI(2,I) )
C              LINTER : -1 SI S3-S4 PARALLELE A S1-S2
C                        0 SI S3-S4 N'INTERSECTE PAS S1-S2 ENTRE CES 2 SOMMETS
C                        1 SI S3-S4   INTERSECTE     S1-S2 ENTRE CES 2 SOMMETS
C              PTI(I) :  2 COORDONNEES DU POINT D'INTERSECTION DU SEGMENT X3Y3-X4Y4
C                          SUR L'ARETE X1Y1-X2Y2
               IMIN = I
            ENDDO
C
C           UNE OU PLUSIEURS ARETES D'INTERSECTION
C           RECHERCHE DU POINT D'INTERSECTION LE PLUS PROCHE DE (XEBL1,YEBL1)
            IMIN = 0
            DMIN = 1D100
            DO I=1,3
               IF( LINTER(I) .EQ. 1 ) THEN
                  D = (XEBL1-PTI(1,I))**2 + (YEBL1-PTI(2,I))**2
                  IF( D .LT. DMIN ) THEN
                     IMIN   = I
                     DMIN = D
                  ENDIF
               ENDIF
            ENDDO

C           LE POINT D'INTERSECTION SUR L'ARETE LA PLUS PROCHE DE XEBL1,YEBL1
            IF( IMIN .NE. 0 ) THEN
C              IL EXISTE UN TEL POINT D'INTERSECTION
               XEBL1 = PTI(1,IMIN)
               YEBL1 = PTI(2,IMIN)
ccc            WRITE(IMPRIM,*) 'convecbf2: Point frontiere Calcule',
ccc     %                       XEBL1, YEBL1,' EF',NEF1
            ELSE
C              IL N'EXISTE PAS DE POINT PROJETE SUR UNE ARETE FRONTIERE.
C              => POINT EXTERIEUR AU DOMAINE
C              LE POINT EST EXTERIEUR AU DOMAINE ET DU TRIANGLE NEF1
C              ------------------------------------------------------
cccC              NO DU DL VITESSE AU SOMMET NOSOMT
ccc               NDL  = NDDLNO( NOSOMT - 1 )
cccC              => VITESSE AU NOEUD NOSOMT NON CONVECTEE
ccc               V(1) = VXYZPN0(NDL+1)
ccc               V(2) = VXYZPN0(NDL+2)
C              VITESSE AU NOEUD NOSOMT AVANT SORTIE DES TRIANGLES
               V(1) = VA(1)
               V(2) = VA(2)
               IERR = 0
               GOTO 9900
            ENDIF

         ENDIF
C
C        CALCUL DES 2 DERNIERES COORDONNEES BARYCENTRIQUES DU
C        POINT (XEBL1,YEBL1) DANS LE TRIANGLE NEF1 QUI LE CONTIENT
C        CE SONT L'ABSCISSE ET ORDONNEE DU POINT DANS L'EF DE REFERENCE
         XRCBL = ( (Xe(1,3)-XEBL1) * (Y1     -YEBL1)
     %           - (X1     -XEBL1) * (Xe(2,3)-YEBL1) ) / DELTAe
         YRCBL = ( (X1     -XEBL1) * (Xe(2,2)-YEBL1)
     %           - (Xe(1,2)-XEBL1) * (Y1     -YEBL1) ) / DELTAe

 30   CONTINUE
ccc
ccc      print*,'XYtn+1=',XEBL1,YEBL1,' XYtn=',XEBL0,YEBL0,
ccc     %'  V=',V,' EF=',NEF1,' CB=',XRCBL,YRCBL
C
C     Calcul de X(tn;tn+1,NOSOMT) fait:
C     Les Pchapeau(X(tn;tn+1,NOSOMT)) POLYNOMES P1 + BULLE
C     DANS LE TRIANGLE NEF1 contenant (XEBL0,YEBL0) sont calcules
C     ----------------------------------------------------------
C     CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
      CB1   = 1.D0 - XRCBL - YRCBL
      BULLE = 9D0 * CB1 * XRCBL * YRCBL
      PRCBL0(1) = CB1   - BULLE
      PRCBL0(2) = XRCBL - BULLE
      PRCBL0(3) = YRCBL - BULLE
      PRCBL0(4) = 3D0 * BULLE
C
C     CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
      V(1) = 0D0
      V(2) = 0D0
      DO I=1,3
C        NO DU DL I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C        AU POINT X(tn;tn+1,NOSOMT)
         NDL  = NDDLNO( NUSOTR(NEF1,I) - 1 )
         V(1) = V(1) +  PRCBL0(I) * VXYZPN0( NDL+1 )
         V(2) = V(2) +  PRCBL0(I) * VXYZPN0( NDL+2 )
      ENDDO
C     LE NO DU DERNIER DL VITESSE DU BARYCENTRE
      NDL  = NTDLHB + 2 * NEF1
      V(1) = V(1) + PRCBL0(4) * VXYZPN0( NDL-1 )
      V(2) = V(2) + PRCBL0(4) * VXYZPN0( NDL   )
C
C     PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
 9900 D = SQRT( V(1)**2 + V(2)**2 ) / VITMAX0
      IF( D .GT. 1.001D0 ) THEN
         V(1) = V(1) / D
         V(2) = V(2) / D
ccc         print*,'convecbf2: nosomt=',nosomt,' d=',d,' V=',V,
ccc     %          ' XRCBL=',XRCBL,' YRCBL=',YRCBL
      ENDIF

      RETURN
      END
