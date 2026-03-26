      SUBROUTINE CONVECBF3( NOSOMT,  NEF00,   NBEF,   NUSOTE,
     %                      MOFACE,  MXFACE,  LFACES,
     %                      DT,      XYZSOM,  NDDLNO,
     %                      NTDLHB,  VXYZPN0, VITMAX0,
     %                      V,       IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU VECTEUR VITESSE CONVECTEE ui(tn,X(tn;tn+1,SOMMETS)
C ----- POUR LE PROBLEME DE NAVIER STOKES AVEC TETRAEDRES BREZZI-FORTIN
C       INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 1 + BULLE P4 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TETRAEDRE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C
C ENTREES:	
C --------
C NOSOMT : NUMERO DU SOMMET DE VITESSE CONVECTEE A CALCULER EN REMONTANT
C          LA CARACTERISTIQUE D'UN TEMPS DT
C NEF00  : NUMERO DE L'EF A TRAITER
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NUSOTE : NUMERO DES 4 SOMMETS DES TETRAEDRES
C MOFACE : NOMBRE DE MOTS DE CHAQUE FACE DU TABLEAU LFACES
C MXFACE : NOMBRE MAXIMAL D'FACES DU TABLEAU LFACES
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
C Rho    : DENSITE VOLUMIQUE DE MASSE DE L'ELEMENT FINI
C XYZSOM : 3 COORDONNEES DES SOMMETS DES TETRAEDRES
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU
C NTDLHB : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET HORS BARYCENTRES
C VXYZPN0: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLHB) au TEMPS tn
C          SUIVI DES 3 VITESSES AU BARYCENTRE DES TETRAEDRES
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
      INTEGER            NOSOMT, NEF00, NBEF, NUSOTE(NBEF,4),
     %                   MOFACE, MXFACE, LFACES(MOFACE,MXFACE), NTDLHB
      INTEGER            NDDLNO(0:*), IERR
      DOUBLE PRECISION   VXYZPN0(1:*), VITMAX0, V(3), VA(3)
C
      INTEGER            I, K, L, ITERT, IMIN, NDL, NS,
     %                   NEF1, NEF0, NBSSPATE
      DOUBLE PRECISION   DETM33, D, PRCBL0(5), VIT0MX
      DOUBLE PRECISION   Xe(3,4), DELTAe,
     %                   CB1, BULLE, X1, Y1, Z1,  XRCBL, YRCBL, ZRCBL
      DOUBLE PRECISION   XEBL0, YEBL0, ZEBL0, XEBL1, YEBL1, ZEBL1,
     %                   XYZ0(3), XYZ1(3)
      EQUIVALENCE      (XYZ0(1),XEBL0), (XYZ0(2),YEBL0), (XYZ0(3),ZEBL0)
      EQUIVALENCE      (XYZ1(1),XEBL1), (XYZ1(2),YEBL1), (XYZ1(3),ZEBL1)
      DOUBLE PRECISION   DDT, CB(3,4), PTI(3,4), DMIN
      INTEGER            LINTER(4), NOFAMX
      INTRINSIC          NINT, SQRT, DBLE, ABS

C     LE NUMERO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER            NOSOFATE(3,4)
      DATA               NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      IERR  = 0

C     NO DU DL VITESSE AU SOMMET NOSOMT
      NDL  = NDDLNO( NOSOMT - 1 )
C     VITESSE AU NOEUD NOSOMT NON CONVECTEE
      V(1) = VXYZPN0(NDL+1)
      V(2) = VXYZPN0(NDL+2)
      V(3) = VXYZPN0(NDL+3)

C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
      DO L = 1, 4
         DO I = 1, 3
            Xe( I, L ) = XYZSOM( I, NUSOTE(NEF00,L) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      Z1 = Xe(3,1)
      DELTAe = ABS( DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                      Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                      Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 ) )

      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'convecbf3: ATTENTION EF',NEF00,
     %                 ' de VOLUME*6=',DELTAe,' NON PRIS EN COMPTE'
         ELSE
            WRITE(IMPRIM,*) 'convecbf3: ATTENTION FE',NEF00,
     %                      ' of VOLUME*6=',DELTAe,' is NOT COMPUTED'
         ENDIF
         WRITE(IMPRIM,*) 'convecbf3: V(',NOSOMT,')=',V
         GOTO 9900
      ENDIF

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP1BP13D( DT,      XYZSOM,  NEF00,   NBEF, NUSOTE, NDDLNO,
     %                 VITMAX0, VXYZPN0, VIT0MX,  NBSSPATE, DDT )
      IF( NBSSPATE .LE. 0 ) GOTO 9900

C     QUEL EST LE NO LOCAL A NEF00 DU SOMMET NOSOMT?
C     ----------------------------------------------
      DO K=1,4
         IF( NOSOMT .EQ. NUSOTE(NEF00,K) ) GOTO 10
      ENDDO

C     COORDONNEES BARYCENTRIQUES 2 et 3          DU SOMMET K
C     COORDONNEES DANS LE TETRAEDRE DE REFERENCE DU SOMMET K
 10   GOTO(11,12,13,14),K
 11   XRCBL = 0D0
      YRCBL = 0D0
      ZRCBL = 0D0
      GOTO 21
 12   XRCBL = 1D0
      YRCBL = 0D0
      ZRCBL = 0D0
      GOTO 21
 13   XRCBL = 0D0
      YRCBL = 1D0
      ZRCBL = 0D0
      GOTO 21
 14   XRCBL = 0D0
      YRCBL = 0D0
      ZRCBL = 1D0
      GOTO 21

C     POINT INITIAL = XYZSOM(NOSOMT)
 21   XEBL1 = XYZSOM(1,NOSOMT)
      YEBL1 = XYZSOM(2,NOSOMT)
      ZEBL1 = XYZSOM(3,NOSOMT)

C     CALCUL DU POINT X(tn;tn+1,NOSOMT) c-a-d du POINT a L'INSTANT tn
C     QUI SERA AU SOMMET NOSOMT a L'INSTANT tn+1
C     ---------------------------------------------------------------
C     INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
      NEF0 = 0
      NEF1 = NEF00

C     NBSSPATE PAS DE TEMPS ENTRE t-DT et t pour calculer
C     X(tn;tn+1,(XEBL1,YEBL1,ZEBL1))
      DO 30 ITERT = 1, NBSSPATE
C
C        SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
         VA(1) = V(1)
         VA(2) = V(2)
         VA(3) = V(3)

C        SAUVEGARDE DU POINT AVANT NOUVELLE REMONTEE DE LA CONVECTION
         XEBL0 = XEBL1
         YEBL0 = YEBL1
         ZEBL0 = ZEBL1

C        CALCUL DES Pchapeau(XYRCBL)   POLYNOMES LAGRANGE de DEGRE 2
         CB1   = 1.D0 - XRCBL - YRCBL - ZRCBL
         BULLE = 64D0 * CB1 * XRCBL * YRCBL * ZRCBL
         PRCBL0(1) = CB1   - BULLE
         PRCBL0(2) = XRCBL - BULLE
         PRCBL0(3) = YRCBL - BULLE
         PRCBL0(4) = ZRCBL - BULLE
         PRCBL0(5) = 4D0 * BULLE

C        CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
         V(1) = 0D0
         V(2) = 0D0
         V(3) = 0D0
         DO I=1,4
C           NO DU DL I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUSOTE(NEF1,I) - 1 )
            V(1) = V(1) + PRCBL0(I) * VXYZPN0( NDL+1 )
            V(2) = V(2) + PRCBL0(I) * VXYZPN0( NDL+2 )
            V(3) = V(3) + PRCBL0(I) * VXYZPN0( NDL+3 )
         ENDDO
C        LE NO DU DERNIER DL VITESSE DU BARYCENTRE
         NDL  = NTDLHB + 3 * NEF1
         V(1) = V(1) + PRCBL0(5) * VXYZPN0( NDL-2 )
         V(2) = V(2) + PRCBL0(5) * VXYZPN0( NDL-1 )
         V(3) = V(3) + PRCBL0(5) * VXYZPN0( NDL   )

C        LE POINT (XEBL0,YEBL0,ZEBL0) A ETE TRANSPORTE de (XEBL1,YEBL1,ZEBL1)
C        AU SOUS-PAS DDT
         XEBL1 = XEBL0 - DDT * V(1)
         YEBL1 = YEBL0 - DDT * V(2)
         ZEBL1 = ZEBL0 - DDT * V(3)

C        SAUVEGARDE DU DERNIER EF PARCOURU
         NEF0 = NEF1

C        DANS QUEL ELEMENT FINI NEF1 du MAILLAGE EST CE POINT XYZEBL1?
         CALL XYZDSTE( NBEF,   4,      NUSOTE,
     %                 MOFACE, MXFACE, LFACES,
     %                 XYZSOM, NEF0,   XEBL1, YEBL1, ZEBL1,
     %                 NEF1,   NOFAMX, IERR )
C        IERR=0 PAS D'ERREUR LE POINT XEBL1, YEBL1, ZEBL1 EST DANS L'EF NEF1
C             1 POINT EXTERIEUR AU MAILLAGE
C               PAS de TETRAEDRE derriere LA FACE) NOFAMX de
C               NEF1 LE DERNIER TETRAEDRE DU PARCOURS POUR ATTEINDRE XYZEBL1
C             2 POINT EXTERIEUR AU MAILLAGE (PAS de TETRAEDRE derriere LA FACE)
C             3 FACE DANS 3 TETRAEDRES
C             4 NOMBRE DE TETRAEDRES PARCOURUS TROP GRAND (VITESSE TROP GRANDE?)
C        NEF1 : SI IERR=0 NUMERO DU TETRAEDRE CONTENANT LE POINT XEBL1,YEBL1,ZEBL1
C              SI IERR>0 NUMERO DU DERNIER TETRAEDRE PARCOURU
C        NOFAMX : NUMERO (1 a 4) DE LA FACE EN CAS DE POINT EXTERIEUR AU MAILLAGE

ccc         IF( IERR .NE. 0 ) THEN
ccc            WRITE(IMPRIM,3333) XEBL1, YEBL1, ZEBL1, V
ccc 3333       FORMAT('convecth3: Point HORS frontiere',T38,'XYZ=',
ccc     %      3G15.6,'  V=',3G15.6)
ccc         ENDIF

         IF( IERR .GE. 2 ) THEN
C           VITESSE AU NOEUD NONOEU AVANT SORTIE DES TRIANGLES
            V(1) = VA(1)
            V(2) = VA(2)
            V(3) = VA(3)
            IERR = 0
            GOTO 9900
         ENDIF

         IF( NEF0 .NE. NEF1 ) THEN
C           RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF1
            DO I=1,4
               NS = NUSOTE(NEF1,I)
               DO K=1,3
                  Xe(K,I) = DBLE( XYZSOM( K, NS ) )
               ENDDO
            ENDDO

C           CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF1
            X1 = Xe(1,1)
            Y1 = Xe(2,1)
            Z1 = Xe(3,1)
            DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                       Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                       Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )
            IF( DELTAe .LE. 0D0 ) GOTO 9900
         ENDIF

         IF( IERR .EQ. 1 ) THEN

C           LE POINT EST EXTERIEUR AU DOMAINE ET DU TETRAEDRE NEF1
C           -----------------------------------------------------
C           CALCUL DU POINT D'INTERSECTION DU SEGMENT
C           XEBL1,YEBL1,ZEBL1-XEBL0,YEBL0,ZEBL0 ET DE LA FACE
C           SUR LA FRONTIERE DU TETRAEDRE NEF1
            DO I=1,4
               CALL INARTR( XYZ0, XYZ1,
     %                      Xe(1,NOSOFATE(1,I)),
     %                      Xe(1,NOSOFATE(2,I)),
     %                      Xe(1,NOSOFATE(3,I)),
     %                      LINTER(I), PTI(1,I), CB(1,I) )
C                 LINTER : -2 SI S1=S2
C                       -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE
C                        0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE SES 3 ST
C                        1 SI S1-S2   INTERSECTE     LE TRIANGLE ENTRE SES 3 ST
C                          OU SUR UNE DE SES 3 ARETES
C                 PTI : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=0
C                 CB  : LES 3 COORDONNEES BARYCENTRIQUES DE PT DANS LE TETRA
            ENDDO

C           UNE OU PLUSIEURS FACES D'INTERSECTION
C           RECHERCHE DU POINT D'INTERSECTION LE PLUS PROCHE DE XYZ1
            IMIN = 0
            DMIN = 1D100
            DO I=1,4
               IF( LINTER(I) .EQ. 1 ) THEN
                  D = (XYZ1(1)-PTI(1,I)) ** 2
     %              + (XYZ1(2)-PTI(2,I)) ** 2
     %              + (XYZ1(3)-PTI(3,I)) ** 2
                  IF( D .LT. DMIN ) THEN
                     IMIN = I
                     DMIN = D
                  ENDIF
               ENDIF
            ENDDO

C           LE POINT D'INTERSECTION SUR LA FACE LE PLUS PROCHE DE XYZ1
            IF( IMIN .GT. 0 ) THEN
C              IL EXISTE UN TEL POINT D'INTERSECTION
               XEBL1 = PTI(1,IMIN)
               YEBL1 = PTI(2,IMIN)
               ZEBL1 = PTI(3,IMIN)
ccc            WRITE(IMPRIM,3333) XEBL1, YEBL1, ZEBL1
            ELSE
C              IL N'EXISTE PAS DE POINT PROJETE SUR UNE FACE FRONTIERE
C              => POINT EXTERIEUR AU DOMAINE
C              LE POINT EST EXTERIEUR AU DOMAINE ET DU TETRAEDRE NEF1
C              ------------------------------------------------------
cccC              NO DU DL VITESSE AU SOMMET NOSOMT
ccc               NDL  = NDDLNO( NOSOMT - 1 )
cccC              => VITESSE AU NOEUD NOSOMT NON CONVECTEE
ccc               V(1) = VXYZPN0(NDL+1)
ccc               V(2) = VXYZPN0(NDL+2)
ccc               V(3) = VXYZPN0(NDL+3)
C              VITESSE AU NOEUD NOSOMT AVANT SORTIE DES TRIANGLES
               V(1) = VA(1)
               V(2) = VA(2)
               V(3) = VA(3)
               IERR = 0
               GOTO 9900
            ENDIF

         ENDIF

C        CALCUL DES 3 DERNIERES COORDONNEES BARYCENTRIQUES DU
C        POINT (XEBL1,YEBL1,ZEBL1) DANS LE TETRAEDRE NEF1 QUI LE CONTIENT
C        CE SONT ABSCISSE,ORDONNEE,COTE DU POINT DANS L'EF DE REFERENCE
C        Lambda 2
         XRCBL = DETM33( Xe(1,3)-XEBL1, Xe(1,1)-XEBL1, Xe(1,4)-XEBL1,
     %                   Xe(2,3)-YEBL1, Xe(2,1)-YEBL1, Xe(2,4)-YEBL1,
     %                   Xe(3,3)-ZEBL1, Xe(3,1)-ZEBL1, Xe(3,4)-ZEBL1)
     %         / DELTAe
C        Lambda 3
         YRCBL = DETM33( Xe(1,4)-XEBL1, Xe(1,1)-XEBL1, Xe(1,2)-XEBL1,
     %                   Xe(2,4)-YEBL1, Xe(2,1)-YEBL1, Xe(2,2)-YEBL1,
     %                   Xe(3,4)-ZEBL1, Xe(3,1)-ZEBL1, Xe(3,2)-ZEBL1)
     %         / DELTAe
C        Lambda 4
         ZRCBL = DETM33( Xe(1,1)-XEBL1, Xe(1,3)-XEBL1, Xe(1,2)-XEBL1,
     %                   Xe(2,1)-YEBL1, Xe(2,3)-YEBL1, Xe(2,2)-YEBL1,
     %                   Xe(3,1)-ZEBL1, Xe(3,3)-ZEBL1, Xe(3,2)-ZEBL1)
     %         / DELTAe

 30   CONTINUE
ccc
ccc      print*,'XYZtn+1=',XEBL1,YEBL1,ZEBL1,' XYZtn=',XEBL0,YEBL0,ZEBL0,
ccc     %'  V=',V,' EF=',NEF1,' CB=',XRCBL,YRCBL,ZRCBL

C     calcul de X(tn;tn+1,NOSOMT) effectue:
C     Les Pchapeau(X(tn;tn+1,NOSOMT)) POLYNOMES BREZZI-FORTIN
C     DANS LE TETRAEDRE NEF1 contenant (XEBL1,YEBL1,ZEBL1) ont ete calcules
C     =====================================================================
C     CALCUL DES Pchapeau(XYRCBL) P1+BULLE
      CB1   = 1.D0 - XRCBL - YRCBL - ZRCBL
      BULLE = 64D0 * CB1 * XRCBL * YRCBL * ZRCBL
      PRCBL0(1) = CB1   - BULLE
      PRCBL0(2) = XRCBL - BULLE
      PRCBL0(3) = YRCBL - BULLE
      PRCBL0(4) = ZRCBL - BULLE
      PRCBL0(5) = 4D0 * BULLE

C     CALCUL DE LA VITESSE V(tn,(XRCBL,YRCBL)) dans l'EF NEF1
      V(1) = 0D0
      V(2) = 0D0
      V(3) = 0D0
      DO I=1,4
C        NO DU DL I DANS L'EF NEF1 DE LA VITESSE AU TEMPS tn
C        AU POINT X(tn;tn+1,NOSOMT)
         NDL  = NDDLNO( NUSOTE(NEF1,I) - 1 )
         V(1) = V(1) +  PRCBL0(I) * VXYZPN0( NDL+1 )
         V(2) = V(2) +  PRCBL0(I) * VXYZPN0( NDL+2 )
         V(3) = V(3) +  PRCBL0(I) * VXYZPN0( NDL+3 )
      ENDDO
C     LE NO DU DERNIER DL VITESSE DU BARYCENTRE
      NDL  = NTDLHB + 3 * NEF1
      V(1) = V(1) +  PRCBL0(5) * VXYZPN0( NDL-2 )
      V(2) = V(2) +  PRCBL0(5) * VXYZPN0( NDL-1 )
      V(3) = V(3) +  PRCBL0(5) * VXYZPN0( NDL   )
C
C     PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
 9900 D = SQRT( V(1)**2 + V(2)**2 + V(3)**2 ) / VITMAX0
      IF( D .GT. 1.001D0 ) THEN
         V(1) = V(1) / D
         V(2) = V(2) / D
         V(3) = V(3) / D
ccc         print*,'convecbf3: nosomt=',nosomt,' d=',d,' V=',V,
ccc     %          ' XRCBL=',XRCBL,' YRCBL=',YRCBL,' ZRCBL=',ZRCBL
      ENDIF

      RETURN
      END
