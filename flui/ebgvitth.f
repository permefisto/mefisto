      SUBROUTINE EBGVITTH( DT,     Rho,    CoGrPr, Omega0, Omega,
     %                     NBSOM,  NBNOVI, NUDDL,  NUNOSO, XYZNOE,
     %                     NBELEM, NUNOEF,
     %                     NBVVEF, NBSFEF, NUVVEF, NUSFEF,
     %                     NUMISU, NUMASU, LTDESU,
     %                     NUMIVO, NUMAVO, LTDEVO,
     %                     MOFACE, MXFACE, LFACES,
     %                     Wtn,    VITMAX, PRESS0,
     %                     NTDLVI, BG,     NBCHTL, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
C ---- VG W* = (Rho - dt Mhu Laplacien ) W* =
C      = Rho W(tn,X(tn+1,x) - dt CoGrPr GRAD P(tn) + dt Force(tn+1) + dt Grad(-Rho g z)
C      + Rho dt [ - 2 Omega(tn+1) x W(tn)
C                 - ( Omega(tn+1)-Omega(tn) )/dt  x r
C                 -   Omega(tn+1) x ( Omega(tn+1) x r ) ]
C
C      CALCUL DU SECOND MEMBRE DU SYSTEME LINEAIRE:
C      Integrale Rho tP2 Ui(tn,X(tn;tn+1,bl) dX
C    - Integrale tP2 dt CoGrPr GRAD P(tn) dX + dt Integrale tP2 Force(tn) dX
C    + Integrale tP2 Rho dt [ - 2 Omega(tn+1) x W(tn)
C                             - ( Omega(tn+1)-Omega(tn) )/dt  x r
C                             -   Omega(tn+1) x ( Omega(tn+1) x r ) ] dX
C      Pour des TETRAEDRES TAYLOR-HOOD
C
C ENTREES:
C --------
C DT     : PAS CONSTANT DU TEMPS
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C Omega0 : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE a tn
C Omega  : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE a tn+1
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
C XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
C NUDDL  : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
C          A UTILISER POUR RETROUVER LES DL PRESSION DES VECTEURS Wtn
C NUNOSO : NUNOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NUNOEF : NUNOEF(NBELEM,10) NO DES NOEUDS DES NBELEM EF

C NBVVEF : NOMBRE DE VOLUME  POUR 1 EF DE CE TYPE
C NBSFEF : NOMBRE DE FACES   POUR 1 EF DE CE TYPE
C NUVVEF : TABLEAU DU NUMERO DE VOLUME DE CHAQUE ELEMENT FINI
C NUSFEF : TABLEAU DU NUMERO DE SURFACE DES FACES DE CHAQUE EF
C
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS VOLUMES

C MOFACE : NOMBRE DE MOTS POUR LES INFO D'UNE FACE
C MXFACE : NOMBRE MAXIMAL D'FACES DECLAREES
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

C Wtn    : VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDL) au TEMPS tn
C VITMAX : NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C NTDLVI : NOMBRE TOTAL DE DL VITESSE (3*NBNOVI)
C
C SORTIES:
C --------
C BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
C NBCHTL : NOMBRE DE REMONTEES DE LA CARACTERISTIQUE TROP LONGUES
C IERR   : 0 PAS DE PROBLEME, >0 PROBLEME RENCONTRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2013
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      include"./incl/p2p23d.inc"
!     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
!     DOUBLE PRECISION  P2P23D(10,10)
!                       P2P23D(i,j) = integrale P2i P2j dX

      REAL              DT
      DOUBLE PRECISION  Rho, CoGrPr, Omega0(3), Omega(3), BG(NTDLVI)
      INTEGER           NUNOEF(NBELEM,10), NUDDL(0:NBNOVI),
     %                  NUNOSO(NBNOVI)
      INTEGER           NUMISU,NUMASU,LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NUMIVO,NUMAVO,LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      INTEGER           NUVVEF(NBVVEF,NBELEM), NUSFEF(NBSFEF,NBELEM)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NTDLVI, NBSOM, NBVVEF, NBSFEF, NBELEM,
     %                  MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                  NBNOVI, NBCHTL, IERR, NUELEM,
     %                  M, K, MEK, NSK
      DOUBLE PRECISION  DELTAT

      INTEGER           NONOTE(10)
      DOUBLE PRECISION  BE1(34), BE2(30)
      DOUBLE PRECISION  PRESS0(NBSOM), Wtn(1:*), VITMAX


C     LE PAS DE TEMPS EN DOUBLE PRECISION
      DELTAT = DT

C     MISE A ZERO DE BG
      CALL AZEROD( NTDLVI, BG )

C     BOUCLE SUR LES EF de TAYLOR-HOOD
      DO NUELEM = 1, NBELEM

C        NO DES NOEUDS DU TETRAEDRE TAYLOR-HOOD NUELEM
         DO M = 1, 10
            NONOTE(M) = NUNOEF(NUELEM,M)
         ENDDO

C        1) Integrale Rho tP2 ui(tn,X(tn;tn+1,x) dx
C        SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
C        avec 10+10+10 COEFFICIENTS
         CALL F3SNVP2P1( NUELEM, NBELEM, NUNOEF,
     %                   MOFACE, MXFACE, LFACES,
     %                   DT,     Rho,    XYZNOE, NUDDL,
     %                   Wtn,    VITMAX,
     %                   BE1,    IERR )
         IF( IERR .GT. 1 ) GOTO 9999

C        2) -Integrale tP2 dt CoGrPr GRAD P(tn) + dt Integrale tP2 Force(tn)
         CALL F3EX2P2P1( NBNOVI, XYZNOE, NONOTE, NUNOSO,
     %                   NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU,
     %                   NUVVEF(1,NUELEM), NUMIVO, NUMAVO, LTDEVO,
     %                   DELTAT*CoGrPr,  DELTAT, P2P23D,
     %                   NBSOM,  PRESS0,
     %                   BE2 )

C        3)Integrale tP2 Rho dt [ - 2 Omega(tn+1) x W(tn)
C                                 - ( Omega(tn+1)-Omega(tn) )/dt  x r
C                                 -   Omega(tn+1) x ( Omega(tn+1) x r ) ] dX
         CALL F3EBVP2P1( Rho,    DELTAT, Omega0, Omega,
     %                   NBNOVI, XYZNOE, NONOTE, NUDDL, Wtn,
     %                   BE2 )

         IF( IERR .EQ. 4 ) THEN
            NBCHTL = NBCHTL + 1
            IERR = 0
         ENDIF

ccc         IF( NUELEM .EQ. 5 ) THEN
ccc            call affvect( 'NEF=5 BE1 X(tn;tn+1,x)=', 3*10, BE1 )
ccc            call affvect( 'NEF=5 BE2=', 3*10, BE2 )
ccc         ENDIF

C        ASSEMBLAGE DE BE1 + BE2 dans BG = VITX1, VITY1, VITZ1
C        Rho U(tn,X(tn+1,x)) - dt CoGrPr GRAD P(tn) + dt FOmega(tn+1)
ccc         CALL ASMEEX2( 3, 10, NONOTE, BE1, BE2, NBNOVI,  BG )
         DO M=1,10

!           NUMERO DU M-EME NOEUD DU TETRAEDRE TAYLOR-HOOD
            NSK = NONOTE(M)
            MEK = M

!           ASSEMBLAGE DES 2 COEFFICIENTS ELEMENTAIRES DANS LE COEFFICIENT GLOBAL
            DO K=1,3
               BG( NSK ) = BG( NSK ) + BE1( MEK ) + BE2( MEK )
               NSK = NSK + NBNOVI
               MEK = MEK + 10
            ENDDO

         ENDDO

C     FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
      ENDDO

ccc      call affvect( 'EBGVITTH VECTEUR BG', 20, BG )
ccc      call afl1ve(  'EBGVITTH VECTEUR BG', 3*NBNOVI, BG )

 9999 RETURN
      END
