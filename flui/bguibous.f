      SUBROUTINE BGUIBOUS( tn, tn1, Rho, G, CoGrPr, CoBOUS,
     %                     NDIM,   NBSOM,  NBNOVI, NTDLVI,
     %                     NDDLNO, NUNOSO, XYZNOE,
     %                     NUTYEL, NBNOEF, NBELEM, NUNOEF,
     %                     NBVVEF, NBSFEF, NBLAEF,
     %                     NUVVEF, NUSFEF, NULAEF,
     %                     NUMILI, NUMALI, LTDELI,
     %                     NUMISU, NUMASU, LTDESU,
     %                     NUMIVO, NUMAVO, LTDEVO,
     %                     MOARET, MXARET, MNLARE,
     %                     MOFACE, MXFACE, LFACES,
     %                     VXYZPNtn, VITMAXtn, PRESStn,
     %                     TEMPERt0, TEMPERtn1,
     %                     P2P22D, P2P23D,
     %                     BG,     NBCHTL, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG du SYSTEME
C ---- VG U(tn+1) = (Rho - dt Mhu Laplacien ) U(tn+1) =
C      = Rho U(tn,X(tn;tn+1,x)) - dt CoGrPr GRAD P(tn) + dt Force(tn+1)
C      CALCUL DU SECOND MEMBRE DU SYSTEME LINEAIRE:
C      Integrale Rho tP2 Ui(tn,X(tn;tn+1,bl)) dX
C    - Integrale tP2 dt CoGrPr GRAD P(tn) dX
C    + Integrale tP2 dt Force(tn+1) dX
C    - Integrale tP2 P2 dx dt Rho G( Z + CoBOUS(Temper(tn+1)-Temper(t0)) )
C      La PESANTEUR TERRESTRE selon l'axe -Z (-Y en 2d) est AINSI PRISE
C      EN COMPTE car cette FORCE FAIT APPARAITRE la CONVECTION NATURELLE
C      et Force(tn+1) = - Rho G Z  est la FORCE POTENTIELLE

C ENTREES:
C --------
C tn, tn1: INTERVALLE DE TEMPS tn tn+1  => DT PAS du TEMPS = tn1-tn
C DT     : PAS CONSTANT DU TEMPS
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX*
C G      : ACCELERATION DE LA PSEANTEUR TERRESTRE
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C CoBOUS : COEFFICIENT dans l'APPROXIMATION DE BOUSSINESQ
C          Rho0 - Rho = Rho CoBOUS (Temper-Temper0)
C          COEFFICIENT DE PROPORTIONNALITE DE LA DIFFERENCE DE
C          TEMPERATURE A LA VARIATION DE LA MASSE VOLUMIQUE
C          DOIT ENCORE le COEFFICIENT DE DILATATION THERMIQUE du FLUIDE
C                    ( the THERMAL EXPANSION COEFFICIENT )

C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
C NTDLVI : NOMBRE TOTAL DE DL de la VITESSE (=NDIM*NBNOVI)
C XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
C NDDLNO : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
C          A UTILISER POUR RETROUVER LES DL PRESSION DES VECTEURS VXYZPNtn
C NUNOSO : NUNOSO(I) = NUMERO DU SOMMET (1 A NBSOM) DU NOEUD GLOBAL I
C NUTYEL : NUMERO DU TYPE D'EF ( 15 TRIANGLE, 20 TETRAEDRE TAYLOR-HOOD )
C NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
C          ( 6 POUR LE TRIANGLE, 10 POUR LE TETRAEDRE)
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF

C NBVVEF : NOMBRE DE VOLUME  POUR 1 EF DE CE TYPE
C NBSFEF : NOMBRE DE FACES   POUR 1 EF DE CE TYPE
C NBLAEF : NOMBRE D'ARETES   POUR 1 EF DE CE TYPE
C NUVVEF : TABLEAU DU NUMERO DE VOLUME DE CHAQUE ELEMENT FINI
C NUSFEF : TABLEAU DU NUMERO DE SURFACE DES FACES   DE CHAQUE EF
C NULAEF : TABLEAU DU NUMERO DE LIGNE   DES ARETES  DE CHAQUE EF
C
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS LIGNES
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

C MOARET : NOMBRE DE MOTS POUR LES INFO D'UNE ARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DECLAREES
C MNLARE : EN 2D SEULEMENT ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE

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

C VXYZPNtn : VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLVP) au TEMPS tn
C VITMAXtn : NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT
C PRESStn  : DL DE LA PRESSION AUX SOMMETS DU MAILLAGE
C TEMPERt0 : DL DE LA TEMPERATURE INITIALE AUX NOEUDS DU MAILLAGE P2
C TEMPERtn1: DL DE LA TEMPERATURE AU TEMPS tn AUX NOEUDS DU MAILLAGE P2

C P2P22D : P2P22D(i,j) = integrale P2i P2j dX SUR LE TRIANGLE P2
C P2P23D : P2P23D(i,j) = integrale P2i P2j dX SUR LE TETRAEDRE P2
C
C SORTIES:
C --------
C BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
C NBCHTL : NOMBRE DE REMONTEES DE LA CARACTERISTIQUE TROP LONGUES
C IERR   : 0 PAS DE PROBLEME, >0 PROBLEME RENCONTRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray               Mai 2022
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2023
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      REAL              DT, tn, tn1, XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  Rho, G, CoGrPr, CoBOUS, BG(NTDLVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NDDLNO(0:NBNOVI),
     %                  NUNOSO(NBNOVI)
      INTEGER           NUMILI,NUMALI,LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           NUMISU,NUMASU,LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NUMIVO,NUMAVO,LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      INTEGER           NUVVEF(NBVVEF,NBELEM), NUSFEF(NBSFEF,NBELEM),
     %                  NULAEF(NBLAEF,NBELEM)
      INTEGER           NDIM,   NTDLVI, NBSOM,
     %                  NBVVEF, NBSFEF, NBLAEF,
     %                  NUTYEL, NBELEM, NBNOEF,
     %                  MOARET, MXARET, MNLARE,
     %                  MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                  NBNOVI, NBCHTL, IERR, NUELEM,
     %                  M, K, MEK, NSK, TempouVit
      INTEGER           NONOEF(10)

      DOUBLE PRECISION  DELTAT
      DOUBLE PRECISION  BE1(34), BE2(30)
      DOUBLE PRECISION  PRESStn(NBSOM), VXYZPNtn(1:*), VITMAXtn
      DOUBLE PRECISION  TEMPERt0(NBNOVI), TEMPERtn1(NBNOVI)

C     INTEGRALES P2 P2 SUR LE TRIANGLE P2 de REFERENCE
      DOUBLE PRECISION  P2P22D(6,6)
C                       P2P22D(i,j) = integrale P2i P2j dX

C     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
      DOUBLE PRECISION  P2P23D(10,10)
C                       P2P23D(i,j) = integrale P2i P2j dX

ccc      DOUBLE PRECISION  t_cpu_0, t_cpu_1
ccc      call cpu_time(t_cpu_0)
ccc      print*,'ENTREE de bguibous.f sans OpenMP'

C     LE PAS DE TEMPS EN DOUBLE PRECISION
      DT     = tn1 - tn
      DELTAT = DT
      TempouVit = 1

C     MISE A ZERO DE BG
      CALL AZEROD( NTDLVI, BG )

C     BOUCLE SUR LES EF de TAYLOR-HOOD
      DO 10 NUELEM = 1, NBELEM

C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO

         IF( NUTYEL .EQ. 15 ) THEN

C           TRIANGLE TAYLOR_HOOD
C           1) Integrale Rho tP2 ui(tn,X(tn;tn+1,x)) dx
C           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
C           avec 6+6 COEFFICIENTS
            CALL F2SBOP2P1( NUELEM,  1,       NBELEM, NUNOEF,
     %                      MOARET,  MXARET,  MNLARE,
     %                      DT,      Rho,     XYZNOE, NDDLNO,
     %                      TempouVit, TEMPERtn1, VXYZPNtn, VITMAXtn,
     %                      BE1,     IERR )
            IF( IERR .GT. 1 ) THEN
C              ABANDON DE L'EF NUELEM
               GOTO 10
            ENDIF

C           2) -Integrale tP2 dt CoGrPr GRAD P(tn) + dt Integrale tP2 Force(tn+1)
            CALL F2EX2P2P1( NBNOVI, XYZNOE, NONOEF, NUNOSO,
     %                      NULAEF(1,NUELEM), NUMILI, NUMALI, LTDELI,
     %                      NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU,
     %                      DELTAT, CoGrPr, P2P22D,   NBSOM,  PRESStn,
     %                      BE2 )

C           3) -Integrale tP2 P2 dx dt Rho (-Y) G CoBOUS(Temper(tn+1)-Temper(t0))
C               POUR LA PESANTEUR SELON l'axe -Y  SEULEMENT
            CALL F2BOUP2P1( XYZNOE, NONOEF, P2P22D,
     %                      DELTAT, Rho, G, CoBOUS, TEMPERtn1, TEMPERt0,
     %                      BE2 )

         ELSE

C           TETRAEDRE TAYLOR-HOOD
C           1) Integrale Rho tP2 ui(tn,X(tn;tn+1,x) dx
C           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
C           avec 10+10+10 COEFFICIENTS
            CALL F3SBOP2P1( NUELEM,  NBELEM,  NUNOEF,
     %                      MOFACE,  MXFACE,  LFACES,
     %                      DT,      Rho,     XYZNOE, NDDLNO,
     %                      TempouVit, TEMPERtn1, VXYZPNtn, VITMAXtn,
     %                      BE1,     IERR )
            IF( IERR .GT. 1 ) THEN
C              ABANDON DE L'EF NUELEM
               GOTO 10
            ENDIF

C           2) -Integrale tP2 dt CoGrPr GRAD P(tn) + dt Integrale tP2 Force(tn+1)
            CALL F3EX2P2P1( NBNOVI, XYZNOE, NONOEF, NUNOSO,
     %                      NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU,
     %                      NUVVEF(1,NUELEM), NUMIVO, NUMAVO, LTDEVO,
     %                      DELTAT*CoGrPr,DELTAT, P2P23D, NBSOM,PRESStn,
     %                      BE2 )

C           3) -Integrale tP2 P2 dx dt Rho (-Z) G CoBOUS (Temper(tn+1)-Temper(t0))
C               POUR LA PESANTEUR SELON l'axe -Z  SEULEMENT
            CALL F3BOUP2P1( XYZNOE, NONOEF, P2P23D,
     %                      DELTAT, Rho, G, CoBOUS, TEMPERtn1, TEMPERt0,
     %                      BE2 )

         ENDIF
         IF( IERR .EQ. 4 ) THEN
            NBCHTL = NBCHTL + 1
            IERR = 0
C           ABANDON DE L'EF NUELEM
            GOTO 10
         ENDIF

ccc         IF( NUELEM .EQ. 5 ) THEN
ccc            call affvect( 'NEF=5 BE1 X(tn;tn+1,x)=', NDIM*NBNOEF, BE1 )
ccc            call affvect( 'NEF=5 BE2=', NDIM*NBNOEF, BE2 )
ccc         ENDIF

C        ASSEMBLAGE DE BE1 + BE2 dans BG = VITX1, VITY1, VITZ1
C        Rho U(tn,X(tn+1,x)) - dt CoGrPr GRAD P(tn) + dt FOmega(tn+1)
ccc         CALL ASMEEX2( NDIM, NBNOEF, NONOEF, BE1, BE2, NBNOVI,  BG )
         DO M=1,NBNOEF

C           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            NSK = NONOEF(M)
            MEK = M

C           ASSEMBLAGE DES 2 COEFFICIENTS ELEMENTAIRES DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               BG( NSK ) = BG( NSK ) + BE1( MEK ) + BE2( MEK )
               NSK = NSK + NBNOVI
               MEK = MEK + NBNOEF
            ENDDO

         ENDDO

C     FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
 10   ENDDO

ccc      call affvect( 'BGUIBOUS VECTEUR BG', 20, BG )
      call afl1ve(  'BGUIBOUS VECTEUR BG', NDIM*NBNOVI, BG )

cccC     Le TEMPS CPU donne par le systeme en fin d'execution de bguibous.f
ccc      call cpu_time(t_cpu_1)
ccc      print*,'SORTIE de bguibous.f sans OpenMP en',
ccc     %        t_cpu_1-t_cpu_0,' CPU secondes'

      IERR = 0
      RETURN
      END
