      SUBROUTINE BGGRATH( DELTAT, CoGrPr, NDIM,   NBSOM, NBNOVI, NTDLVI,
     %                    NUNOSO, XYZNOE, NUTYEL, NBNOEF, NBELEM,NUNOEF,
     %                    PRESP1,   BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
C ---- (Rho - dt Mhu Laplacien){u(tn+1)-u*} = -dt CoGrPr GRAD(p(tn+1)-p(tn))
C      BG=Integrale sur e  -dt CoGrPr tV(x) Grad PRESP1(x) dx
C
C ENTREES:
C --------
C DELTAT : PAS CONSTANT DU TEMPS DOUBLE PRECISION
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
C NTDLVI : NOMBRE TOTAL DE DL VITESSE (NDIM*NBNOVI)
C XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE P2
C NUNOSO : NUNOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C NUTYEL : NUMERO DU TYPE D'EF ( 15 TRIANGLE, 20 TETRAEDRE TAYLOR-HOOD )
C NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
C          ( 6 POUR LE TRIANGLE, 10 POUR LE TETRAEDRE)
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF
C PRESP1 : DL DE LA PRESSION AUX SOMMETS DU MAILLAGE
C
C SORTIE :
C --------
C BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
C MODIF  : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2022
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      include"./incl/a___npef.inc"

      DOUBLE PRECISION  DELTAT, CoGrPr, BG(NTDLVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NUNOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NUTYEL, NDIM, NBNOEF, NBELEM, NTDLVI,
     %                  NBSOM, NBNOVI
      INTEGER           M, NUELEM, K, NSK, MEK, NONOEF(10)
      DOUBLE PRECISION  BE(30)
      DOUBLE PRECISION  PRESP1(NBSOM)

C     MISE A ZERO DE BG
      CALL AZEROD( NTDLVI, BG )

C     BOUCLE SUR LES EF de TAYLOR-HOOD P2 P1
      DO NUELEM = 1, NBELEM

C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM P2
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO

C       (Integrale sur e  - dt CoGrPr     tV(x) Grad (p(tn+1,x)- p(tn,x)) dx=)
C        Integrale sur e  - dt CoGrPr tDiv V(x)      (p(tn+1,x)- p(tn,x)) dx
         IF( NUTYEL .EQ. 15 ) THEN

C           TRIANGLE TAYLOR-HOOD

Cccc           Integrale sur e  - dt CoGrPr tDiv V(x) PRESP1(x) dx
ccc            CALL F2EX4P2P1( DELTAT*CoGrPr, NBNOVI,XYZNOE,NONOEF, NUNOSO,
ccc     %                      NBSOM, PRESP1,  BE )

C           Integrale sur e  - dt CoGrPr  tV(x) Grad PRESP1(x) dx=
            CALL F2EX4P2P1GRAD( XYZNOE, NONOEF, NUNOSO, DELTAT, CoGrPr,
     %                          NBNOVI, NBSOM,  PRESP1,   BE )

         ELSE

C           TETRAEDRE TAYLOR-HOOD

Cccc           Integrale sur e  - dt CoGrPr tDiv V(x) PRESP1(x) dx
ccc            CALL F3EX4P2P1( DELTAT*CoGrPr, NBNOVI,XYZNOE,NONOEF, NUNOSO,
ccc     %                      NBSOM, PRESP1,  BE )

C           Integrale sur e  - dt CoGrPr  tV(x) Grad PRESP1(x) dx=
            CALL F3EX4P2P1GRAD( XYZNOE, NONOEF, NUNOSO, DELTAT, CoGrPr,
     %                          NBNOVI, NBSOM,  PRESP1,   BE )

         ENDIF

ccc         IF( NUELEM .EQ. 5 ) THEN
ccc            call affvect(
ccc     %      'BE: - dt CoGrPr tV(x) Grad (p(tn+1,x)- p(tn,x))',
ccc     %       NDIM*NBNOEF, BE )
ccc         ENDIF

C        ASSEMBLAGE DE BE dans BG
ccc         CALL ASMEEX1( NDIM, NBNOEF, NONOEF, BE, NBNOVI, BG )

         DO M=1,NBNOEF
C           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            NSK = NONOEF(M)
            MEK = M

C           ASSEMBLAGE DU COEFFICIENT ELEMENTAIRE DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               BG( NSK ) = BG( NSK ) + BE( MEK )
               NSK = NSK + NBNOVI
               MEK = MEK + NBNOEF
            ENDDO
         ENDDO


C     FIN DE LA BOUCLE SUR LES EF DE LA CONSTRUCTION DE BG
      ENDDO

ccc      call affvect( 'BGGRATH BG=', 20, BG )
C     FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG

      RETURN
      END
