SUBROUTINE BGVITTH( tn, tn1, Rho, CoGrPr, NDIM, NBSOM, NBNOVI, &
                    NDDLNO, NUNOSO, XYZNOE, &
                    NUTYEL, NBNOEF, NBELEM, NUNOEF, &
                    NBVVEF, NBSFEF, NBLAEF,  &
                    NUVVEF, NUSFEF, NULAEF,  &
                    NUMILI, NUMALI, LTDELI,  &
                    NUMISU, NUMASU, LTDESU,  &
                    NUMIVO, NUMAVO, LTDEVO,  &
                    MOARET, MXARET, MNLARE,  &
                    MOFACE, MXFACE, LFACES,  &
                    VXYZPNtn, VITMXtn, PRESStn,&
                    P2P22D,   P2P23D, &
                    NTDLVI, BG, NBCHTL, IERR )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
! ---- VG U(tn+1) = (Rho - dt Mhu Laplacien ) U(tn+1) =
!      = Rho U(tn,X(tn+1,x) - dt CoGrPr GRAD P(tn) + dt Force(tn+1)

!      CALCUL DU SECOND MEMBRE DU SYSTEME LINEAIRE:
!      Integrale Rho tP2 Ui(tn,X(tn;tn+1,bl) dX
!    - Integrale tP2 dt CoGrPr GRAD P(tn) dX + dt Integrale tP2 Force(tn) dX
!
! ENTREES:
! --------
! tn, tn1: INTERVALLE DE TEMPS tn tn+1  => DT PAS du TEMPS = tn1-tn
! DT     : PAS CONSTANT DU TEMPS
! Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
! CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
! NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NDDLNO : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
!          A UTILISER POUR RETROUVER LES DL PRESSION DES VECTEURS VXYZPNtn
! NUNOSO : NUNOSO(I) = NUMERO DU SOMMET (1 A NBSOM) DU NOEUD GLOBAL I
! NUTYEL : NUMERO DU TYPE D'EF ( 15 TRIANGLE, 20 TETRAEDRE TAYLOR-HOOD )
! NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
!          ( 6 POUR LE TRIANGLE, 10 POUR LE TETRAEDRE)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF

! NBVVEF : NOMBRE DE VOLUME  POUR 1 EF DE CE TYPE
! NBSFEF : NOMBRE DE FACES   POUR 1 EF DE CE TYPE
! NBLAEF : NOMBRE D'ARETES   POUR 1 EF DE CE TYPE
! NUVVEF : TABLEAU DU NUMERO DE VOLUME DE CHAQUE ELEMENT FINI
! NUSFEF : TABLEAU DU NUMERO DE SURFACE DES FACES   DE CHAQUE EF
! NULAEF : TABLEAU DU NUMERO DE LIGNE   DES ARETES  DE CHAQUE EF
!
! NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
! NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
! LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS LIGNES
!
! NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
! NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
! LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS SURFACES
!
! NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
! NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
! LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS VOLUMES

! MOARET : NOMBRE DE MOTS POUR LES INFO D'UNE ARETE
! MXARET : NOMBRE MAXIMAL D'ARETES DECLAREES
! MNLARE : EN 2D SEULEMENT ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE

! MOFACE : NOMBRE DE MOTS POUR LES INFO D'UNE FACE
! MXFACE : NOMBRE MAXIMAL D'FACES DECLAREES
! LFACES : TABLEAU DES FACES DU MAILLAGE  (cf hachag.f)
!          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
!          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
!          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
!          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
!                       0 SI TRIANGLE
!          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
!                       1 CUBE EST UN TETRA ou PENTA ou HEXAEDRE
!          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
!                       0 SI PAS DE 1-ER  CUBE
!          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
!                       ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
!          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
!          LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE

! VXYZPNtn: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDL) au TEMPS tn
! VITMXtn : NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT
! PRESStn : DL DE LA PRESSION AUX SOMMETS DU MAILLAGE
! P2P22D : P2P22D(i,j) = integrale P2i P2j dX SUR LE TRIANGLE
! P2P23D : P2P23D(i,j) = integrale P2i P2j dX SUR LE TETRAEDRE
! NTDLVI : NOMBRE TOTAL DE DL VITESSE (NDIM*NBNOVI)
!
! SORTIES:
! --------
! BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
! NBCHTL : NOMBRE DE REMONTEES DE LA CARACTERISTIQUE TROP LONGUES
! IERR   : 0 PAS DE PROBLEME, >0 PROBLEME RENCONTRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
   IMPLICIT NONE
      include"./incl/langue.inc"
      include "./incl/threads.inc"
      include "./incl/donflu.inc95"
      REAL              DT, tn, tn1, XYZNOE(3,NBNOVI)

      DOUBLE PRECISION  Rho, CoGrPr, BG(NTDLVI)
      DOUBLE PRECISION  BE1(34), BE2(30), D
      DOUBLE PRECISION  PRESStn(NBSOM), VXYZPNtn(1:*), VITMXtn
      DOUBLE PRECISION  DELTAT

      INTEGER           NUNOEF(NBELEM,NBNOEF), NDDLNO(0:NBNOVI)
      INTEGER           NUNOSO(NBNOVI)
      INTEGER           NUMILI,NUMALI,LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           NUMISU,NUMASU,LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NUMIVO,NUMAVO,LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      INTEGER           NUVVEF(NBVVEF,NBELEM)
      INTEGER           NUSFEF(NBSFEF,NBELEM)
      INTEGER           NULAEF(NBLAEF,NBELEM)
      INTEGER           NDIM,   NTDLVI, NBSOM,  &
                        NBVVEF, NBSFEF, NBLAEF, &
                        NUTYEL, NBELEM, NBNOEF, &
                        MOARET, MXARET, MNLARE, &
                        MOFACE, MXFACE, LFACES(MOFACE,MXFACE), &
                        NBNOVI, NBCHTL, IERR, IER
      INTEGER           NUELEM, NONOEF(10), M, K, MEK, NSK

!     INTEGRALES  P2 P2 SUR LE TRIANGLE P2 de REFERENCE
      DOUBLE PRECISION  P2P22D(6,6)
!                       P2P22D(i,j) = integrale P2i P2j dX

!     INTEGRALES  P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
      DOUBLE PRECISION  P2P23D(10,10)
!                       P2P23D(i,j) = integrale P2i P2j dX
!!!      DOUBLE PRECISION  t_cpu_0, t_cpu_1

!!!      call cpu_time(t_cpu_0)
!!!      print*,'ENTREE de bgvitth.f95 avec 1 ATOMIC sur BG'

!     LE PAS DE TEMPS EN DOUBLE PRECISION
      DT     = tn1 - tn
      DELTAT = DT

!     MISE A ZERO DE BG
      CALL AZEROD( NTDLVI, BG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NUELEM,M,K,MEK,NSK,NONOEF,BE1,BE2,D,IER)
!     BOUCLE SUR LES EF de TAYLOR-HOOD
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO

         IF( NUTYEL .EQ. 15 ) THEN

!           TRIANGLE TAYLOR_HOOD
!           1) Integrale Rho tP2 ui(tn,X(tn;tn+1,x) dx
!           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
!           avec 6+6 COEFFICIENTS
            CALL F2SNVP2P1( NUELEM,  1,      NBELEM, NUNOEF,&
                            MOARET,  MXARET, MNLARE,&
                            DT,      Rho,    XYZNOE, NDDLNO,&
                            VXYZPNtn, VITMXtn,&
                            BE1,     IER )
!!!            IF( IER .GT. 1 ) EXIT

!           2) -Integrale tP2 dt CoGrPr GRAD P(tn) + dt Integrale tP2 Force(tn)
            CALL F2EX2P2P1( NBNOVI, XYZNOE,   NONOEF, NUNOSO,&
                            NULAEF(1,NUELEM), NUMILI, NUMALI, LTDELI,&
                            NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU,&
                            DELTAT, CoGrPr, P2P22D,&
                            NBSOM,  PRESStn,&
                            BE2 )

         ELSE

!           TETRAEDRE TAYLOR-HOOD
!           1) Integrale Rho tP2 ui(tn,X(tn;tn+1,x) dx
!           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
!           avec 10+10+10 COEFFICIENTS
            CALL F3SNVP2P1( NUELEM,   NBELEM, NUNOEF,&
                            MOFACE,   MXFACE, LFACES,&
                            DT,       Rho,    XYZNOE, NDDLNO,&
                            VXYZPNtn, VITMXtn,&
                            BE1,     IER )
!!!            IF( IER .GT. 1 ) EXIT

!           2) -Integrale tP2 dt CoGrPr GRAD P(tn) + dt Integrale tP2 Force(tn)
            CALL F3EX2P2P1( NBNOVI, XYZNOE,   NONOEF, NUNOSO,&
                            NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU,&
                            NUVVEF(1,NUELEM), NUMIVO, NUMAVO, LTDEVO,&
                            DELTAT*CoGrPr,  DELTAT, P2P23D,&
                            NBSOM,  PRESStn,&
                            BE2 )

         ENDIF


         IF( IER .EQ. 4 ) THEN
!           SOMME DES CARACTERISTIQUES TROP LONGUES
!$OMP ATOMIC
            NBCHTL = NBCHTL + 1
         ENDIF

         DO M=1,NBNOEF

!           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            NSK = NONOEF( M )
            MEK = M

!           ASSEMBLAGE DES 2 COEFFICIENTS ELEMENTAIRES DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               D = BE1( MEK ) + BE2( MEK )
!$OMP ATOMIC
               BG( NSK ) = BG( NSK ) + D

               NSK = NSK + NBNOVI
               MEK = MEK + NBNOEF
            ENDDO

         ENDDO

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!      call affvect( 'BGVITTH VECTEUR BG', 20,     BG )
!!!      call afl1ve(  'BGVITTH VECTEUR BG', NTDLVI, BG )

!!!!     Le TEMPS CPU donne par le systeme en fin d'execution de bgvitth.f95
!!!      call cpu_time(t_cpu_1)
 
!!!      IF( LANGAG .EQ. 0 ) THEN
!!!         print*,'SORTIE de bgvitth.f95 avec 1 ATOMIC sur BG en',&
!!!                 t_cpu_1-t_cpu_0,' CPU secondes'
!!!      ELSE
!!!         print*,'EXIT bgvitth.f95 with 1 ATOMIC on BG. CPU seconds=',&
!!!                 t_cpu_1-t_cpu_0
!!!      ENDIF

      IERR = 0
      RETURN
END SUBROUTINE BGVITTH
