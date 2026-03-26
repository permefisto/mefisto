SUBROUTINE BGGRATH( DELTAT, CoGrPr, NDIM,   NBSOM,  NBNOVI, NTDLVI, &
                    NUNOSO, XYZNOE, NUTYEL, NBNOEF, NBELEM, NUNOEF, &
                    PRESP1,   BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
! ---- ( Rho - dt Mhu Laplacien ){u(tn+1)-u*} = -dt CoGrPr GRAD(p(tn+1)-p(tn))
!  BG=(Integrale sur e  - dt CoGrPr     tV(x) Grad (p(tn+1,x)- p(tn,x)) dx=)

! ENTREES:
! --------
! DELTAT : PAS CONSTANT DU TEMPS DOUBLE PRECISION
! CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
! NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! NTDLVI : NOMBRE TOTAL DE DL VITESSE (NDIM*NBNOVI)
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NUNOSO : NUNOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! NUTYEL : NUMERO DU TYPE D'EF ( 15 TRIANGLE, 20 TETRAEDRE TAYLOR-HOOD )
! NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
!          ( 6 POUR LE TRIANGLE, 10 POUR LE TETRAEDRE)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF
! PRESP1 : DL DE LA PRESSION AUX SOMMETS DU MAILLAGE

! SORTIE :
! --------
! BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
! MODIF  : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2022
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      include"./incl/donflu.inc95"
      include"./incl/a___npef.inc95"

      DOUBLE PRECISION  D, DELTAT, CoGrPr, BG(NTDLVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NUNOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NUTYEL, NDIM, NBNOEF, NBELEM, NTDLVI,&
                        NBSOM, NBNOVI
      INTEGER           M, K, NSK, MEK, NUELEM, NONOEF(10)
      DOUBLE PRECISION  BE(30)
      DOUBLE PRECISION  PRESP1(NBSOM)

!     MISE A ZERO DE BG
      CALL AZEROD( NTDLVI, BG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NUELEM,M,K,NONOEF,NSK,MEK,BE,D )
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
!     BOUCLE SUR LES EF de TAYLOR-HOOD
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO

!       (Integrale sur e  - dt CoGrPr     tV(x) Grad (p(tn+1,x)- p(tn,x)) dx=)
!        Integrale sur e  - dt !oGrPr tDiv V(x)      (p(tn+1,x)- p(tn,x)) dx
         IF( NUTYEL .EQ. 15 ) THEN

!           TRIANGLE TAYLOR-HOOD

!!!!           Integrale sur e  - dt CoGrPr tDiv V(x)  (p(tn+1,x)- p(tn,x)) dx
!!!            CALL F2EX4P2P1( DELTAT*CoGrPr, NBNOVI,XYZNOE,NONOEF, NUNOSO, &
!!!                            NBSOM, PRESP1,  BE )

!           Integrale sur e  - dt CoGrPr  tV(x) Grad (p(tn+1,x)- p(tn,x)) dx=
            CALL F2EX4P2P1GRAD( XYZNOE, NONOEF, NUNOSO, DELTAT, CoGrPr,  &
                                NBNOVI, NBSOM,  PRESP1, BE )


         ELSE

!           TETRAEDRE TAYLOR-HOOD

!!!!           Integrale sur e  - dt CoGrPr tDiv V(x) PRESP1(x) dx
!!!            CALL F3EX4P2P1( DELTAT*CoGrPr, NBNOVI,XYZNOE,NONOEF, NUNOSO, &
!!!                            NBSOM, PRESP1,  BE )

!           Integrale sur e  - dt CoGrPr  tV(x) Grad PRESP1(x) dx=
            CALL F3EX4P2P1GRAD( XYZNOE, NONOEF, NUNOSO, DELTAT, CoGrPr, &
                                NBNOVI, NBSOM,  PRESP1,   BE )

         ENDIF

!!!         IF( NUELEM .EQ. 5 ) THEN
!!!            call affvect( &
!!!            'BE: - dt CoGrPr tV(x) Grad (p(tn+1,x)- p(tn,x))', &
!!!             NDIM*NBNOEF, BE )
!!!         ENDIF

!        ASSEMBLAGE DE BE dans BG
!!!        CALL ASMEEX1( NDIM, NBNOEF, NONOEF, BE, NBNOVI, BG )

         DO M=1,NBNOEF

!           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            NSK = NONOEF(M)
            MEK = M

!           ASSEMBLAGE DU COEFFICIENT ELEMENTAIRE DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               D = BE( MEK )
!$OMP ATOMIC
               BG( NSK ) = BG( NSK ) + D

               NSK = NSK + NBNOVI
               MEK = MEK + NBNOEF
            ENDDO

         ENDDO

      ENDDO     !FIN DE LA BOUCLE SUR LES EF DE LA CONSTRUCTION DE BG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!      call affvect( 'BGGRATH BG=', 20,     BG )
!      call afl1ve(  'BGGRATH BG=', NTDLVI, BG )

      RETURN
END SUBROUTINE BGGRATH
