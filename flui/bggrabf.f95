SUBROUTINE BGGRABF( DELTAT, CoGrPr, NDIM,   NBSOM,  NBNOVI, &
                    XYZNOE, NUTYEL, NBNOEF, NBELEM, NUNOEF, &
                    PRESP1, NOBARY, AGGAUSS,   BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
! ---- (Rho - dt Mhu Laplacien){u(tn+1)-u*} = -dt CoGrPr GRAD(p(tn+1)-p(tn))
!      BG=(Integrale sur e  -dt CoGrPr tV(x) Grad PRESP1(x) dx

! ENTREES:
! --------
! DELTAT : PAS CONSTANT DU TEMPS DOUBLE PRECISION
! CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
! NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NUTYEL : NUMERO DU TYPE D'EF ( 13 TRIANGLE, 19 TETRAEDRE BREZZI-FORTIN )
! NBNOEF : NOMBRE DE SOMMETS D'UN EF DE CE TYPE
!          ( 3 POUR LE TRIANGLE, 4 POUR LE TETRAEDRE)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF
! NOBARY : NUMERO DU DL BARYCENTRE DANS UN EF
! PRESP1 : DL DE LA PRESSION P1 AUX SOMMETS DU MAILLAGE

! SORTIE :
! --------
! BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      DOUBLE PRECISION  PIVOT, D, DELTAT, CoGrPr, &
                        BG(NBNOVI,NDIM), AGGAUSS(NOBARY,NBELEM)
      INTEGER           NUNOEF(NBELEM,NBNOEF)
      REAL              XYZNOE(3,NBSOM)
      INTEGER           NUTYEL, NDIM, NBNOEF, NBELEM, NOBARY, &
                        NBSOM, NBNOVI
      INTEGER           M, K, N, MEK, NUELEM, NONOEF(5)
      DOUBLE PRECISION  BE(19)
      DOUBLE PRECISION  PRESP1(NBSOM)

!     MISE A ZERO DE BG
      CALL AZEROD( NDIM*NBNOVI, BG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NUELEM,N,M,K,NONOEF,MEK,BE,D,PIVOT )
!     BOUCLE SUR LES EF de BREZZI-FORTIN
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO
!        NO DU BARYCENTRE
         NONOEF(NOBARY) = NBSOM + NUELEM

!        Integrale sur ef  - dt CoGrPr tV(x) Grad (p(tn+1,x)- p(tn,x)) dx
         IF( NUTYEL .EQ. 13 ) THEN

!           TRIANGLE BREZZI_FORTIN
            CALL F2EX4P1BP1GRAD( XYZNOE, NONOEF, DELTAT, CoGrPr, &
                                 NBSOM,  PRESP1, BE )

!           ASSEMBLAGE DES DEGRES DE LIBERTE 4 ET 8 DE BE DANS BG
!           COMPRESSION DE 8 A 6 DEGRES DE LIBERTE DE BE
            CALL GAS2P1BP1( NBSOM, NUELEM, BE, NBNOVI, BG )

         ELSE
!
!           TETRAEDRE BREZZI_FORTIN
            CALL F3EX4P1BP1GRAD( XYZNOE, NONOEF, DELTAT, CoGrPr, &
                                 NBSOM,  PRESP1, BE )

!           ASSEMBLAGE  DES DEGRES DE LIBERTE 5 10 15 DE BE DANS BG
!           COMPRESSION DES DEGRES DE LIBERTE 5 10 15 DE BE
            CALL GAS3P1BP1( NBSOM, NUELEM, BE, NBNOVI, BG )

         ENDIF

!!!         IF( NUELEM .EQ. 5 ) THEN
!!!            call affvect( &
!!!           'BGGRABF: BE=-dt CoGrPr tV(x) Grad (p(tn+1,x)- p(tn,x))', &
!!!            NDIM*NBNOEF, BE )
!!!         ENDIF

!        ASSEMBLAGE de BE dans BG
!         CALL ASMEEX1( NDIM, NBNOEF, NONOEF, BE, NBNOVI, BG )

         DO M=1,NBNOEF

!           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            N   = NONOEF(M)
            MEK = M
!           ASSEMBLAGE DU COEFFICIENT ELEMENTAIRE DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               D = BE( MEK )
!$OMP ATOMIC
               BG( N,K ) = BG( N,K ) + D

               MEK = MEK + NBNOEF
            ENDDO

         ENDDO

      ENDDO    ! FIN DE LA BOUCLE SUR LES EF DE LA CONSTRUCTION DE BG
!$OMP END DO


!     LA BOUCLE SUR LES ELEMENTS FINIS BREZZI-FORTIN POUR FAIRE GAUSS
!     SUR LE SECOND MEMBRE DU SYSTEME LINEAIRE A RESOUDRE
!     ---------------------------------------------------------------
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        TRIANGLE  BREZZI-FORTIN P1BULLEP3 EN VITESSE et P1 EN PRESSION
!        TETRAEDRE BREZZI-FORTIN P1BULLEP4 EN VITESSE et P1 EN PRESSION
!        --------------------------------------------------------------
!        METHODE DE GAUSS SUR LES NDIM DL DU BARYCENTRE SUR LE SECOND
!        MEMBRE DU SYSTEME A RESOUDRE
!        CALL GAR3P1BP1( NBSOM, NUELEM, NONOEF, AGGAUSS(1,NUELEM), NBNOVI, BG )

!        TRIANGULATION DE GAUSS SUR LA LIGNE NOBARY DU VECTEUR
         DO K = 1, NDIM
            PIVOT = BG( NBSOM+NUELEM, K ) / AGGAUSS(NOBARY,NUELEM)
            DO M = 1, NBNOEF
               D = AGGAUSS(M,NUELEM) * PIVOT
               N = NUNOEF(NUELEM,M)
!$OMP ATOMIC
               BG(N,K) = BG(N,K) - D

            ENDDO
         ENDDO

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!      call affvect( 'BGGRABF BG=', 20,     BG )
!      call afl1ve(  'BGGRABF BG=', NTDLVI, BG )

      RETURN
END SUBROUTINE BGGRABF
