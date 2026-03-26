SUBROUTINE P12DAGGC( COEF,   NBSOM,  XYZSOM, &
                     NBNOEF, NBELEM, NUNOEF, NONOSO, &
                     NBCOPG, LPLIGN, LPCOLO, PG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE MORSE GLOBALE LES
! ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
!          POUR UNE INTERPOLATION P1 EN 2D DE FONCTIONS DE BASE LAMBDA(1:3)
! ENTREES:
! --------
! COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
!         (3 POUR BREZZI-FORTIN & 6 POUR TAYLOR-HOOD)
! NBELEM : NOMBRE DE TRIANGLES DU MAILLAGE
! NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
! NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM

! NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
! LPLIGN : POINTEUR SUR LES COEFFICIENTS DIAGONAUX      DE PG
! LPCOLO : NUMERO DES COLONNES DES COEFFICIENTS STOCKES DE PG

! MODIFIES:
! ---------
! PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT  NONE
      include"./incl/threads.inc"
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,*)
      INTEGER           NBSOM, NBNOEF, NBELEM, NBCOPG
      INTEGER           NUNOEF(NBELEM,NBNOEF), &
                        NONOSO(1:*), &
                        LPLIGN(0:NBSOM), &
                        LPCOLO(1:NBCOPG)
      DOUBLE PRECISION  COEF, PG(1:NBCOPG)

      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, AE(6), &
                        C, EPSii, EPSij
      INTEGER           I, M, NUELEM, NOSOTR(3)

      INTRINSIC         ABS

!     MISE A ZERO DE LA MATRICE MORSE GLOBALE PG
      CALL AZEROD( NBCOPG, PG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NOSOTR, NUELEM, I, M ) &
!$OMP PRIVATE( X21, Y21, X31, Y31, X32, Y32, DELTA, AE, C, EPSii, EPSij )
!     BOUCLE SUR LES EF de TAYLOR-HOOD
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO 100 NUELEM = 1, NBELEM

!        NUMERO XYZSOM DES 3 SOMMETS DU TRIANGLE NUELEM
!        XYZSOM CONTIENT LES XYZ DES SOMMETS P1 SI BREZZI-FORTIN
!        XYZSOM CONTIENT LES XYZ DES NOEUDS  P2 SI TAYLOR-YOUNG
         DO I=1,3
            NOSOTR(I) = NUNOEF(NUELEM,I)
         ENDDO

!        [DFM1] [DLa]
         X21 = XYZSOM(1,NOSOTR(2)) - XYZSOM(1,NOSOTR(1))
         X31 = XYZSOM(1,NOSOTR(3)) - XYZSOM(1,NOSOTR(1))
         X32 = XYZSOM(1,NOSOTR(3)) - XYZSOM(1,NOSOTR(2))

         Y21 = XYZSOM(2,NOSOTR(2)) - XYZSOM(2,NOSOTR(1))
         Y31 = XYZSOM(2,NOSOTR(3)) - XYZSOM(2,NOSOTR(1))
         Y32 = XYZSOM(2,NOSOTR(3)) - XYZSOM(2,NOSOTR(2))

!        JACOBIEN DE F
         DELTA = ABS( X21 * Y31 - X31 * Y21 )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'p12daggc.f95: ATTENTION EF',NUELEM, &
                      ' de SURFACE*2=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*,'p12daggc.f95: ATTENTION FE',NUELEM, &
                      ' of SURFACE*2=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

!        COEFFICIENT MULTIPLICATEUR / ( 2 * JACOBIEN de DF )
         C = COEF / ( 2D0 * DELTA )

!        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 24D0
         EPSii = EPSij * 2

!        AE = Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) dX sur l'EF NUELEM
!           + EPSILON Integrale tLa La dX DELTA
         AE(1) = (  Y32 * Y32 + X32 * X32 ) * C + EPSii
         AE(2) = (- Y32 * Y31 - X32 * X31 ) * C + EPSij
         AE(3) = (  Y31 * Y31 + X31 * X31 ) * C + EPSii
         AE(4) = (  Y32 * Y21 + X32 * X21 ) * C + EPSij
         AE(5) = (- Y31 * Y21 - X31 * X21 ) * C + EPSij
         AE(6) = (  Y21 * Y21 + X21 * X21 ) * C + EPSii

!!!         print *,'p12daggc.f95: NOSOTR=',NOSOTR
!!!         print *,'p12daggc.f95: AE=',(AE(m),m=1,6)

         IF( NBNOEF .GT. 3 ) THEN
!           EF TAYLOR-YOUNG
!           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
!           NUNOEF SONT LES NUMEROS DES NOEUDS P2
!           NONOSO SONT LES NUMEROS DE SOMMET DES NOEUDS P2
!           MAIS LPLIGN, LPCOLO SONT DEFINIS SUR UN MAILLAGE P1
            DO I=1,3
               NOSOTR(I) = NONOSO( NUNOEF(NUELEM,I) )
            ENDDO
         ENDIF

!        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE MORSE GLOBALE
!        LES MATRICES ELEMENTAIRE ET GLOBALE SONT ICI SYMETRIQUES
         CALL ASMEGC( 3, NOSOTR, 1, AE, AE, 1, LPLIGN, LPCOLO, PG )

!!!         print *,'p12daggc.f95: PG GC=',(PG(m),m=1,20)

 100  ENDDO   ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER PG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!      call affvect( 'p12daggc.f95: PG=', 5,      PG )
!!!      call afl1ve(  'p12daggc.f95: PG=', NBCOPG, PG )


      RETURN
END SUBROUTINE P12DAGGC
