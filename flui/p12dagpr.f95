SUBROUTINE P12DAGPR( COEF,   NBSOM,  XYZSOM, &
                     NBNOEF, NBEF,   NONOEF, NONOSO, &
                     NBCOPG, LPDIAG, PG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE PROFIL GLOBALE LES
! ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
!          POUR UNE INTERPOLATION P1 EN 2D AVEC LES 3 FONCTIONS DE BASE LAMBDAi
! ENTREES:
! --------
! COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
!         (3 POUR BREZZI-FORTIN & 6 POUR TAYLOR-HOOD)
! NBEF   : NOMBRE D'EF DU MAILLAGE
! NONOEF : NUMERO DES NBNOEF NOEUDS DES NBEF EF
! NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM

! NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
! LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX   DE PG

! MODIFIES:
! ---------
! PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    Avril 2013
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT  NONE
      include"./incl/threads.inc"
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,*)
      INTEGER           NBSOM, NBNOEF, NBEF, NONOEF(NBEF,NBNOEF), &
                        NONOSO(1:*), LPDIAG(0:NBSOM), NBCOPG
      DOUBLE PRECISION  PG(1:NBCOPG), COEF

      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, AE(6), &
                        C, EPSii, EPSij, S
      INTEGER           I, J, NEF, NS1, NS2, NS3, IEG, JEG, KG, KE
      INTRINSIC         ABS

!     MISE A ZERO DE LA MATRICE PROFIL GLOBALE PG
      CALL AZEROD( NBCOPG, PG )

!     BOUCLE SUR LES EF
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( X21, Y21, X31, Y31, X32, Y32, DELTA, AE ) &
!$OMP PRIVATE( C, EPSii, EPSij, S, I, J, NEF, NS1, NS2, NS3, IEG, JEG, KG, KE )
!$OMP DO SCHEDULE( STATIC, NBEF/NBTHREADS )
      DO 100 NEF = 1, NBEF

!        NUMERO DES 3 SOMMETS DU TRIANGLE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)

!        [DFM1] [DLa]
         X21 = XYZSOM(1,NS2) - XYZSOM(1,NS1)
         X31 = XYZSOM(1,NS3) - XYZSOM(1,NS1)
         X32 = XYZSOM(1,NS3) - XYZSOM(1,NS2)

         Y21 = XYZSOM(2,NS2) - XYZSOM(2,NS1)
         Y31 = XYZSOM(2,NS3) - XYZSOM(2,NS1)
         Y32 = XYZSOM(2,NS3) - XYZSOM(2,NS2)

!        JACOBIEN de DF
         DELTA = ABS( X21 * Y31 - X31 * Y21 )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*, 'p12dagpr.f95: ATTENTION EF',NEF, &
                       ' de SURFACE*2=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*, 'p12dagpr.f95: ATTENTION FE',NEF, &
                       ' of SURFACE*2=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

!        COEFFICIENT MULTIPLICATEUR / ( 2 * JACOBIEN de DF )
         C = COEF / ( 2D0 * DELTA )

!        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 24D0
         EPSii = EPSij * 2

!        AE = Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) dX sur l'EF NEF
!           + EPSILON Integrale tLa La dX DELTA
         AE(1) = (  Y32 * Y32 + X32 * X32 ) * C + EPSii
         AE(2) = (- Y32 * Y31 - X32 * X31 ) * C + EPSij
         AE(3) = (  Y31 * Y31 + X31 * X31 ) * C + EPSii
         AE(4) = (  Y32 * Y21 + X32 * X21 ) * C + EPSij
         AE(5) = (- Y31 * Y21 - X31 * X21 ) * C + EPSij
         AE(6) = (  Y21 * Y21 + X21 * X21 ) * C + EPSii

!        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE PROFIL GLOBALE
         KE = 0

         IF( NBNOEF .LE. 3 ) THEN

!           TRIANGLE BREZZI-FORTIN (NOEUDS=SOMMETS)
!           XYZSOM CONTIENT LES XYZ DES SOMMETS P1
!           NONOEF SONT LES NUMEROS DES SOMMETS P1
            DO I = 1, 3

!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOEF(NEF,I)

               DO J = 1, I

!                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOEF(NEF,J)

!                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

!                 LE NUMERO DU COEFFICIENT DE LA MATRICE ELEMENTAIRE
                  KE = KE + 1

!                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  S = AE(KE)
!$OMP ATOMIC
                  PG( KG ) = PG( KG ) + S

               ENDDO
            ENDDO

         ELSE

!           TRIANGLE TAYLOR-HOOD (NOEUDS=SOMMETS+MILIEUX ARETES)
!           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
!           NONOEF SONT LES NUMEROS DES NOEUDS P2
            DO I = 1, 3

!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOSO( NONOEF(NEF,I) )

               DO J = 1, I

!                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOSO( NONOEF(NEF,J) )

!                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

!                 LE NUMERO DU COEFFICIENT DE LA MATRICE ELEMENTAIRE
                  KE = KE + 1

!                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  S = AE(KE)
!$OMP ATOMIC
                  PG( KG ) = PG( KG ) + S

               ENDDO
            ENDDO

         ENDIF

 100  ENDDO   ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER PG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

     call affvect( 'p12dagpr.f95: PG=', 20,     PG )
     call afl1ve(  'p12dagpr.f95: PG=', NBCOPG, PG )

      RETURN
END SUBROUTINE P12DAGPR
