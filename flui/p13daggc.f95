SUBROUTINE P13DAGGC( COEF,   NBSOM,  XYZSOM, &
                     NBNOEF, NBELEM, NUNOEF, NONOSO, &
                     NBCOPG, LPLIGN, LPCOLO, PG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE MORSE GLOBALE LES
! ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + Epsilon (p,q)
!          POUR UNE INTERPOLATION P1 EN 3D DE FONCTIONS DE BASE LAMBDA(1:4)

! ENTREES:
! --------
! COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
!         (4 POUR BREZZI-FORTIN & 10 POUR TAYLOR-HOOD)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
! NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM

! NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
! LPLIGN : POINTEUR SUR LES COEFFICIENTS DIAGONAUX      DE PG
! LPCOLO : NUMERO DES COLONNES DES COEFFICIENTS STOCKES DE PG

! SORTIE :
! ---------
! PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT  NONE
      include"./incl/threads.inc"
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-6)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NBSOM, NBNOEF, NBELEM, NBCOPG
      INTEGER           NUNOEF(NBELEM,NBNOEF), &
                        NONOSO(1:*), &
                        LPLIGN(0:NBSOM), &
                        LPCOLO(1:NBCOPG)
      DOUBLE PRECISION  COEF, PG(1:NBCOPG)

      DOUBLE PRECISION  AE(10), DF(3,3), DFM1(3,3), DELTA, DFM1DLa(3,4), &
                        S, X, Y, Z, C, EPSii, EPSij
      INTEGER           NUELEM, NOSOTE(4), I, J, K, M
      INTRINSIC         ABS

!     MISE A ZERO DE LA MATRICE MORSE GLOBALE PG
      CALL AZEROD( NBCOPG, PG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NUELEM, NOSOTE, I, J, K, M, AE, DF, DFM1, DELTA, DFM1DLa ) &
!$OMP PRIVATE( C, S, X, Y, Z, EPSii, EPSij )
!     BOUCLE SUR LES EF
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO 100 NUELEM = 1, NBELEM

!        NUMERO DE NOEUD DES 4 SOMMETS DU TETRAEDRE NUELEM
         DO I=1,4
            NOSOTE(I) = NUNOEF(NUELEM,I)
         ENDDO
 
!        CONSTRUCTION DE LA MATRICE DF
         I = NOSOTE(1)
         X = XYZSOM(1,I)
         Y = XYZSOM(2,I)
         Z = XYZSOM(3,I)

         J = NOSOTE(2)
         DF(1,1) = XYZSOM(1,J) - X
         DF(1,2) = XYZSOM(2,J) - Y
         DF(1,3) = XYZSOM(3,J) - Z

         K = NOSOTE(3)
         DF(2,1) = XYZSOM(1,K) - X
         DF(2,2) = XYZSOM(2,K) - Y
         DF(2,3) = XYZSOM(3,K) - Z

         M = NOSOTE(4)
         DF(3,1) = XYZSOM(1,M) - X
         DF(3,2) = XYZSOM(2,M) - Y
         DF(3,3) = XYZSOM(3,M) - Z

!        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)) &
                   + DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3)) &
                   + DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'p13daggc.f95: ATTENTION EF',NUELEM, &
                      ' de VOLUME*6=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*,'p13daggc.f95: ATTENTION FE',NUELEM, &
                      ' of VOLUME*6=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

!        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
!        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
         DFM1(1,1) = DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)
         DFM1(2,1) = DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1)
         DFM1(3,1) = DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2)

         DFM1(1,2) = DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3)
         DFM1(2,2) = DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1)
         DFM1(3,2) = DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2)

         DFM1(1,3) = DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)
         DFM1(2,3) = DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1)
         DFM1(3,3) = DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2)

!        [DFM1] [DLa]
         DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
         DFM1DLa(1,2) =  DFM1(1,1)
         DFM1DLa(1,3) =  DFM1(1,2)
         DFM1DLa(1,4) =  DFM1(1,3)

         DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
         DFM1DLa(2,2) =  DFM1(2,1)
         DFM1DLa(2,3) =  DFM1(2,2)
         DFM1DLa(2,4) =  DFM1(2,3)

         DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
         DFM1DLa(3,2) =  DFM1(3,1)
         DFM1DLa(3,3) =  DFM1(3,2)
         DFM1DLa(3,4) =  DFM1(3,3)

!        COEF / ( 6D0 * JACOBIEN de DF )
         C = COEF / ( 6D0 * DELTA )

!        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 120D0
         EPSii = EPSij * 2

!        AE = COEF  Integrale  t([DFM1] [DLa]) ([DFM1] [DLa]) dx
         M = 0
         DO I = 1, 4
            DO J = 1, I
               S = 0D0
               DO K = 1, 3
                  S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
               ENDDO
               IF( I .EQ. J ) THEN
                  S = S * C + EPSii
               ELSE
                  S = S * C + EPSij
               ENDIF
               M = M + 1
               AE(M) = S
            ENDDO
         ENDDO

         IF( NBNOEF .GT. 4 ) THEN
!           NUMERO DE SOMMET DES 4 SOMMETS DU TETRAEDRE NUELEM INTERPOLE P2
!           EF TAYLOR-YOUNG
!           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
!           NUNOEF SONT LES NUMEROS DES NOEUDS P2
!           NONOSO SONT LES NUMEROS DE SOMMET DES NOEUDS P2
!           MAIS LPLIGN, LPCOLO SONT DEFINIS SUR LE MAILLAGE P1
            DO I=1,4
               NOSOTE(I) = NONOSO( NUNOEF(NUELEM,I) )
            ENDDO
         ENDIF

!        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE MORSE GLOBALE
!        LES MATRICES ELEMENTAIRE ET GLOBALE SONT ICI SYMETRIQUES
         CALL ASMEGC( 4, NOSOTE, 1, AE, AE, 1, LPLIGN, LPCOLO, PG )

 100  CONTINUE   ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER PG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!      call affvect( 'p13daggc.f95: PG=', 20, PG )
!      call afl1ve(  'p13daggc.f95: PG=', NBCOPG, PG )

      RETURN
END SUBROUTINE P13DAGGC
