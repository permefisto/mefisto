SUBROUTINE P13DAGPR( COEF,   NBSOM,  XYZSOM, &
                     NBNOEF, NBEF,   NONOEF, NONOSO, &
                     NBCOPG, LPDIAG, PG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: CONSTRUIRE ET ASSEMBLER DANS LA MATRICE PROFIL GLOBALE LES
! ---- MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
!      POUR UNE INTERPOLATION P1 EN 3D DE FONCTIONS DE BASE COORDONNEES
!      BARYCENTRIQUES LAMBDA(1:4)

! ENTREES:
! --------
! COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
! NBNOEF : NOMBRE DE SOMMETS (4) D'UN EF
! NBEF   : NOMBRE D'EF DU MAILLAGE
! NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
! NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM
! NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
! LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX   DE PG

! MODIFIES:
! ---------
! PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    AVRIL 2013
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT  NONE
      include"./incl/threads.inc"
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NBSOM, NBNOEF, NBEF, NONOEF(NBEF,NBNOEF), &
                        NONOSO(1:*), LPDIAG(0:NBSOM), NBCOPG
      DOUBLE PRECISION  PG(1:NBCOPG), COEF

      DOUBLE PRECISION  DF(3,3), DFM1(3,3), DELTA, DFM1DLa(3,4), S, D, &
                        X, Y, Z, EPSii, EPSij
      INTEGER           I, J, K, NEF, NS1, NS2, NS3, NS4, IEG, JEG, KG
      INTRINSIC         ABS

!     MISE A ZERO DE LA MATRICE PROFIL GLOBALE PG
      CALL AZEROD( NBCOPG, PG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NEF, NS1,NS2,NS3,NS4, I, J, K, DF, DFM1, DELTA, DFM1DLa ) &
!$OMP PRIVATE( D, S, X, Y, Z, EPSii, EPSij, IEG, JEG, KG )
!     BOUCLE SUR LES EF
!$OMP DO SCHEDULE( STATIC, NBEF/NBTHREADS )
      DO 100 NEF = 1, NBEF

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
         NS4 = NONOEF(NEF,4)

!        CONSTRUCTION DE LA MATRICE DF
         X = XYZSOM(1,NS1)
         DF(1,1) = XYZSOM(1,NS2) - X
         DF(2,1) = XYZSOM(1,NS3) - X
         DF(3,1) = XYZSOM(1,NS4) - X

         Y = XYZSOM(2,NS1)
         DF(1,2) = XYZSOM(2,NS2) - Y
         DF(2,2) = XYZSOM(2,NS3) - Y
         DF(3,2) = XYZSOM(2,NS4) - Y

         Z = XYZSOM(3,NS1)
         DF(1,3) = XYZSOM(3,NS2) - Z
         DF(2,3) = XYZSOM(3,NS3) - Z
         DF(3,3) = XYZSOM(3,NS4) - Z

!        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)) &
                    +DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3)) &
                    +DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'p13dagpr.f95: ATTENTION EF',NEF, &
                      ' de VOLUME*6=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*,'p13dagpr.f95: ATTENTION FE',NEF, &
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

!        COEFFICIENT MULTIPLICATEUR
         D = COEF / ( 6D0 * DELTA )

!        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 120D0
         EPSii = EPSij * 2

         IF( NBNOEF .LE. 4 ) THEN

!           TETRAEDRE BREZZI-FORTIN (NOEUDS=SOMMETS)
            DO I = 1, 4

!              ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS MATRICE PROFIL GLOBALE
!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOEF(NEF,I)

!              Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) sur l'EF NEF
               DO J = 1, I

                  S = 0D0
                  DO K = 1, 3
                     S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
                  ENDDO

!                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOEF(NEF,J)

!                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

!                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  IF( I .EQ. J ) THEN
                     S = S * D + EPSii
                  ELSE
                     S = S * D + EPSij
                  ENDIF
!$OMP ATOMIC
                  PG( KG ) = PG( KG ) + S

               ENDDO
            ENDDO

         ELSE

!           TETRAEDRE TAYLOR-HOOD (NOEUDS=SOMMETS+MILIEUX des ARETES)
            DO I = 1, 4

!              ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS MATRICE PROFIL GLOBALE
!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOSO( NONOEF(NEF,I) )

!              Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) sur l'EF NEF
               DO J = 1, I

                  S = 0D0
                  DO K = 1, 3
                     S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
                  ENDDO

!                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOSO( NONOEF(NEF,J) )

!                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

!                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  IF( I .EQ. J ) THEN
                     S = S * D + EPSii
                  ELSE
                     S = S * D + EPSij
                  ENDIF
!$OMP ATOMIC
                  PG( KG ) = PG( KG ) + S

               ENDDO
            ENDDO

         ENDIF

 100  CONTINUE   ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER PG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!     call affvect( 'p13dagpr.f95: PG=', 20, PG )
!     call afl1ve(  'p13dagpr.f95: PG=', NBCOPG, PG )

      RETURN
END SUBROUTINE P13DAGPR
