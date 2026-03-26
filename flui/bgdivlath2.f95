SUBROUTINE BGDIVLATH2( Rho,    dtMhu,  VITEXY, NBNOVI, XYZNOE, &
                       NONOSO, NBELEM, NUNOEF, NULAEF, &
                       NBSOM,  BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -LAPLACIEN P = -(Rho - dtMhu Laplacien) Div  VITXYZ*
!                -dP/dn = dtMhu n . Laplacien VITXYZ*
!          POUR L'ELEMENT FINI TRIANGLE TAYLOR-HOOD
!          Div VITXYZ* EST INTERPOLEE P1 EXACTEMENT

! ENTREES:
! --------
! Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
! dtMhu  : Coefficient du Laplacien  Integrale dtMhu DP2 DP2 dX
! VITEXY : 2 COMPOSANTES X Y DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
!          DE LA TRIANGULATION
! NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
! XYZNOE : (3,NBNOVI) TABLEAU DES 3 COORDONNEES DES NBNOVI NOEUDS
! NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,6) NO DES 6 NOEUDS DES NBELEM EF
! NULAEF : TABLEAU DU NUMERO DE LIGNE DES 3 ARETES DE CHAQUE TRIANGLE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE ou NOMBRE DE COMPOSANTES DE BG

! SORTIE :
! --------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE

! REMARQUE: SI UN EF EST DE SURFACE NULLE, SA CONTRIBUTION A BG EST NULLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      include"./incl/langue.inc"
      include"./incl/donflu.inc95"
      INTEGER         NBSOM,  NBELEM, NBNOVI, &
                      NUNOEF(NBELEM,6), NULAEF(3,NBELEM), NONOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  Rho, dtMhu, VITEXY(NBNOVI,2), BG(NBSOM)

      INTEGER           I, J, K, L, M, NS1, NS2, NS3, NUELEM, NS, N, &
                        K1, K2
      DOUBLE PRECISION  DIVUP2(3), TDLDL(3,3), X21, X31, Y21, Y31, &
                        DELTA, DFM1(2,2), DFM1DLA(2,3), VN(2), &
                        X1, Y1, X2, Y2, V1, V2, V3, S, &
                        Rho24, dtMhu2

!     [DP2(3Sommets)]   SUR LE TRIANGLE de REFERENCE
      DOUBLE PRECISION  DP2S(2,6,3)
      DATA              DP2S /         &
                   -0.3D+01, -0.3D+01, &
                   -0.1D+01,  0.0D+00, &
                    0.0D+00, -0.1D+01, &
                    0.4D+01,  0.0D+00, &
                    0.0D+00,  0.0D+00, &
                    0.0D+00,  0.4D+01, &
                    0.1D+01,  0.1D+01, &
                    0.3D+01,  0.0D+00, &
                    0.0D+00, -0.1D+01, &
                   -0.4D+01, -0.4D+01, &
                    0.0D+00,  0.4D+01, &
                    0.0D+00,  0.0D+00, &
                    0.1D+01,  0.1D+01, &
                   -0.1D+01,  0.0D+00, &
                    0.0D+00,  0.3D+01, &
                    0.0D+00,  0.0D+00, &
                    0.4D+01,  0.0D+00, &
                   -0.4D+01, -0.4D+01  /

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     CONTRIBUTION DE LA SURFACE AU SECOND MEMBRE
!     ===========================================
      Rho24  = Rho / 24D0
      dtMhu2 = dtMhu / 2D0

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NUELEM,I, J, K, L, M, NS1, NS2, NS3 ) &
!$OMP PRIVATE( NS, N, K1, K2, DIVUP2, TDLDL, X21, X31, Y21, Y31 ) &
!$OMP PRIVATE( DELTA, DFM1, DFM1DLA, VN, X1, Y1, X2, Y2, V1, V2, V3, S )
!     BOUCLE SUR LES EF de TAYLOR-HOOD
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO 100 NUELEM = 1, NBELEM

!        NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
         NS1 = NUNOEF(NUELEM,1)
         NS2 = NUNOEF(NUELEM,2)
         NS3 = NUNOEF(NUELEM,3)

!        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X1  = XYZNOE(1,NS1)
         X21 = XYZNOE(1,NS2) - X1
         X31 = XYZNOE(1,NS3) - X1

         Y1  = XYZNOE(2,NS1)
         Y21 = XYZNOE(2,NS2) - Y1
         Y31 = XYZNOE(2,NS3) - Y1

!        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DELTA = ABS( X21*Y31 - X31*Y21 )
         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT *,'BGDIVLATH2: ATTENTION EF',NUELEM, &
                       ' de SURFACE*2=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT *,'BGDIVLATH2: ATTENTION FE',NUELEM, &
                       ' of SURFACE*2=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

!        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
!        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1
         DFM1(1,1) =  Y31 / DELTA
         DFM1(2,1) = -X31 / DELTA

         DFM1(1,2) = -Y21 / DELTA
         DFM1(2,2) =  X21 / DELTA

!        CALCUL de -Div [P2] {U(3 Sommets)}
!        ----------------------------------
         DO M = 1, 3
!           AU SOMMET M
            S = 0D0
            DO I = 1, 2
!              dUi/dxi
               DO J = 1, 2
                  DO L = 1, 6
!                    NO GLOBAL DU NOEUD J DU TRIANGLE NUELEM
                     NS = NUNOEF(NUELEM,L)
                     S = S - DFM1(I,J) * DP2S(J,L,M) * VITEXY(NS,I)
                  ENDDO
               ENDDO
            ENDDO
            DIVUP2( M ) = S
         ENDDO

!        DFM1 DLambda(2,3)
         DFM1DLa(1,1) = Y21 - Y31
         DFM1DLa(2,1) = X31 - X21

         DFM1DLa(1,2) = Y31
         DFM1DLa(2,2) =-X31

         DFM1DLa(1,3) =-Y21
         DFM1DLa(2,3) = X21

!        tDLambda tDFM1  DFM1 DLambda
         DO I=1,3
            DO J=1,3
               S = 0D0
               DO K=1,2
                  S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
               ENDDO
               TDLDL(I,J) = S
            ENDDO
         ENDDO

!        PREMIERE COMPOSANTE DU SECOND MEMBRE BE(1)
         V1 = Rho24 * 2D0 + dtMhu2 * TDLDL(1,1)
         V2 = Rho24       + dtMhu2 * TDLDL(1,2)
         V3 = Rho24       + dtMhu2 * TDLDL(1,3)
!        ASSEMBLAGE DE BE(1) DANS BG( NONOSO( NOSOMMET1 ) )
         NS = NONOSO( NS1 )
         S  = ( V1*DIVUP2(1) + V2*DIVUP2(2) + V3*DIVUP2(3) ) * DELTA
!$OMP ATOMIC
         BG(NS) = BG(NS) + S

!        SECONDE COMPOSANTE DU SECOND MEMBRE BE(2)
         V1 = Rho24       + dtMhu2 * TDLDL(2,1)
         V2 = Rho24 * 2D0 + dtMhu2 * TDLDL(2,2)
         V3 = Rho24       + dtMhu2 * TDLDL(2,3)
!        ASSEMBLAGE DE BE(2) DANS BG( NONOSO( NOSOMMET2 ) )
         NS = NONOSO( NS2 )
         S  = ( V1*DIVUP2(1) + V2*DIVUP2(2) + V3*DIVUP2(3) ) * DELTA
!$OMP ATOMIC
         BG(NS) = BG(NS) + S

!        TROISIEME COMPOSANTE DU SECOND MEMBRE BE(3)
         V1 = Rho24       + dtMhu2 * TDLDL(3,1)
         V2 = Rho24       + dtMhu2 * TDLDL(3,2)
         V3 = Rho24 * 2D0 + dtMhu2 * TDLDL(3,3)
!        ASSEMBLAGE DE BE(3) DANS BG( NONOSO( NOSOMMET3 ) )
         NS = NONOSO( NS3 )
         S  = ( V1*DIVUP2(1) + V2*DIVUP2(2) + V3*DIVUP2(3) ) * DELTA
!$OMP ATOMIC
         BG(NS) = BG(NS) + S

!        CONTRIBUTION DES ARETES FRONTALIERES AU SECOND MEMBRE
!        =====================================================
         DO K1=1,3

!           NO DE LA LIGNE DE L'ARETE K1
            N = NULAEF(K1,NUELEM)
            IF( N .GT. 0 ) THEN

!              TOUTE LIGNE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
!              LE NUMERO DES 2 SOMMETS DE L'ARETE K1 DU TRIANGLE
               IF( K1 .NE. 3 ) THEN
                  K2 = K1+1
               ELSE
                  K2 = 1
               ENDIF
               NS1 = NUNOEF(NUELEM,K1)
               X1  = XYZNOE(1,NS1)
               Y1  = XYZNOE(2,NS1)

               NS2 = NUNOEF(NUELEM,K2)
               X2  = XYZNOE(1,NS2)
               Y2  = XYZNOE(2,NS2)

!              LE VECTEUR ORTHOGONAL AU VECTEUR TANGENT DE L'ARETE K
               VN(1) = Y2 - Y1
               VN(2) = X1 - X2

!!!!           LE JACOBIEN SUR L'ARETE K1
!!!            DELTA = SQRT( VN(1)**2 + VN(2)**2 )

!              LE VECTEUR NORMAL UNITAIRE
!!!            ( NON /DELTA COMPENSEE PAR LE JACOBIEN )
!!!            VN(1) = VN(1) / DELTA
!!!            VN(2) = VN(2) / DELTA

!              COMPOSANTE K1=K2 DU SECOND MEMBRE BE
               V1 = dtMhu2 * &
                  ( VN(1) * ( DFM1DLA(1,1) * DIVUP2(1) &
                            + DFM1DLA(1,2) * DIVUP2(2) &
                            + DFM1DLA(1,3) * DIVUP2(3) ) &
                  + VN(2) * ( DFM1DLA(2,1) * DIVUP2(1) &
                            + DFM1DLA(2,2) * DIVUP2(2) &
                            + DFM1DLA(2,3) * DIVUP2(3) ) )

!              ASSEMBLAGE DE Be(K1) DANS BG( NONOSO( NOSOMMET1 ) )
               NS = NONOSO( NS1 )
!$OMP ATOMIC
               BG(NS) = BG(NS) - V1

!              BG = BG + Be(K2)
!              ASSEMBLAGE DE Be(K2) DANS BG( NONOSO( NOSOMMET2 ) )
               NS = NONOSO( NS2 )
!$OMP ATOMIC
               BG(NS) = BG(NS) - V1

            ENDIF

!           FIN TRAITEMENT ARETE K1
         ENDDO

!        FIN TRAITEMENT DE L'EF NUELEM
 100  CONTINUE

!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!      call affvect( 'bgdivlath2.f95: BG=', 20, BG )
!      call afl1ve(  'bgdivlath2.f95: BG=', NBSOM, BG )

      RETURN
END SUBROUTINE BGDIVLATH2
