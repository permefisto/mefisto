SUBROUTINE BGDIVTH2( NBSOM,  XYZNOE, NBNOEF, NBELEM, NUNOEF,  NONOSO, &
                     COEFBG, NBNOVI, VITEXY, SURFAC, INTDIVV, BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoGrPr LAPLACIEN P = COEFBG * ( Div VITEXY
!                                            - Integrale Div VITEXY/SURFAC)
!          -dt CoGrPr       dP/dn = 0
!          pour des TRIANGLES TAYLOR-HOOD

! BG=Integrale COEFBG * tP1 ( Div VITEXY - Integrale Div VITXYZ / SURFAC ) dx

! ENTREES:
! --------
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUMERO DES NBNOEF NOEUDS DES NBELEM EF
! NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! COEFBG : COEFFICIENT DU SECOND MEMBRE
! VITEXY : 2 COMPOSANTES X Y DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
!          DE LA TRIANGULATION
! SURFAC : SURFACE DU MAILLAGE 2D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITEXY SUR LE MAILLAGE

! SORTIE :
! ---------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      INTEGER           NBSOM, NBNOEF, NBELEM, NBNOVI
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NONOSO(NBNOVI)
      DOUBLE PRECISION  COEFBG, VITEXY(NBNOVI,2), BG(NBSOM), &
                        SURFAC, INTDIVV, MOYDIVV

      DOUBLE PRECISION  X, Y, X21, X31, Y21, Y31, DFM1(2,2), DET
      INTEGER           I, J, K, L, NUELEM, NS
      INTEGER           NOSOTR(3)

!     INTEGRALES P1 DP2 SUR LE TRIANGLE de REFERENCE
!     P1DP2(i,k,j) = integrale P1i DPj/dxk dx dy sur le TRIANGLE UNITE
!     DOUBLE PRECISION  P1DP22D(3,2,6)
      include"./incl/p1dp22d.inc95"

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     MOYENNE SURFACIQUE DE DIV VITEXY / Integrale P1 dx
      MOYDIVV = INTDIVV / SURFAC / 6D0

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NUELEM,I,J,K,L,NS,NOSOTR,X,Y,X21,X31,Y21,Y31,DFM1,DET)
!     BOUCLE SUR LES EF de TAYLOR-HOOD
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
         NOSOTR(1) = NUNOEF(NUELEM,1)
         NOSOTR(2) = NUNOEF(NUELEM,2)
         NOSOTR(3) = NUNOEF(NUELEM,3)

!        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         I   = NOSOTR(1)
         X   = XYZNOE(1,I)
         X21 = XYZNOE(1,NOSOTR(2)) - X
         X31 = XYZNOE(1,NOSOTR(3)) - X

         Y   = XYZNOE(2,I)
         Y21 = XYZNOE(2,NOSOTR(2)) - Y
         Y31 = XYZNOE(2,NOSOTR(3)) - Y

!        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DET = ABS( X21*Y31 - X31*Y21 )

!        L'INTEGRALE tPression DIV Vitesse dx

!        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
!        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DET
!        CAR COMPENSE PAR LA MULTIPLICATION PAR DELTA DE L'INTEGRATION
         DFM1(1,1) =  Y31
         DFM1(2,1) = -X31

         DFM1(1,2) = -Y21
         DFM1(2,2) =  X21

!        DET integrale tLambda ( DF-1 ligne1  DP ue1
!                               +DF-1 ligne2  DP ue2 
!                              - INTDIVV/SURFAC ) dX
         Y = MOYDIVV * DET
         DO I=1,3

!           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            X = 0D0
            DO J=1,6
!              NO GLOBAL DU NOEUD J DU TRIANGLE TAYLOR-HOOD
               NS = NUNOEF(NUELEM,J)
               DO K=1,2
                  DO L=1,2
!                    P1DP22D(i,l,j) = integrale Lambdai dPj/dxl dx dy dz
                     X = X + DFM1(K,L) * P1DP22D(I,L,J) * VITEXY(NS,K)
                  ENDDO
               ENDDO
            ENDDO

!           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTR(I) ) )
            NS = NONOSO( NOSOTR(I) )
            X  = COEFBG * ( X - Y )
!$OMP ATOMIC
            BG(NS) = BG(NS) + X

         ENDDO

      ENDDO

!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!      call affvect( 'bgdivth2.f95: BG=', 20,    BG )
!      call afl1ve(  'bgdivth2.f95: BG=', NBSOM, BG )

      RETURN
END SUBROUTINE BGDIVTH2
