SUBROUTINE VORTTH2( NBNOVI, XYZNOE, NBELEM, NUNOEF, NONOSO, &
                    NCAS0, NCAS1, vitx, vity, NBSOM, AG,   ROTV )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CALCULER le TOURBILLON ou VORTICITE ou ROT Vitesse
! ----     PAR PROJECTION SUR LES POLYNOMES P1 SUR CHAQUE TETRAEDRE
!          Integrale tP1 P1 dX = Integrale tP1 Rot U dX
!          Triangle de TAYLOR-HOOD

! ENTREES:
! --------
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBELEM : NOMBRE DE TETRAEDRES P2 DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,10) NO DES NOEUDS DES NBELEM EF
! NONOSO : NO NOEUD  => NO SOMMET
! NCAS0:NCAS1 : LES VECTEURS VITESSE STOCKES
! vitx   : COMPOSANTE X des VECTEURS VITESSE NCAS0:NCAS1
! vity   : COMPOSANTE Y des VECTEURS VITESSE NCAS0:NCAS1
! NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION ET DE LIGNES DE AG
! AG     : MATRICE DIAGONALE AUXILIAIRE DE NBSOM VALEURS

! SORTIES:
! --------
! ROTV   : ROTV(NBSOM,NCAS0:NCAS1) le TOURBILLON ou VORTICITE (dV/dx-dU/dy)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2013
! MODIFS : ALAIN PERRONNET             St PIERRE du PERRAY  Mai     2021
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      INTEGER           NBNOVI, NBELEM, NUNOEF(NBELEM,6),NONOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  AG(NBSOM), ROTV(NBSOM,NCAS0:NCAS1),  &
                        DFM1(2,2), DET, X1, Y1, X21, X31, Y21, Y31
      INTEGER           NUELEM, J, L, M, N, NS, NCAS, NCAS0, NCAS1, NBSOM
      INTRINSIC         SQRT

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity

!     INTEGRALES P1 DP2 SUR LE TETRAEDRE de REFERENCE
!     P1DP22D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
      include"./incl/p1dp22d.inc95"

!     MISE A ZERO DE AG et ROTV
      CALL AZEROD( NBSOM, AG )
      CALL AZEROD( NBSOM*(NCAS1-NCAS0+1), ROTV )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NUELEM, J, L, M, N, NS, NCAS, DFM1, DET, X1, Y1, X21, X31, Y21, Y31 )
!     BOUCLE SUR LES TETRAEDRES P2 POUR CONSTRUIRE AG
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
         J = NUNOEF(NUELEM,1)
         L = NUNOEF(NUELEM,2)
         M = NUNOEF(NUELEM,3)

!        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X1  = XYZNOE(1,J)
         X21 = XYZNOE(1,L) - X1
         X31 = XYZNOE(1,M) - X1

         Y1  = XYZNOE(2,J)
         Y21 = XYZNOE(2,L) - Y1
         Y31 = XYZNOE(2,M) - Y1

!        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe / 6D0
         DET = ABS( X21*Y31 - X31*Y21 ) / 6D0

!        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
!        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 NON / DET
         DFM1(1,1) =  Y31
         DFM1(2,1) = -X31

         DFM1(1,2) = -Y21
         DFM1(2,2) =  X21

         DO M=1,3

!           NO DU SOMMET M DU TRIANGLE P2
            NS = NONOSO( NUNOEF(NUELEM,M) )

!           ASSEMBLAGE DE LA DIAGONALE AG
!$OMP       ATOMIC
            AG( NS ) = AG( NS ) + DET

            DO NCAS = NCAS0, NCAS1

               Y1 = 0D0
               DO L=1,2
                  DO J=1,6
                     N = NUNOEF(NUELEM,J)
                     X1 = P1DP22D(M,L,J)
                     Y1 = Y1 + DFM1(1,L) * X1 * vity(NCAS)%dptab(N) &
                             - DFM1(2,L) * X1 * vitx(NCAS)%dptab(N)
                  ENDDO
               ENDDO

!              ASSEMBLAGE DES SECONDS MEMBRES DU SYSTEME LINEAIRE
!$OMP          ATOMIC
               ROTV( NS, NCAS ) = ROTV( NS, NCAS ) + Y1

            ENDDO

         ENDDO
         
      ENDDO
! FIN DE LA BOUCLE DES TRIANGLES DU MAILLAGE
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////


! RESOLUTION DU SYSTEME LINEAIRE AVEC AG MATRICE DIAGONALE
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NCAS, NS )
!$OMP DO SCHEDULE( STATIC, NBSOM/NBTHREADS )
      DO NS = 1, NBSOM
         DO NCAS = NCAS0, NCAS1
            ROTV( NS, NCAS ) = ROTV( NS, NCAS ) / AG(NS)
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

RETURN
END SUBROUTINE VORTTH2
