 SUBROUTINE VORTBF3( NBNOVI, XYZNOE, NBELEM, NUNOEF, &
                     NCAS0, NCAS1, vitx, vity, vitz, NBSOM, AG,  ROTV )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CALCULER le TOURBILLON ou VORTICITE ou ROT Vitesse
! ----     PAR PROJECTION SUR LES POLYNOMES P1 SUR CHAQUE TETRAEDRE
!          Integrale tP1 P1 dX = Integrale tP1 Rot U dX
!          Tetraedre de BREZZI-FORTIN

! ENTREES:
! --------
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRES DES EF
! XYZNOE : XYZNOE(3,NBNOVI) 3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBELEM : NOMBRE DE TETRAEDRES P1B DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,4) NO DES SOMMETS DES NBELEM EF
! NCAS0:NCAS1 : LES VECTEURS VITESSE STOCKES
! vitx   : COMPOSANTE X des VECTEURS VITESSE NCAS0:NCAS1
! vity   : COMPOSANTE Y des VECTEURS VITESSE NCAS0:NCAS1
! vitz   : COMPOSANTE Z des VECTEURS VITESSE NCAS0:NCAS1
! NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION ET DE LIGNES DE AG
! AG     : MATRICE DIAGONALE AUXILIAIRE DE NBSOM VALEURS

! SORTIES:
! --------
! ROTV   : ROTV(NBSOM, NCAS0:NCAS1, 4) le TOURBILLON (Composante X Y Z et Module)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2013
! MODIFS : ALAIN PERRONNET             St PIERRE du PERRAY  Mai     2021
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      INTEGER           NBNOVI, NBELEM, NUNOEF(NBELEM,4)
      REAL              XYZNOE( 3, NBNOVI )
      DOUBLE PRECISION  AG(NBSOM), ROTV( NBSOM, NCAS0:NCAS1, 4 ), &
                        DF(3,3), DFM1(3,3), DET, X, Y, Z, V1, V2
      INTEGER           NUELEM, I, J, K, L, M, N, NS, &
                        NCAS, NCAS0, NCAS1, K1, K2, &
                        NBSOM, NONOEF(5)
      INTRINSIC         SQRT

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity, vitz

!     INTEGRALES P1 DP1B SUR LE TETRAEDRE de REFERENCE
!     DP1BLa3d(i,k,j) = integrale dP1Bi/dxk Lambdaj dX
      include"./incl/dp1bla3d.inc95"

!     MISE A ZERO DE AG et ROTV
      CALL AZEROD( NBSOM, AG )
      CALL AZEROD( NBSOM*(NCAS1-NCAS0+1)*4, ROTV )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NUELEM, I, J, K, K1, K2, L, M, N, NS, NCAS, NONOEF ) &
!$OMP PRIVATE( DF, DFM1, DET, X, Y, Z, V1, V2 )

!     BOUCLE SUR LES TETRAEDRES P1B POUR CONSTRUIRE AG
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         I = NUNOEF(NUELEM,1)
         NONOEF(1) = I
         J = NUNOEF(NUELEM,2)
         NONOEF(2) = J
         K = NUNOEF(NUELEM,3)
         NONOEF(3) = K
         L = NUNOEF(NUELEM,4)
         NONOEF(4) = L
!        NUMERO DE NOEUD VITESSE DU BARYCENTRE
         NONOEF(5) = NBSOM + NUELEM

!        CONSTRUCTION DE LA MATRICE DF
         X = XYZNOE(1,I)
         Y = XYZNOE(2,I)
         Z = XYZNOE(3,I)

         DF(1,1) = XYZNOE(1,J) - X
         DF(2,1) = XYZNOE(1,K) - X
         DF(3,1) = XYZNOE(1,L) - X

         DF(1,2) = XYZNOE(2,J) - Y
         DF(2,2) = XYZNOE(2,K) - Y
         DF(3,2) = XYZNOE(2,L) - Y

         DF(1,3) = XYZNOE(3,J) - Z
         DF(2,3) = XYZNOE(3,K) - Z
         DF(3,3) = XYZNOE(3,L) - Z

!        LE DETERMINANT DE DF / 24
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) &
             + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) ) &
             + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
         DET = DET / 24D0
!
!        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
!        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DET
!        CAR COMPENSE PAR LA MULTIPLICATION PAR DET DE L'INTEGRATION
         DFM1(1,1) = DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)
         DFM1(2,1) = DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1)
         DFM1(3,1) = DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2)

         DFM1(1,2) = DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3)
         DFM1(2,2) = DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1)
         DFM1(3,2) = DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2)

         DFM1(1,3) = DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)
         DFM1(2,3) = DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1)
         DFM1(3,3) = DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2)

         DO M=1,4

!           NO DU SOMMET M DU TETRAEDRE
            NS = NONOEF(M)

!           ASSEMBLAGE DE LA DIAGONALE AG
!$OMP       ATOMIC
            AG( NS ) = AG( NS ) + DET

            DO K=1,3
               IF( K .LE. 2 ) THEN
                  K1 = K+1
               ELSE
                  K1 = 1
               ENDIF
               IF( K1 .LE. 2 ) THEN
                  K2 = K1+1
               ELSE
                  K2 = 1
               ENDIF

               DO NCAS = NCAS0, NCAS1
                  Z = 0D0
                  DO L=1,3
                     DO J=1,5
                        N = NONOEF(J)
                        X = DP1BLA3D(J,L,M)

                        GOTO( 11, 12, 13 ), K1
11                      V1 = vitx(NCAS)%dptab(N)
                        GOTO 20
12                      V1 = vity(NCAS)%dptab(N)
                        GOTO 20
13                      V1 = vitz(NCAS)%dptab(N)

20                      GOTO( 21, 22, 23 ), K2
21                      V2 = vitx(NCAS)%dptab(N)
                        GOTO 30
22                      V2 = vity(NCAS)%dptab(N)
                        GOTO 30
23                      V2 = vitz(NCAS)%dptab(N)

30                      Z = Z + DFM1(K1,L) * X * V2 &
                              - DFM1(K2,L) * X * V1
                     ENDDO
                  ENDDO
!                 ASSEMBLAGE DES SECONDS MEMBRES DU SYSTEME LINEAIRE
!$OMP             ATOMIC
                  ROTV( NS, NCAS, K ) = ROTV( NS, NCAS, K ) + Z
               ENDDO

            ENDDO

         ENDDO
         
      ENDDO
! FIN DE LA BOUCLE DES TETRAEDRES DU MAILLAGE
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////


! RESOLUTION DU SYSTEME LINEAIRE AVEC AG MATRICE DIAGONALE
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NCAS, K, NS )
!$OMP DO SCHEDULE( STATIC, NBSOM/NBTHREADS )
      DO NS = 1, NBSOM
         DO NCAS = NCAS0, NCAS1
            DO K = 1, 3
               ROTV( NS, NCAS, K ) = ROTV( NS, NCAS, K ) / AG(NS)
            ENDDO
!           LE MODULE
            ROTV( NS, NCAS, 4 ) = SQRT( ROTV( NS, NCAS, 1 ) ** 2 &
                                      + ROTV( NS, NCAS, 2 ) ** 2 &
                                      + ROTV( NS, NCAS, 3 ) ** 2 )
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

RETURN
END SUBROUTINE VORTBF3
