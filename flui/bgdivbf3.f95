SUBROUTINE BGDIVBF3( NBSOM,  XYZSOM, NBSTEF, NBELEM, NOSTEF, &
                     COEFBG, NBNOVI, VITXYZ, VOLUME, INTDIVV,  BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoPres LAPLACIEN P = COEFBG * Div VITXYZ(tn+1)
!                                            - Integrale Div VITXYZ/VOLUME)
!          -dt CoGrPr       dP/dn = 0
!          pour des TETRAEDREs BREZZI-FORTIN

! BG=Integrale COEFBG * tP1 ( Div VITXYZ - Integrale Div VITXYZ / VOLUME ) dx

! ENTREES:
! --------
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
! NBSTEF : NOMBRE DE SOMMETS D'UN EF = 4
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NOSTEF : NUMERO DES NBSTEF SOMMETS DES NBELEM EF
! COEFBG : COEFFICIENT DU SECOND MEMBRE
! NBNOVI : NOMBRE DE NOEUDS (SOMMETS+BARYCENTRES) DU MAILLAGE
! VITXYZ : 3 COMPOSANTES DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
! VOLUME : VOLUME DU MAILLAGE 3D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITXYZ SUR LE MAILLAGE

! SORTIE :
! ---------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      REAL              XYZSOM(3,NBNOVI)
      INTEGER           NBSOM, NBSTEF, NBELEM, NOSTEF(NBELEM,NBSTEF), &
                        NBNOVI
      DOUBLE PRECISION  COEFBG, VITXYZ(NBNOVI,3), BG(NBSOM), &
                        VOLUME, INTDIVV, MOYDIVV

      DOUBLE PRECISION  DET, DF(3,3), DFM1(3,3), X, Y, Z
      INTEGER           NOSOTE(5), NUELEM, I, J, K, L, NS

!     DOUBLE PRECISION  DP1BLa3d(5,3,4)
!     DP1BLa3d(i,k,j)   = integrale dP1Bi/dxk Lambdaj dX
      include"./incl/dp1bla3d.inc95"

!     MISE A ZERO DE LA MATRICE MORSE GLOBALE
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     MOYENNE VOLUMIQUE DE DIV VITXYZ / Integrale P1 dx
      MOYDIVV = INTDIVV / VOLUME / 24D0

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( DET, DF, DFM1, X,Y,Z, NOSOTE, NUELEM, I,J,K,L, NS )
!     BOUCLE SUR LES EF de BREZZI-FORTIN
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         NOSOTE(1) = NOSTEF(NUELEM,1)
         NOSOTE(2) = NOSTEF(NUELEM,2)
         NOSOTE(3) = NOSTEF(NUELEM,3)
         NOSOTE(4) = NOSTEF(NUELEM,4)
         NOSOTE(5) = NBSOM + NUELEM

!        CONSTRUCTION DE LA MATRICE DF
         X = XYZSOM(1,NOSOTE(1))
         DF(1,1) = XYZSOM(1,NOSOTE(2)) - X
         DF(2,1) = XYZSOM(1,NOSOTE(3)) - X
         DF(3,1) = XYZSOM(1,NOSOTE(4)) - X

         Y = XYZSOM(2,NOSOTE(1))
         DF(1,2) = XYZSOM(2,NOSOTE(2)) - Y
         DF(2,2) = XYZSOM(2,NOSOTE(3)) - Y
         DF(3,2) = XYZSOM(2,NOSOTE(4)) - Y

         Z = XYZSOM(3,NOSOTE(1))
         DF(1,3) = XYZSOM(3,NOSOTE(2)) - Z
         DF(2,3) = XYZSOM(3,NOSOTE(3)) - Z
         DF(3,3) = XYZSOM(3,NOSOTE(4)) - Z

!        LE DETERMINANT DE DF = VOLUME DU TETRAEDRE * 6
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) &
             + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) ) &
             + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )

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

!        DET integrale tLambda ( DF-1 ligne1  DP ue1
!                               +DF-1 ligne2  DP ue2
!                               +DF-1 ligne3  DP ue3
!                              - INTDIVV/VOLUME ) dX
         Y = MOYDIVV * DET
         DO I=1,4

!           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            X = 0D0
            DO J=1,5
!              NO GLOBAL DU NOEUD J DU TETRAEDRE BREZZI-FORTIN
               NS = NOSOTE(J)
               DO K=1,3
!                 COMPOSANTE K DE LA VITESSE AU NOEUD J 
                  Z = VITXYZ(NS,K)
                  DO L=1,3
!                    DP1BLa3d(j,l,i) = integrale dP1Bj/dxl Lambdai dX
                     X = X + DFM1(K,L) * DP1BLa3d(J,L,I) * Z
                  ENDDO
               ENDDO
            ENDDO

!           ASSEMBLAGE DE BE(I) DANS BG(NOSOTE(I))
            X  = COEFBG * ( X - Y )
            NS = NOSOTE(I)
!$OMP ATOMIC
            BG(NS) = BG(NS) + X

         ENDDO

      ENDDO

!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!      call affvect( 'bgdivbf3.f95: BG=', 20,    BG )
!!!      call afl1ve(  'bgdivbf3.f95: BG=', NBSOM, BG )

      RETURN
END SUBROUTINE BGDIVBF3
