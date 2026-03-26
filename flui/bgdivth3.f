      SUBROUTINE BGDIVTH3( NBSOM,  XYZNOE,
     %                     NBNOEF, NBELEM, NUNOEF, NONOSO,
     %                     COEFBG, NBNOVI, VITXYZ, VOLUME, INTDIVV, BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoGrPr LAPLACIEN P = COEFBG * ( Div VITXYZ
!                                            - Integrale Div VITXYZ/VOLUME)
!          -dt CoGrPr       dP/dn = 0
!          pour des TETRAEDREs TAYLOR-HOOD

! BG=Integrale COEFBG * tP1 ( Div VITXYZ - Integrale Div VITXYZ / VOLUME ) dx

! ENTREES:
! --------
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
! NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! COEFBG : COEFFICIENT DU SECOND MEMBRE
! NBNOVI : NOMBRE DE NOEUDS VITESSE
! VITXYZ : 3 COMPOSANTES DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
!          DU MAILLAGE
! VOLUME : VOLUME DU MAILLAGE 3D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITXYZ SUR LE MAILLAGE

! SORTIES:
! ---------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOEF, NBELEM, NBNOVI,
     %                  NUNOEF(NBELEM,NBNOEF), NONOSO(NBNOVI)
      DOUBLE PRECISION  COEFBG, VITXYZ(NBNOVI,3), BG(NBSOM),
     %                  VOLUME, INTDIVV, MOYDIVV
!
      DOUBLE PRECISION  DET, DF(3,3), DFM1(3,3), X, Y, Z
      INTEGER           NOSOTE(4), I, J, K, L, NS, NUELEM

!     DOUBLE PRECISION  P1DP23D(4,3,10)
!     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
      include"./incl/p1dp23d.inc"

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     MOYENNE VOLUMIQUE DE DIV VITXYZ / Integrale P1 dx
      MOYDIVV = INTDIVV / VOLUME / 24D0

!     BOUCLE SUR LES EF de TAYLOR-HOOD
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         NOSOTE(1) = NUNOEF(NUELEM,1)
         NOSOTE(2) = NUNOEF(NUELEM,2)
         NOSOTE(3) = NUNOEF(NUELEM,3)
         NOSOTE(4) = NUNOEF(NUELEM,4)

!        CONSTRUCTION DE LA MATRICE DF
         X = XYZNOE(1,NOSOTE(1))
         DF(1,1) = XYZNOE(1,NOSOTE(2)) - X
         DF(2,1) = XYZNOE(1,NOSOTE(3)) - X
         DF(3,1) = XYZNOE(1,NOSOTE(4)) - X

         Y = XYZNOE(2,NOSOTE(1))
         DF(1,2) = XYZNOE(2,NOSOTE(2)) - Y
         DF(2,2) = XYZNOE(2,NOSOTE(3)) - Y
         DF(3,2) = XYZNOE(2,NOSOTE(4)) - Y

         Z = XYZNOE(3,NOSOTE(1))
         DF(1,3) = XYZNOE(3,NOSOTE(2)) - Z
         DF(2,3) = XYZNOE(3,NOSOTE(3)) - Z
         DF(3,3) = XYZNOE(3,NOSOTE(4)) - Z

!        LE DETERMINANT DE DF = VOLUME DU TETRAEDRE * 6
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %       + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %       + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
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

!        DET integrale tLambda ( DF-1 ligne1  DP ue1
!                              + DF-1 ligne2  DP ue2
!                              + DF-1 ligne3  DP ue3
!                              - INTDIVV/VOLUME ) dX
         Y = MOYDIVV * DET
         DO I=1,4

!           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            X = 0D0
            DO J=1,10
!              NO GLOBAL DU NOEUD J DU TETRAEDRE TAYLOR-HOOD
               NS = NUNOEF(NUELEM,J)
               DO K=1,3
!                 COMPOSANTE K DE LA VITESSE AU NOEUD J 
                  Z = VITXYZ(NS,K)
                  DO L=1,3
!                    P1DP23D(i,l,j) = integrale Lambdai dP2j/dxl dx dy dz
                     X = X + DFM1(K,L) * P1DP23D(I,L,J) * Z
                  ENDDO
               ENDDO
            ENDDO

!           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTE(I) ) )
            X  = COEFBG * ( X - Y )
            NS = NONOSO( NOSOTE(I) )
            BG(NS) = BG(NS) + X

         ENDDO

      ENDDO

!!!      call affvect( 'bgdivth3.f: BG=', 20,    BG )
!!!      call afl1ve(  'bgdivth3.f: BG=', NBSOM, BG )

      RETURN
      END SUBROUTINE BGDIVTH3
