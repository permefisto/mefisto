      SUBROUTINE F3EBPP2P1( Rho, DELTAT, Omega, NBSOM, NBNOVI, XYZNOE,
     %                      NBNOEF, NBELEM, NUNOEF, NONOSO,
     %                      WSTAR, VOLUME, INTDIVWS, Wtn, Wm,  BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoGrPr LAPLACIEN P = -Rho * ( Div WSTAR
!                                          - Integrale Div WSTAR/VOLUME )
!                      + 2 Rho Dt Div( Omega(tn+1) x ( Wn+1,m - W(tn) ) )
!          -dt CoGrPr dP/dn = 0
!          pour des TETRAEDREs TAYLOR-HOOD (issu de bodivlath3.f)

! BG=Integrale tP1 * ( -Rho * ( Div WSTAR - Integrale Div WSTAR / VOLUME)
!                    + 2 Rho Dt Div( Omega(tn+1) x (Wn+1,m-W(tn) )  dx

! ENTREES:
! --------
! Rho    : DENSITE VOLUMIQUE DE MASSE DU FLUIDE
! DELTAT : PAS DE TEMPS DU SCHEMA EN TEMPS
! Omega  : 3 COMPOSANTES DE LA VITESSE ANGULAIRE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBNOVI : NOMBRE DE NOEUDS VITESSE DU MAILLAGE
! XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN EF
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
! NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! WSTAR  : 3 COMPOSANTES DE LA VITESSE W*     EN LES NBNOVI NOEUDS VITESSE
! Wtn    : 3 COMPOSANTES DE LA VITESSE W(tn)  EN LES NBNOVI NOEUDS VITESSE
! Wm     : 3 COMPOSANTES DE LA VITESSE Wn+1,m EN LES NBNOVI NOEUDS VITESSE
! VOLUME : VOLUME DU MAILLAGE 3D
! INTDIVWS: INTEGRALE DE LA DIVERGENCE DE WSTAR SUR LE MAILLAGE

! SORTIE :
! ---------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Fevrier 2013
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOEF, NBELEM, NBNOVI,
     %                  NUNOEF(NBELEM,NBNOEF), NONOSO(NBNOVI)
      DOUBLE PRECISION  Rho, DELTAT, Omega(3),
     %                  WSTAR(NBNOVI,3), Wtn(NBNOVI,3), Wm(NBNOVI,3),
     %                  BG(NBSOM),
     %                  VOLUME, INTDIVWS, MOYDIVWS
!
      DOUBLE PRECISION  DET, DF(3,3), DFM1(3,3), X, Y, Z, DEUXDT
      INTEGER           NOSOTE(4), I, J, K, K1, K2, L, NS, NUELEM
      INTRINSIC         ABS

!     DOUBLE PRECISION  P1DP23D(4,3,10)
!     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
      include"./incl/p1dp23d.inc"

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     MOYENNE VOLUMIQUE DE DIV WSTAR / Integrale P1 dx
      MOYDIVWS = INTDIVWS / VOLUME / 24D0

      DEUXDT = 2D0 * DELTAT

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
         DET= ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %           + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %           + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) )
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
!                              - INTDIVWS/VOLUME ) dX
         Y = MOYDIVWS * DET
         DO I=1,4

!           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            X = 0D0
            DO J=1,10
!              NO GLOBAL DU NOEUD J DU TETRAEDRE TAYLOR-HOOD
               NS = NUNOEF(NUELEM,J)
               DO K=1,3
                  IF( K .LT. 3 ) THEN
                     K1 = K+1
                  ELSE
                     K1 = 1
                  ENDIF
                  IF( K1 .LT. 3 ) THEN
                     K2 = K1+1
                  ELSE
                     K2 = 1
                  ENDIF
!                 COMPOSANTE K DE LA VITESSE AU NOEUD J 
                  Z = DEUXDT * ( Omega(K1) * (Wm(NS,K2)-Wtn(NS,K2) )
     %                         - Omega(K2) * (Wm(NS,K1)-Wtn(NS,K1) ) )
     %              - WSTAR(NS,K) 
                  DO L=1,3
!                    P1DP23D(i,l,j) = integrale Lambdai dP2j/dxl dx dy dz
                     X = X + DFM1(K,L) * P1DP23D(I,L,J) * Z
                  ENDDO
               ENDDO
            ENDDO

!           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTE(I) ) )
            NS = NONOSO( NOSOTE(I) )
            BG(NS) = BG(NS) + Rho * ( X - Y )

         ENDDO

      ENDDO

!!!      call affvect( 'f3ebpp2p1.f: BG=', 20,    BG )
!!!      call afl1ve(  'f3ebpp2p1.f: BG=', NBSOM, BG )

      RETURN
      END SUBROUTINE F3EBPP2P1
