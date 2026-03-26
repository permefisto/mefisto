      SUBROUTINE BGDIVBF2( NBSOM,  XYZSOM, NBSTEF, NBELEM, NOSTEF,
     %                     COEFBG, NBNOVI, VITEXY, SURFAC, INTDIVV, BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoGrPr LAPLACIEN P = COEFBG * (Div VITEXY(tn+1)
!                                            - Integrale Div VITEXY/SURFAC)
!          -dt CoGrPr       dP/dn = 0
!          pour des TRIANGLES TRIANGLE BREZZI-FORTIN
!
! ENTREES:
! --------
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
! NBSTEF : NOMBRE DE SOMMETS D'UN EF
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NOSTEF : NUMERO DES NBSTEF SOMMETS DES NBELEM EF
! COEFBG : COEFFICIENT DU SECOND MEMBRE
! VITEXY : 2 COMPOSANTES X Y DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
!          DE LA TRIANGULATION
! SURFAC : SURFACE DU MAILLAGE 2D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITEXY SUR LE MAILLAGE
!
! SORTIE :
! ---------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NBSOM, NBSTEF, NBELEM, NBNOVI
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NOSTEF(NBELEM,NBSTEF)
      DOUBLE PRECISION  COEFBG, VITEXY(NBNOVI,2), BG(NBSOM),
     %                  SURFAC, INTDIVV, MOYDIVV

      DOUBLE PRECISION  DET, DFM1(2,2), X, Y, X21, X31, Y21, Y31
      INTEGER           I, J, NUELEM, NS, NOSOTR(4)

!     DOUBLE PRECISION  DP1BLa2d(4,2,3)
!     DP1BLA(i,k,j) = integrale DP1Bi/dxk Lambdaj dX
      include"./incl/dp1bla2d.inc"

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     MOYENNE SURFACIQUE DE DIV VITEXY / Integrale P1 dx
      MOYDIVV = INTDIVV / SURFAC / 6D0

!     BOUCLE SUR LES EF de BREZZI-FORTIN
      DO NUELEM = 1, NBELEM

!        NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
         NOSOTR(1) = NOSTEF(NUELEM,1)
         NOSOTR(2) = NOSTEF(NUELEM,2)
         NOSOTR(3) = NOSTEF(NUELEM,3)
!        NUMERO DU BARYCENTRE
         NOSOTR(4) = NBSOM + NUELEM

!        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X   = XYZSOM(1,NOSOTR(1))
         X21 = XYZSOM(1,NOSOTR(2)) - X
         X31 = XYZSOM(1,NOSOTR(3)) - X

         Y   = XYZSOM(2,NOSOTR(1))
         Y21 = XYZSOM(2,NOSOTR(2)) - Y
         Y31 = XYZSOM(2,NOSOTR(3)) - Y

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

!        L'INTEGRALE tPression Div Vitesse dX
!        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
!        CAR COMPENSE PAR LA MULTIPLICATION PAR DELTA DE L'INTEGRATION
!        COEFBG Integrale tLambda ( DF-1ligne1  DP1B Ue1
!                                  +DF-1ligne2  DP1B Ue2 
!                                 - INTDIVV/SURFAC ) dX
!        -------------------------------------------------------------
         Y = MOYDIVV * DET
         DO I=1,3

!           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            X = 0D0
            DO J=1,4
!              NO GLOBAL DU NOEUD J DU TRIANGLE
               NS = NOSOTR(J)
!              DP1BLa2d(j,k,i) = integrale dPj/dxk Lambdai dX
               X = X + ( DFM1(1,1) * DP1BLa2d(J,1,I) +
     %                   DFM1(1,2) * DP1BLa2d(J,2,I) ) * VITEXY(NS,1)
     %               + ( DFM1(2,1) * DP1BLa2d(J,1,I) +
     %                   DFM1(2,2) * DP1BLa2d(J,2,I) ) * VITEXY(NS,2)
            ENDDO

!           ASSEMBLAGE DE BE(I) DANS BG(NOSOTR(I))
            NS = NOSOTR(I)
            X  = COEFBG * ( X - Y )
            BG(NS) = BG(NS) + X

         ENDDO

      ENDDO

!      call affvect( 'bgdivbf2.f: BG=', 20,    BG )
!      call afl1ve(  'bgdivbf2.f: BG=', NBSOM, BG )

      RETURN
      END SUBROUTINE BGDIVBF2
