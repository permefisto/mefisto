      SUBROUTINE LAPLACB( NBSOM,  XYZSOM, NBNOEF, NBEF, NONOEF, NONOSO,
     %                    NBNOVI, NCAS0, NCAS1, NCAS, vitx, vity,
     %                    BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     ASSEMBLER LE SECOND MEMBRE ROT VITESSE INTERPOLEE PAR
C ----     soit BREZZI-FORTIN  P1+BULLE
C          soit TAYLOR-HOOD    P2
C          ET POUR UNE INTERPOLATION P1 DES FONCTIONS TEST
C          Calcul de  Integrale e  tLambda ( dv/dx-du/dy ) dX
C          Integration exacte sur le triangle e
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF NOEUDS DES NBEF EF
C NONOSO : NONOSO(I)=NO DU SOMMET DE 1 A NBSOM DU NOEUD I
C NBNOVI : NOMBRE DE NOEUDS VITESSE DU MAILLAGE
C NCAS0:NCAS1 : LES VECTEURS VITESSE STOCKES
C NCAS   : NUMERO DU CAS A TRAITE
C vitx   : COMPOSANTE X des VECTEURS VITESSE NCAS0:NCAS1
C vitx   : COMPOSANTE Y des VECTEURS VITESSE NCAS0:NCAS1

C MODIFIES:
C ---------
C BG     : COEFFICIENTS DU SECOND MEMBRE GLOBAL ASSEMBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Avril   2011
C MODIFS : ALAIN PERRONNET             St PIERRE du PERRAY  Mai     2021
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NONOEF(NBEF,NBNOEF), NONOSO(NBNOVI),
     %                  NBSOM, NBNOEF, NBEF,
     %                  NBNOVI, NCAS0, NCAS1, NCAS

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity

      DOUBLE PRECISION  BG(NBSOM)
C
      DOUBLE PRECISION  DFM1(2,2), S
      INTEGER           NS, NS1, NS2, NS3, NOSOTR(4), NEF, I, J
      EQUIVALENCE      (NS1,NOSOTR(1)),(NS2,NOSOTR(2)),(NS3,NOSOTR(3))
C
C     DOUBLE PRECISION  DP1BLa2d(4,2,3)
C     DP1BLA2D(i,k,j) = integrale DP1Bi/dxk Lambdaj dX
      include"./incl/dp1bla2d.inc"
C
C     INTEGRALES P1 DP2 SUR LE TRIANGLE de REFERENCE
C     P1DP22D(i,k,j) = integrale P1i DP2j/dxk dx dy sur TRIANGLE UNITE
C     DOUBLE PRECISION  P1DP22D(3,2,6)
      include"./incl/p1dp22d.inc"
C
      DO NEF = 1, NBEF
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
C
C        [DF]-1 SANS DIVISION PAR DELTA COMPENSEE ENSUITE
         DFM1(1,1) = XYZSOM(2,NS3) - XYZSOM(2,NS1)
         DFM1(1,2) = XYZSOM(2,NS1) - XYZSOM(2,NS2)
         DFM1(2,1) = XYZSOM(1,NS1) - XYZSOM(1,NS3)
         DFM1(2,2) = XYZSOM(1,NS2) - XYZSOM(1,NS1)
C
C        Integrale e  tLambda ( dv/dx-du/dy ) dX
C
         IF( NBNOEF .LE. 4 ) THEN
C
C           EF BREZZI-FORTIN avec P1+Bulle
C           avec DP1BLA2D(j,k,i) = integrale DP1Bj/dxk Lambdai dX
            NOSOTR(4) = NBSOM + NEF
            DO I = 1, 3
C
               S = 0D0
               DO J = 1, 4
                  NS = NOSOTR(J)
                  S = S + ( DFM1(1,1) *DP1BLA2D(J,1,I)
     %                     +DFM1(1,2) *DP1BLA2D(J,2,I) )
     %                   * vity(NCAS)%dptab(NS)
     %                  - ( DFM1(2,1) *DP1BLA2D(J,1,I)
     %                     +DFM1(2,2) *DP1BLA2D(J,2,I) )
     %                   * vitx(NCAS)%dptab(NS)
               ENDDO
C
C              ASSEMBLAGE DANS LE VECTEUR GLOBAL
               NS = NOSOTR(I)
               BG(NS) = BG(NS) + S
C
            ENDDO
C
         ELSE
C
C           EF TAYLOR-HOOD avec P2. P1i=Lambdai
C           avec P1DP22D(i,k,j) = integrale P1i DP2j/dxk  dX
            DO I = 1, 3
C
               S = 0D0
               DO J = 1, 6
                  NS = NONOEF(NEF,J)
                  S = S + ( DFM1(1,1) * P1DP22D(I,1,J)
     %                     +DFM1(1,2) * P1DP22D(I,2,J) )
     %                    * vity(NCAS)%dptab(NS)
     %                  - ( DFM1(2,1) * P1DP22D(I,1,J)
     %                     +DFM1(2,2) * P1DP22D(I,2,J) )
     %                    * vitx(NCAS)%dptab(NS)
               ENDDO
C
C              ASSEMBLAGE DANS LE VECTEUR GLOBAL AU SOMMET NS
               NS = NONOSO( NONOEF(NEF,I) )
               BG(NS) = BG(NS) + S
C
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
