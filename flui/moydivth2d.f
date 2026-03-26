      SUBROUTINE MOYDIVTH2D( NBNOVI, XYZNOE, NBELEM, NUNOEF, VITXY,
     %                       SURFACE, INTDIVV )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CALCULER LA SURFACE DU MAILLAGE 2D ET
! ----     L'INTEGRALE DE LA DIVERGENCE DE LA VITESSE SUR LE MAILLAGE
!          Integrale Div P2 dX VITXYe

! ENTREES:
! --------
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBELEM : NOMBRE DE TRIANGLES P2 DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,6) NO DES NOEUDS DES NBELEM EF
! VITXY  : VITESSE EN X, Y EN LES NBNOVI NOEUDS (INTERPOLATION P2)

! SORTIES:
! --------
! SURFACE: SURFACE DU MAILLAGE 2D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITXY SUR LE MAILLAGE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NBNOVI, NBELEM, NUNOEF(NBELEM,6)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  VITXY(NBNOVI,3), SURFACE, INTDIVV
      DOUBLE PRECISION  DFM1(2,2), DET, X, Y, X21, Y21, X31, Y31
      INTEGER           NUELEM, I, J, K

!     DOUBLE PRECISION  DP22D(2,6)
!     DP22D(k,j) = integrale DP2j/dxk dx dy sur le TRIANGLE UNITE
      include"./incl/dp22d.inc"

      SURFACE = 0D0
      INTDIVV = 0D0

!     BOUCLE SUR LES TRIANGLES P2
      DO NUELEM = 1, NBELEM

!        NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
         I = NUNOEF(NUELEM,1)
         J = NUNOEF(NUELEM,2)
         K = NUNOEF(NUELEM,3)

!        CONSTRUCTION DE LA MATRICE DF-1
         X = XYZNOE(1,I)
         Y = XYZNOE(2,I)
         X21 = (XYZNOE(1,J) - X)
         Y21 = (XYZNOE(2,J) - Y)
         X31 = (XYZNOE(1,K) - X)
         Y31 = (XYZNOE(2,K) - Y)

!        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DET = ABS( X21 * Y31 - X31 * Y21 )

!        SURFACE DU TRIANGLE * 2
         SURFACE = SURFACE + DET
!
!        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
!        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DET
!        CAR COMPENSE PAR LA MULTIPLICATION PAR DET DE L'INTEGRATION
         DFM1(1,1) =  Y31
         DFM1(2,1) = -X31
         DFM1(1,2) = -Y21
         DFM1(2,2) =  X21

         Y = 0D0
         DO K=1,2
            DO I=1,6
               X = VITXY( NUNOEF(NUELEM,I), K )
               DO J=1,2
                  Y = Y + DFM1(K,J) * DP22D(J,I) * X
               ENDDO
            ENDDO
         ENDDO
         INTDIVV = INTDIVV + Y

      ENDDO     ! FIN DE LA BOUCLE DES TRIANGLES DU MAILLAGE

      SURFACE = SURFACE / 2D0

      print *,'moydivth2d: SURFACE=',SURFACE,' INTEGRALE div Vxy=',
     % INTDIVV,' INTEGRALE div Vxy/Surface=',INTDIVV/SURFACE

      RETURN
      END SUBROUTINE MOYDIVTH2D
