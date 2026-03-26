      SUBROUTINE MOYDIVTH3D( NBNOVI, XYZNOE, NBELEM, NUNOEF, VITXYZ,
     %                       VOLUME, INTDIVV )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CALCULER LE VOLUME DU MAILLAGE ET
! ----     L'INTEGRALE DE DIVERGENCE VITESSE SUR LE MAILLAGE
!          Integrale Div P2 dX VITXYZe

! ENTREES:
! --------
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBELEM : NOMBRE DE TETRAEDRES P2 DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,10) NO DES NOEUDS DES NBELEM EF
! VITXYZ : VITESSE EN X, Y, Z EN LES NBNOVI NOEUDS (INTERPOLATION P2)

! SORTIES:
! --------
! VOLUME : VOLUME DU MAILLAGE 3D
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITXYZ SUR LE MAILLAGE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NBNOVI, NBELEM, NUNOEF(NBELEM,10)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  VITXYZ(NBNOVI,3), VOLUME, INTDIVV
      DOUBLE PRECISION  DF(3,3), DFM1(3,3), DET, X, Y, Z
      INTEGER           NUELEM, I, J, K, L

!     DOUBLE PRECISION  DP23D(3,10)
!     DP23D(k,j) = integrale DP2j/dxk dx dy dz sur le TETRAEDRE UNITE
      include"./incl/dp23d.inc"

      VOLUME  = 0D0
      INTDIVV = 0D0

!     BOUCLE SUR LES TETRAEDRES P2
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         I = NUNOEF(NUELEM,1)
         J = NUNOEF(NUELEM,2)
         K = NUNOEF(NUELEM,3)
         L = NUNOEF(NUELEM,4)

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

!        LE DETERMINANT DE DF
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %       + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %       + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )

!        VOLUME DU TETRAEDRE * 6
         VOLUME = VOLUME + DET
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

         Z = 0D0
         DO K=1,3
            DO I=1,10
               X = VITXYZ( NUNOEF(NUELEM,I), K )
               DO J=1,3
                  Z = Z + DFM1(K,J) * DP23D(J,I) * X
               ENDDO
            ENDDO
         ENDDO
         INTDIVV = INTDIVV + Z

      ENDDO     ! FIN DE LA BOUCLE DES TETRAEDRES DU MAILLAGE

      VOLUME = VOLUME / 6D0

      print *,'moydivth3d.f95: VOLUME=',VOLUME,' INTEGRALE div Vxyz=',
     %  INTDIVV,' INTEGRALE div Vxyz/VOLUME=',INTDIVV/VOLUME

      RETURN
      END SUBROUTINE MOYDIVTH3D
