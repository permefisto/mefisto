SUBROUTINE MOYDIVBF3D( NBSOM,  XYZSOM, NBELEM, NUSOTE, NBNOVI, VITXYZ, & 
                       VOLUME, INTDIVV )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CALCULER LE VOLUME DU MAILLAGE ET
! ----     L'INTEGRALE DE LA DIVERGENCE DE LA VITESSE SUR LE MAILLAGE
!          Integrale Div P1B dX VITXYZe

! ENTREES:
! --------
! NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
! XYZSOM : XYZSOM(3,NBSOM)  3 COORDONNEES DES SOMMETS DU MAILLAGE
! NBELEM : NOMBRE DE TETRAEDRES DU MAILLAGE
! NUSOTE : NUSOTE(NBELEM,4) NO DES 4 SOMMETS DES NBELEM TETRAEDRES
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRE DES TETRAEDRES
!          NBSOM + NBELEM
! VITXYZ : VITESSE EN X, Y, Z EN LES NBNOVI NOEUDS (INTERPOLATION P1+BULLE)

! SORTIES:
! --------
! VOLUME : VOLUME DU MAILLAGE 3D DES TETRAEDRES
! INTDIVV: INTEGRALE DE LA DIVERGENCE DE VITXYZ SUR LE MAILLAGE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
   IMPLICIT NONE
      INTEGER           NBSOM, NBNOVI, NBELEM, NUSOTE(NBELEM,4)
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  VITXYZ(NBNOVI,3), VOLUME, INTDIVV
      DOUBLE PRECISION  DF(3,3), DFM1(3,3), DET, X, Y, Z
      INTEGER           NUELEM, I, J, K, NUNO1T(5)

!     DOUBLE PRECISION  DP1B3D(3,5)
!     DP1B3D(k,j) = integrale DP1Bj/dxk dx dy dz sur le TETRAEDRE UNITE
      include"./incl/dp1b3d.inc95"

      VOLUME  = 0D0
      INTDIVV = 0D0

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NUELEM, I, J, K, NUNO1T, DF, DFM1, DET, X, Y, Z)
!     BOUCLE SUR LES TETRAEDRES P2
!$OMP DO REDUCTION(+:VOLUME) REDUCTION(+:INTDIVV)
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         DO I=1,4
            NUNO1T(I) = NUSOTE(NUELEM,I)
         ENDDO
!        NUMERO DU BARYCENTRE
         NUNO1T(5) = NBSOM + NUELEM

!        CONSTRUCTION DE LA MATRICE DF
         I = NUNO1T(1)
         X = XYZSOM(1,I)
         Y = XYZSOM(2,I)
         Z = XYZSOM(3,I)

         DF(1,1) = XYZSOM(1,NUNO1T(2)) - X
         DF(2,1) = XYZSOM(1,NUNO1T(3)) - X
         DF(3,1) = XYZSOM(1,NUNO1T(4)) - X

         DF(1,2) = XYZSOM(2,NUNO1T(2)) - Y
         DF(2,2) = XYZSOM(2,NUNO1T(3)) - Y
         DF(3,2) = XYZSOM(2,NUNO1T(4)) - Y

         DF(1,3) = XYZSOM(3,NUNO1T(2)) - Z
         DF(2,3) = XYZSOM(3,NUNO1T(3)) - Z
         DF(3,3) = XYZSOM(3,NUNO1T(4)) - Z

!        LE DETERMINANT DE DF
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) &
             + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) ) &
             + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )

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
            DO I=1,5
               X = VITXYZ( NUNO1T(I), K )
               DO J=1,3
                  Z = Z + DFM1(K,J) * DP1B3D(J,I) * X
               ENDDO
            ENDDO
         ENDDO
         INTDIVV = INTDIVV + Z

      ENDDO     ! FIN DE LA BOUCLE DES TETRAEDRES DU MAILLAGE
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

      VOLUME = VOLUME / 6D0

  print *,'moydivbf3d.f95: VOLUME=',VOLUME,' INTEGRALE div Vxyz=',INTDIVV,&
      ' INTEGRALE div Vxyz/VOLUME=',INTDIVV/VOLUME


      RETURN
END SUBROUTINE MOYDIVBF3D
