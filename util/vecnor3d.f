      SUBROUTINE VECNOR3D( P1, P2, P3, NORMVEC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C GOAL   :  CALCULATE THE THREE COMPONENTS OF THE NORMAL VECTOR
C --------

C INPUTS :
C --------
C P1 P2 P3 : THE 3 COORDINATES OF THE 3 VERTICES

C OUTPUTS:
C --------
C NORMVEC: THE THREE COMPONENTS OF THE NORMAL VECTOR
C          ( P2 - P1 ) x ( P3 - P1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTHOR : Dai-Ni Hsieh National Taiwan University TAIPEI TAIWAN 01/2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY Decembre 2015
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY Octobre  2017
C23456---------------------------------------------------------------012
      DOUBLE PRECISION P1(3), P2(3), P3(3), NORMVEC(3)

      NORMVEC(1) = ( P2(2) - P1(2) ) * ( P3(3) - P1(3) )
     %           - ( P2(3) - P1(3) ) * ( P3(2) - P1(2) )

      NORMVEC(2) = ( P2(3) - P1(3) ) * ( P3(1) - P1(1) )
     %           - ( P2(1) - P1(1) ) * ( P3(3) - P1(3) )

      NORMVEC(3) = ( P2(1) - P1(1) ) * ( P3(2) - P1(2) )
     %           - ( P2(2) - P1(2) ) * ( P3(1) - P1(1) )

      RETURN
      END
