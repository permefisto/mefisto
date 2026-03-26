      SUBROUTINE VECNOR3( P1, P2, P3, NORMVEC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C GOAL   : CALCULATE THE THREE COMPONENTS OF THE NORMAL VECTOR
C -------- XYZ of POINTS REAL SIMPLE PRECISION
C
C INPUTS :
C --------
C P1 P2 P3 : THE 3 COORDINATES OF THE 3 VERTICES
C
C OUTPUTS:
C --------
C NORMVEC: THE THREE COMPONENTS OF THE NORMAL VECTOR
C          ( P2 - P1 ) x ( P3 - P1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTHOR : Dai-Ni Hsieh National Taiwan University TAIPEI TAIWAN 01/2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Octobre 2017
C23456---------------------------------------------------------------012
      REAL              P1(3), P2(3), P3(3)
      DOUBLE PRECISION  NORMVEC(3), DBLE
C
      NORMVEC(1) = ( DBLE(P2(2)) - P1(2) ) * ( DBLE(P3(3)) - P1(3) ) -
     %             ( DBLE(P2(3)) - P1(3) ) * ( DBLE(P3(2)) - P1(2) )
      NORMVEC(2) = ( DBLE(P2(3)) - P1(3) ) * ( DBLE(P3(1)) - P1(1) ) -
     %             ( DBLE(P2(1)) - P1(1) ) * ( DBLE(P3(3)) - P1(3) )
      NORMVEC(3) = ( DBLE(P2(1)) - P1(1) ) * ( DBLE(P3(2)) - P1(2) ) -
     %             ( DBLE(P2(2)) - P1(2) ) * ( DBLE(P3(1)) - P1(1) )

      RETURN
      END
