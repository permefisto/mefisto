      SUBROUTINE RENVMO( MOT, M )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RENVERSER SUR LUI MEME LE TABLEAU M EN INVERSANT LE SIGNE
C ----- M(1) EST PERMUTE AVEC -M(MOT)
C       M(2) EST PERMUTE AVEC -M(MOT-1) ...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      INTEGER  M(MOT)
C
      DO 10 I = 1,MOT/2
         J = MOT + 1 - I
         K    =  M(I)
         M(I) = -M(J)
         M(J) = -K
 10   CONTINUE
      END
