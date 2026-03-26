      SUBROUTINE MATAXI
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CONSTRUCTION DE LA MATRICE IDENTITE DANS AXOMAT de
C -----   ~/incl/trvari.inc           VERSION xvue
C SORTIE :
C --------
C AXOMAT(3,3) DANS LE COMMON TRVAR6 DE ./incl/trvari.inc"
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C2345X789............................................................012
      include"./incl/trvari.inc"
C
      DO 20 J=1,3
         DO 10 I=1,3
            AXOMAT(I,J) = 0.0
 10      CONTINUE
         AXOMAT(J,J) = 1.0
 20   CONTINUE
      END
