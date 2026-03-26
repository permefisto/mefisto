      SUBROUTINE MATAXI6
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CONSTRUCTION DE LA MATRICE IDENTITE DANS AXOMAT de
C -----   ~/incl/trvari.inc           VERSION xvue
C SORTIE :
C --------
C AXOMAT(6,6) DANS LE COMMON TRVAR4 DE ./incl/trvari.inc"
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C2345X789............................................................012
      include"./incl/trvari.inc"
C
      DO 20 J=1,6
         DO 10 I=1,6
            AXOMAT(I,J) = 0.0
 10      CONTINUE
         AXOMAT(J,J) = 1.0
 20   CONTINUE
      END
