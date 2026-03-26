      SUBROUTINE AZEROR ( L , RTAB )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISATION A ZERO D UN TABLEAU RTAB DE L VARIABLES REELLES
C -----
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     FEVRIER 1980
C ......................................................................
      REAL    RTAB(L)
      DO 1 I = 1 , L
         RTAB( I ) = 0.0
    1 CONTINUE
      END
