      INTEGER FUNCTION NOPRE3( I )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NUMERO PRECEDENT I DANS LE SENS CIRCULAIRE  1 2 3 1 ...
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C2345X7..............................................................012
      IF( I .EQ. 1 ) THEN
         NOPRE3 = 3
      ELSE
         NOPRE3 = I - 1
      ENDIF
      END
