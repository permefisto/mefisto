      INTEGER FUNCTION NOSUI3( I )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NUMERO SUIVANT I DANS LE SENS CIRCULAIRE  1 2 3 1 ...
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C2345X7..............................................................012
      IF( I .EQ. 3 ) THEN
         NOSUI3 = 1
      ELSE
         NOSUI3 = I + 1
      ENDIF
      END
