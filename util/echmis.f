      REAL FUNCTION ECHMIS( X, XMIN, XMXMI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORME X EN L'ENTIER LE PLUS PROCHE DANS [0,2**20]
C -----  XMIN=>0 , XMXMI+XMIN=2**20
C
C ENTREES:
C --------
C XMIN   : VAUT ZERO APRES TRANSFORMATION
C XMXMI  : XMIN + XMXMI VAUT 2**20 APRES TRANSFORMATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1994
C2345X7..............................................................012
      PARAMETER (MAXINT=2**20)
C
      ECHMIS = XMIN + X * XMXMI / MAXINT
      END
