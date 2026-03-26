      DOUBLE PRECISION FUNCTION DIS2ST( DXYZ1, DXYZ2 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA DISTANCE ENTRE 2 POINTS DXYZ1 ET DXYZ2
C -----
C
C ENTREES :
C ---------
C DXYZ1 DXYZ2 : LES 2 POINTS DE R ** 3
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  JANVIER 1986
C....................................................................
      DOUBLE PRECISION  DXYZ1(3),DXYZ2(3)
      INTRINSIC         REAL
C
      DIS2ST = SQRT( ( DXYZ2(1) - DXYZ1(1) ) ** 2 +
     %               ( DXYZ2(2) - DXYZ1(2) ) ** 2 +
     %               ( DXYZ2(3) - DXYZ1(3) ) ** 2 )
C
      RETURN
      END
