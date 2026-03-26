      REAL FUNCTION DIST2P( PT1 , PT2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA DISTANCE ENTRE 2 POINTS PT1 ET PT2
C -----
C
C ENTREES :
C ---------
C PT1 PT2 : LES 2 POINTS DE R ** 3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1986
C23456---------------------------------------------------------------012
      REAL    PT1(3),PT2(3)
C
      DIST2P = SQRT( (PT2(1)-PT1(1))**2 +
     %               (PT2(2)-PT1(2))**2 +
     %               (PT2(3)-PT1(3))**2 )
      END
