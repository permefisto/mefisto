      DOUBLE PRECISION  FUNCTION DIPTDR( PT , P1DR , P2DR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : CALCULER LA DISTANCE ENTRE UN POINT ET UNE DROITE
C ----- DEFINIE PAR 2 POINTS P1DR ET P2DR
C
C       VERSION DOUBLE PRECISION DE DIS2PD
C
C ENTREES :
C ---------
C PT        : LE POINT DE R ** 2
C P1DR P2DR : LES 2 POINTS DE R ** 2  DE LA DROITE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMRIQUE PARIS  JANVIER 1986
C....................................................................012
      DOUBLE PRECISION  PT(2),P1DR(2),P2DR(2), A, B, C
C
C     LES COEFFICIENTS DE LA DROITE A X + BY + C =0
      A = P2DR(2) - P1DR(2)
      B = P1DR(1) - P2DR(1)
      C = - A * P1DR(1) - B * P1DR(2)
C
C     LA DISTANCE = | A * X + B * Y + C | / SQRT( A*A + B*B )
      DIPTDR = ABS( A * PT(1) + B * PT(2) + C ) / SQRT( A*A + B*B )
      END
