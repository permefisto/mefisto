      REAL FUNCTION DISTPD( PT , P1DR , P2DR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :       CALCULER LA DISTANCE ENTRE UN POINT ET UNE DROITE
C -----       DEFINIE PAR 2 POINTS P1DR ET P2DR
C             (VOIR AUSSI ptprdr.f)
C
C ENTREES :
C ---------
C PT        : LE POINT DE R ** 3
C P1DR P2DR : LES 2 POINTS DE R ** 3  DE LA DROITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMRIQUE PARIS  JANVIER 1986
C.......................................................................
      REAL    PT(3),P1DR(3),P2DR(3),P0(3)
C
C     LA PROJECTION ORTHOGONALE P0 DE PT SUR LA DROITE VERIFIE
C
C     (X2-X1)(X0-XPT) + (Y2-Y1)(Y0-YPT) + (Z2-Z1)(Z0-ZPT) = 0
C
C      X0-X1   Y0-Y1   Z0-Z1
C      ----- = ----- = -----
C      X2-X1   Y2-Y1   Z2-Z1
C
      A = P2DR(1) - P1DR(1)
      B = P2DR(2) - P1DR(2)
      C = P2DR(3) - P1DR(3)
C
      D = PT(1)   - P1DR(1)
      E = PT(2)   - P1DR(2)
      F = PT(3)   - P1DR(3)
C
      DELTA = 1. / ( A*A + B*B + C*C )
C
C     LE VECTEUR P0-PT AVEC P0 LA PROJECTION DE PT SUR LA DROITE
C     A POUR COMPOSANTES
      P0(1) = ( B * ( E * A - B * D ) + C * ( F * A - C * D ) ) * DELTA
      P0(2) = ( C * ( F * B - C * E ) + A * ( D * B - A * E ) ) * DELTA
      P0(3) = ( A * ( D * C - A * F ) + B * ( E * C - B * F ) ) * DELTA
C
C     LA DISTANCE PT P0
      DISTPD = SQRT( P0(1)**2 + P0(2)**2 + P0(3)**2 )

      RETURN
      END
