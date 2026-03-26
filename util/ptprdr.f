      SUBROUTINE PTPRDR( PT, P1DR, P2DR, P0 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :       CALCULER LE POINT P0 PROJECTION DU POINT PT SUR LA DROITE
C -----       DEFINIE PAR LES 2 POINTS P1DR ET P2DR
C             (VOIR AUSSI distpd.f)

C ENTREES :
C ---------
C PT       : LE POINT DE R ** 3
C P1DR P2DR: LES 2 POINTS DE R ** 3 DE DEFINITION DE LA DROITE

C SORTIES :
C ---------
C P0      : LE POINT PROJETE SUR LA DROITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1998
C2345X7..............................................................012
      REAL   PT(3), P1DR(3), P2DR(3), P0(3)

C     LA PROJECTION ORTHOGONALE P0 DE PT SUR LA DROITE VERIFIE

C     (X2-X1)(X0-XPT) + (Y2-Y1)(Y0-YPT) + (Z2-Z1)(Z0-ZPT) = 0

C      X0-X1   Y0-Y1   Z0-Z1
C      ----- = ----- = -----
C      X2-X1   Y2-Y1   Z2-Z1

      A = P2DR(1) - P1DR(1)
      B = P2DR(2) - P1DR(2)
      C = P2DR(3) - P1DR(3)

      D = PT(1)   - P1DR(1)
      E = PT(2)   - P1DR(2)
      F = PT(3)   - P1DR(3)

      DELTA = 1. / ( A*A + B*B + C*C )

C     LE VECTEUR P0-PT AVEC P0 LA PROJECTION DE PT SUR LA DROITE
C     A POUR COMPOSANTES
      P0(1) = ( B * ( E * A - B * D ) + C * ( F * A - C * D ) ) * DELTA
      P0(2) = ( C * ( F * B - C * E ) + A * ( D * B - A * E ) ) * DELTA
      P0(3) = ( A * ( D * C - A * F ) + B * ( E * C - B * F ) ) * DELTA

CCCC     LA DISTANCE PT P0
CCC      DISTPD = SQRT( P0(1)**2 + P0(2)**2 + P0(3)**2 )

C     LE POINT P0 PROJETE DE PT SUR LA DROITE P1DR-P2DR
      P0(1) = P0(1) + PT(1)
      P0(2) = P0(2) + PT(2)
      P0(3) = P0(3) + PT(3)

      RETURN
      END
