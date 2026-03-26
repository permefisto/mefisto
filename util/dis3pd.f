      REAL FUNCTION DIS3PD( PT , P1DR , P2DR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    CALCULER LA DISTANCE ENTRE UN POINT ET UNE DROITE
C -----    DEFINIE PAR 2 POINTS P1DR ET P2DR  DE R**3
C
C ENTREES :
C ---------
C PT        : LE POINT DE R ** 3
C P1DR P2DR : LES 2 POINTS DE R ** 3  DE LA DROITE
C
C SORTIE :
C --------
C DIS3PD : DISTANCE DU POINT PT A LA DROITE P1DR P2DR
C           0 SI PT SUR LA DROITE
C          -1 SI DROITE REDUITE A UN POINT P1DR=P2DR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  PT(3), P1DR(3), P2DR(3)
      DOUBLE PRECISION  X, Y, Z, LAMBDA
      INTRINSIC         REAL, SQRT
C
C     LE POINT SUR P1DR P2DR  EST TEL QUE PQ ORTHOGONAL A P1DR P2DR
C     RECHERCHE DE LA PLUS GRANDE DIFFERENCE X2-X1 Y2-Y1 Z2-Z1
      IF( ABS(P2DR(1)-P1DR(1)) .GE. ABS(P2DR(2)-P1DR(2)) ) THEN
         I1 = 1
      ELSE
         I1 = 2
      ENDIF
      IF( ABS(P2DR(3)-P1DR(3)) .GT. ABS(P2DR(I1)-P1DR(I1)) ) THEN
         I1 = 3
      ENDIF
      IF( ABS(P2DR(I1)-P1DR(I1)) .LE. 0.0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'DIS3PD: DROITE REDUITE A UN POINT'
         ELSE
            WRITE(IMPRIM,*) 'DIS3PD: LINE REDUCED TO ONE POINT'
         ENDIF
         DIS3PD = -1.0
         RETURN
      ENDIF
C
C     CALCUL DU POINT ORTHOGONAL SUR LA DROITE
      IF( I1 .LT. 3 ) THEN
         I2 = I1 + 1
      ELSE
         I2 = 1
      ENDIF
      IF( I2 .LT. 3 ) THEN
         I3 = I2 + 1
      ELSE
         I3 = 1
      ENDIF
C
      X = P2DR(I1) - P1DR(I1)
      Y = P2DR(I2) - P1DR(I2)
      Z = P2DR(I3) - P1DR(I3)
C
      X = ( PT(I1) * (X**2) + P1DR(I1) * ( Y**2 + Z**2 )
     %      + X * ( (PT(I2)-P1DR(I2))*Y + (PT(I3)-P1DR(I3))*Z ) )
     %  / ( X**2 + Y**2 + Z**2 )
C
      LAMBDA = ( X - P1DR(I1) ) / ( P2DR(I1) - P1DR(I1) )
C
      Y = P1DR(I2) + LAMBDA * Y
C
      Z = P1DR(I3) + LAMBDA * Z
C
C     LA DISTANCE DE PT A LA DROITE P1DR-P2DR
      DIS3PD = REAL( SQRT((X-PT(I1))**2 +(Y-PT(I2))**2 +(Z-PT(I3))**2) )
C
      RETURN
      END
