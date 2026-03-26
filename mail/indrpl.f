      SUBROUTINE INDRPL( S1, S2, P1, P2, P3, PT, NOCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES DU POINT PT, INTERSECTION DE LA
C -----     DROITE S1-S2 ET DU PLAN DEFINI PAR P1, P2, P3
C
C ENTREES :
C ---------
C S1,S2   : LES 2 POINTS QUI DEFINISSENT LA DROITE
C P1,P2,P3: LES 3 POINTS QUI DEFINISSENT LE PLAN
C
C SORTIES :
C ---------
C PT      : LES 3 COORDONNEES DU POINT D'INTERSECTION SI NOCODE=0
C NOCODE  : 1 SI LA DROITE EST PARALLELE AU PLAN
C           2 SI S1=S2
C           0 SI LE POINT D'INTERSECTION PT A ETE CALCULE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1993
C23456...............................................................012
      DOUBLE PRECISION  S1(3), S2(3), P1(3), P2(3), P3(3), PT(3)
      DOUBLE PRECISION  A, B, C, D, E, COSA, F
C
C     LE CARRE DE LA DISTANCE S1-S2
      E = (S2(1)-S1(1))**2+(S2(2)-S1(2))**2+(S2(3)-S1(3))**2
      IF( E .LE. 0D0 ) THEN
         NOCODE = 2
         RETURN
      ENDIF
C
C     LE VECTEUR NORMAL AU PLAN DES 3 POINTS
      A = ( P2(2) - P1(2) ) * ( P3(3) - P1(3) )
     %  - ( P2(3) - P1(3) ) * ( P3(2) - P1(2) )
      B = ( P2(3) - P1(3) ) * ( P3(1) - P1(1) )
     %  - ( P2(1) - P1(1) ) * ( P3(3) - P1(3) )
      C = ( P2(1) - P1(1) ) * ( P3(2) - P1(2) )
     %  - ( P2(2) - P1(2) ) * ( P3(1) - P1(1) )
      D = A * A + B * B + C * C
C
C     LE COSINUS DE L'ANGLE ENTRE LE VECTEUR NORMAL ET S1-S2
      F = (S2(1)-S1(1))*A + (S2(2)-S1(2))*B + (S2(3)-S1(3))*C
C
C     COSINUS( S1-S2 , NORMALE AU PLAN DU TRIANGLE )
      COSA = F / SQRT(E*D)
C
      IF( ABS(COSA) .LT. 1D-4 ) THEN
C        COSA=1D-5 => A=89.999427 degres
C        COSA=1D-4 => A=89.994270 degres
C        COSA=1D-3 => A=89.942704 degres
C        S1-S2 PARALLELE AU PLAN
         NOCODE = 1
         RETURN
      ENDIF
C
C     LES 3 COORDONNEES DU POINT D'INTERSECTION DROITE-PLAN
      A = ( (P1(1)-S1(1))*A + (P1(2)-S1(2))*B + (P1(3)-S1(3))*C ) / F
      PT(1) = S1(1) + A * ( S2(1) - S1(1) )
      PT(2) = S1(2) + A * ( S2(2) - S1(2) )
      PT(3) = S1(3) + A * ( S2(3) - S1(3) )
      NOCODE= 0
C
      RETURN
      END
