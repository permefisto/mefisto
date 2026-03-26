      SUBROUTINE CEBOCI( P1, P2, P3, P4,  CENTRE, VOLUTE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES 3 COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C -----   AU TETRAEDRE P1 P2 P3 P4 ET DU CARRE DE SON RAYON
C         CALCUL DU VOLUME ORIENTE DU TETRAEDRE P1P2P3P4
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C CENTRE : 3 COORDONNEES DU CENTRE ET CARRE DU RAYON DE LA BOULE CIRCONSCRITE
C          0,0,0,-1 SI TETRAEDRE INCORRECT
C VOLUTE : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE
C                                  =0 SI TETRAEDRE DEGENERE
C                                  <0 SI TETRAEDRE MAL ORIENTE)
C IERR   :  0 SI CALCUL CORRECT SANS PROBLEME RENCONTRE
C          +1 SI TETRAEDRE DEGENERE
C          -1 SI TETRAEDRE MAL ORIENTE DE VOLUME NEGATIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS  MAI 1991
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  P1(3),P2(3),P3(3),P4(3),CENTRE(4),VOLUTE
      DOUBLE PRECISION  X1,Y1,Z1,X12,Y12,Z12,X13,Y13,Z13,X14,Y14,Z14,
     %                  A,B,C,D,E,F,G,VOLU6T

      IERR = 0

      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)

      X12 = X1 - P2(1)
      Y12 = Y1 - P2(2)
      Z12 = Z1 - P2(3)

      X13 = X1 - P3(1)
      Y13 = Y1 - P3(2)
      Z13 = Z1 - P3(3)

      X14 = X1 - P4(1)
      Y14 = Y1 - P4(2)
      Z14 = Z1 - P4(3)

      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13

C     VOLU6T EST 6 FOIS LE VOLUME DU TETRAEDRE
C     VOLU6T EST LE DETERMINANT DE LA MATRICE [P2-P1, P3-P1, P4-P1]
      VOLU6T = X12 * D + X13 * E + X14 * F

C     LE VOLUME DU TETRAEDRE
      VOLUTE = VOLU6T / 6D0

      IF( VOLU6T .LT. 0.D0 ) THEN

         NBLGRC(NRERR) = 1
         KERR(1) = 'ceboci: TETRAEDRE MAL ORIENTE de VOLUME<0'
         CALL LEREUR
         WRITE(IMPRIM,10020) VOLUTE, P1, P2, P3, P4
10020    FORMAT('ceboci: VOLUME DU TETRAEDRE =',G25.16,'<0'/
     %         ('  X=',G25.16,'   Y=',G25.16,'   Z=',G25.16))
         IERR = -1
         GOTO 9000

ccc      ELSE IF( VOLU6T .LE. SQRT(A*B*C)*1D-3 ) THEN   19/8/2014

      ELSE IF( VOLU6T .LE. 0D0 ) THEN

         NBLGRC(NRERR) = 1
         KERR(1) = 'ceboci: TETRAEDRE DEGENERE'
         CALL LEREUR
         WRITE(IMPRIM,10020) VOLUTE, P1, P2, P3, P4
         IERR = 1
         GOTO 9000

      ENDIF

      A = X12 * (P1(1)+P2(1)) + Y12*(P1(2)+P2(2)) + Z12 * (P1(3)+P2(3))
      B = X13 * (P1(1)+P3(1)) + Y13*(P1(2)+P3(2)) + Z13 * (P1(3)+P3(3))
      C = X14 * (P1(1)+P4(1)) + Y14*(P1(2)+P4(2)) + Z14 * (P1(3)+P4(3))

      G = -0.5D0 / VOLU6T

C     XYZ DU CENTRE DE LA SPHERE CIRCONSCRITE
      CENTRE(1) = ( A * D + B * E + C * F ) * G

      CENTRE(2) = ( A * ( X14 * Z13 - X13 * Z14 )
     %            + B * ( X12 * Z14 - X14 * Z12 )
     %            + C * ( X13 * Z12 - X12 * Z13 ) ) * G

      CENTRE(3) = ( A * ( X13 * Y14 - X14 * Y13 )
     %            + B * ( X14 * Y12 - X12 * Y14 )
     %            + C * ( X12 * Y13 - X13 * Y12 ) ) * G

C     LE CARRE DU RAYON DE LA BOULE CIRCONSCRITE
      CENTRE(4) = ( X1 - CENTRE(1) ) ** 2
     %          + ( Y1 - CENTRE(2) ) ** 2
     %          + ( Z1 - CENTRE(3) ) ** 2

C     VERIFICATION: A SUPPRIMER ENSUITE
      A = ( P2(1) - CENTRE(1) ) ** 2
     %  + ( P2(2) - CENTRE(2) ) ** 2
     %  + ( P2(3) - CENTRE(3) ) ** 2

      B = ( P3(1) - CENTRE(1) ) ** 2
     %  + ( P3(2) - CENTRE(2) ) ** 2
     %  + ( P3(3) - CENTRE(3) ) ** 2

      C = ( P4(1) - CENTRE(1) ) ** 2
     %  + ( P4(2) - CENTRE(2) ) ** 2
     %  + ( P4(3) - CENTRE(3) ) ** 2

      PRINT *,'ceboci: CENTRE=',(CENTRE(I),I=1,3),
     %        ' RAYON**2=',CENTRE(4),A,B,C

      RETURN

C     TETRAEDRE DE VOLUME <= 0
C     CENTRE XYZ NULS et CARRE DU RAYON =-1 POUR EVITER DIVISION/R=0D0
 9000 CENTRE(1) = 0D0
      CENTRE(2) = 0D0
      CENTRE(3) = 0D0
      CENTRE(4) =-1D0

      RETURN
      END
