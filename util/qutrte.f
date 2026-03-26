      SUBROUTINE QUTRTE( P1, P2, P3, QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      CALCUL DE LA QUALITE DU TRIANGLE P1 P2 P3
C -----
C ENTREES:
C --------
C P1,P2,P3 : LES 3 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C
C SORTIE :
C --------
C QUALIT :   QUALITE DU TRIANGLE  2 RAC(3) r / PLUS LONGUE ARETE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     JANVIER 1995
C2345X7..............................................................012
      DOUBLE PRECISION  D2UXR3
      PARAMETER        (D2UXR3 = 3.4641016151377544D0)
C                       D2UXR3 = 2 * SQRT(3)
      DOUBLE PRECISION  P1(3),P2(3),P3(3), A,B,C,P
      INTRINSIC         SQRT
C
      A = (P2(1) - P1(1)) ** 2
     %  + (P2(2) - P1(2)) ** 2
     %  + (P2(3) - P1(3)) ** 2
C
      B = (P3(1) - P2(1)) ** 2
     %  + (P3(2) - P2(2)) ** 2
     %  + (P3(3) - P2(3)) ** 2
C
      C = (P1(1) - P3(1)) ** 2
     %  + (P1(2) - P3(2)) ** 2
     %  + (P1(3) - P3(3)) ** 2
C
      A = SQRT ( A )
      B = SQRT ( B )
      C = SQRT ( C )
C
      P = ( A + B + C ) * 0.5D0
C
      IF ( A*B*C .NE. 0D0 ) THEN
C        CRITERE : 2 RACINE(3) * RAYON_INSCRIT / PLUS LONGUE ARETE
         QUALIT = REAL( D2UXR3 * SQRT( ABS( (P-A) / P * (P-B) * (P-C) ))
     %                  / MAX(A,B,C) )
      ELSE
         QUALIT  = 0.0
      ENDIF
      RETURN
      END
