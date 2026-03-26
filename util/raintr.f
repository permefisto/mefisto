      SUBROUTINE RAINTR( P1, P2, P3, R, S )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      CALCULER LE RAYON DU CERCLE INSCRIT DU TRIANGLE PLAN
C -----
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C            SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C R       : RAYON DU CERCLE INSCRIT AU TRIANGLE SI SON ORIENTATION
C           EST DIRECTE  ET 0 SINON
C S       : 2 FOIS LA SURFACE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
C     LA LONGUEUR DES 3 COTES
      A = SQRT( (P2(1)-P1(1))**2 + (P2(2)-P1(2))**2 )
      B = SQRT( (P3(1)-P2(1))**2 + (P3(2)-P2(2))**2 )
      C = SQRT( (P1(1)-P3(1))**2 + (P1(2)-P3(2))**2 )
C
C     2 FOIS LA SURFACE DU TRIANGLE = 2S
      S = (P2(1)-P1(1)) * (P3(2)-P1(2)) - (P2(2)-P1(2)) * (P3(1)-P1(1))
      IF( S .LE. 0 ) THEN
C
C        TRIANGLE DEGENERE OU INDIRECT
         R = 0.0
      ELSE
C
C        RAYON DU CERCLE INSCRIT AU TRIANGLE  R = S / P ET P=(A+B+C)/2
         R = S / ( A + B + C )
      ENDIF
      END
