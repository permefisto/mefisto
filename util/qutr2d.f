      SUBROUTINE QUTR2D( P1, P2, P3, QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LA QUALITE D'UN TRIANGLE DE R**2
C -----     2 COORDONNEES DES 3 SOMMETS EN DOUBLE PRECISION
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 3 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C            SENS DIRECT POUR UNE SURFACE ET QUALITE >0
C SORTIES :
C ---------
C QUALIT : VALEUR DE LA QUALITE DU TRIANGLE ENTRE 0 ET 1 (EQUILATERAL)
C          1 ETANT LA QUALITE OPTIMALE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1995
C2345X7..............................................................012
      DOUBLE PRECISION  D2UXR3
      PARAMETER        (D2UXR3 = 3.4641016151377544D0 )
C                       D2UXR3 = 2 * SQRT(3)
      DOUBLE PRECISION  P1(2), P2(2), P3(2), A, B, C, P
      INTRINSIC         SQRT
C
C     LA LONGUEUR DES 3 COTES
      A = SQRT( (P2(1)-P1(1))**2 + (P2(2)-P1(2))**2 )
      B = SQRT( (P3(1)-P2(1))**2 + (P3(2)-P2(2))**2 )
      C = SQRT( (P1(1)-P3(1))**2 + (P1(2)-P3(2))**2 )
C
C     DEMI PERIMETRE
      P = (A+B+C) * 0.5D0
C
      IF ( (A*B*C) .NE. 0D0 ) THEN
C        CRITERE : 2 RACINE(3) * RAYON_INSCRIT / PLUS LONGUE ARETE
         QUALIT = REAL( D2UXR3 * SQRT( ABS( (P-A) / P * (P-B) * (P-C) ))
     %          / MAX(A,B,C) )
      ELSE
         QUALIT  = 0.
      ENDIF
C
C
C     AUTRES CRITERES POSSIBLES
C     CRITERE : 2 * RAYON_INSCRIT / RAYON_CIRCONSCRIT
C     QUALIT  = 8D0 * (P-A) * (P-B) * (P-C) / (A * B * C)
C     CRITERE : 3*SQRT(3.) * RAY_INSCRIT / DEMI PERIMETRE
C     QUALIT  = 3*SQRT(3.) * SQRT ((P-A)*(P-B)*(P-C) / P**3)
C     CRITERE : 2*SQRT(3.) * RAY_INSCRIT / MAX( DES ARETES )
C     QUALIT  = 2*SQRT(3.) * SQRT( (P-A)*(P-B)*(P-C) / P ) / MAX(A,B,C)
      RETURN
      END
