      SUBROUTINE QUATRID3S( XYZS1, XYZS2, XYZS3,  QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA QUALITE D'UN TRIANGLE DANS R**3 (DOUBLE PRECISION)
C -----   QUALITE: 2 RACINE(3) x RAYON CERCLE INSCRIT
C                / H LONGUEUR DE LA PLUS LONGUE ARETE

C         VERSION AVEC REELS DOUBLE PRECISION

C ENTREES:
C --------
C XYZS1  : 3 COORDONNEES DU SOMMET 1 DU TRIANGLE
C XYZS2  : 3 COORDONNEES DU SOMMET 2 DU TRIANGLE
C XYZS3  : 3 COORDONNEES DU SOMMET 3 DU TRIANGLE

C SORTIE :
C --------
C QUALIT : QUALITE DU TRIANGLE (VALEUR COMPRISE ENTRE 0 ET 1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC       AVRIL 1998
C2345X7..............................................................012
      DOUBLE PRECISION  D2UXR3
      PARAMETER        (D2UXR3 = 3.4641016151377544D0)
C                       D2UXR3 = 2 * SQRT(3)
      DOUBLE PRECISION  XYZS1(3), XYZS2(3), XYZS3(3), QUALIT, A, B, C, P
C
C     LA LONGUEUR DES 3 COTES DU TRIANGLE
      A = SQRT( (XYZS2(1) - XYZS1(1)) ** 2
     %        + (XYZS2(2) - XYZS1(2)) ** 2
     %        + (XYZS2(3) - XYZS1(3)) ** 2 )
C
      B = SQRT( (XYZS3(1) - XYZS2(1)) ** 2
     %        + (XYZS3(2) - XYZS2(2)) ** 2
     %        + (XYZS3(3) - XYZS2(3)) ** 2 )
C
      C = SQRT( (XYZS1(1) - XYZS3(1)) ** 2
     %        + (XYZS1(2) - XYZS3(2)) ** 2
     %        + (XYZS1(3) - XYZS3(3)) ** 2 )
C
C     LE DEMI-PERIMETRE DU TRIANGLE
      P = (A+B+C) / 2
C
      IF ( A*B*C .NE. 0 ) THEN
C
C        CRITERE: 2*SQRT(3.) * RAYON_INSCRIT / MAX( LONGUEUR DES ARETES )
C                 RAYON_INSCRIT = SQRT( (P-A) (P-B) (P-C) / P )
C                 P = PERIMETRE/2
         QUALIT = D2UXR3 * SQRT( ABS( (P-A) / P * (P-B) * (P-C) ) )
     %          / MAX(A,B,C)
C
      ELSE
C
C        TRIANGLE DEGENERE
         QUALIT = 0.0
C
      ENDIF
C
C     AUTRES CRITERES POSSIBLES DE QUALITE DU TRIANGLE
C     ================================================
C     CRITERE : 2 * RAYON_INSCRIT / RAYON_CIRCONSCRIT
C     QUALIT  = 8 * (XP-A)*(XP-B)*(XP-C) / P
C
C     CRITERE : 3*SQRT(3.) * RAY_INSCRIT / DEMI PERIMETRE
C     QUALIT  = 3*SQRT(3.) * SQRT ((XP-A)*(XP-B)*(XP-C) / XP**3)
C
C     CRITERE : 2*SQRT(3.) * RAY_INSCRIT / MAX( DES ARETES )
C     QUALIT  = 2*SQRT(3.) * SQRT( (XP-A)*(XP-B)*(XP-C) / XP ) / MAX( A, B, C )
C
      RETURN
      END
