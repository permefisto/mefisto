      SUBROUTINE QUATRI( NOSOTR, XYZSOM,  QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA QUALITE D'UN TRIANGLE DANS R**3
C -----   QUALITE: 2 RACINE CARREE(3) x RAYON CERCLE INSCRIT
C                / H LONGUEUR DE LA PLUS LONGUE ARETE
C         L'ORIENTION DES SOMMETS N'INTERVIENT PAS CAR CALCUL DANS R3

C         VERSION AVEC REELS IMPLE PRECISION

C ENTREES:
C --------
C NOSOTR : NUMERO DANS LE TABLEAU XYZSOM DES 3 SOMMETS DU TRIANGLE
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C QUALIT : QUALITE DU TRIANGLE (VALEUR COMPRISE ENTRE 0 ET 1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC       AVRIL 1998
C2345X7..............................................................012
      DOUBLE PRECISION  D2UXR3
      PARAMETER        (D2UXR3 = 3.4641016151377544D0)
C                       D2UXR3 = 2 * SQRT(3)
      INTEGER           NOSOTR(1:3)
      REAL              XYZSOM(3,*), QUALIT
      DOUBLE PRECISION  A, B, C, P
C
C     LA LONGUEUR DES 3 COTES DU TRIANGLE
      NS1 = NOSOTR(1)
      NS2 = NOSOTR(2)
      NS3 = NOSOTR(3)
C
      A = SQRT( (XYZSOM(1,NS2) - XYZSOM(1,NS1)) ** 2
     %        + (XYZSOM(2,NS2) - XYZSOM(2,NS1)) ** 2
     %        + (XYZSOM(3,NS2) - XYZSOM(3,NS1)) ** 2 )
C
      B = SQRT( (XYZSOM(1,NS3) - XYZSOM(1,NS2)) ** 2
     %        + (XYZSOM(2,NS3) - XYZSOM(2,NS2)) ** 2
     %        + (XYZSOM(3,NS3) - XYZSOM(3,NS2)) ** 2 )
C
      C = SQRT( (XYZSOM(1,NS1) - XYZSOM(1,NS3)) ** 2
     %        + (XYZSOM(2,NS1) - XYZSOM(2,NS3)) ** 2
     %        + (XYZSOM(3,NS1) - XYZSOM(3,NS3)) ** 2 )
C
C     LE DEMI-PERIMETRE DU TRIANGLE
      P = (A+B+C) / 2
C
      IF ( A*B*C .NE. 0 ) THEN
C
C        CRITERE: 2*SQRT(3.) * RAYON_INSCRIT / MAX( LONGUEUR DES ARETES )
C                 RAYON_INSCRIT = SQRT( (P-A) (P-B) (P-C) / P )
C                 P = PERIMETRE/2
         QUALIT = REAL( D2UXR3 * SQRT( ABS( (P-A) / P * (P-B) * (P-C) ))
     %                         / MAX(A,B,C) )
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
