      SUBROUTINE QUATR2D( NOSOTR, PXYD,  QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA QUALITE D'UN TRIANGLE DANS R**2
C -----   QUALITE: 2 RACINE(3) x RAYON CERCLE INSCRIT
C                / H LONGUEUR DE LA PLUS LONGUE ARETE
C
C ENTREES:
C --------
C NOSOTR : NUMERO DES 3 SOMMETS DANS PXYD DES SOMMETS DU TRIANGLE
C PXYD   : 2 COORDONNEES DES SOMMETS DU MAILLAGE + DISTANCE SOUHAITEE
C
C SORTIE :
C --------
C QUALIT : QUALITE DU TRIANGLE (VALEUR COMPRISE ENTRE 0 ET 1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY    AVRIL 2008
C2345X7..............................................................012
      DOUBLE PRECISION  D2UXR3
      PARAMETER        (D2UXR3 = 3.4641016151377544D0)
C                       D2UXR3 = 2 * SQRT(3)
      INTEGER           NOSOTR(1:3)
      DOUBLE PRECISION  PXYD(3,*), A, B, C, P
      REAL              QUALIT
C
C     LA LONGUEUR DES 3 COTES DU TRIANGLE
      NS1 = NOSOTR(1)
      NS2 = NOSOTR(2)
      NS3 = NOSOTR(3)
C
      A = SQRT( (PXYD(1,NS2) - PXYD(1,NS1)) ** 2
     %        + (PXYD(2,NS2) - PXYD(2,NS1)) ** 2 )
C
      B = SQRT( (PXYD(1,NS3) - PXYD(1,NS2)) ** 2
     %        + (PXYD(2,NS3) - PXYD(2,NS2)) ** 2 )
C
      C = SQRT( (PXYD(1,NS1) - PXYD(1,NS3)) ** 2
     %        + (PXYD(2,NS1) - PXYD(2,NS3)) ** 2 )
C
C     LE DEMI-PERIMETRE DU TRIANGLE
      P = (A+B+C) / 2
C
      IF ( A*B*C .NE. 0 ) THEN
C
C        CRITERE: 2*SQRT(3.) * RAYON_INSCRIT / MAX( LONGUEUR DES ARETES )
C                 RAYON_INSCRIT = SQRT( (P-A) (P-B)(P-C) / P )
C                 P = PERIMETRE/2
         QUALIT = REAL( D2UXR3 * SQRT( ABS( (P-A) / P * (P-B) * (P-C) ))
     %                  / MAX(A,B,C) )
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
