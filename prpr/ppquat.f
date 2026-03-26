      PROGRAM PPQUAT
      REAL    P1(3),P2(3),P3(3),P4(3),SURFTR(4),SQRT3
C
C     LE TETRAEDRE EQUILATERAL
      SQRT3 = SQRT(3.)

      P1(1) = SQRT3/6.
      P1(2) = -1./2.
      P1(3) = 0.
C
      P2(1) = SQRT3/6.
      P2(2) = 1./2.
      P2(3) = 0.
C
      P3(1) = -SQRT3/3.
      P3(2) = 0.
      P3(3) = 0.
C
      P4(1) = 0.
      P4(2) = 0.
      P4(3) = SIN( ACOS( 1.0 / SQRT3 ) )
C
      CALL QUATET( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C
      CALL QUATET0( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C
      PRINT *, ARMIN, ARMAX, SURFTR, VOLUME, QUALIT
C
C     LE TRIANGLE RECTANGLE UNITE
      P1(1) = 0.
      P1(2) = 0.
      P1(3) = 0.
C
      P2(1) = 1.
      P2(2) = 0.
      P2(3) = 0.
C
      P3(1) = 0.
      P3(2) = 1.
      P3(3) = 0.
C
      P4(1) = 0.
      P4(2) = 0.
      P4(3) = 1.
C
      CALL QUATET( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C
      CALL QUATET0( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C
C     LE TRIANGLE QUELCONQUE
      P1(1) = 0.
      P1(2) = 0.
      P1(3) = 0.
C
      P2(1) = 1.
      P2(2) = 0.
      P2(3) = 0.
C
      P3(1) = 0.
      P3(2) = 1.
      P3(3) = 0.
C
      P4(1) = 0.
      P4(2) = 0.2
      P4(3) = 1.5
C
      CALL QUATET( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C
      CALL QUATET0( P1, P2, P3, P4,
     %             ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
      END
      SUBROUTINE QUATET( P1, P2, P3, P4,
     %                   ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES LONGUEURS DES ARETES MIN ET MAX
C -----   SURFACE DES 4 FACES, VOLUME
C         QUALITE: ALFA RAYON INSCRIT / H LONGUEUR DE LA PLUS LONGUE ARETE
C            OU    ALFA = 4/(3*SIN(ARCOS(1/SQRT(3)))
C         (RAPPEL : R = 3 VOLUME TETRA / SURFACE DES 4 FACES)
C         DU TETRAEDRE DE SOMMETS P1 P2 P3 P4
C
C ENTREES:
C --------
C P1,P2,P3,P4 : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C
C SORTIE :
C --------
C ARMIN  : LONGUEUR DE LA PLUS PETITE ARETE DU TETRAEDRE
C ARMAX  : LONGUEUR DE LA PLUS GRANDE ARETE DU TETRAEDRE
C SURFTR : SURFACE DES 4 FACES TRIANGULAIRES DU TETRAEDRE
C          123 234 341 412
C VOLUME : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE)
C                                  <0 SI TETRAEDRE DEGENERE
C                                        OU MAL ORIENTE
C QUALIT : QUALITE DU TETRAEDRE (VALEUR COMPRISE ENTRE 0 ET 1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     JANVIER 1995
C2345X7..............................................................012
      REAL              P1(3), P2(3), P3(3), P4(3),
     %                  ARMIN, ARMAX, SURFTR(4), VOLUME, QUALIT
      DOUBLE PRECISION  ALFA
C
C     VARIABLES AUXILIAIRES
      DOUBLE PRECISION  X1,Y1,Z1,X12,Y12,Z12,X13,Y13,Z13,X14,Y14,Z14
      DOUBLE PRECISION  A,B,C,D,E,F,SURF,SURFTE,VOLU6T
      DOUBLE PRECISION  V1(3),V2(3),V3(3)
C
      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)
C
      X12 = X1 - P2(1)
      Y12 = Y1 - P2(2)
      Z12 = Z1 - P2(3)
      A   = X12 ** 2 + Y12 ** 2 + Z12 ** 2
C
      X13 = X1 - P3(1)
      Y13 = Y1 - P3(2)
      Z13 = Z1 - P3(3)
      B   = X13 ** 2 + Y13 ** 2 + Z13 ** 2
C
      X14 = X1 - P4(1)
      Y14 = Y1 - P4(2)
      Z14 = Z1 - P4(3)
      C   = X14 ** 2 + Y14 ** 2 + Z14 ** 2
C
      ARMIN = MIN( A, B, C )
      ARMAX = MAX( A, B, C )
C
      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13
C
      VOLU6T = X12 * D + X13 * E + X14 * F
C
C     VOLU6T EST 6 FOIS LE VOLUME DU TETRAEDRE
C     ========================================
      IF( VOLU6T .LE. 0.D0 ) THEN
         QUALIT = 0
         GOTO 9999
      ENDIF
C
C     LA SURFACE DES 4 TRIANGLES FACES DU TETRAEDRE
C     =============================================
      SURFTE = 0D0
      DO 90 I=1,4
C
         GOTO( 20, 30, 40, 50 ),I
C
C        LA FACE 1 DU TETRAEDRE
 20      DO 25 J=1,3
            V1 (J) = P3(J) - P1(J)
            V2 (J) = P2(J) - P1(J)
 25      CONTINUE
         GOTO 60
C
C        LA FACE 2 DU TETRAEDRE
 30      A = 0D0
         B = 0D0
         DO 35 J=1,3
            V1(J) = P3(J) - P2(J)
            V2(J) = P4(J) - P2(J)
            A = A + V1(J) ** 2
            B = B + V2(J) ** 2
 35      CONTINUE
         GOTO 60
C
C        LA FACE 3 DU TETRAEDRE
 40      C = 0D0
         DO 45 J=1,3
            V1(J) = P1(J) - P3(J)
            V2(J) = P4(J) - P3(J)
            C = C + V2(J) ** 2
 45      CONTINUE
         GOTO 60
C
C        LA FACE 4 DU TETRAEDRE
 50      DO 55 J=1,3
            V1(J) = P1(J) - P4(J)
            V2(J) = P2(J) - P4(J)
 55      CONTINUE
C
C        PRODUIT VECTORIEL V1 * V2 => V3
 60      DO 70 J=1,3
            IF( J .EQ. 3 ) THEN
               J1 = 1
            ELSE
               J1 = J + 1
            ENDIF
            IF( J1 .EQ. 3 ) THEN
               J2 = 1
            ELSE
               J2 = J1 + 1
            ENDIF
            V3( J ) = V1( J1 ) * V2( J2 ) - V1( J2 ) * V2( J1 )
 70      CONTINUE
C
C        2*SURFACE DE LA FACE I DU TETRAEDRE
         SURF = 0D0
         DO 80 J=1,3
            SURF = SURF + V3(J) ** 2
 80      CONTINUE
         SURF = SQRT( SURF )
C
C        LA SURFACE DU TRIANGLE I DU TETRAEDRE
         SURFTR( I ) = SURF * 0.5D0
C
         SURFTE = SURFTE + SURF
 90   CONTINUE
C
C     ARETES MINIMALE ET MAXIMALE
C     ===========================
      ARMIN = SQRT( MIN( ARMIN, A, B, C ) )
      ARMAX = SQRT( MAX( ARMAX, A, B, C ) )
C
C     CRITERE DE QUALITE : ALFA * RAYON INSCRIT / LONGUEUR DE LA PLUS LONGUE ARE
C     ==========================================================================
      ALFA = 4D0 / SIN( ACOS( 1D0 / SQRT(3.D0) ) )
      PRINT * , 'ALFA = ',ALFA
      QUALIT = ALFA * VOLU6T / ( SURFTE * ARMAX )
C
C     VOLUME EXACT DU TETRAEDRE
C     =========================
 9999 VOLUME = VOLU6T / 6D0
      PRINT *, 'CRITERE R/H'
      PRINT *, ARMIN, ARMAX, SURFTR, VOLUME, ' QUALITE R/H=',QUALIT
      PRINT *
      END
      SUBROUTINE QUATET0( P1, P2, P3, P4,
     %                   ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES LONGUEURS DES ARETES MIN ET MAX
C -----   SURFACE DES 4 FACES, VOLUME
C         QUALITE 3 RAYON INSCRIT / RAYON CIRCONSCRIT
C         DU TETRAEDRE DE SOMMETS P1 P2 P3 P4
C
C ENTREES:
C --------
C P1,P2,P3,P4 : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C
C SORTIE :
C --------
C ARMIN  : LONGUEUR DE LA PLUS PETITE ARETE DU TETRAEDRE
C ARMAX  : LONGUEUR DE LA PLUS GRANDE ARETE DU TETRAEDRE
C SURFTR : SURFACE DES 4 FACES TRIANGULAIRES DU TETRAEDRE
C          123 234 341 412
C VOLUME : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE)
C                                  <0 SI TETRAEDRE DEGENERE
C                                        OU MAL ORIENTE
C QUALIT : QUALITE DU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC    DECEMBRE 1991
C2345X7..............................................................012
      REAL              P1(3),P2(3),P3(3),P4(3),
     %                  ARMIN,ARMAX,SURFTR(4),VOLUME,QUALIT
C
      DOUBLE PRECISION  CENTRE(4)
      DOUBLE PRECISION  X1,Y1,Z1,X12,Y12,Z12,X13,Y13,Z13,X14,Y14,Z14
      DOUBLE PRECISION  A,B,C,D,E,F,G,RAYON,SURF,SURFTE,VOLU6T
      DOUBLE PRECISION  V1(3),V2(3),V3(3)
C
      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)
C
      X12 = X1 - P2(1)
      Y12 = Y1 - P2(2)
      Z12 = Z1 - P2(3)
      A   = X12 ** 2 + Y12 ** 2 + Z12 ** 2
C
      X13 = X1 - P3(1)
      Y13 = Y1 - P3(2)
      Z13 = Z1 - P3(3)
      B   = X13 ** 2 + Y13 ** 2 + Z13 ** 2
C
      X14 = X1 - P4(1)
      Y14 = Y1 - P4(2)
      Z14 = Z1 - P4(3)
      C   = X14 ** 2 + Y14 ** 2 + Z14 ** 2
      ARMIN = MIN( A, B, C )
      ARMAX = MAX( A, B, C )
C
      D = Y13 * Z14 - Y14 * Z13
      E = Y14 * Z12 - Y12 * Z14
      F = Y12 * Z13 - Y13 * Z12
C
      G = X12 * D + X13 * E + X14 * F
C
C     G EST -6 FOIS LE VOLUME DU TETRAEDRE
C     ====================================
      VOLU6T = -G
      IF( VOLU6T .LE. 0.D0 ) THEN
         QUALIT = 0
         GOTO 9999
      ENDIF
C
      G = 0.5D0 / G
C
      A = X12 * ( X1 + P2(1) ) +
     %    Y12 * ( Y1 + P2(2) ) +
     %    Z12 * ( Z1 + P2(3) )
C
      B = X13 * ( X1 + P3(1) ) +
     %    Y13 * ( Y1 + P3(2) ) +
     %    Z13 * ( Z1 + P3(3) )
C
      C = X14 * ( X1 + P4(1) ) +
     %    Y14 * ( Y1 + P4(2) ) +
     %    Z14 * ( Z1 + P4(3) )
C
C     LE CENTRE DE LA BOULE CIRCONSCRITE
C     ==================================
      CENTRE(1) = ( A * D + B * E + C * F ) * G
      CENTRE(2) = ( A * ( X14 * Z13 - X13 * Z14 )
     %            + B * ( X12 * Z14 - X14 * Z12 )
     %            + C * ( X13 * Z12 - X12 * Z13 ) ) * G
      CENTRE(3) = ( A * ( X13 * Y14 - X14 * Y13 )
     %            + B * ( X14 * Y12 - X12 * Y14 )
     %            + C * ( X12 * Y13 - X13 * Y12 ) ) * G
C
C     LE CARRE DU RAYON DE LA BOULE CIRCONSCRITE
C     ==========================================
      CENTRE(4) = ( X1 - CENTRE(1) ) ** 2
     %          + ( Y1 - CENTRE(2) ) ** 2
     %          + ( Z1 - CENTRE(3) ) ** 2
      RAYON = SQRT( CENTRE(4) )
C
C     LA SURFACE DES 4 TRIANGLES FACES DU TETRAEDRE
C     =============================================
      SURFTE = 0D0
      DO 90 I=1,4
         GOTO( 20, 30, 40, 50 ),I
C
C        LA FACE 1 DU TETRAEDRE
 20      DO 25 J=1,3
            V1 (J) = P3(J) - P1(J)
            V2 (J) = P2(J) - P1(J)
 25      CONTINUE
         GOTO 60
C
C        LA FACE 2 DU TETRAEDRE
 30      A = 0D0
         B = 0D0
         DO 35 J=1,3
            V1(J) = P3(J) - P2(J)
            V2(J) = P4(J) - P2(J)
            A = A + V1(J) ** 2
            B = B + V2(J) ** 2
 35      CONTINUE
         GOTO 60
C
C        LA FACE 3 DU TETRAEDRE
 40      C = 0D0
         DO 45 J=1,3
            V1(J) = P1(J) - P3(J)
            V2(J) = P4(J) - P3(J)
            C = C + V2(J) ** 2
 45      CONTINUE
         GOTO 60
C
C        LA FACE 4 DU TETRAEDRE
 50      DO 55 J=1,3
            V1(J) = P1(J) - P4(J)
            V2(J) = P2(J) - P4(J)
 55      CONTINUE
C
C        PRODUIT VECTORIEL V1 * V2 => V3
 60      DO 70 J=1,3
            IF( J .EQ. 3 ) THEN
               J1 = 1
            ELSE
               J1 = J + 1
            ENDIF
            IF( J1 .EQ. 3 ) THEN
               J2 = 1
            ELSE
               J2 = J1 + 1
            ENDIF
            V3( J ) = V1( J1 ) * V2( J2 ) - V1( J2 ) * V2( J1 )
 70      CONTINUE
C
C        2*SURFACE DE LA FACE I DU TETRAEDRE
         SURF = 0D0
         DO 80 J=1,3
            SURF = SURF + V3(J) ** 2
 80      CONTINUE
         SURF = SQRT( SURF )
C
C        LA SURFACE DU TRIANGLE I DU TETRAEDRE
         SURFTR( I ) = SURF * 0.5D0
C
         SURFTE = SURFTE + SURF
 90   CONTINUE
C
C     ARETES MINIMALE ET MAXIMALE
C     ===========================
      ARMIN = SQRT( MIN( ARMIN, A, B, C ) )
      ARMAX = SQRT( MAX( ARMAX, A, B, C ) )
C
C     CRITERE DE QUALITE : 3 * RAYON INSCRIT / RAYON CIRCONSCRIT
C     ==========================================================
      QUALIT  = 3D0 * VOLU6T / ( SURFTE * RAYON )
C
C     VOLUME EXACT DU TETRAEDRE
 9999 VOLUME = VOLU6T / 6D0
      PRINT *, 'CRITERE R/R'
      PRINT *, ARMIN, ARMAX, SURFTR, VOLUME, ' QUALITE R/R=',QUALIT
      PRINT *, 'RAYON_CIR=', RAYON
      PRINT *
      END
