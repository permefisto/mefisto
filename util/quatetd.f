      SUBROUTINE QUATETD( P1, P2, P3, P4,
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES LONGUEURS DES ARETES MIN ET MAX
C -----   SURFACE DES 4 FACES, VOLUME
C         QUALITE: ALFA RAYON INSCRIT / H LONGUEUR DE LA PLUS LONGUE ARETE
C            OU    ALFA = 4/(SIN(ARCOS(1/SQRT(3)))
C         (RAPPEL : r = 3 VOLUME TETRA / SURFACE DES 4 FACES)
C         DU TETRAEDRE DE SOMMETS P1 P2 P3 P4

C ENTREES:
C --------
C P1,P2,P3,P4 : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE

C SORTIES:
C --------
C ARMIN  : LONGUEUR DE LA PLUS PETITE ARETE DU TETRAEDRE
C ARMAX  : LONGUEUR DE LA PLUS GRANDE ARETE DU TETRAEDRE
C SURFTR : SURFACE DES 4 FACES TRIANGULAIRES DU TETRAEDRE
C          123 234 341 412
C VOLUTE : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE)
C                                  <0 SI TETRAEDRE DEGENERE
C                                        OU MAL ORIENTE
C QUALTE : QUALITE DU TETRAEDRE 
C          VOLUME<0    Q=-1  SOMMET DE L'AUTRE COTE DE LA FACE OPPOSEE
C          VOLUME=0    Q= 0  4 SOMMETS DANS UN MEME PLAN
C          VOLUME>0  0<Q<=1  TETRAEDRE CORRECT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     Janvier 1995
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2016
C2345X7..............................................................012
      DOUBLE PRECISION  ALFA
      PARAMETER        (ALFA=4.89897948556636D0)
C     PARAMETER        (ALFA= 4D0 / SIN( ACOS( 1D0 / SQRT(3D0) ) )

      REAL              QUALTE
      DOUBLE PRECISION  P1(3), P2(3), P3(3), P4(3),
     %                  ARMIN, ARMAX, SURFTR(4), VOLUTE

C     VARIABLES AUXILIAIRES
      DOUBLE PRECISION  X1,Y1,Z1,X12,Y12,Z12,X13,Y13,Z13,X14,Y14,Z14,
     %                  A,B,C,D,E,F,SURF,SURFTE,VOLU6T,
     %                  V1(3),V2(3),V3(3)
      INTRINSIC         REAL, MIN, MAX

      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)

      X12 = X1 - P2(1)
      Y12 = Y1 - P2(2)
      Z12 = Z1 - P2(3)
C     LE CARRE DE LA LONGUEUR DE L'ARETE 12
      A   = X12 ** 2 + Y12 ** 2 + Z12 ** 2

      X13 = X1 - P3(1)
      Y13 = Y1 - P3(2)
      Z13 = Z1 - P3(3)
C     LE CARRE DE LA LONGUEUR DE L'ARETE 13
      B   = X13 ** 2 + Y13 ** 2 + Z13 ** 2

      X14 = X1 - P4(1)
      Y14 = Y1 - P4(2)
      Z14 = Z1 - P4(3)
C     LE CARRE DE LA LONGUEUR DE L'ARETE 14
      C   = X14 ** 2 + Y14 ** 2 + Z14 ** 2

      ARMIN = MIN( A, B, C )
      ARMAX = MAX( A, B, C )

      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13

C     VOLU6T EST 6 FOIS LE VOLUME DU TETRAEDRE
      VOLU6T = X12 * D + X13 * E + X14 * F

C     VOLUME EXACT DU TETRAEDRE
      VOLUTE = VOLU6T / 6D0

C     LA SURFACE DES 4 TRIANGLES FACES DU TETRAEDRE
C     =============================================
      SURFTE = 0D0
      DO 90 I=1,4

         GOTO( 20, 30, 40, 50 ),I

C        LA FACE 1 DU TETRAEDRE
 20      DO J=1,3
            V1(J) = P3(J) - P1(J)
            V2(J) = P2(J) - P1(J)
         ENDDO
         GOTO 60

C        LA FACE 2 DU TETRAEDRE
 30      A = 0D0
         B = 0D0
         DO J=1,3
            V1(J) = P3(J) - P2(J)
            V2(J) = P4(J) - P2(J)
            A = A + V1(J) ** 2
            B = B + V2(J) ** 2
         ENDDO
         GOTO 60

C        LA FACE 3 DU TETRAEDRE
 40      C = 0D0
         DO J=1,3
            V1(J) = P1(J) - P3(J)
            V2(J) = P4(J) - P3(J)
            C = C + V2(J) ** 2
         ENDDO
         GOTO 60

C        LA FACE 4 DU TETRAEDRE
 50      DO J=1,3
            V1(J) = P1(J) - P4(J)
            V2(J) = P2(J) - P4(J)
         ENDDO

C        PRODUIT VECTORIEL V1 * V2 => V3
 60      DO J=1,3
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
         ENDDO

C        2*SURFACE DE LA FACE I DU TETRAEDRE
         SURF = 0D0
         DO J=1,3
            SURF = SURF + V3(J) ** 2
         ENDDO
         SURF = SQRT( SURF )

C        LA SURFACE DU TRIANGLE I DU TETRAEDRE
         SURFTR( I ) = SURF / 2

         SURFTE = SURFTE + SURF

 90   CONTINUE

C     ARETES MINIMALE ET MAXIMALE
C     ===========================
      ARMIN = SQRT( MIN( ARMIN, A, B, C ) )
      ARMAX = SQRT( MAX( ARMAX, A, B, C ) )

C     CRITERE DE QUALITE : ALFA * RAYON INSCRIT/LONGUEUR DE LA PLUS LONGUE ARETE
C     ==========================================================================
      IF( VOLU6T .LT. 0.D0 ) THEN
         QUALTE = -1
      ELSE
         QUALTE = REAL( ALFA * VOLU6T / ( SURFTE * ARMAX ) )

cccC     AUTRE CRITERE : 3 * RAYON INSCRIT / RAYON CIRCONSCRIT
ccc      QUATET  = REAL( 3D0 * VOLU6T / ( SURFTE * SQRT(CENTRE(4)) )

      ENDIF

      RETURN
      END
