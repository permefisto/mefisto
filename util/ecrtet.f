      SUBROUTINE ECRTET( P1 , P2 , P3 , P4 , ECRASE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE L'ECRASEMENT DU TETRAEDRE P1 P2 P3 P4
C ----- C-A-D ECRASE=MIN( HAUTEUR **2 / SURFACE FACE OPPOSEE )
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C
C SORTIE :
C --------
C ECRASE : MIN( HAUTEUR **2 / SURFACE FACE OPPOSEE )
C          -I SI LA FACE TRIANGULAIRE I EST DEGENEREE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS JANVIER 1987
C...............................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/gsmenu.inc"
      REAL              P1(3),P2(3),P3(3),P4(3)
      REAL              ABC(3),A,B,C,A2,B2,C2,P1P2(3)
      EQUIVALENCE       (ABC(1),A),(ABC(2),B),(ABC(3),C),
     %                  (P1P2(1),X21),(P1P2(2),Y21),(P1P2(3),Z21)
      DATA              RMAX / 1E38 /
C
      ECRASE = RMAX
C
C     BOUCLE SUR LES 4 FACES DU TETRAEDRE
      DO 1000 I = 1,4
C        LA HAUTEUR EST ISSUE DU SOMMET I
         GOTO(1,2,3,4) , I
C
C        SOMMET 1 FACE 2 3 4
 1       S   = P2(1)
         X21 = S - P1(1)
         X32 = P3(1) - S
         X42 = P4(1) - S
C
         S   = P2(2)
         Y21 = S - P1(2)
         Y32 = P3(2) - S
         Y42 = P4(2) - S
C
         S   = P2(3)
         Z21 = S - P1(3)
         Z32 = P3(3) - S
         Z42 = P4(3) - S
         GOTO 5
C
C        SOMMET 2 FACE 3 4 1
 2       S   = P3(1)
         X21 = S - P2(1)
         X32 = P4(1) - S
         X42 = P1(1) - S
C
         S   = P3(2)
         Y21 = S - P2(2)
         Y32 = P4(2) - S
         Y42 = P1(2) - S
C
         S   = P3(3)
         Z21 = S - P2(3)
         Z32 = P4(3) - S
         Z42 = P1(3) - S
         GOTO 5
C
C        SOMMET 3 FACE 4 1 2
 3       S   = P4(1)
         X21 = S - P3(1)
         X32 = P1(1) - S
         X42 = P2(1) - S
C
         S   = P4(2)
         Y21 = S - P3(2)
         Y32 = P1(2) - S
         Y42 = P2(2) - S
C
         S   = P4(3)
         Z21 = S - P3(3)
         Z32 = P1(3) - S
         Z42 = P2(3) - S
         GOTO 5
C
C        SOMMET 4 FACE 1 2 3
 4       S   = P1(1)
         X21 = S - P4(1)
         X32 = P2(1) - S
         X42 = P3(1) - S
C
         S   = P1(2)
         Y21 = S - P4(2)
         Y32 = P2(2) - S
         Y42 = P3(2) - S
C
         S   = P1(3)
         Z21 = S - P4(3)
         Z32 = P2(3) - S
         Z42 = P3(3) - S
         GOTO 5
C
C        COMPOSANTES DE LA NORMALE
 5       A = Y32 * Z42 - Y42 * Z32
         B = Z32 * X42 - Z42 * X32
         C = X32 * Y42 - X42 * Y32
C
C        SURFACE*2 DU TRIANGLE S2 S3 S4
         A2 = A * A
         B2 = B * B
         C2 = C * C
         S2 = A2 + B2 + C2
         IF( S2 .LE. 0. ) GOTO 9999
C
C        REPERAGE DU PLUS GRAND DES A , B , C
         IF( A2 .GE. B2 ) THEN
            IF( A2 .GE. C2 ) THEN
               J1 = 1
               J2 = 2
               J3 = 3
            ELSE
               J1 = 3
               J2 = 1
               J3 = 2
            ENDIF
         ELSE
            IF( B2 .GE. C2 ) THEN
               J1 = 2
               J2 = 3
               J3 = 1
            ELSE
               J1 = 3
               J2 = 1
               J3 = 2
            ENDIF
         ENDIF
C
C        ATTENTION : ABC  EN EQUIVALENCE AVEC A , B , C
C                    P1P2 EN EQUIVALENCE AVEC X21 , Y21 , Z21
C
         HHS2S = ABC(J1) * ( ABC(J1) * P1P2(J1) + ABC(J2) * P1P2(J2) +
     %                       ABC(J3) * P1P2(J3) ) / S2
         HHS2S = HHS2S * HHS2S / ( ABC(J1) ** 2 )
         IF( HHS2S .LT. ECRASE ) ECRASE = HHS2S
 1000 CONTINUE
C
C     NORMALISATION EVENTUELLE A FAIRE
      RETURN
C
C     LE TRIANGLE I EST DEGENERE
 9999 ECRASE = -I
      NBLGRC(NRERR) = 1
      KERR(1) = 'FACE TRIANGULAIRE DEGENEREE'
      CALL LEREUR
      WRITE(IMPRIM,19999) I,(P1(J1),P2(J1),P3(J1),P4(J1),J1=1,3)
19999 FORMAT(' ERREUR ECRTET:LA FACE',I2,' TRIANGULAIRE EST DEGENEREE'/
     %' LES 4 SOMMETS SONT'/' X=',4G14.6/' Y=',4G14.6/' Z=',4G14.6/)
      END
