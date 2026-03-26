         SUBROUTINE VOEXH5( NF, HEXYZ, CUXYZ, NBS1, NBS2, NBS1S2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES COORDONNEES DES SOMMETS
C -----    DU MAILLAGE DES FACES DU CUBE UNITE
C ENTREES :
C ---------
C NF      : NUMERO DE LA FACE
C HEXYZ   : COORDONNEES DES SOMMETS DE LA FACE NF DE L'HEXAEDRE
C NBS1    : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C NBS2    : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C NBS1S2  : NOMBRE DE POINTS DANS CHAQUE FACE (MAJORATION)
C
C SORTIES :
C ---------
C CUXYZ   : COORDONNEES DES SOMMETS DE LA FACE NF DU CUBE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         JANVIER 1989
C23456---------------------------------------------------------------012
C
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              XYZ1(3),XYZ2(3),XYZ3(3),XYZ4(3)
      REAL              HEXYZ(3,6,NBS1S2),CUXYZ(3,6,NBS1S2)
C
C-----------------------------------------------------------------------
C     INITIALISATION
C-----------------------------------------------------------------------
C
      IF (NF .EQ. 1) THEN
         J1 = 1
         J2 = 2
         J3 = 3
         XYZ = 0.
      ELSE IF (NF .EQ. 2) THEN
         J1 = 1
         J2 = 3
         J3 = 2
         XYZ = 0.
      ELSE IF (NF .EQ. 3) THEN
         J1 = 2
         J2 = 3
         J3 = 1
         XYZ = 0.
      ELSE IF (NF .EQ. 4) THEN
         J1 = 1
         J2 = 2
         J3 = 3
         XYZ = 1.
      ELSE IF (NF .EQ. 5) THEN
         J1 = 1
         J2 = 3
         J3 = 2
         XYZ = 1.
      ELSE IF (NF .EQ. 6) THEN
         J1 = 2
         J2 = 3
         J3 = 1
         XYZ = 1.
      ELSE
         J1 = 0
         J2 = 0
         J3 = 0
         XYZ=-1.
      ENDIF
C
C-----------------------------------------------------------------------
C        CALCUL DES COORDONNEES DES SOMMETS DU BORD
C-----------------------------------------------------------------------
C
         NBS3=NBS1*(NBS2-1)
C
         RL1 = 0.
         RL3 = 0.
         DO 10 J=1,3
            XYZ1(J)=HEXYZ(J,NF,1)
            XYZ3(J)=HEXYZ(J,NF,NBS3+1)
 10      CONTINUE
         CUXYZ(J1,NF,1) = 0.
         CUXYZ(J2,NF,1) = 0.
         CUXYZ(J1,NF,NBS3+1) = 0.
         CUXYZ(J2,NF,NBS3+1) = 1.
C
         DO 1 I = 2 , NBS1
C
C        ** LE COTE 1
C
            DO 11 J=1,3
               XYZ2(J)=HEXYZ(J,NF,I)
 11         CONTINUE
            XL = XYZ2(1) - XYZ1(1)
            YL = XYZ2(2) - XYZ1(2)
            ZL = XYZ2(3) - XYZ1(3)
            RL1 = RL1 + SQRT( XL*XL + YL*YL + ZL*ZL )
            CUXYZ(J1,NF,I) = RL1
            DO 12 J=1,3
               XYZ1(J)=XYZ2(J)
 12         CONTINUE
C
C        ** LE COTE 3
C
            DO 13 J=1,3
               XYZ4(J)=HEXYZ(J,NF,NBS3+I)
 13         CONTINUE
            XL = XYZ4(1) - XYZ3(1)
            YL = XYZ4(2) - XYZ3(2)
            ZL = XYZ4(3) - XYZ3(3)
            RL3 = RL3 + SQRT( XL*XL + YL*YL + ZL*ZL )
            CUXYZ(J1,NF,NBS3+I) = RL3
            DO 14 J=1,3
               XYZ3(J)=XYZ4(J)
 14         CONTINUE
1        CONTINUE
C
C        CALCUL DE LA LONGUEUR DES COTES 1 ET 3
C
         RR1 = RL1
         RR3 = RL3
C
C        NORMALISATION
C
         DO 2 I = 1 , NBS1
C        ** LE COTE 1
            CUXYZ(J1,NF,I) = CUXYZ(J1,NF,I) / RR1
            CUXYZ(J2,NF,I) = 0.
C        ** LE COTE 3
            CUXYZ(J1,NF,NBS3+I) = CUXYZ(J1,NF,NBS3+I) / RR3
            CUXYZ(J2,NF,NBS3+I) = 1.
2        CONTINUE
C
         RL2 = 0.
         RL4 = 0.
         DO 20 J=1,3
            XYZ1(J)=HEXYZ(J,NF,NBS1)
            XYZ3(J)=HEXYZ(J,NF,1)
 20      CONTINUE
         CUXYZ(J1,NF,NBS1) = 1.
         CUXYZ(J2,NF,NBS1) = 0.
         CUXYZ(J1,NF,1) = 0.
         CUXYZ(J2,NF,1) = 0.
C
         DO 3 I = 2 , NBS2
C
C        ** LE COTE 2
C
            DO 21 J=1,3
               XYZ2(J)=HEXYZ(J,NF,NBS1*I)
 21         CONTINUE
            XL = XYZ2(1) - XYZ1(1)
            YL = XYZ2(2) - XYZ1(2)
            ZL = XYZ2(3) - XYZ1(3)
            RL2 = RL2 + SQRT( XL*XL + YL*YL + ZL*ZL )
            CUXYZ(J2,NF,NBS1*I) = RL2
            DO 22 J=1,3
               XYZ1(J)=XYZ2(J)
 22         CONTINUE
C
C        ** LE COTE 4
C
            DO 23 J=1,3
               XYZ4(J)=HEXYZ(J,NF,NBS1*(I-1)+1)
 23         CONTINUE
            XL = XYZ4(1) - XYZ3(1)
            YL = XYZ4(2) - XYZ3(2)
            ZL = XYZ4(3) - XYZ3(3)
            RL4 = RL4 + SQRT( XL*XL + YL*YL + ZL*ZL )
            CUXYZ(J2,NF,NBS1*(I-1)+1) = RL4
            DO 24 J=1,3
               XYZ3(J)=XYZ4(J)
 24         CONTINUE
3        CONTINUE
C
C        CALCUL DE LA LONGUEUR DES COTES 2 ET 4
C
         RR2 = RL2
         RR4 = RL4
C
C        NORMALISATION
C
         DO 4 I = 1 , NBS2
C        ** LE COTE 2
            CUXYZ(J1,NF,NBS1*I) = 1.
            CUXYZ(J2,NF,NBS1*I) = CUXYZ(J2,NF,NBS1*I) / RR2
C        ** LE COTE 4
            CUXYZ(J1,NF,NBS1*(I-1)+1) = 0.
            CUXYZ(J2,NF,NBS1*(I-1)+1) = CUXYZ(J2,NF,NBS1*(I-1)+1) / RR4
4        CONTINUE
C
C-----------------------------------------------------------------------
C        CALCUL DES COORDONNEES DES SOMMETS INTERNES DE LA FACE
C-----------------------------------------------------------------------
C
         DO 5 J = 2 , NBS2 - 1
C           LES COORDONNEES DES POINTS  M2 ET M4 DES COTES 2 ET 4
            Y2 = CUXYZ(J2,NF,NBS1*J)
            Y4 = CUXYZ(J2,NF,NBS1*(J-1)+1)
            DO 6 I = 2 , NBS1 - 1
C              LES COORDONNEES DES POINTS  M1 ET M3 DES COTES 1 ET 3
               X1 = CUXYZ(J1,NF,I)
               X3 = CUXYZ(J1,NF,NBS3+I)
C              NUMERO DU POINT M A CREER DANS LE DOMAINE
               NUSO = NBS1 * (J-1) + I
C              INTERSECTION DES DROITES M1M3 ET M2M4
               X3X1 = X3 - X1
               Y2Y4 = Y2 - Y4
               IF(ABS(X3X1) .LT. EPSXYZ) THEN
                  XM = X1
                  YM = Y4 + Y2Y4 * XM
               ELSE
                  IF(ABS(Y2Y4) .LT. EPSXYZ) THEN
                    YM = Y2
                    XM = X1 + X3X1 * YM
                  ELSE
                    YM = ( Y4 + Y2Y4 * X1 ) / ( 1. - Y2Y4 * X3X1 )
                    XM = X1 + X3X1 * YM
                  END IF
               END IF
C
C              LES COORDONNEES DU POINT M
C
               CUXYZ(J1,NF,NUSO) = XM
               CUXYZ(J2,NF,NUSO) = YM
 6          CONTINUE
 5       CONTINUE
C
C-----------------------------------------------------------------------
C        RETOUR A LA DIMENSION 3
C-----------------------------------------------------------------------
      DO 7 NS=1,NBS1*NBS2
         CUXYZ(J3,NF,NS) = XYZ
 7    CONTINUE
C
      RETURN
      END
