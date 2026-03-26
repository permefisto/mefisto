         SUBROUTINE  TRCTR1( NBSA, COSO, XYSTR1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER TOUS LES SOMMETS DU TRIANGLE RECTANGLE UNITE
C -----    A PARTIR DES SOMMETS DES 3 COTES COURBES D'UN TRIANGLE
C
C ENTREES:
C --------
C NBSA   : NOMBRE DE SOMMETS SUR CHAQUE ARETE DU TRIANGLE UNITE
C COSO   : COORDONNEES DES SOMMETS DES 3 COTES DU TRIANGLE COURBE
C
C SORTIES:
C --------
C XYSTR1 : COORDONNEES XY DES NBSA*(NBSA+1)/2 SOMMETS DU TRIANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C234567..............................................................012
      REAL    COSO(3,*)
      REAL    XYSTR1(2,*)
      REAL    XYZ(3,0:6)
C
C     LA FONCTION FORMULE DES NUMEROS DES SOMMETS DU TRIANGLE
      NUSOTR(I,J) = ( I * I - I ) / 2 + J
C
C     LE NOMBRE TOTAL DE SOMMETS DU MAILLAGE REGULIER DU TRIANGLE
      NBTST = ( NBSA * NBSA + NBSA ) / 2
C
      RL1 = 0
      RL2 = 0
      RL3 = 0
      NUM1 = 1
      NUM2 = NBTST + 1 - NBSA
      NUM3 = NBTST
      XYSTR1(1,NUM1) = 0
      XYSTR1(2,NUM1) = 0
      XYSTR1(1,NUM2) = 0
      XYSTR1(2,NUM2) = 0
      XYSTR1(1,NUM3) = 0
      XYSTR1(2,NUM3) = 0
C
      DO 20 I = 2 , NBSA
C
C        L'ARETE I-1 DU COTE S1 S2
         NUM10 = NUM1
         NUM1  = NUM1 + I - 1
         RL1 = RL1 + SQRT( ( COSO(1,NUM1) - COSO(1,NUM10) ) ** 2
     %                   + ( COSO(2,NUM1) - COSO(2,NUM10) ) ** 2
     %                   + ( COSO(3,NUM1) - COSO(3,NUM10) ) ** 2 )
         XYSTR1(1,NUM1) = RL1
C
C        L'ARETE I-1 DU COTE S2 S3
         NUM20 = NUM2
         NUM2  = NUM2 + 1
         RL2 = RL2 + SQRT( ( COSO(1,NUM2) - COSO(1,NUM20) ) ** 2
     %                   + ( COSO(2,NUM2) - COSO(2,NUM20) ) ** 2
     %                   + ( COSO(3,NUM2) - COSO(3,NUM20) ) ** 2 )
         XYSTR1(2,NUM2) = RL2
C
C        L'ARETE I-1 DU COTE S3 S1  DANS CET ORDRE!
         NUM30 = NUM3
         NUM3  = NUM3 - NBSA - 2 + I
         RL3 = RL3 + SQRT( ( COSO(1,NUM3) - COSO(1,NUM30) ) ** 2
     %                   + ( COSO(2,NUM3) - COSO(2,NUM30) ) ** 2
     %                   + ( COSO(3,NUM3) - COSO(3,NUM30) ) ** 2 )
         XYSTR1(2,NUM3) = RL3
C
 20   CONTINUE
C
C     NORMALISATION A 1 DES COORDONNES CURVILIGNES
      NUM1 = 1
      NUM2 = NBTST + 1 - NBSA
      NUM3 = NBTST
      DO 30 I = 2, NBSA
C
C        LE COTE 1
         NUM1 = NUM1 + I - 1
         XYSTR1(1,NUM1) = XYSTR1(1,NUM1) / RL1
         XYSTR1(2,NUM1) = 0
C
C        LE COTE 2
         NUM2 = NUM2 + 1
         XYSTR1(2,NUM2) = XYSTR1(2,NUM2) / RL2
         XYSTR1(1,NUM2) = 1 - XYSTR1(2,NUM2)
C
C        LE COTE 3
         NUM3 = NUM3 - NBSA - 2 + I
         XYSTR1(1,NUM3) = 0
         XYSTR1(2,NUM3) = 1 - ( XYSTR1(2,NUM3) / RL3 )
C
 30   CONTINUE
C
C     ---------------------------------------------------------------
C     CALCUL DES 2 COORDONNEES DES SOMMETS INTERNES AU TRIANGLE UNITE
C     A DISTANCE MINIMALE DE 3 DROITES
C     ---------------------------------------------------------------
      DO 70 I = 3 , NBSA - 1
         DO 50 J = 2 , I - 1
C
C          LA DROITE D1 DE SOMMETS S I-J+1,1 S M,M-I+J
C          LE NUMERO DES 2 SOMMETS DANS LA NUMEROTATION DES SOMMETS DU TRIANGLE
C          SOMMET I-J+1,1 DU TRIANGLE
           I1 = NUSOTR(I-J+1,1)
           XYZ(1,1) = XYSTR1(1,I1)
           XYZ(2,1) = 0
           XYZ(3,1) = 0
C          SOMMET NBSA,NBSA-I+J DU TRIANGLE
           I2 = NUSOTR(NBSA,NBSA-I+J)
           XYZ(1,4) = 1 - XYSTR1(2,I2)
           XYZ(2,4) = XYSTR1(2,I2)
           XYZ(3,4) = 0
C
C          LA DROITE D2 DE SOMMETS S NBSA,J S J,J
C          LE NUMERO DES 2 SOMMETS DANS LA NUMEROTATION DES SOMMETS DU TRIANGLE
C          SOMMET M,J DU TRIANGLE
           J1 = NUSOTR(NBSA,J)
           XYZ(1,2) = 1 - XYSTR1(2,J1)
           XYZ(2,2) = XYSTR1(2,J1)
           XYZ(3,2) = 0
C          SOMMET J,J DU TRIANGLE
           J2 = NUSOTR(J,J)
           XYZ(1,5) = 0
           XYZ(2,5) = XYSTR1(2,J2)
           XYZ(3,5) = 0
C
C          LA DROITE D3 DE SOMMETS S I,I S I,1 (I DESIGNE LA DIAGONALE, J Y)
C          LE NUMERO DES 2 SOMMETS DANS LA NUMEROTATION DES SOMMETS DU TRIANGLE
C          SOMMET I,I DU TRIANGLE
           K1 = NUSOTR(I,I)
           XYZ(1,3) = 0
           XYZ(2,3) = XYSTR1(2,K1)
           XYZ(3,3) = 0
C          SOMMET I,1 DU TRIANGLE
           K2 = NUSOTR(I,1)
           XYZ(1,6) = XYSTR1(1,K2)
           XYZ(2,6) = 0
           XYZ(3,6) = 0
C
C          CALCUL DU POINT LE PLUS PROCHE DES DROITES D1, D2 ET D3
C          =======================================================
           CALL MINDIS( 3, XYZ )
C
C          NUMERO DU POINT (I,J) A CREER DANS LE TRIANGLE UNITE
           NUSO = NUSOTR(I,J)
C
C          LES 3 COORDONNEES DU POINT I,J DANS LE TRIANGLE UNITE
           XYSTR1(1,NUSO) = XYZ(1,0)
           XYSTR1(2,NUSO) = XYZ(2,0)
C
 50      CONTINUE
 70   CONTINUE
      END
