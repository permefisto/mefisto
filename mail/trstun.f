         SUBROUTINE TRSTUN( NB, POXY1, POXY2, POXY3, COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER LE TRIANGLE UNITE A PARTIR D'UN TRIANGLE
C -----    DEFINI PAR SES BORDS MAILLES SELON DES LIGNES
C          LES 3 COTES ONT LE MEME NOMBRE DE SOMMETS
C          LE TRIANGLE UNITE EST UN TRIANGLE STRUCTURE
C
C ENTREES :
C ---------
C NB     : NOMBRE DE POINTS SUR CHAQUE COTE
C COSO   : COORDONNEES DES SOMMETS DES 3 COTES DU TRIANGLE COURBE
C          POUR LA NUMEROTATION DES SOMMETS D'UN TRIANGLE STRUCTURE
C
C ENTREES ET SORTIES :
C --------------------
C POXY   : COORDONNEES DES SOMMETS DU BORD DU TRIANGLE UNITE
C COSO   : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C   ATTENTION: LA 3-EME COORDONNEE DE COSO EST EN FAIT QUELCONQUE!
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: P. JOLY A. PERRONNET  ANALYSE NUMERIQUE PARIS    OCTOBRE 1996
C234567..............................................................012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              POXY1(NB), POXY2(NB), POXY3(NB)
      REAL              COSO(3,*)
C
      REAL              XYZ(3,0:6)
      REAL              XYZS1(3), XYZS2(3), XYZS3(3)
C
C     NUMERO DU SOMMET (I,J)  SELON L'ORDRE ANTI-DIAGONALE
      NOSOMT(I,J) = I * ( I-1) / 2 + J
C
C-----------------------------------------------------------------------
C     CALCUL DES COORDONNEES DES POINTS SUR LE BORD DU TRIANGLE UNITE
C-----------------------------------------------------------------------
      NBS=(NB*NB+NB)/2
C
C     RLi EST LA LONGUEUR DU COTE COURBE i DU TRIANGLE COURBE
      RL1 = 0.
      RL2 = 0.
      RL3 = 0.
C
C     LA COORDONNEE CURVILIGNE SUR [0,1] DE CHACUN DES 3 COTES CORBES
      POXY1(1) = 0.
      POXY2(1) = 0.
      POXY3(1) = 0.
C
      NUM1 = 1
      NUM2 = NBS + 1 - NB
      NUM3 = NBS
      DO 1 J=1,3
         XYZ(J,1)=COSO(J,NUM1)
         XYZ(J,3)=COSO(J,NUM2)
         XYZ(J,5)=COSO(J,NUM3)
 1    CONTINUE
C
      DO 2 I = 2 , NB
         NUM1 = NUM1 + I - 1
         NUM2 = NUM2 + 1
         NUM3 = NUM3 - NB - 2 + I
         DO 3 J=1,3
            XYZ(J,2)=COSO(J,NUM1)
            XYZ(J,4)=COSO(J,NUM2)
            XYZ(J,6)=COSO(J,NUM3)
 3       CONTINUE
C
         XL = XYZ(1,2) - XYZ(1,1)
         YL = XYZ(2,2) - XYZ(2,1)
         ZL = XYZ(3,2) - XYZ(3,1)
         RL1 = RL1 + SQRT( XL*XL + YL*YL + ZL*ZL )
         POXY1(I) = RL1
C
         XL = XYZ(1,4) - XYZ(1,3)
         YL = XYZ(2,4) - XYZ(2,3)
         ZL = XYZ(3,4) - XYZ(3,3)
         RL2 = RL2 + SQRT( XL*XL + YL*YL + ZL*ZL )
         POXY2(I) = RL2
C
         XL = XYZ(1,6) - XYZ(1,5)
         YL = XYZ(2,6) - XYZ(2,5)
         ZL = XYZ(3,6) - XYZ(3,5)
         RL3 = RL3 + SQRT( XL*XL + YL*YL + ZL*ZL )
         POXY3(I) = RL3
C
         DO 4 J=1,3
            XYZ(J,1)=XYZ(J,2)
            XYZ(J,3)=XYZ(J,4)
            XYZ(J,5)=XYZ(J,6)
 4       CONTINUE
 2    CONTINUE
C
C     NORMALISATION A 1 DE LA LONGUEUR => PARAMETRE DANS [0,1]
      DO 5 I = 1 , NB
C        ** LE COTE 1
         POXY1(I) = POXY1(I) / RL1
C        L'ORDONNEE VAUT 0.
C        ** LE COTE 2
         POXY2(I) = POXY2(I) / RL2
C        L'ABSCISSE VAUT 1. - Y
C        ** LE COTE 3
         POXY3(I) = POXY3(I) / RL3
C        L'ABSCISSE VAUT 1.
5     CONTINUE
C
C-----------------------------------------------------------------------
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES AU TRIANGLE UNITE
C-----------------------------------------------------------------------
      DO 6 J=1,3
         XYZS1(J)=COSO(J,1)
         XYZS2(J)=COSO(J,NBS+1-NB)
         XYZS3(J)=COSO(J,NBS)
 6    CONTINUE
C
C     BOUCLE SUR LES POINTS INTERNES DU DOMAINE
      DO 7 I = 3 , NB - 1
         DO 8 J = 2 , I - 1
C
C           LA DROITE D1
            I1 = I + 1 - J
            I2 = NB - I + J
C           LA DROITE D2
            J1 = NB + 1 - J
            J2 = J
C           LA DROITE D3
            K1 = NB + 1 - I
            K2 = I
C
            XYZ(1,1) = POXY1(I1)
            XYZ(2,1) = 0.
            XYZ(3,1) = 0.
C
            XYZ(1,4) = 1. - POXY2(I2)
            XYZ(2,4) = POXY2(I2)
            XYZ(3,4) = 0.
C
            XYZ(1,2) = 0.
            XYZ(2,2) = 1. - POXY3(J1)
            XYZ(3,2) = 0.
C
            XYZ(1,5) = 1. - POXY2(J2)
            XYZ(2,5) = POXY2(J2)
            XYZ(3,5) = 0.
C
            XYZ(1,3) = 0.
            XYZ(2,3) = 1. - POXY3(K1)
            XYZ(3,3) = 0.
C
            XYZ(1,6) = POXY1(K2)
            XYZ(2,6) = 0.
            XYZ(3,6) = 0.
C
C           CALCUL DU POINT LE PLUS PROCHE DES DROITES D1, D2 ET D3
C           =======================================================
            CALL MINDIS( 3, XYZ )
C
C           NUMERO DU SOMMET INTERNE A CREER DANS LE DOMAINE
            NUSO = NOSOMT(I,J)
            COSO(1,NUSO) = XYZ(1,0)
            COSO(2,NUSO) = XYZ(2,0)
C
 8       CONTINUE
 7    CONTINUE
C
C     SENS C1:S1S2 POUR POXY1   C2:S2S3 POUR POXY2  C3:S3S1 POUR POXY3
C
C             S3
C             |  \
C             |     \
C       C3   \/       /\ C2
C             |          \
C             |            \
C             S1---->-------S2
C                C1
C
C     CALCUL DES 2 COORDONNEES DES SOMMETS DES COTES DU TRIANGLE UNITE
C     ----------------------------------------------------------------
C     LE SOMMET (0,0)
      COSO(1,1) = 0
      COSO(2,1) = 0
C
C     LE COTE 1
C     POXY1 L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 1
      DO 16 I=2,NB-1
         NUSO = NOSOMT(I,1)
         COSO(1,NUSO) = POXY1(I)
         COSO(2,NUSO) = 0
 16   CONTINUE
C
C     LE SOMMET (1,0)
      NUSO = NOSOMT(NB,1)
      COSO(1,NUSO) = 1
      COSO(2,NUSO) = 0
C
C     LE COTE 2
C     POXY2 L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 2
      DO 17 I=2,NB-1
         NUSO = NOSOMT(NB,I)
         COSO(1,NUSO) = 1 - POXY2(I)
         COSO(2,NUSO) = POXY2(I)
 17   CONTINUE
C
C     LE SOMMET (0,1)
      NUSO = NOSOMT(NB,NB)
      COSO(1,NUSO) = 0
      COSO(2,NUSO) = 1
C
C     LE COTE 3
C     POXY3 L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 3
      DO 18 I=2,NB-1
         NUSO = NOSOMT(NB+1-I,NB+1-I)
         COSO(1,NUSO) = 0
         COSO(2,NUSO) = 1 - POXY3(I)
 18   CONTINUE
C
      RETURN
      END
