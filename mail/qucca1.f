         SUBROUTINE  QUCCA1( NX, NY, COSO, XYSCA1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER TOUS LES SOMMETS DU CARRE UNITE
C -----    A PARTIR DES SOMMETS DES 4 COTES COURBES D'UN QUADRANGLE
C
C ENTREES:
C --------
C NX     : NOMBRE DE POINTS SUR CHAQUE COTE 1 3 DU QUADRANGLE
C NY     : NOMBRE DE POINTS SUR CHAQUE COTE 2 4 DU QUADRANGLE
C COSO   : COORDONNEES DES SOMMETS DES 4 COTES DU QUADRANGLE COURBE
C
C SORTIES:
C --------
C XYSCA1 : COORDONNEES XY DES NX*NY SOMMETS DU CARRE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C234567..............................................................012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL     COSO(3,NX,NY)
      REAL     XYSCA1(2,NX,NY)
C
C     LE NOMBRE TOTAL DE SOMMETS DU MAILLAGE REGULIER DU QUADRANGLE
      NBS = NX * NY
C
      RL1 = 0.
      RL2 = 0.
      NUM1 = 1
      NUM2 = NBS - NX + 1
      XYSCA1(1, 1,1) = 0.
      XYSCA1(1,NX,1) = 0.
C
      DO 10 I = 2 , NX
C
C        L'ARETE I-1 DU COTE S1 S2
         I0 = I - 1
         RL1 = RL1 + SQRT( ( COSO(1,I,1) - COSO(1,I0,1) ) ** 2
     %                   + ( COSO(2,I,1) - COSO(2,I0,1) ) ** 2
     %                   + ( COSO(3,I,1) - COSO(3,I0,1) ) ** 2 )
         XYSCA1(1,I,1) = RL1
C
C        L'ARETE I-1 DU COTE S4 S3
         RL2 = RL2 + SQRT( ( COSO(1,I,NY) - COSO(1,I0,NY) ) ** 2
     %                   + ( COSO(2,I,NY) - COSO(2,I0,NY) ) ** 2
     %                   + ( COSO(3,I,NY) - COSO(3,I0,NY) ) ** 2 )
         XYSCA1(1,I,NY) = RL2
C
 10   CONTINUE
C
C     NORMALISATION A 1 DES COORDONNES CURVILIGNES DES COTES 1 3
      DO 20 I = 2, NX
C
C        LE COTE 1
         XYSCA1(1,I,1) = XYSCA1(1,I,1) / RL1
         XYSCA1(2,I,1) = 0
C
C        LE COTE 3
         XYSCA1(1,I,NY) = XYSCA1(1,I,NY) / RL2
         XYSCA1(2,I,NY) = 1
C
 20   CONTINUE
C
C     LES COTES 2 ET 4
      RL1 = 0.
      RL2 = 0.
      XYSCA1(2,NX,1) = 0.
      XYSCA1(2, 1,1) = 0.
C
      DO 30 J = 2 , NY
C
C        L'ARETE I-1 DU COTE S2 S3
         J0 = J - 1
         RL1 = RL1 + SQRT( ( COSO(1,NX,J) - COSO(1,NX,J0) ) ** 2
     %                   + ( COSO(2,NX,J) - COSO(2,NX,J0) ) ** 2
     %                   + ( COSO(3,NX,J) - COSO(3,NX,J0) ) ** 2 )
         XYSCA1(2,NX,J) = RL1
C
C        L'ARETE I-1 DU COTE S1 S4
         RL2 = RL2 + SQRT( ( COSO(1,1,J) - COSO(1,1,J0) ) ** 2
     %                   + ( COSO(2,1,J) - COSO(2,1,J0) ) ** 2
     %                   + ( COSO(3,1,J) - COSO(3,1,J0) ) ** 2 )
         XYSCA1(2,1,J) = RL2
C
 30   CONTINUE
C
C     NORMALISATION A 1 DES COORDONNES CURVILIGNES DES COTES 2 4
      DO 40 J = 2, NY
C
C        LE COTE 2
         XYSCA1(1,NX,J) = 1
         XYSCA1(2,NX,J) = XYSCA1(2,NX,J) / RL1
C
C        LE COTE 4
         XYSCA1(1,1,J) = 0
         XYSCA1(2,1,J) = XYSCA1(2,1,J) / RL2
C
 40   CONTINUE
C
C     ---------------------------------------------------------------
C     CALCUL DES 2 COORDONNEES DES SOMMETS INTERNES AU CARRE UNITE
C     INTERSECTION DES DROITES S I,1 S I,NY ET S 1,J NX,J
C     ---------------------------------------------------------------
      DO 70 J = 2, NY-1
C
C        L'ORDONNEE DU SOMMET M2: S 1,J DU COTE 2
         Y2 = XYSCA1(2,NX,J)
C        L'ORDONNEE DU SOMMET M4: S NX,J DU COTE 4
         Y4 = XYSCA1(2,1,J)
C
         DO 50 I = 2, NX-1
C
C              L'ABSCISSE DU SOMMET M1: S I,1 SUR LE COTE 1
               X1 = XYSCA1(1,I,1)
C              L'ABSCISSE DU SOMMET M3: S I,NY SUR LE COTE 3
               X3 = XYSCA1(1,I,NY)
C
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
                    YM = ( Y4 + Y2Y4 * X1 ) / ( 1 - Y2Y4 * X3X1 )
                    XM = X1 + X3X1 * YM
                  ENDIF
               ENDIF
C
C              LES COORDONNEES DU POINT INTERNE I,J
               XYSCA1(1,I,J) = XM
               XYSCA1(2,I,J) = YM
C
 50      CONTINUE
 70   CONTINUE
      END
