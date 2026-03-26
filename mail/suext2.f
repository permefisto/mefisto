         SUBROUTINE  SUEXT2( NB, POXY1, POXY2, POXY3, COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN TRIANGLE DEFINI PAR SES BORDS
C -----    PAR INTERPOLATION TRANSFINIE DE DEGRE 1 (CF CRAS A. PERRONNET)
C
C ENTREES:
C --------
C NB     : NOMBRE DE POINTS SUR CHAQUE COTE
C COSO   : COORDONNEES DES SOMMETS DES 3 COTES DU TRIANGLE COURBE
C
C SORTIES:
C --------
C POXYi  : COORDONNEES DES SOMMETS DU COTE i DU TRIANGLE UNITE
C COSO   : COORDONNEES DE TOUS LES SOMMETS DU TRIANGLE COURBE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: P. JOLY & A. PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1997
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL              XYZ(3,0:6), XYZ1(3),   XYZ2(3), F(3,6)
      REAL              POXY1(NB),  POXY2(NB), POXY3(NB)
      REAL              COSO(3,*)
      EQUIVALENCE      (XYZ(1,0),CBA2), (XYZ(2,0),CBA3), (XYZ(3,0),CBA1)
C
C-----------------------------------------------------------------------
C     CALCUL DES COORDONNEES DES POINTS SUR LE BORD DU TRIANGLE UNITE
C-----------------------------------------------------------------------
      NBS=(NB*NB+NB)/2
C
      RL1 = 0.
      RL2 = 0.
      RL3 = 0.
      POXY1(1) = 0.
      POXY2(1) = 0.
      POXY3(1) = 0.
      NUM1 = 1
      NUM2 = NBS + 1 - NB
      NUM3 = NBS
      DO 1 J=1,3
         XYZ(J,1)=COSO(J,NUM1)
         XYZ(J,3)=COSO(J,NUM2)
         XYZ(J,5)=COSO(J,NUM3)
 1    CONTINUE
C
      DO 7 I = 2 , NB
C
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
C
 7    CONTINUE
C
C     NORMALISATION A 1
      DO 8 I = 1 , NB
C
C        LE COTE 1
         POXY1(I) = POXY1(I) / RL1
C        L'ORDONNEE VAUT 0.
C
C        LE COTE 2
         POXY2(I) = POXY2(I) / RL2
C        L'ABSCISSE VAUT 1. - Y
C
C        LE COTE 3
         POXY3(I) = POXY3(I) / RL3
C        L'ABSCISSE VAUT 0.
C
 8    CONTINUE
C
C-----------------------------------------------------------------------
C     CALCUL DES 2 COORDONNEES DES SOMMETS INTERNES AU TRIANGLE UNITE
C-----------------------------------------------------------------------
C     BOUCLE SUR LES POINTS INTERNES DU DOMAINE
      DO 70 I = 3 , NB - 1
         DO 60 J = 2 , I - 1
C
C          LA DROITE D1
           I1 = I + 1 - J
           I2 = NB - I + J
C          LA DROITE D2
           J1 = NB + 1 - J
           J2 = J
C          LA DROITE D3
           K1 = NB + 1 - I
           K2 = I
C
C          NUMERO DU POINT M A CREER DANS LE DOMAINE
           NUSO = (I*I-I)/2 + J
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
C          CALCUL DU POINT LE PLUS PROCHE DES DROITES D1, D2 ET D3
C          =======================================================
           CALL MINDIS( 3, XYZ )
C          ABSCISSE ET ORDONNEE DU POINT A DISTANCE MINIMALE
C          CF L'EQUIVALENCE!
C          CBA2 = XYZ(1,0)
C          CBA3 = XYZ(2,0)
           CBA1 = 1. - CBA2 - CBA3
C
C          INTERPOLATION SUR CHACUNE DES TROIS COURBES
C          ===========================================
C          LE COTE 1
C          ---------
C          LE POINT S(1-L2,L2,0)
           IM1=1
           DO 10 I1=2,NB
               IF( CBA2 .LE. POXY1(I1)+EPSXYZ ) THEN
                  IM1=I1
                  GO TO 11
               END IF
 10        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS CBA2 SUR LE COTE 1 UNITE
 11        NUM1=IM1*(IM1-1)/2+1
           NUM2=NUM1-IM1+1
           DO 12 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 12        CONTINUE
           RAP = ( CBA2       - POXY1(IM1-1) )
     %         / ( POXY1(IM1) - POXY1(IM1-1) )
C          LE POINT SUR LE COTE 1 COURBE
           DO 13 K=1,3
              F(K,1) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 13        CONTINUE
C
C          LE POINT S(L1,1-L1,0)
           IM1=1
           DO 16 I1=2,NB
               IF( 1.-CBA1 .LE. POXY1(I1)+EPSXYZ ) THEN
                  IM1=I1
                  GO TO 17
               END IF
 16        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS CBA2 SUR LE COTE 1 UNITE
 17        NUM1=IM1*(IM1-1)/2+1
           NUM2=NUM1-IM1+1
           DO 18 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 18        CONTINUE
           RAP = ( 1.-CBA1    - POXY1(IM1-1) )
     %         / ( POXY1(IM1) - POXY1(IM1-1) )
C          LE POINT SUR LE COTE 1 COURBE
           DO 19 K=1,3
              F(K,2) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 19        CONTINUE
C
C          LE COTE 2
C          ---------
C          LE POINT S(0,1-L3,L3)
           IM2=1
           DO 20 I2=2,NB
              IF( CBA3 .LE. POXY2(I2)+EPSXYZ ) THEN
                 IM2=I2
                 GOTO 21
              END IF
 20        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS (1-CBA3,CBA3)
C          SUR LE COTE 2 DU TRIANGLE UNITE
 21        NUM1=NBS-NB+IM2
           NUM2=NUM1-1
           DO 22 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 22        CONTINUE
           RAP = ( CBA3       - POXY2(IM2-1) )
     S         / ( POXY2(IM2) - POXY2(IM2-1) )
C          LE POINT SUR LE COTE 2 COURBE
           DO 23 K=1,3
              F(K,3) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 23        CONTINUE
C
C          LE POINT S(0,L2,1-L2)
           IM2=1
           DO 26 I2=2,NB
              IF( 1.-CBA2 .LE. POXY2(I2)+EPSXYZ ) THEN
                 IM2=I2
                 GOTO 27
              ENDIF
 26        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS (1-CBA3,CBA3)
C          SUR LE COTE 2 DU TRIANGLE UNITE
 27        NUM1=NBS-NB+IM2
           NUM2=NUM1-1
           DO 28 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 28        CONTINUE
           RAP = ( 1.-CBA2    - POXY2(IM2-1) )
     S         / ( POXY2(IM2) - POXY2(IM2-1) )
C          LE POINT SUR LE COTE 2 COURBE
           DO 29 K=1,3
              F(K,4) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 29        CONTINUE
C
C          LE COTE 3
C          LE POINT S(L1,0,1-L1)
           IM3=1
           DO 30 I3=2,NB
              IF( CBA1 .LE. POXY3(I3)+EPSXYZ ) THEN
                 IM3=I3
                 GO TO 31
              END IF
 30        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS SUR LE COTE 3 UNITE
 31        I3=NB+1-IM3
           NUM1=I3*(I3+1)/2
           NUM2=NUM1+I3+1
           DO 32 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 32        CONTINUE
           RAP = ( CBA1       - POXY3(IM3-1) )
     %         / ( POXY3(IM3) - POXY3(IM3-1) )
C          LE POINT SUR LE COTE 3 COURBE
           DO 33 K=1,3
              F(K,5) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 33        CONTINUE
C
C          LE POINT S(1-L3,0,L3)
           IM3=1
           DO 36 I3=2,NB
              IF( 1.-CBA3 .LE. POXY3(I3)+EPSXYZ ) THEN
                 IM3=I3
                 GO TO 37
              END IF
 36        CONTINUE
C          LES NUMEROS DES SOMMETS ENCADRANTS SUR LE COTE 3 UNITE
 37        I3=NB+1-IM3
           NUM1=I3*(I3+1)/2
           NUM2=NUM1+I3+1
           DO 38 J1=1,3
              XYZ1(J1) = COSO(J1,NUM1)
              XYZ2(J1) = COSO(J1,NUM2)
 38        CONTINUE
           RAP = ( 1.-CBA3    - POXY3(IM3-1) )
     %         / ( POXY3(IM3) - POXY3(IM3-1) )
C          LE POINT SUR LE COTE 3 COURBE
           DO 39 K=1,3
              F(K,6) = XYZ1(K) * RAP + XYZ2(K) * (1. - RAP)
 39        CONTINUE
C
C          LA FORMULE DE TRANSFORMATION: TRIANGLE UNITE -> TRIANGLE COURBE
C          ---------------------------------------------------------------
           DO 50 K=1,3
              COSO(K,NUSO)= CBA1 * ( F(K,6) + F(K,1) - COSO(K,1) )
     %                    + CBA2 * ( F(K,3) + F(K,2) - COSO(K,NBS+1-NB))
     %                    + CBA3 * ( F(K,4) + F(K,5) - COSO(K,NBS) )
 50        CONTINUE
C
 60      CONTINUE
 70   CONTINUE
C
      RETURN
      END
