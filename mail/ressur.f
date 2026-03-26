      SUBROUTINE RESSUR( NBS1, NBS2, COSO, FCI, FCJ, 
     S                   MAT, A, B, TRAV, OMEGA, TEST, TESTL, NCHT)
C***********************************************************************
C BUT :    RESOLUTION SURRELAXATION
C***********************************************************************
C
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C           FCI    : FONCTIONS DE CONTROLE EN I
C           FCJ    : FONCTIONS DE CONTROLE EN J
C
C TRAVAIL:
C           COSO   : COORDONNEES DES NOEUDS DU MAILLAGE
C           MAT    : MATRICE DU SYSTEME
C           A      : MATRICE DU SYSTEME
C           B      : MATRICE DU SYSTEME
C           TRAV   : MATRICE DU SYSTEME
C           VEC    : MATRICE DU SYSTEME
C
C SORTIES:
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C           TEST  : TEST D'ERREUR
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DIMENSION COSO(3,NBS1*NBS2)
      REAL      MAT(NBS1*NBS2,3)
      DIMENSION A(NBS1*NBS2,-1:1),B(NBS1*NBS2,2)
      DIMENSION FCI(NBS1*NBS2),FCJ(NBS1*NBS2)
      DIMENSION TRAV(NBS1*NBS2,2)
C
      TESTL = 0.00001
      TESTM = 0
      TEST1 = 0
      TEST2 = 0
C***********************************************************************
C             REMPLISSAGE ET CALCUL DES MATRICES MAT ET A
C***********************************************************************
      DO 30 J=2,NBS2-1
        DO 10 I=2,NBS1-1
          NE = (J-1)*NBS1+I
          X1 = COSO(1,NE+1)-COSO(1,NE-1)
          X2 = COSO(1,NE+NBS1)-COSO(1,NE-NBS1)
          Y1 = COSO(2,NE+1)-COSO(2,NE-1)
          Y2 = COSO(2,NE+NBS1)-COSO(2,NE-NBS1)
          CX2 = 0.25*(X2**2+Y2**2)
          CY2 = 0.25*(X1**2+Y1**2)
          CXY = 0.125*(X1*X2+Y1*Y2)
          MAT(NE,1) = CX2
          MAT(NE,2) = CY2
          MAT(NE,3) = CXY
C-----------------------------------------------------------------------
C                REMPLISSAGE DE LA MATRICE TRIDIAGONALE A
C-----------------------------------------------------------------------
          A(NE,-1) =  CX2+0.5*CX2*FCI(NE)
          A(NE,0)  = -2*(CX2+CY2)
          A(NE,1)  =  CX2-0.5*CX2*FCI(NE)
   10   CONTINUE
C-----------------------------------------------------------------------
C              FACTORISATION L x U DE LA MATRICE TRIDIAGONALE A
C-----------------------------------------------------------------------
        DO 20 I=3,NBS1-1
          NE = (J-1)*NBS1+I
          A(NE,-1) = A(NE,-1)/A(NE-1,0)
          A(NE,0)  = A(NE,0)-A(NE,-1)*A(NE-1,1)
   20   CONTINUE
   30 CONTINUE
C***********************************************************************
C               BOUCLE POUR LA RESOLUTION DU SYSTEME LINEAIRE
C                        PENDANT 5 ITERATIONS
C***********************************************************************
      DELTAX1 = 0.
      DELTAY1 = 0.
      DENX1   = 0.
      DENY1   = 0.
      NSORT1  = 0
      NSORT2  = 0
      DO 260 KK=1,5
        IF (NCHT.EQ.1) THEN
          DELTAX1 = 0.
          DELTAY1 = 0.
          DENX1   = 0.
          DENY1   = 0.
        ELSE
          TEST1   = 0.
          TEST2   = 0.
        ENDIF
        IF (NSORT1.EQ.0) THEN
          IF (NSORT2.EQ.0) THEN
C***********************************************************************
C   RESOLUTION EN PARALLELE DES SYSTEMES LINEAIRES A X = B ET A Y = B
C***********************************************************************
            DO 90 J=2,NBS2-1
              DO 40 I=2,NBS1-1
                NE  = (J-1)*NBS1+I
                NEG = NE-NBS1
                NED = NE+NBS1
                CX2 = MAT(NE,1)
                CY2 = MAT(NE,2)
                CXY = MAT(NE,3)
C-----------------------------------------------------------------------
C   REMPLISSAGE DES VECTEURS B POUR LA RESOLUTION DE : A TRAV = B
C-----------------------------------------------------------------------
                B(NE,1) = CXY*COSO(1,NEG-1)  -CY2*COSO(1,NEG)
     S                   -CXY*COSO(1,NEG+1)  -CX2*COSO(1,NE-1)
     S                   +2*(CX2+CY2)*COSO(1,NE)
     S                   -CX2*COSO(1,NE+1)   -CXY*COSO(1,NED-1)
     S                   -CY2*COSO(1,NED)    +CXY*COSO(1,NED+1)
     S                   +CX2*(COSO(1,NE+1)-COSO(1,NE-1))*FCI(NE)/2
     S                   +CY2*(COSO(1,NED) -COSO(1,NEG) )*FCJ(NE)/2
                B(NE,1) = OMEGA*B(NE,1)
                B(NE,2) = CXY*COSO(2,NEG-1)  -CY2*COSO(2,NEG)
     S                   -CXY*COSO(2,NEG+1)  -CX2*COSO(2,NE-1)
     S                   +2*(CX2+CY2)*COSO(2,NE)
     S                   -CX2*COSO(2,NE+1)   -CXY*COSO(2,NED-1)
     S                   -CY2*COSO(2,NED)    +CXY*COSO(2,NED+1)
     S                   +CX2*(COSO(2,NE+1)-COSO(2,NE-1))*FCI(NE)/2
     S                   +CY2*(COSO(2,NED) -COSO(2,NEG) )*FCJ(NE)/2
                B(NE,2) = OMEGA*B(NE,2)
   40         CONTINUE
C-----------------------------------------------------------------------
C         RESOLUTION EN DELTA(X) ET EN DELTA(Y) STOCKES DANS TRAV
C-----------------------------------------------------------------------
              DO 70 K=1,2
                NE = (J-1)*NBS1+2
                TRAV(NE,K) = B(NE,K)
                DO 50 I=3,NBS1-1
                  NE = (J-1)*NBS1+I
                  TRAV(NE,K) = B(NE,K)-A(NE,-1)*TRAV(NE-1,K)
   50           CONTINUE
                NE = J*NBS1-1
                TRAV(NE,K) = TRAV(NE,K)/A(NE,0)
                DO 60 I=NBS1-2,2,-1
                  NE = (J-1)*NBS1+I
                  TRAV(NE,K)=(TRAV(NE,K)-A(NE,1)*TRAV(NE+1,K))/A(NE,0)
   60           CONTINUE
   70         CONTINUE
C-----------------------------------------------------------------------
C              CALCUL DES NOUVELLES COORDONNEES X ET Y
C-----------------------------------------------------------------------
              DO 80 I=2,NBS1-1
                NE = (J-1)*NBS1+I
                IF (NCHT.EQ.1) THEN
                  DELTAX1 = DELTAX1 + ABS(TRAV(NE,1))
                  DENX1 = DENX1 + ABS(COSO(1,NE))
                  DELTAY1 = DELTAY1 + ABS(TRAV(NE,2))
                  DENY1 = DENY1 + ABS(COSO(2,NE))
                ELSE
                  IF (ABS(COSO(1,NE)).GT.TESTL) THEN
                    TEST1 = AMAX1(TEST1,ABS( TRAV(NE,1)/COSO(1,NE)) )
                  ELSE
                    TEST1 = AMAX1(TEST1,ABS(TRAV(NE,1)))
                  ENDIF
                  IF (ABS(COSO(2,NE)).GT.TESTL) THEN
                    TEST2 = AMAX1(TEST2,ABS( TRAV(NE,2)/COSO(2,NE)) )
                  ELSE
                    TEST2 = AMAX1(TEST2,ABS(TRAV(NE,2)))
                  ENDIF
                ENDIF
                COSO(1,NE) = COSO(1,NE)+TRAV(NE,1)
                COSO(2,NE) = COSO(2,NE)+TRAV(NE,2)
   80         CONTINUE
   90       CONTINUE
          ELSE
C***********************************************************************
C            RESOLUTION SEULE DU SYSTEME LINEAIRE A X = B
C***********************************************************************
            DO 140 J=2,NBS2-1
              DO 100 I=2,NBS1-1
                NE  = (J-1)*NBS1+I
                NEG = NE-NBS1
                NED = NE+NBS1
                CX2 = MAT(NE,1)
                CY2 = MAT(NE,2)
                CXY = MAT(NE,3)
C-----------------------------------------------------------------------
C      REMPLISSAGE DU VECTEUR B POUR LA RESOLUTION DE : A TRAV = B
C-----------------------------------------------------------------------
                B(NE,1) = CXY*COSO(1,NEG-1)  -CY2*COSO(1,NEG)
     S                   -CXY*COSO(1,NEG+1)  -CX2*COSO(1,NE-1)
     S                   +2*(CX2+CY2)*COSO(1,NE)
     S                   -CX2*COSO(1,NE+1)   -CXY*COSO(1,NED-1)
     S                   -CY2*COSO(1,NED)    +CXY*COSO(1,NED+1)
     S                   +CX2*(COSO(1,NE+1)-COSO(1,NE-1))*FCI(NE)/2
     S                   +CY2*(COSO(1,NED) -COSO(1,NEG) )*FCJ(NE)/2
                B(NE,1) = OMEGA*B(NE,1)
  100         CONTINUE
C-----------------------------------------------------------------------
C                 RESOLUTION EN DELTA(X) STOCKE DANS TRAV
C-----------------------------------------------------------------------
              NE = (J-1)*NBS1+2
              TRAV(NE,1) = B(NE,1)
              DO 110 I=3,NBS1-1
                NE = (J-1)*NBS1+I
                TRAV(NE,1) = B(NE,1)-A(NE,-1)*TRAV(NE-1,1)
  110         CONTINUE
              NE = J*NBS1-1
              TRAV(NE,1) = TRAV(NE,1)/A(NE,0)
              DO 120 I=NBS1-2,2,-1
                NE = (J-1)*NBS1+I
                TRAV(NE,1)=(TRAV(NE,1)-A(NE,1)*TRAV(NE+1,1))/A(NE,0)
  120         CONTINUE
C-----------------------------------------------------------------------
C                  CALCUL DES NOUVELLES COORDONNEES X
C-----------------------------------------------------------------------
              DO 130 I=2,NBS1-1
                NE = (J-1)*NBS1+I
                IF (NCHT.EQ.1) THEN
                  DELTAX1    = DELTAX1    + ABS( TRAV(NE,1) )
                  DENX1      = DENX1      + ABS( COSO(1,NE) )
                ELSE
                  IF (ABS(COSO(1,NE)).GT.TESTL) THEN
                    TEST1 = AMAX1(TEST1,ABS( TRAV(NE,1)/COSO(1,NE)) )
                  ELSE
                    TEST1 = AMAX1(TEST1,ABS(TRAV(NE,1)))
                  ENDIF
                ENDIF
                COSO(1,NE) = COSO(1,NE) + TRAV(NE,1)
  130         CONTINUE
  140       CONTINUE
          ENDIF
        ELSE
C***********************************************************************
C            RESOLUTION SEULE DU SYSTEME LINEAIRE A Y = B
C***********************************************************************
          DO 190 J=2,NBS2-1
            DO 150 I=2,NBS1-1
              NE  = (J-1)*NBS1+I
              NEG = NE-NBS1
              NED = NE+NBS1
              CX2 = MAT(NE,1)
              CY2 = MAT(NE,2)
              CXY = MAT(NE,3)
C-----------------------------------------------------------------------
C      REMPLISSAGE DU VECTEUR B POUR LA RESOLUTION DE : A TRAV = B
C-----------------------------------------------------------------------
              B(NE,2) = CXY*COSO(2,NEG-1)  -CY2*COSO(2,NEG)
     S                 -CXY*COSO(2,NEG+1)  -CX2*COSO(2,NE-1)
     S                 +2*(CX2+CY2)*COSO(2,NE)
     S                 -CX2*COSO(2,NE+1)   -CXY*COSO(2,NED-1)
     S                 -CY2*COSO(2,NED)    +CXY*COSO(2,NED+1)
     S                 +CX2*(COSO(2,NE+1)-COSO(2,NE-1))*FCI(NE)/2
     S                 +CY2*(COSO(2,NED) -COSO(2,NEG) )*FCJ(NE)/2
              B(NE,2) = OMEGA*B(NE,2)
  150       CONTINUE
C-----------------------------------------------------------------------
C                 RESOLUTION EN DELTA(Y) STOCKE DANS TRAV
C-----------------------------------------------------------------------
            NE = (J-1)*NBS1+2
            TRAV(NE,2) = B(NE,2)
            DO 160 I=3,NBS1-1
              NE = (J-1)*NBS1+I
              TRAV(NE,2) = B(NE,2)-A(NE,-1)*TRAV(NE-1,2)
  160       CONTINUE
            NE = J*NBS1-1
            TRAV(NE,2) = TRAV(NE,2)/A(NE,0)
            DO 170 I=NBS1-2,2,-1
              NE = (J-1)*NBS1+I
              TRAV(NE,2)=(TRAV(NE,2)-A(NE,1)*TRAV(NE+1,2))/A(NE,0)
  170       CONTINUE
C-----------------------------------------------------------------------
C                  CALCUL DES NOUVELLES COORDONNEES Y
C-----------------------------------------------------------------------
            DO 180 I=2,NBS1-1
              NE = (J-1)*NBS1+I
              IF (NCHT.EQ.1) THEN
                DELTAY1    = DELTAY1    + ABS( TRAV(NE,2) )
                DENY1      = DENY1      + ABS( COSO(2,NE) )
              ELSE
                IF (ABS(COSO(2,NE)).GT.TESTL) THEN
                  TEST2 = AMAX1(TEST2,ABS( TRAV(NE,2)/COSO(2,NE)) )
                ELSE
                  TEST2 = AMAX1(TEST2,ABS(TRAV(NE,2)))
                ENDIF
              ENDIF
              COSO(2,NE) = COSO(2,NE) + TRAV(NE,2)
  180       CONTINUE
  190     CONTINUE
        ENDIF
C-----------------------------------------------------------------------
C             TEST DE SORTIE DES ITERATIONS DE RESOLUTION
C-----------------------------------------------------------------------
        IF (NCHT.EQ.1) THEN
          IF (NSORT1.EQ.0) THEN
            TEST1 = DELTAX1/DENX1
          ENDIF
          IF (NSORT2.EQ.0) THEN
            TEST2 = DELTAY1/DENY1
          ENDIF
        ENDIF
CCC        IF (NSORT1.EQ.0) THEN
CCC          IF (NSORT2.EQ.0) THEN
CCC            WRITE (*,1020) KK,TEST1,TEST2
CCC          ELSE
CCC            WRITE (*,1021) KK,TEST1
CCC          ENDIF
CCC        ELSE
CCC          WRITE (*,1022) KK,TEST2
CCC        ENDIF
CCC 1020   FORMAT('   ITERATION : ',I4,' TEST X : ',1PE13.6,
CCC     S                              ' TEST Y : ',1PE13.6)
CCC 1021   FORMAT('   ITERATION : ',I4,' TEST X : ',1PE13.6,
CCC     S                              ' TEST Y : ',13X)
CCC 1022   FORMAT('   ITERATION : ',I4,' TEST X : ',13X,
CCC     S                              ' TEST Y : ',1PE13.6)
        IF (KK.EQ.1) THEN
          TESTM =  AMAX1(TEST1,TEST2)
        ENDIF
        IF ((NSORT1.EQ.0).AND.(TEST1.LT.TESTL)) THEN
          NSORT1 = KK
        ENDIF
        IF ((NSORT2.EQ.0).AND.(TEST2.LT.TESTL)) THEN
          NSORT2 = KK
        ENDIF
        IF ((NSORT1.NE.0).AND.(NSORT2.NE.0)) THEN
          TEST = TESTM
          RETURN
        ENDIF
  260 CONTINUE
      TEST = TESTM
      RETURN
      END
