      SUBROUTINE TA1STG( CB,     DIRARE,
     S                   XMM1,   YMM1,
     S                   NBSOCT,
     S                   ABCUC1, ABCUC2, ABCUC3,
     S                   MNXYS1, MNXYS2, MNXYS3,
     S                   MNXYT1, MNXYT2, MNXYT3,
     S                   MNNTG1, MNNTG2, MNNTG3,
     S                   XYZST,  XYZTGA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LES 3 COORDONNEES ET LES 3 COMPOSANTES DES 2
C -----  DERIVEES TANGENTIELLES AU POINT DE COORDONNEES BARYCENTRIQUES
C        (CB(1),CB(2),CB(3)) DU TRIANGLE k PLAN DE LONGUEURS DES COTES
C        EGALES A CELLES DU TRIANGLE K COURBE DE R**3
C        LA FORMULE DE LA TRANSFORMATION TRANSFINIE EST DANS LE
C        CRAS A. PERRONNET 1/1998 SOIT ICI
C
C        s(CB1,CB2,CB3) = CB1 * ( L1(CB2) + L3(1-CB3) - ST1 )
C                       + CB2 * ( L2(CB3) + L1(1-CB1) - ST2 )
C                       + CB3 * ( L3(CB1) + L2(1-CB2) - ST3 )
C
C        AVEC STi   LES 3 COORDONNEES DU SOMMET i DU TRIANGLE COURBE K
C             Li(u) LE PARAMETRAGE [0,1] --> COTE COURBE i DU TRIANGLE K
C
C ENTREES :
C ---------
C CB      : LES 3 COORDONNEES DU POINT
C           (UN SOMMET D'UN SOUS-TRIANGLE DU MAILLAGE DU TRIANGLE k)
C DIRARE  : LES 2 COMPOSANTES DES 2 DIRECTIONS AU POINT
C           (SELON LES 2 COTES DU SOUS-TRIANGLE DU TRIANGLE k AU POINT)
C XMM1    : COMPOSANTES "HOMOGENES" DES VECTEURS=COTES DU TRIANGLE k
C           ( X m - X  m mod 3 + 1 ) / DELTA    m=1,2,3
C YMM1    : ( Y m - Y  m mod 3 + 1 ) / DELTA    m=1,2,3
C
C NBSOCT  : NOMBRE DE POINTS DES 3 LIGNES=COTES DU TRIANGLE COURBE K
C ABCUC123: ABSCISSE CURVILIGNE HOMOGENE DES SOMMETS DES 3 COTES
C           SELON LE SENS   S1=(0,0), S2=(1,0), S3(0,1)
C           COTE 1: S1->S2   COTE 2: S2->S3   COTE 3: S3->S1
C MNXYS123: ADRESSE MCN DE LA PREMIERE COORDONNEE DES SOMMETS   DU COTE 123
C MNXYT123: ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TANGENTES DU COTE 123
C           0 SI PAS DE TANGENTES
C MNNTG123: ADRESSE MCN DES 2 NUMEROS DES TANGENTES DES ARETES  DU COTE 123
C           0 SI PAS DE TANGENTES
C
C SORTIES :
C ---------
C XYZST   : 3 COORDONNEES DU POINT DE K IMAGE DU POINT CB PAR LA
C           TRANSFORMATION TRANSFINIE DU TRIANGLE k SUR LE TRIANGLE K
C XYZTGA  : 3 COMPOSANTES DES 2 DERIVEES DANS LES DIRECTIONS DIRARE DE k
C           C-A-D 3 COMPOSANTES DES 2 TANGENTES AUX COTES AU POINT IMAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1998
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      REAL              CB(3), DIRARE(2,2), XMM1(3), YMM1(3)
      INTEGER           NBSOCT(3)
      REAL              ABCUC1(1:*), ABCUC2(1:*), ABCUC3(1:*)
      REAL              XYZST(3), XYZTGA(3,2)
C
      REAL              F(3,6), DF(3,6), TG2AR(3,2)
C
C     =========================================================
C     LE POINT (XM,YM) EST IL L'UN DES 3 SOMMETS DU TRIANGLE k?
C     =========================================================
      IF( CB(1) .GT. 1.0-EPSXYZ ) THEN
C
C        CB EST LE SOMMET (1,0,0) DU TRIANGLE k
C        SOMMET 1 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 1
         XYZST(1) = RMCN(MNXYS1)
         XYZST(2) = RMCN(MNXYS1+1)
         XYZST(3) = RMCN(MNXYS1+2)
C
C        SOMMET 1 : DS/DU  COTE 1  ARETE P3 HERMITE => C1
         CALL TG2ARL( 1, RMCN(MNXYS1), RMCN(MNXYS1+3),
     S                MNXYT1, MNNTG1, TG2AR )
C        LA LONGUEUR DE L'ARETE 1 DU COTE 1 DU TRIANGLE k
         S1 = ABCUC1(2) - ABCUC1(1)
         XYZTGA(1,1) = TG2AR(1,1) / S1
         XYZTGA(2,1) = TG2AR(2,1) / S1
         XYZTGA(3,1) = TG2AR(3,1) / S1
C
C        SOMMET 1 : DS/DV  COTE 3  DERNIERE ARETE P3 HERMITE => C1
         NBS = NBSOCT(3)
         NBA = NBSOCT(3) - 1
         MN  = MNXYS3 + 3 * NBA
         CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                MNXYT3, MNNTG3, TG2AR )
C        LA LONGUEUR DE LA DERNIERE ARETE DU COTE 3 DU TRIANGLE k
         S1 = ABCUC3(NBS) - ABCUC3(NBA)
         XYZTGA(1,2) = TG2AR(1,2) / S1
         XYZTGA(2,2) = TG2AR(2,2) / S1
         XYZTGA(3,2) = TG2AR(3,2) / S1
         RETURN
      ENDIF
C
      IF( CB(2) .GT. 1.0-EPSXYZ ) THEN
C
C        CB EST LE SOMMET (0,1,0) DU TRIANGLE k
C        SOMMET 2 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 2
         XYZST(1) = RMCN(MNXYS2)
         XYZST(2) = RMCN(MNXYS2+1)
         XYZST(3) = RMCN(MNXYS2+2)
C
C        SOMMET 2 : DS/DU  COTE 1  DERNIERE ARETE P3 HERMITE => C1
         NBS = NBSOCT(1)
         NBA = NBS - 1
         MN  = MNXYS1 + 3 * NBA
         CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                MNXYT1, MNNTG1, TG2AR )
C        LA LONGUEUR DE LA DERNIERE ARETE DU COTE 1 DU TRIANGLE k
         S1 = ABCUC1(NBS) - ABCUC1(NBA)
         XYZTGA(1,1) = -TG2AR(1,2) / S1
         XYZTGA(2,1) = -TG2AR(2,2) / S1
         XYZTGA(3,1) = -TG2AR(3,2) / S1
C
C        SOMMET 2 : DS/DV  COTE 2  PREMIERE ARETE P3 HERMITE => C1
C        Ds(1,0)(S2S3) = ds/du(1,0)(-1) + ds/dv(1,0)(1)  =>
C        ds/dv(1,0)    = Ds(1,0)(S2S3)  + ds/du(1,0)
         CALL TG2ARL( 1, RMCN(MNXYS2), RMCN(MNXYS2+3),
     S                MNXYT2, MNNTG2, TG2AR )
C        LA LONGUEUR DE LA PREMIERE ARETE DU COTE 2 DU TRIANGLE k
         S1 = ABCUC2(2) - ABCUC2(1)
         XYZTGA(1,2) = TG2AR(1,1) / S1 + XYZTGA(1,1)
         XYZTGA(2,2) = TG2AR(2,1) / S1 + XYZTGA(2,1)
         XYZTGA(3,2) = TG2AR(3,1) / S1 + XYZTGA(3,1)
         RETURN
      ENDIF
C
      IF( CB(3) .GT. 1.0-EPSXYZ ) THEN
C
C        CB EST LE SOMMET (0,0,1) DU TRIANGLE k
C        SOMMET 3 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 3
         XYZST(1) = RMCN(MNXYS3  )
         XYZST(2) = RMCN(MNXYS3+1)
         XYZST(3) = RMCN(MNXYS3+2)
C
C        SOMMET 3 : DS/DV  COTE 3  -PREMIERE ARETE P3 HERMITE => C1
         CALL TG2ARL( 1, RMCN(MNXYS3), RMCN(MNXYS3+3),
     S                MNXYT3, MNNTG3, TG2AR )
C        LA LONGUEUR DE DE LA PREMIERE ARETE DU COTE 2 DU TRIANGLE k
         S1 = ABCUC3(2) - ABCUC3(1)
         XYZTGA(1,2) = -TG2AR(1,1) / S1
         XYZTGA(2,2) = -TG2AR(2,1) / S1
         XYZTGA(3,2) = -TG2AR(3,1) / S1
C
C        SOMMET 3: DS/DU  COTE 2  DERNIERE ARETE P3 HERMITE => C1
C        Ds(0,1)(S3S2) = ds/du(1,0)(1) + ds/dv(0,1)(-1)  =>
C        ds/du(0,1)    = Ds(0,1)(S3S2) + ds/dv(0,1)
         NBS = NBSOCT(2)
         NBA = NBSOCT(2) - 1
         MN  = MNXYS2 + 3 * NBA
         CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                MNXYT2, MNNTG2, TG2AR )
C        LA LONGUEUR DE L'ARETE NBA DU COTE 2 DU TRIANGLE k
         S1 = ABCUC2(NBS) - ABCUC2(NBA)
         XYZTGA(1,1) = TG2AR(1,2) / S1 + XYZTGA(1,2)
         XYZTGA(2,1) = TG2AR(2,2) / S1 + XYZTGA(2,2)
         XYZTGA(3,1) = TG2AR(3,2) / S1 + XYZTGA(3,2)
         RETURN
      ENDIF
C
C     ==============================================================
C     ICI LE POINT EST INTERNE OU SUR L'UN DES 3 COTES DU TRIANGLE k
C     ==============================================================
C
C     RECHERCHE DU POINT PROJETE L1(CB2) SUR LE COTE LIGNE 1
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(2), NBSOCT(1), ABCUC1, MNXYS1, MNXYT1, MNNTG1,
     %             F(1,1), DF(1,1) )
C     RECHERCHE DU POINT PROJETE L1(1-CB1) SUR LE COTE LIGNE 1
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.0-CB(1), NBSOCT(1), ABCUC1, MNXYS1, MNXYT1, MNNTG1,
     %             F(1,2), DF(1,2) )
C
C     RECHERCHE DU POINT PROJETE L2(CB3) SUR LE COTE LIGNE 2
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(3), NBSOCT(2), ABCUC2, MNXYS2, MNXYT2, MNNTG2,
     %             F(1,3), DF(1,3) )
C     RECHERCHE DU POINT PROJETE L2(1-CB2) SUR LE COTE LIGNE 2
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.0-CB(2), NBSOCT(2), ABCUC2, MNXYS2, MNXYT2, MNNTG2,
     %             F(1,4), DF(1,4) )
C
C     RECHERCHE DU POINT PROJETE L3(CB1) SUR LE COTE LIGNE 3
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(1), NBSOCT(3), ABCUC3, MNXYS3, MNXYT3, MNNTG3,
     %             F(1,5), DF(1,5) )
C     RECHERCHE DU POINT PROJETE L3(1-CB3) SUR LE COTE LIGNE 3
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.0-CB(3), NBSOCT(3), ABCUC3, MNXYS3, MNXYT3, MNNTG3,
     %             F(1,6), DF(1,6) )
C
C     LE POINT DE K IMAGE DU POINT CB DE k ET SES 2 DERIVEES TANGENTES
C     ----------------------------------------------------------------
C     ADRESSE-1 DES 3 SOMMETS DU TRIANGLE COURBE K
      MN1 = MNXYS1 - 1
      MN2 = MNXYS2 - 1
      MN3 = MNXYS3 - 1
C
      DO 60 K=1,3
C
         A1 = F(K,1) + F(K,6) - RMCN(MN1+K)
         A2 = F(K,3) + F(K,2) - RMCN(MN2+K)
         A3 = F(K,5) + F(K,4) - RMCN(MN3+K)
C
C        (CB) SUR LE TRIANGLE k --> s(CB) SUR LE TRIANGLE COURBE
         XYZST(K) = CB(1) * A1 + CB(2) * A2 + CB(3) * A3
C
         DO 50 L=1,2
C
C           LE POINT (CB) SUR LE TRIANGLE k --> Ds(CB)(DIRARE L)
C                                               SUR LE TRIANGLE COURBE K
            DTGX = DIRARE(1,L)
            DTGY = DIRARE(2,L)
C
            XYZTGA(K,L) = ( YMM1(2) * DTGX - XMM1(2) * DTGY ) * A1
     %      + CB(1) * ( DF(K,1) * ( YMM1(3)*DTGX - XMM1(3)*DTGY )
     %                - DF(K,6) * ( YMM1(1)*DTGX - XMM1(1)*DTGY ) )
     %                  + ( YMM1(3) * DTGX - XMM1(3) * DTGY ) * A2
     %      + CB(2) * ( DF(K,3) * ( YMM1(1)*DTGX - XMM1(1)*DTGY )
     %                - DF(K,2) * ( YMM1(2)*DTGX - XMM1(2)*DTGY ) )
     %                  + ( YMM1(1) * DTGX - XMM1(1) * DTGY ) * A3
     %      + CB(3) * ( DF(K,5) * ( YMM1(2)*DTGX - XMM1(2)*DTGY )
     %                - DF(K,4) * ( YMM1(3)*DTGX - XMM1(3)*DTGY ) )
C
 50       CONTINUE
C
 60   CONTINUE
      END
