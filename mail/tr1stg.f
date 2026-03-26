      SUBROUTINE TR1STG( XM,     YM,     NBSOCT,
     S                   ABCUC1, ABCUC2, ABCUC3,
     S                   MNXYS1, MNXYS2, MNXYS3,
     S                   MNXYT1, MNXYT2, MNXYT3,
     S                   MNNTG1, MNNTG2, MNNTG3,
     S                   XYZST,  XYZDER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES ET LES 3 COMPOSANTES DES 2
C -----     DERIVEES PURES AU POINT (XM,YM) DE L'EF UNITE POUR
C           LA TRANSFORMATION TRANSFINIE CF CRAS A. PERRONNET 1/1998
C
C ENTREES :
C ---------
C XM,YM   : 2 COORDONNEES DU POINT DANS L'EF UNITE
C NBSOCT  : NOMBRE DE POINTS LES 3 LIGNES COTES DU QUADRANGLE COURBE
C ABCUC123: ABSCISSE CURVILIGNE HOMOGENE DES SOMMETS DES 3 COTES
C            SELON LE SENS   S1=(0,0), S2=(1,0), S3(0,1)
C            COTE 1: S1->S2   COTE 2: S2->S3   COTE 3: S3->S1
C MNXYS123: ADRESSE MCN DE LA PREMIERE COORDONNEE DES SOMMETS   DU COTE 123
C MNXYT123: ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TANGENTES DU COTE 123
C MNNTG123: ADRESSE MCN DES 2 NUMEROS DES TANGENTES DES ARETES  DU COTE 123
C
C SORTIES :
C ---------
C XYZST   : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
C XYZDER  : 3 COMPOSANTES DES 2 DERIVEES DS/DU ET DS/DV DU POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           NBSOCT(3)
      REAL              ABCUC1(1:*), ABCUC2(1:*), ABCUC3(1:*)
      REAL              XYZST(3), XYZDER(3,2)
      REAL              CB(3), F(3,6), DF(3,6), TG2AR(3,2)
C
C     LE POINT (XM,YM) EST IL SUR L'UN DES 3 SOMMETS?
C     ===============================================
      IF( XM .LT. EPSXYZ .AND. YM .LT. EPSXYZ ) THEN
C
C        (XM,YM) EST LE SOMMET (0,0) DE L'EF UNITE
C        SOMMET 1 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 1
         XYZST(1) = RMCN(MNXYS1)
         XYZST(2) = RMCN(MNXYS1+1)
         XYZST(3) = RMCN(MNXYS1+2)
C
C        SOMMET 1 : DS/DU  COTE 1  ARETE P3 HERMITE => C1
         CALL TG2ARL( 1, RMCN(MNXYS1), RMCN(MNXYS1+3),
     S                MNXYT1, MNNTG1, TG2AR )
C        LA LONGUEUR DE L'ARETE 1 DU COTE 1 DU TRIANGLE UNITE
         S1 = ABCUC1(2) - ABCUC1(1)
         XYZDER(1,1) = TG2AR(1,1) / S1
         XYZDER(2,1) = TG2AR(2,1) / S1
         XYZDER(3,1) = TG2AR(3,1) / S1
C
C        SOMMET 1 : DS/DV  COTE 3  DERNIERE ARETE P3 HERMITE => C1
         NBS = NBSOCT(3)
         NBA = NBSOCT(3) - 1
         MN  = MNXYS3 + 3 * NBA
         CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                MNXYT3, MNNTG3, TG2AR )
C        LA LONGUEUR DE LA DERNIERE ARETE DU COTE 3 DU TRIANGLE UNITE
         S1 = ABCUC3(NBS) - ABCUC3(NBA)
         XYZDER(1,2) = TG2AR(1,2) / S1
         XYZDER(2,2) = TG2AR(2,2) / S1
         XYZDER(3,2) = TG2AR(3,2) / S1
         RETURN
      ENDIF
C
      IF( XM .GT. 1.0-EPSXYZ .AND. YM .LT. EPSXYZ ) THEN
C
C        (XM,YM) EST LE SOMMET (1,0) DU TRIANGLE UNITE
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
C        LA LONGUEUR DE LA DERNIERE ARETE DU COTE 1 DE L'EF UNITE
         S1 = ABCUC1(NBS) - ABCUC1(NBA)
         XYZDER(1,1) = -TG2AR(1,2) / S1
         XYZDER(2,1) = -TG2AR(2,2) / S1
         XYZDER(3,1) = -TG2AR(3,2) / S1
C
C        SOMMET 2 : DS/DV  COTE 2  PREMIERE ARETE P3 HERMITE => C1
C        Ds(1,0)(S2S3) = ds/du(1,0)(-1) + ds/dv(1,0)(1)  =>
C        ds/dv(1,0)    = Ds(1,0)(S2S3)  + ds/du(1,0)
         CALL TG2ARL( 1, RMCN(MNXYS2), RMCN(MNXYS2+3),
     S                MNXYT2, MNNTG2, TG2AR )
C        LA LONGUEUR DE LA PREMIERE ARETE DU COTE 2 DE L'EF UNITE
         S1 = ABCUC2(2) - ABCUC2(1)
         XYZDER(1,2) = TG2AR(1,1) / S1 + XYZDER(1,1)
         XYZDER(2,2) = TG2AR(2,1) / S1 + XYZDER(2,1)
         XYZDER(3,2) = TG2AR(3,1) / S1 + XYZDER(3,1)
         RETURN
      ENDIF
C
      IF( XM .LT. EPSXYZ . AND. YM .GT. 1.0-EPSXYZ ) THEN
C
C        SOMMET 3 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 3
         XYZST(1) = RMCN(MNXYS3  )
         XYZST(2) = RMCN(MNXYS3+1)
         XYZST(3) = RMCN(MNXYS3+2)
C
C        SOMMET 3 : DS/DV  COTE 3  -PREMIERE ARETE P3 HERMITE => C1
         CALL TG2ARL( 1, RMCN(MNXYS3), RMCN(MNXYS3+3),
     S                MNXYT3, MNNTG3, TG2AR )
C        LA LONGUEUR DE DE LA PREMIERE ARETE DU COTE 2 DE L'EF UNITE
         S1 = ABCUC3(2) - ABCUC3(1)
         XYZDER(1,2) = -TG2AR(1,1) / S1
         XYZDER(2,2) = -TG2AR(2,1) / S1
         XYZDER(3,2) = -TG2AR(3,1) / S1
C
C        SOMMET 3: DS/DU  COTE 2  DERNIERE ARETE P3 HERMITE => C1
C        Ds(0,1)(S3S2) = ds/du(1,0)(1) + ds/dv(0,1)(-1)  =>
C        ds/du(0,1)    = Ds(0,1)(S3S2) + ds/dv(0,1)
         NBS = NBSOCT(2)
         NBA = NBSOCT(2) - 1
         MN  = MNXYS2 + 3 * NBA
         CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                MNXYT2, MNNTG2, TG2AR )
C        LA LONGUEUR DE L'ARETE NBA DU COTE 2 DE L'EF UNITE
         S1 = ABCUC2(NBS) - ABCUC2(NBA)
         XYZDER(1,1) = TG2AR(1,2) / S1 + XYZDER(1,2)
         XYZDER(2,1) = TG2AR(2,2) / S1 + XYZDER(2,2)
         XYZDER(3,1) = TG2AR(3,2) / S1 + XYZDER(3,2)
         RETURN
      ENDIF
C
C     LE POINT EST INTERNE OU SUR L'UN DES 3 COTES DU TRIANGLE UNITE
C     ==============================================================
C     LES 3 COORDONNEES BARYCENTRIQUES DU POINT (XM,YM) DANS LE TRIANGLE UNITE
      CB(1) = 1.0 - XM - YM
      CB(2) = XM
      CB(3) = YM
C
C     LA FORMULE (CRAS 1998 A.P.) DE L'INTERPOLATION TRANSFINIE DU
C     TRIANGLE UNITE AU TRIANGLE COURBE DE SOMMETS ST
C     s(u,v) = CB1 * ( L1(CB2) + L3(1-CB3) - ST(1) )
C            + CB2 * ( L2(CB3) + L1(1-CB1) - ST(2) )
C            + CB3 * ( L3(CB1) + L2(1-CB2) - ST(3) )
C
C     RECHERCHE DU POINT PROJETE L1(CB2) SUR LE COTE LIGNE 1
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(2), NBSOCT(1), ABCUC1, MNXYS1, MNXYT1, MNNTG1,
     %             F(1,1), DF(1,1) )
C     RECHERCHE DU POINT PROJETE L1(1-CB1) SUR LE COTE LIGNE 1
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.-CB(1), NBSOCT(1), ABCUC1, MNXYS1, MNXYT1, MNNTG1,
     %             F(1,2), DF(1,2) )
C
C     RECHERCHE DU POINT PROJETE L2(CB3) SUR LE COTE LIGNE 2
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(3), NBSOCT(2), ABCUC2, MNXYS2, MNXYT2, MNNTG2,
     %             F(1,3), DF(1,3) )
C     RECHERCHE DU POINT PROJETE L2(1-CB2) SUR LE COTE LIGNE 2
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.-CB(2), NBSOCT(2), ABCUC2, MNXYS2, MNXYT2, MNNTG2,
     %             F(1,4), DF(1,4) )
C
C     RECHERCHE DU POINT PROJETE L3(CB1) SUR LE COTE LIGNE 3
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( CB(1), NBSOCT(3), ABCUC3, MNXYS3, MNXYT3, MNNTG3,
     %             F(1,5), DF(1,5) )
C     RECHERCHE DU POINT PROJETE L3(1-CB3) SUR LE COTE LIGNE 3
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( 1.-CB(3), NBSOCT(3), ABCUC3, MNXYS3, MNXYT3, MNNTG3,
     %             F(1,6), DF(1,6) )
C
C     LE POINT (XM,YM) DU MAILLAGE ET SES 2 DERIVEES
C     ----------------------------------------------
C     ADRESSE-1 DES 3 SOMMETS DU TRIANGLE COURBE
      MN1 = MNXYS1 - 1
      MN2 = MNXYS2 - 1
      MN3 = MNXYS3 - 1
C
      DO 50 K=1,3
C
C        (u,v) SUR LE TRIANGLE UNITE --> s(u,v) SUR LE TRIANGLE COURBE
         XYZST(K) = CB(1) * ( F(K,1) + F(K,6) - RMCN(MN1+K) )
     %            + CB(2) * ( F(K,3) + F(K,2) - RMCN(MN2+K) )
     %            + CB(3) * ( F(K,5) + F(K,4) - RMCN(MN3+K) )
C
 50   CONTINUE
C
      DO 60 K=1,3
C
C        (u,v) SUR LE TRIANGLE UNITE --> ds(u,v)/du SUR LE  TRIANGLE COURBE
         XYZDER(K,1) = - F(K,1) - F(K,6) + RMCN(MN1+K)
     %                 + F(K,3) + F(K,2) - RMCN(MN2+K)
     %                 + CB(1) * DF(K,1)
     %                 + CB(2) * DF(K,2)
     %                 - CB(3) * ( DF(K,5) + DF(K,4) )
C
C        (u,v) SUR LE TRIANGLE UNITE --> ds(u,v)/dv SUR LE  TRIANGLE COURBE
         XYZDER(K,2) = - F(K,1) - F(K,6) + RMCN(MN1+K)
     %                 + F(K,5) + F(K,4) - RMCN(MN3+K)
     %                 - CB(1) *   DF(K,6)
     %                 + CB(2) * ( DF(K,3) + DF(K,2) )
     %                 - CB(3) *   DF(K,5)
C
 60   CONTINUE
      END
