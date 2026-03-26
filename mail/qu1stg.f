      SUBROUTINE QU1STG( U,      V,      NBSOCT,
     S                   ABCUC1, ABCUC2, ABCUC3, ABCUC4,
     S                   MNXYS1, MNXYS2, MNXYS3, MNXYS4,
     S                   MNXYT1, MNXYT2, MNXYT3, MNXYT4,
     S                   MNNTG1, MNNTG2, MNNTG3, MNNTG4,
     S                   XYZST,  XYZDER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES ET LES 3 COMPOSANTES DES 2
C -----     DERIVEES PURES AU POINT (U,V) DU CARRE UNITE POUR
C           LA TRANSFORMATION TRANSFINIE DE COOK GORDON ET HALL
C
C ENTREES :
C ---------
C U,V     : 2 COORDONNEES DU POINT DANS LE CARRE UNITE
C NBSOCT  : NOMBRE DE POINTS LES 4 LIGNES COTES DU QUADRANGLE COURBE
C ABCUC1234: ABSCISSE CURVILIGNE HOMOGENE DES SOMMETS DES 4 COTES
C            SELON LE SENS   S1=(0,0), S2=(1,0), S3(1,1), S4=(0,1)
C                            COTE 1: S1->S2   COTE 2: S2->S3
C                            COTE 3: S4->S3   COTE 4: S1->S4
C MNXYS1234: ADRESSE MCN DE LA PREMIERE COORDONNEE DES SOMMETS   DU COTE 1234
C MNXYT1234: ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TANGENTES DU COTE 1234
C MNNTG1234: ADRESSE MCN DES 2 NUMEROS DES TANGENTES DES ARETES  DU COTE 1234
C
C SORTIES :
C ---------
C XYZST   : 3 COORDONNEES DU POINT DE COORDONNEES PARAMETRIQUES  (U,V)
C XYZDER  : 3 COMPOSANTES DES 2 DERIVEES DS/DU ET DS/DV AU POINT (U,V)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1996
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           NBSOCT(4)
      REAL              ABCUC1(1:*), ABCUC2(1:*),
     %                  ABCUC3(1:*), ABCUC4(1:*)
      REAL              XYZST(3), XYZDER(3,2)
      REAL              F(3,4), DF(3,4), TG2AR(3,2)
C
C     LE POINT (U,V) EST IL SUR L'UN DES 4 SOMMETS?
C     ===============================================
      IF( U .LT. EPSXYZ ) THEN
C
C        U EST ASSIMILE A ZERO => IL EST SUR LE COTE 4
C
         IF( V .LT. EPSXYZ ) THEN
C           V EST ASSIMILE A ZERO => (U,V) EST LE SOMMET (0,0) DU CARRE
C
C           SOMMET 1 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 1
            XYZST(1) = RMCN(MNXYS1)
            XYZST(2) = RMCN(MNXYS1+1)
            XYZST(3) = RMCN(MNXYS1+2)
C
C           SOMMET 1 : DS/DU  COTE 1  ARETE P3 HERMITE => C1
            CALL TG2ARL( 1, RMCN(MNXYS1), RMCN(MNXYS1+3),
     S                   MNXYT1, MNNTG1, TG2AR )
C           LA LONGUEUR DE L'ARETE 1 DU COTE 1 DU CARRE
            S1 = ABCUC1(2) - ABCUC1(1)
            XYZDER(1,1) = TG2AR(1,1) / S1
            XYZDER(2,1) = TG2AR(2,1) / S1
            XYZDER(3,1) = TG2AR(3,1) / S1
C
C           SOMMET 1 : DS/DV  COTE 4  ARETE P3 HERMITE => C1
            CALL TG2ARL( 1, RMCN(MNXYS4), RMCN(MNXYS4+3),
     S                   MNXYT4, MNNTG4, TG2AR )
C           LA LONGUEUR DE L'ARETE 1 DU COTE 4 DU CARRE
            S1 = ABCUC4(2) - ABCUC4(1)
            XYZDER(1,2) = TG2AR(1,1) / S1
            XYZDER(2,2) = TG2AR(2,1) / S1
            XYZDER(3,2) = TG2AR(3,1) / S1
            RETURN
         ENDIF
C
         IF( V .GT. 1.0-EPSXYZ ) THEN
C           V EST ASSIMILE A 1 => (U,V) EST LE SOMMET (0,1) DU CARRE
C
C           SOMMET 4 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 3
            XYZST(1) = RMCN(MNXYS3)
            XYZST(2) = RMCN(MNXYS3+1)
            XYZST(3) = RMCN(MNXYS3+2)
C
C           SOMMET 4 : DS/DU  COTE 3  ARETE P3 HERMITE => C1
            CALL TG2ARL( 1, RMCN(MNXYS3), RMCN(MNXYS3+3),
     S                   MNXYT3, MNNTG3, TG2AR )
C           LA LONGUEUR DE LA PREMIERE ARETE DU COTE 3 DU CARRE
            S1 = ABCUC3(2) - ABCUC3(1)
            XYZDER(1,1) = TG2AR(1,1) / S1
            XYZDER(2,1) = TG2AR(2,1) / S1
            XYZDER(3,1) = TG2AR(3,1) / S1
C
C           SOMMET 4 : DS/DV  COTE 4  DERNIERE ARETE P3 HERMITE => C1
            NBS = NBSOCT(4)
            NBA = NBSOCT(4) - 1
            MN  = MNXYS4 + 3 * NBA
            CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                   MNXYT4, MNNTG4, TG2AR )
C           LA LONGUEUR DE LA DERNIERE ARETE DU COTE 4 DU CARRE
            S1 = ABCUC4(NBS) - ABCUC4(NBA)
            XYZDER(1,2) = -TG2AR(1,2) / S1
            XYZDER(2,2) = -TG2AR(2,2) / S1
            XYZDER(3,2) = -TG2AR(3,2) / S1
            RETURN
         ENDIF
C
C        LE POINT (U,V) EST SUR LE COTE 4 ENTRE LES 2 EXTREMITES
C        IL EST TRAITE COMME UN POINT INTERNE
      ENDIF
C
      IF( U .GT. 1.0-EPSXYZ ) THEN
C
C        U EST ASSIMILE A 1 => IL EST SUR LE COTE 2
C
         IF( V .LT. EPSXYZ ) THEN
C
C           SOMMET 2 : LES 3 COORDONNEES = CELLES DU PREMIER SOMMET DU COTE 2
            XYZST(1) = RMCN(MNXYS2)
            XYZST(2) = RMCN(MNXYS2+1)
            XYZST(3) = RMCN(MNXYS2+2)
C
C           SOMMET 2 : DS/DU  COTE 1  DERNIERE ARETE P3 HERMITE => C1
            NBS = NBSOCT(1)
            NBA = NBSOCT(1) - 1
            MN  = MNXYS1 + 3 * NBA
            CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                   MNXYT1, MNNTG1, TG2AR )
C           LA LONGUEUR DE LA DERNIERE ARETE NBA DU COTE 1 DU CARRE
            S1 = ABCUC1(NBS) - ABCUC1(NBA)
            XYZDER(1,1) = -TG2AR(1,2) / S1
            XYZDER(2,1) = -TG2AR(2,2) / S1
            XYZDER(3,1) = -TG2AR(3,2) / S1
C
C           SOMMET 2 : DS/DV  COTE 2  PREMIERE ARETE P3 HERMITE => C1
            CALL TG2ARL( 1, RMCN(MNXYS2), RMCN(MNXYS2+3),
     S                   MNXYT2, MNNTG2, TG2AR )
C           LA LONGUEUR DE L'ARETE 1 DU COTE 2 DU CARRE
            S1 = ABCUC4(2) - ABCUC4(1)
            XYZDER(1,2) = TG2AR(1,1) / S1
            XYZDER(2,2) = TG2AR(2,1) / S1
            XYZDER(3,2) = TG2AR(3,1) / S1
            RETURN
         ENDIF
C
         IF( V .GT. 1.0-EPSXYZ ) THEN
C
C           SOMMET 3 : LES 3 COORDONNEES = CELLES DU DERNIER SOMMET DU COTE 3
            NBS = NBSOCT(3)
            NBA = NBSOCT(3) - 1
            MN  = MNXYS3 + 3 * NBA
            XYZST(1) = RMCN(MN  )
            XYZST(2) = RMCN(MN+1)
            XYZST(3) = RMCN(MN+2)
C
C           SOMMET 3: DS/DU  COTE 3  DERNIERE ARETE P3 HERMITE => C1
            CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                   MNXYT3, MNNTG3, TG2AR )
C           LA LONGUEUR DE L'ARETE NBA DU COTE 3 DU CARRE
            S1 = ABCUC3(NBS) - ABCUC3(NBA)
            XYZDER(1,1) = -TG2AR(1,2) / S1
            XYZDER(2,1) = -TG2AR(2,2) / S1
            XYZDER(3,1) = -TG2AR(3,2) / S1
C
C           SOMMET 3 : DS/DV  COTE 2  DERNIERE ARETE P3 HERMITE => C1
            NBS = NBSOCT(2)
            NBA = NBSOCT(2) - 1
            MN  = MNXYS2 + 3 * NBA
            CALL TG2ARL( NBA, RMCN(MN-3), RMCN(MN),
     S                   MNXYT2, MNNTG2, TG2AR )
C           LA LONGUEUR DE L'ARETE DU COTE 2 DU CARRE
            S1 = ABCUC2(NBS) - ABCUC2(NBA)
            XYZDER(1,2) = -TG2AR(1,2) / S1
            XYZDER(2,2) = -TG2AR(2,2) / S1
            XYZDER(3,2) = -TG2AR(3,2) / S1
            RETURN
         ENDIF
C
C        LE POINT (U,V) EST SUR LE COTE 2 ENTRE LES 2 EXTREMITES
C        IL EST TRAITE COMME UN POINT INTERNE
      ENDIF
C
C     LE POINT EST INTERNE OU SUR L'UN DES 4 COTES DU CARRE UNITE
C     ===========================================================
C     RECHERCHE DU POINT PROJETE SUR LA LIGNE 1
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( U, NBSOCT(1), ABCUC1, MNXYS1, MNXYT1, MNNTG1,
     %             F(1,1), DF(1,1) )
C
C     RECHERCHE DU POINT PROJETE SUR LA LIGNE 2
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( V, NBSOCT(2), ABCUC2, MNXYS2, MNXYT2, MNNTG2,
     %             F(1,2), DF(1,2) )
C
C     RECHERCHE DU POINT PROJETE SUR LA LIGNE 3
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( U, NBSOCT(3), ABCUC3, MNXYS3, MNXYT3, MNNTG3,
     %             F(1,3), DF(1,3) )
C
C     RECHERCHE DU POINT PROJETE SUR LA LIGNE 4
C     ET DE LA DERIVEE TANGENTIELLE EN CE POINT
      CALL PTPRCT( V, NBSOCT(4), ABCUC4, MNXYS4, MNXYT4, MNNTG4,
     %             F(1,4), DF(1,4) )
C
C     LE POINT (U,V) DU MAILLAGE ET SES 2 DERIVEES
C     D'APRES LA FORMULE DE COOK / GORDON ET HALL
C     INTERPOLATION TRANSFINIE
C     ----------------------------------------------
      MN1 = MNXYS1 - 1
      MN2 = MNXYS2 - 1
      MN3 = MNXYS2 + 3 * NBSOCT(2) - 4
      MN4 = MNXYS3 - 1
      U1 = 1 - U
      V1 = 1 - V
C
      DO 50 K=1,3
C
C        (u,v) SUR LE CARRE UNITE --> s(u,v) SUR LE QUADRANGLE
         XYZST(K) = V1 * ( F(K,1) - U1 * RMCN(MN1+K) )
     S            + U  * ( F(K,2) - V1 * RMCN(MN2+K) )
     S            + V  * ( F(K,3) - U  * RMCN(MN3+K) )
     S            + U1 * ( F(K,4) - V  * RMCN(MN4+K) )
C
 50   CONTINUE
C
      DO 60 K=1,3
C
C        (u,v) SUR LE CARRE UNITE --> ds(u,v)/du SUR LE QUADRANGLE
         XYZDER(K,1) = F(K,2) - F(K,4)
     S          + V  * (DF(K,3)-RMCN(MN3+K)+RMCN(MN4+K))
     S          + V1 * (DF(K,1)-RMCN(MN2+K)+RMCN(MN1+K))
C
C        (u,v) SUR LE CARRE UNITE --> ds(u,v)/dv SUR LE QUADRANGLE
         XYZDER(K,2) = F(K,3) - F(K,1)
     S          + U  * (DF(K,2)-RMCN(MN3+K)+RMCN(MN2+K))
     S          + U1 * (DF(K,4)-RMCN(MN4+K)+RMCN(MN1+K))
C
 60   CONTINUE
      END
