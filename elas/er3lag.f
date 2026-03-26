      SUBROUTINE ER3LAG( PENALI,
     &                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     &                   NOOBLA, NUMILI, NUMALI, LTDELI,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     &                   NUTYEL, NBPOLY, NPI,    POLY,
     &                   G1,     G2,     G3,
     &                   IP,     F,      POIDEL, DP,
     &                   ELAS,   AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE RAIDEUR DES ELEMENTS FINIS 3D
C -----    LAGRANGE ISOPARAMETRIQUE SAUF TETR 3P1D
C
C ENTREES:
C --------
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE FIXATION
C          CETTE VALEUR 1/EPSILON DOIT ETRE PLUS GRANDE QUE LES
C          COEFFICIENTS DU VECTEUR ELEMENTAIRE ET GLOBAL
C
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'EF
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE LA SURFACE DE CHAQUE FACE DE L'EF
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES
C          DES OBJETS SURFACES
C
C NOOBVC : NUMERO DU VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES VOLUMES
C NUMAVO : NUMERO MAXIMAL DES VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES VOLUMES DE L'OBJET
C
C NUTYEL : NUMERO DU TYPE DE L'EF
C NBPOLY : NOMBRE DE POLYNOMES DE CET EF
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C
C G1,G2,G3 : TABLEAUX AUXILIAIRES
C IP     : IP(J) = POSITION DU J-EME D.L. COMPOSANTE PAR COMPOSANTE
C          DANS LA NUMEROTATION PAR NOEUDS
C F      : 3 COORDONNEES DES NPI POINTS D INTEGRATION DE L EF
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L EF
C
C SORTIES:
C --------
C ELAS   : TENSEUR SYMETRIQUE DE L ELASTICITE LINEAIRE SYMETRIQUE
C AE     : MATRICE DE RAIDEUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1998
C23456---------------------------------------------------------------012
      include"./incl/donela.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      INTEGER          LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER          LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER          LTDESU(1:MXDOEL,NUMISU:NUMASU)
      INTEGER          LTDEVO(1:MXDOEL,NUMIVO:NUMAVO)
      INTEGER          NOOBPS(1:*),
     %                 NOOBLA(1:*),
     %                 NOOBSF(1:*)
      INTEGER          NONOFK(1:8)
C
      DOUBLE PRECISION PENALI
      DOUBLE PRECISION POLY(NBPOLY,NPI)
      DOUBLE PRECISION F(NPI,3), POIDEL(NPI),
     %                 DP(3,NBPOLY,NPI),
     %                 AE(*)
      DOUBLE PRECISION G1(NBPOLY*(NBPOLY+1)/2),
     %                 G2(NBPOLY,3),
     %                 G3(NBPOLY,NBPOLY)
      DOUBLE PRECISION AS(6), ANS(3,3), COEFDE(3,27)
      DOUBLE PRECISION ELAS(21), AUX
      INTEGER          IP(3*NBPOLY)
C
C     INITIALISATION DE AE
C     --------------------
      L = 3 * NBPOLY
      CALL AZEROD( L*(L+1)/2, AE )
C
C     CONSTRUCTION DES MATRICES   (TDED11  TDED12  TDED13)
C     -------------------------   (        TDED22  TDED23)
C                                 (                TDED33)
C
C     PUIS, CALCUL DE LA MATRICE AE : ( AE11  AE12  AE12 )
C     -----------------------------   (       AE22  AE12 )
C                                     (             AE33 )
      DO 10 L=1,NPI
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE L'ELASTICITE AU POINT
C        D'INTEGRATION L DE L EF
         CALL REELAS( 4, NOOBVC, 3, F(L,1), F(L,2), F(L,3),
     %                LTDEVO(LPYOUN,NOOBVC), ELAS)
C
C        TDED11
         AUX = POIDEL(L)
         AS(1) = AUX * ELAS( 1)
         AS(2) = AUX * ELAS( 7)
         AS(3) = AUX * ELAS(10)
         AS(4) = AUX * ELAS(16)
         AS(5) = AUX * ELAS(19)
         AS(6) = AUX * ELAS(21)
C
C        CALCUL DE LA MATRICE AE11 SYMETRIQUE
C        calcul du produit   g1(nbpoly,nbpoly) = Tdp(nbpoly,3) * a11(3,3) * dp(3
         CALL TABA8D( NBPOLY, 3, DP(1,1,L), AS, G1, G2 )
C        transfert dans la matrice AE11
         CALL PLONAD( 1, NBPOLY, NBPOLY, 1, 1, IP, G1, G1, AE )
C
C        TDED22
         AS(1) = AUX * ELAS(10)
         AS(2) = AUX * ELAS( 8)
         AS(3) = AUX * ELAS( 3)
         AS(4) = AUX * ELAS(14)
         AS(5) = AUX * ELAS(12)
         AS(6) = AUX * ELAS(15)
C
C        CALCUL DE LA MATRICE AE22 SYMETRIQUE
C        calcul du produit   g1(nbpoly,nbpoly) = Tdp(nbpoly,3) * a22(3,3) * dp(3
         CALL TABA8D( NBPOLY, 3, DP(1,1,L), AS, G1, G2 )
C        transfert dans la matrice AE22
         CALL PLONAD( 1, NBPOLY, NBPOLY, 2, 2, IP, G1, G1, AE )
C
C        TDED33
         AS(1) = AUX * ELAS(21)
         AS(2) = AUX * ELAS(20)
         AS(3) = AUX * ELAS(15)
         AS(4) = AUX * ELAS(18)
         AS(5) = AUX * ELAS(13)
         AS(6) = AUX * ELAS( 6)
C
C        CALCUL DE LA MATRICE AE33 SYMETRIQUE
C        calcul du produit   g1(nbpoly,nbpoly) = Tdp(nbpoly,3) * a33(3,3) * dp(3
         CALL TABA8D( NBPOLY, 3, DP(1,1,L), AS, G1, G2 )
C        transfert dans la matrice AE33
         CALL PLONAD( 1, NBPOLY, NBPOLY, 3, 3, IP, G1, G1, AE )
C
C
C        TDED12 NON SYMETRIQUE
         ANS(1,1) = AUX * ELAS( 7)
         ANS(1,2) = AUX * ELAS( 2)
         ANS(1,3) = AUX * ELAS(11)
C
         ANS(2,1) = AUX * ELAS(10)
         ANS(2,2) = AUX * ELAS( 8)
         ANS(2,3) = AUX * ELAS(14)
C
         ANS(3,1) = AUX * ELAS(19)
         ANS(3,2) = AUX * ELAS(17)
         ANS(3,3) = AUX * ELAS(20)
C
C        CALCUL DE LA MATRICE AE12 NON SYMETRIQUE
C        calcul du produit   g2(nbpoly,3) = Tdp(3,nbpoly) * a12(3,3)
         CALL TAB0D( NBPOLY, 3, 3, DP(1,1,L), ANS, G2 )
C        calcul du produit   g3(nbpoly,nbpoly) =  g2(nbpoly,3) * dp(3,nbpoly)
         CALL AB0D( NBPOLY, 3, NBPOLY, G2, DP(1,1,L), G3 )
C
C        ON PLONGE G3 DANS AE12
         CALL PLONAD( -1, NBPOLY, NBPOLY, 1, 2, IP, G3, G3, AE )
C
C        TDED13
         ANS(1,1) = AUX * ELAS(16)
         ANS(1,2) = AUX * ELAS(11)
         ANS(1,3) = AUX * ELAS( 4)
C
         ANS(2,1) = AUX * ELAS(19)
         ANS(2,2) = AUX * ELAS(14)
         ANS(2,3) = AUX * ELAS( 9)
C
         ANS(3,1) = AUX * ELAS(21)
         ANS(3,2) = AUX * ELAS(20)
         ANS(3,3) = AUX * ELAS(18)
C
C        CALCUL DE LA MATRICE AE13 NON SYMETRIQUE
C        calcul du produit   g2(nbpoly,3) = Tdp(3,nbpoly) * a13(3,3)
         CALL TAB0D( NBPOLY, 3, 3, DP(1,1,L), ANS, G2 )
C        calcul du produit   g3(nbpoly,nbpoly) =  g2(nbpoly,3) * dp(3,nbpoly)
         CALL AB0D( NBPOLY, 3, NBPOLY, G2, DP(1,1,L), G3 )
C
C        ON PLONGE G3 DANS AE13
         CALL PLONAD( -1, NBPOLY, NBPOLY, 1, 3, IP, G3, G3, AE )
C
C        TDED23
         ANS(1,1) = AUX * ELAS(19)
         ANS(1,2) = AUX * ELAS(14)
         ANS(1,3) = AUX * ELAS( 9)
C
         ANS(2,1) = AUX * ELAS(17)
         ANS(2,2) = AUX * ELAS(12)
         ANS(2,3) = AUX * ELAS( 5)
C
         ANS(3,1) = AUX * ELAS(20)
         ANS(3,2) = AUX * ELAS(15)
         ANS(3,3) = AUX * ELAS(13)
C
C        CALCUL DE LA MATRICE AE23 NON SYMETRIQUE
C        calcul du produit   g2(nbpoly,3) = Tdp(3,nbpoly) * a23(3,3)
         CALL TAB0D( NBPOLY, 3, 3, DP(1,1,L), ANS, G2 )
C        calcul du produit   g3(nbpoly,nbpoly) =  g2(nbpoly,3) * dp(3,nbpoly)
         CALL AB0D( NBPOLY, 3, NBPOLY, G2, DP(1,1,L), G3 )
C
C        ON PLONGE G3 DANS AE23
         CALL PLONAD( -1, NBPOLY, NBPOLY, 2, 3, IP, G3, G3, AE )
C
 10   CONTINUE
C
C     CONTRIBUTION DU COEFFICIENT DU DEPLACEMENT SUR LA SURFACE
C     SOUS LA FORME DE COEFFICIENTS DIAGONAUX DU TYPE DE LA MASSE
C     ===========================================================
C     LA SURFACE SUPPORTE T ELLE UN COEFFICIENT DU DEPLACEMENT?
C     ADRESSE MCN DU COEFFICIENT DU DEPLACEMENT
      MN = LTDEVO(LPCOED,NOOBVC)
      IF( MN .GT. 0 ) THEN
         CALL EN3LAG( NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     &                NBPOLY, NPI,    POLY,
     &                F,      POIDEL, COEFDE,
     &                AE  )
      ENDIF
C
      IF( PENALI .EQ. 0D0 ) RETURN
C
C     =======================================================
C     CONTRIBUTION DES FACES A LA PENALISATION DE LA FIXATION
C     =======================================================
C
      DO 160 K=1,NFACE
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FIXATION PENALISEE?
            MN = LTDESU( LPFIXA, NOOB )
            IF( MN .LE. 0 ) GOTO 160
C
C           UN TABLEAU FIXATION EXISTE POUR CETTE SURFACE
C           CALCUL DE LA CONTRIBUTION AU SECOND MEMBRE ELEMENTAIRE
C           ------------------------------------------------------
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
C           LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
            NBCOFI = MCN( MN + WBCOFI )
C
            DO 140 I=1,NBNOFK
C
C              LE NUMERO ELEMENTAIRE DU NOEUD I DE LA FACE K
               NI = NONOFK(I)
C
               DO 120 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                  M = MCN( MN + WUCOFI - 1 + J )
C
C                 LE NUMERO DU DL DU SOMMET NI
                  II = NI * 3 - 3 + M
C
C                 CALCUL DE LA CONTRIBUTION A LA FIXATION PENALISEE
                  II = II * ( II + 1 ) / 2
                  AE(II) = AE(II) + PENALI
C
 120           CONTINUE
C
 140        CONTINUE
         ENDIF
 160  CONTINUE
C
C     ========================================================
C     CONTRIBUTION DES ARETES A LA PENALISATION DE LA FIXATION
C     ========================================================
      DO 250 K=1,NARET
C
C        NO DE LIGNE DE L'ARETE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. EST IL SUPPORT D'UN CONTACT PENALISE?
            MN = LTDELI(LPFIXA,NOOB)
            IF( MN .GT. 0 ) THEN
C
C              OUI: UN TABLEAU FIXATION PENALISE EXISTE POUR CE COTE
C              RECHERCHE DES NUMEROS LOCAUX DES NOEUDS DE L'ARETE K
               NONOFK(1) = NOSOAR(1,K)
               NONOFK(2) = NOSOAR(2,K)
               IF( NBNOAR(K) .GT. 0 ) NONOFK(3) = NONOAR(1,K)
C
               DO 220 I=1,2+NBNOAR(K)
C
C                 LE NUMERO DU I-EME NOEUD DE L'ARETE K
                  NI = NONOFK(I)
C
C                 LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
                  NBCOFI = MCN( MN + WBCOFI )
C
                  DO 210 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL DU SOMMET K
                     II = NI * 3 - 3 + M
C
C                    CALCUL DE LA CONTRIBUTION A LA FIXATION PENALISEE
                     II = II * ( II + 1 ) / 2
                     AE(II) = AE(II) + PENALI
C
 210              CONTINUE
C
 220           CONTINUE
C
            ENDIF
         ENDIF
 250  CONTINUE
C
C     =========================================================
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DE LA FIXATION
C     =========================================================
      DO 300 K=1,NBNSOM
C
C        NO DE POINT DU SOMMET K
         NOOB = NOOBPS(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
            MN = LTDEPO( LPFIXA, NOOB )
            IF( MN .GT. 0 ) THEN
C
C              OUI: UN TABLEAU FIXATION PENALISEE EXISTE POUR CE POINT
C              LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
               NBCOFI = MCN( MN + WBCOFI )
C
               DO 270 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                  M = MCN( MN + WUCOFI - 1 + J )
C
C                 LE NUMERO DU DL DU SOMMET K
                  II = K * 3 - 3 + M
C
C                 CALCUL DE LA CONTRIBUTION A LA FIXATION PENALISEE
                  II = II * ( II + 1 ) / 2
                  AE(II) = AE(II) + PENALI
C
 270           CONTINUE
            ENDIF
         ENDIF
 300  CONTINUE
C
      RETURN
      END
