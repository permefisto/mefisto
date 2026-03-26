      SUBROUTINE ER3P1D( X,      PENALI,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                   ELAS,   DELTA,  DP,     IP,
     %                   AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR D'UN TETRAEDRE 3P1D
C -----
C
C ENTREES:
C --------
C X      : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
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
C NOOBVC : NUMERO DE L'OBJET VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES OBJETS VOLUMES
C
C ELAS   : TENSEUR SYMETRIQUE DE L'ELASTICITE
C DELTA  : JACOBIEN DE LA TRANSFORMATION EF REFERENCE -> EF
C DP     : GRADIENT DES POLYNOMES DE BASE AUX SOMMETS DU TETRAEDRE
C IP     : IP(J) = POSITION DI J-EME D.L. COMPOSANTE PAR COMPOSANTE
C          DANS LA RENUMEROTATION PAR NOEUDS
C
C SORTIES:
C --------
C AE     : MATRICE ELEMENTAIRE SYMETRIQUE DE L'ELASTICITE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : KYNDT Jerome - WALTER David A.PERRONNET UPMC PARIS  JUIN 1996
C MODIF  : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      include "./incl/donela.inc"
      include "./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOEL,NUMIVO:NUMAVO)
      INTEGER           NOOBPS(1:4),
     %                  NOOBLA(1:6),
     %                  NOOBSF(1:4)
      INTEGER           NONOFK(3)
      INTEGER           IP(12)
      REAL              X(4,3)
      DOUBLE PRECISION  PENALI
      DOUBLE PRECISION  DELTA, DP(3,4), ELAS(21), AE(78)
C
      DOUBLE PRECISION  XD, YD, ZD, DELPOI
      DOUBLE PRECISION  A11( 6 ), A12(3,3), A13(3,3),
     %                            A22( 6 ), A23(3,3),
     %                                      A33( 6 ),
     %                  G1(10),
     %                  G2(4,3),
     %                  G3(4,4)
C
C     INITIALISATION DES MATRICES
C     ===========================
      CALL AZEROD( 78 , AE )
C
      A11(2)=0D0
      A11(4)=0D0
      A11(5)=0D0
C
      A22(2)=0D0
      A22(4)=0D0
      A22(5)=0D0
C
      A33(2)=0D0
      A33(4)=0D0
      A33(5)=0D0
C
      CALL AZEROD( 9 , A12 )
      CALL AZEROD( 9 , A13 )
      CALL AZEROD( 9 , A23 )
C
C     CONSTRUCTION DES MATRICES   (A11   A12  A13)
C     =========================   (      A22  A23)
C                                 (           A33)
      DO 5 K=1,4
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE L'ELASTICITE AU POINT
C        D'INTEGRATION K
         XD = X(K,1)
         YD = X(K,2)
         ZD = X(K,3)
         CALL REELAS( 4, NOOBVC, 3, XD, YD, ZD,
     %                LTDEVO(LPYOUN,NOOBVC), ELAS )
C
         A11(1)=ELAS(1)
         A11(3)=ELAS(10)
         A11(6)=ELAS(21)
C
         A22(1)=ELAS(10)
         A22(3)=ELAS(3)
         A22(6)=ELAS(15)
C
         A33(1)=ELAS(21)
         A33(3)=ELAS(15)
         A33(6)=ELAS(6)
C
         A12(1,2)=ELAS(2)
         A12(2,1)=ELAS(10)
C
         A13(1,3)=ELAS(4)
         A13(3,1)=ELAS(21)
C
         A23(2,3)=ELAS(5)
         A23(3,2)=ELAS(15)
C
C        CALCUL DE LA MATRICE AE11
C        -------------------------
C        [G1]=T[DP]*[A11]*[DP]
         CALL TABA8D(4, 3, DP, A11, G1, G2)
C        TRANSFERT DANS LA MATRICE AE11 PAR BLOCS 4x4
         CALL PLONAD(1, 4, 4, 1, 1, IP, G1, G1, AE)
C
C        CALCUL DE LA MATRICE AE22
C        -------------------------
C        [G1]=T[DP]*[A22]*[DP]
         CALL TABA8D(4, 3, DP, A22, G1, G2)
C        TRANSFERT DANS LA MATRICE AE22 PAR BLOCS 4x4
         CALL PLONAD(1, 4, 4, 2, 2, IP, G1, G1, AE)
C
C        CALCUL DE LA MATRICE AE33
C        -------------------------
C        [G1]=T[DP]*[A33]*[DP]
         CALL TABA8D(4, 3, DP, A33, G1, G2)
C        TRANSFERT DANS LA MATRICE AE33 PAR BLOCS 4x4
         CALL PLONAD(1, 4, 4, 3, 3, IP, G1, G1, AE)
C
C        CALCUL DE LA MATRICE AE12
C        -------------------------
C        [G2]=T[DP]*[A12]
         CALL TAB0D(4,3,3,DP,A12,G2)
C        [G3]=[G2]*[DP]
         CALL AB0D(4,3,4,G2,DP,G3) 
C        ON PLONGE G3 DANS AE12 PAR BLOCS 4x4
         CALL PLONAD(-1, 4, 4, 1, 2, IP, G3, G3, AE)
C
C        CALCUL DE LA MATRICE AE13
C        -------------------------
C        [G2]=T[DP]*[A13]
         CALL TAB0D(4,3,3,DP,A13,G2)
C        [G3]=[G2]*[DP]
         CALL AB0D(4,3,4,G2,DP,G3)
C        ON PLONGE G3 DANS AE13 PAR BLOCS 4x4
         CALL PLONAD(-1, 4, 4, 1, 3, IP, G3, G3, AE)
C
C        CALCUL DE LA MATRICE AE23
C        -------------------------
C        [G2]=T[DP]*[A23]
         CALL TAB0D(4,3,3,DP,A23,G2)
C        [G3]=[G2]*[DP]
         CALL AB0D(4,3,4,G2,DP,G3)
C        ON PLONGE G3 DANS AE23 PAR BLOCS 4x4
         CALL PLONAD(-1, 4, 4, 2, 3, IP, G3, G3, AE)
C
    5 CONTINUE
C
C     AE = AE * DELTA * POIDS
      DELPOI = DELTA / 24D0
      DO 30 I=1,78
         AE( I ) = AE( I ) * DELPOI
   30 CONTINUE
C
C     CONTRIBUTION DU COEFFICIENT DU DEPLACEMENT SUR LE VOLUME
C     SOUS LA FORME DE COEFFICIENTS DIAGONAUX DU TYPE DE LA MASSE
C     -----------------------------------------------------------
C     LE VOLUME SUPPORTE T IL UN COEFFICIENT DU DEPLACEMENT?
C     ADRESSE MCN DU COEFFICIENT DU DEPLACEMENT
      MN = LTDEVO(LPCOED,NOOBVC)
      IF( MN .GT. 0 ) THEN
C
C        JACOBIEN * POIDS  CALCULE AU DESSUS
C        DELPOI = DELTA / 24D0
         K = 0
         N = 0
         DO 28 L=1,4
C
C           RECHERCHE DES COMPOSANTES DU COEFFICIENT DU DEPLACEMENT
C           AU POINT D'INTEGRATION L = SOMMET L DU TRIANGLE
            XD = X(L,1)
            YD = X(L,2)
            ZD = X(L,3)
            CALL RECOED( 4, NOOBVC, 3, XD,YD,ZD, MN, ELAS )
C
C           RAIDEUR = RAIDEUR + COEFFICIENT * POIDS * DELTA
C           RANGEMENT DIRECT DES DL PAR NOEUDS
            K = K + 1
            N = N + K
            AE(N) = AE(N) + ELAS(1) * DELPOI
            K = K + 1
            N = N + K
            AE(N) = AE(N) + ELAS(2) * DELPOI
            K = K + 1
            N = N + K
            AE(N) = AE(N) + ELAS(3) * DELPOI
C
 28      CONTINUE
C
      ENDIF
C
C     ============================================================
C     CONTRIBUTION DES FACES A LA CONDITION DE DIRICHLET PENALISEE
C     ============================================================
      IF( PENALI .EQ. 0D0 ) RETURN
C
      DO 100 K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE UTILISATEUR
            MN = LTDESU(LPFIXA,NOOB)
            IF( MN .LE. 0 ) GOTO 100
C           UN TABLEAU FIXATION EXISTE POUR CETTE SURFACE
C
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
C           LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
            NBCOFI = MCN( MN + WBCOFI )
C
            DO 60 J=1,NBCOFI
C
C              LE NUMERO DE LA COMPOSANTE FIXEE
               M = MCN( MN + WUCOFI - 1 + J )
C
               DO 40 L=1,3
C
C                 LE NUMERO DANS LE TETRAEDRE DU SOMMET L DE LA FACE K
                  N = NONOFK( L )
C
C                 LE NUMERO DU DL
                  N = N * 3 - 3 + M
C
C                 SOMMATION AVEC LA MATRICE ELEMENTAIRE
                  N = N * ( N + 1 ) / 2
                  AE( N ) = AE( N ) + PENALI
C
 40            CONTINUE
 60         CONTINUE
         ENDIF
 100  CONTINUE
C
C     =============================================================
C     CONTRIBUTION DES ARETES A LA CONDITION DE DIRICHLET PENALISEE
C     =============================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 140 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION PENALISEE
               MN = LTDELI(LPFIXA,NOOB)
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU FIXATION EXISTE POUR CETTE LIGNE
C                 FIXATION PENALISEE = PENALI
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 111, 112, 113, 114, 115, 116 ) , K
 111              N1 = 1
                  N2 = 2
                  GOTO 118
 112              N1 = 2
                  N2 = 3
                  GOTO 118
 113              N1 = 3
                  N2 = 1
                  GOTO 118
 114              N1 = 1
                  N2 = 4
                  GOTO 118
 115              N1 = 2
                  N2 = 4
                  GOTO 118
 116              N1 = 3
                  N2 = 4
C
C                 LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
 118              NBCOFI = MCN( MN + WBCOFI )
C
                  DO 120 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL DU 1-ER SOMMET
                     N = N1 * 3 - 3 + M
C                    SOMMATION AVEC LA MATRICE ELEMENTAIRE
                     N = N * ( N + 1 ) / 2
                     AE( N ) = AE( N ) + PENALI
C
C                    LE NUMERO DU DL DU 2-EME SOMMET
                     N = N2 * 3 - 3 + M
C                    SOMMATION AVEC LA MATRICE ELEMENTAIRE
                     N = N * ( N + 1 ) / 2
                     AE( N ) = AE( N ) + PENALI
C
 120              CONTINUE
               ENDIF
            ENDIF
 140     CONTINUE
C
C     ==============================================================
C     CONTRIBUTION DES SOMMETS A LA CONDITION DE DIRICHLET PENALISEE
C     ==============================================================
         DO 200 K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION PENALISEE
               MN = LTDEPO( LPFIXA, NOOB )
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU FIXATION EXISTE POUR CE POINT
C                 FIXATION PENALISEE = PENALI
C                 LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
                  NBCOFI = MCN( MN + WBCOFI )
C
                  DO 160 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL DU SOMMET
                     N = K * 3 - 3 + M
C
C                    SOMMATION AVEC LA MATRICE ELEMENTAIRE
                     N = N * ( N + 1 ) / 2
                     AE( N ) = AE( N ) + PENALI
C
 160              CONTINUE
C
               ENDIF
            ENDIF
 200     CONTINUE
      ENDIF
C
      RETURN
      END
