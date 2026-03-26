      SUBROUTINE ER2P1D( X,      PENALI,
     &                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     &                   NOOBLA, NUMILI, NUMALI, LTDELI,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   ELAS,   RAIDEU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE RAIDEUR DU TRIANGLE 2P1D
C -----    LAGRANGE DE DEGRE 1. INTEGRATION AUX 3 SOMMETS DU TRIANGLE
C
C ENTREES:	
C --------
C X      : 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C PENALI : 1/EPSILON DE LA PENALISATION DU DEPLACEMENT IMPOSE
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL  DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL  DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES POINTS DE L'OBJET
C
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES LIGNES DE L'OBJET
C NUMALI : NUMERO MAXIMAL DES LIGNES DE L'OBJET
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES LIGNES DE L'OBJET
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES SURFACES DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C
C SORTIES:
C --------
C ELAS   : TENSEUR SYMETRIQUE DE L ELASTICITE SYMETRIQUE
C          PUIS COEFFICIENT DU DEPLACEMENT
C RAIDEU : MATRICE DE RAIDEUR 6X6 STOCKEE SYMETRIQUE PLEINE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1999
C23456---------------------------------------------------------------012
      include "./incl/donela.inc"
      include "./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      REAL              X(3,2)
      DOUBLE PRECISION  PENALI,
     &                  ELAS(6),
     &                  RAIDEU(21)
      INTEGER           NOOBPS(3), NOOBLA(3), NOOBSF
      INTEGER           LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
C
      DOUBLE PRECISION  DELTA, DELPOI, XD, YD
      DOUBLE PRECISION  X21, Y12, X13, Y31, X32, Y23
      DOUBLE PRECISION  X21X21,X21X13,X21X32,
     %                         X13X13,X13X32,
     %                                X32X32
      DOUBLE PRECISION  Y12Y12,Y12Y31,Y12Y23,
     %                         Y31Y31,Y31Y23,
     %                                Y23Y23
      DOUBLE PRECISION  X21Y12,X21Y31,X21Y23,
     %                  X13Y12,X13Y31,X13Y23,
     %                  X32Y12,X32Y31,X32Y23
C
      DOUBLE PRECISION  E(6),E1,E2,E3,E4,E5,E6
      EQUIVALENCE      (E(1),E1),(E(2),E2),(E(3),E3),
     %                 (E(4),E4),(E(5),E5),(E(6),E6)
C
C     CONTRIBUTION DE LA SURFACE
C     ==========================
      X21 = X(2,1) - X(1,1)
      X13 = X(1,1) - X(3,1)
      X32 = X(3,1) - X(2,1)
C
      Y12 = X(1,2) - X(2,2)
      Y31 = X(3,2) - X(1,2)
      Y23 = X(2,2) - X(3,2)
C
      X21X21 = X21 * X21
      X21X13 = X21 * X13
      X21X32 = X21 * X32
C
      X13X13 = X13 * X13
      X13X32 = X13 * X32
C
      X32X32 = X32 * X32
C
      Y12Y12 = Y12 * Y12
      Y12Y31 = Y12 * Y31
      Y12Y23 = Y12 * Y23
C
      Y31Y31 = Y31 * Y31
      Y31Y23 = Y31 * Y23
C
      Y23Y23 = Y23 * Y23
C
      X21Y12 = X21 * Y12
      X21Y31 = X21 * Y31
      X21Y23 = X21 * Y23
C
      X13Y12 = X13 * Y12
      X13Y31 = X13 * Y31
      X13Y23 = X13 * Y23
C
      X32Y12 = X32 * Y12
      X32Y31 = X32 * Y31
      X32Y23 = X32 * Y23
C
C     POIDS / JACOBIEN
      DELTA  = ABS( X21Y31 - X13Y12 )
      DELPOI = 1D0 / ( 6D0 * DELTA )
C
      DO 5 L=1,6
         E(L) = 0D0
 5    CONTINUE
C
      DO 20 K=1,3
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE L'ELASTICITE AU SOMMET K
         XD = X(K,1)
         YD = X(K,2)
         CALL REELAS( 3, NOOBSF, 2, XD,YD,0D0,
     %                LTDESU(LPYOUN,NOOBSF), ELAS )
C
C        LA SOMME DES COEFFICIENTS DE L'ELASTICITE AUX 3 SOMMETS
         DO 15 L=1,6
            E(L) = E(L) + ELAS(L)
 15      CONTINUE
C
 20   CONTINUE
C
C     E = ELASTICITE * POIDS / JACOBIEN
      DO 25 L=1,6
         E(L) = E(L) * DELPOI
 25   CONTINUE
C
C     CONSTRUCTION DES MATRICES   ( A11   A12 )
C     -------------------------   (       A22 )
C     LA SOUS-MATRICE SYMETRIQUE A11
      IF( E4 .EQ. 0D0 ) THEN
C
         RAIDEU( 1) =   E1 * Y23Y23
     %                + E6 * X32X32
C
         RAIDEU( 4) =   E1 * Y31Y23
     %                + E6 * X13X32
C
         RAIDEU( 6) =   E1 * Y31Y31
     %                + E6 * X13X13
C
         RAIDEU(11) =   E1 * Y12Y23
     %                + E6 * X21X32
C
         RAIDEU(13) =   E1 * Y12Y31
     %                + E6 * X21X13
C
         RAIDEU(15) =   E1 * Y12Y12
     %                + E6 * X21X21
C
      ELSE
C
         RAIDEU( 1) =   E1 * Y23Y23
     %                + E4 * X32Y23 * 2
     %                + E6 * X32X32
C
         RAIDEU( 4) =   E1 * Y31Y23
     %                + E4 * ( X32Y31 + X13Y23 )
     %                + E6 * X13X32
C
         RAIDEU( 6) =   E1 * Y31Y31
     %                + E4 * X13Y31 * 2
     %                + E6 * X13X13
C
         RAIDEU(11) =   E1 * Y12Y23
     %                + E4 * ( X32Y12 + X21Y23 )
     %                + E6 * X21X32
C
         RAIDEU(13) =   E1 * Y12Y31
     %                + E4 * ( X13Y12 + X21Y31 )
     %                + E6 * X21X13
C
         RAIDEU(15) =   E1 * Y12Y12
     %                + E4 * X21Y12 * 2
     %                + E6 * X21X21
      ENDIF
C
C     LA SOUS+MATRICE SYMETRIQUE A22
      IF( E5 .EQ. 0D0 ) THEN
C
         RAIDEU( 3) =   E6 * Y23Y23
     %                + E3 * X32X32
C
         RAIDEU( 8) =   E6 * Y31Y23
     %                + E3 * X13X32
C
         RAIDEU(10) =   E6 * Y31Y31
     %                + E3 * X13X13
C
         RAIDEU(17) =   E6 * Y12Y23
     %                + E3 * X21X32
C
         RAIDEU(19) =   E6 * Y12Y31
     %                + E3 * X21X13
C
         RAIDEU(21) =   E6 * Y12Y12
     %                + E3 * X21X21
C
      ELSE
C
         RAIDEU( 3) =   E6 * Y23Y23
     %                + E5 * X32Y23 * 2
     %                + E3 * X32X32
C
         RAIDEU( 8) =   E6 * Y31Y23
     %                + E5 * ( X32Y31 + X13Y23 )
     %                + E3 * X13X32
C
         RAIDEU(10) =   E6 * Y31Y31
     %                + E5 * X13Y31 * 2
     %                + E3 * X13X13
C
         RAIDEU(17) =   E6 * Y12Y23
     %                + E5 * ( X32Y12 + X21Y23 )
     %                + E3 * X21X32
C
         RAIDEU(19) =   E6 * Y12Y31
     %                + E5 * ( X13Y12 + X21Y31 )
     %                + E3 * X21X13
C
         RAIDEU(21) =   E6 * Y12Y12
     %                + E5 * X21Y12 * 2
     %                + E3 * X21X21
      ENDIF
C
C     LA SOUS+MATRICE 3 NON SYMETRIQUE A12
      IF( E4 .EQ. 0D0 .AND. E5 .EQ. 0D0 ) THEN
C
         RAIDEU( 2) =  (E2+E6) * X32Y23
C
         RAIDEU( 5) =   E2 * X32Y31
     %                + E6 * X13Y23
C
         RAIDEU(12) =   E2 * X32Y12
     %                + E6 * X21Y23
C
         RAIDEU( 7) =   E2 * X13Y23
     %                + E6 * X32Y31
C
         RAIDEU( 9) =  (E2+E6) * X13Y31
C
         RAIDEU(14) =   E2 * X13Y12
     %                + E6 * X21Y31
C
         RAIDEU(16) =   E2 * X21Y23
     %                + E6 * X32Y12
C
         RAIDEU(18) =   E2 * X21Y31
     %                + E6 * X13Y12
C
         RAIDEU(20) =  (E2+E6) * X21Y12
C
      ELSE
C
         RAIDEU( 2) =     E4   * Y23Y23
     %                +(E2+E6) * X32Y23
     %                +   E5   * X32X32
C
         RAIDEU( 5) =   E4 * Y31Y23
     %                + E2 * X32Y31
     %                + E6 * X13Y23
     %                + E5 * X13X32
C
         RAIDEU(12) =   E4 * Y12Y23
     %                + E2 * X32Y12
     %                + E6 * X21Y23
     %                + E5 * X21X32
C
         RAIDEU( 7) =   E4 * Y31Y23
     %                + E2 * X13Y23
     %                + E6 * X32Y31
     %                + E5 * X13X32
C
         RAIDEU( 9) =     E4   * Y31Y31
     %                +(E2+E6) * X13Y31
     %                +   E5   * X13X13
C
         RAIDEU(14) =   E4 * Y12Y31
     %                + E2 * X13Y12
     %                + E6 * X21Y31
     %                + E5 * X21X13
C
         RAIDEU(16) =   E4 * Y12Y23
     %                + E2 * X21Y23
     %                + E6 * X32Y12
     %                + E5 * X21X32
C
         RAIDEU(18) =   E4 * Y12Y31
     %                + E2 * X21Y31
     %                + E6 * X13Y12
     %                + E5 * X21X13
C
         RAIDEU(20) =     E4   * Y12Y12
     %                +(E2+E6) * X21Y12
     %                +   E5   * X21X21
      ENDIF
C
C     CONTRIBUTION DU COEFFICIENT DU DEPLACEMENT SUR LA SURFACE
C     SOUS LA FORME DE COEFFICIENTS DIAGONAUX DU TYPE DE LA MASSE
C     -----------------------------------------------------------
C     LA SURFACE SUPPORTE T ELLE UN COEFFICIENT DU DEPLACEMENT?
C     ADRESSE MCN DU COEFFICIENT DU DEPLACEMENT
      MN = LTDESU(LPCOED,NOOBSF)
      IF( MN .GT. 0 ) THEN
C
C        JACOBIEN * POIDS
         DELPOI = DELTA / 6D0
         K = 0
         N = 0
         DO 28 L=1,3
C
C           RECHERCHE DES COMPOSANTES DU COEFFICIENT DU DEPLACEMENT
C           AU POINT D'INTEGRATION L = SOMMET L DU TRIANGLE
            XD = X(L,1)
            YD = X(L,2)
            CALL RECOED( 3, NOOBSF, 2, XD,YD,0D0, MN, ELAS )
C
C           RAIDEUR = RAIDEUR + COEFFICIENT * POIDS * DELTA
C           RANGEMENT DIRECT DES DL PAR NOEUDS
            K = K + 1
            N = N + K
            RAIDEU(N) = RAIDEU(N) + ELAS(1) * DELPOI
            K = K + 1
            N = N + K
            RAIDEU(N) = RAIDEU(N) + ELAS(2) * DELPOI
C
 28      CONTINUE
C
      ENDIF
C
      IF( PENALI .EQ. 0D0 ) RETURN
C
C     CONTRIBUTION DES ARETES A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     ------------------------------------------------------------------
      DO 35 K = 1,3
C
C        LA LIGNE DU COTE K SUPPORTE T ELLE UNE FIXATION ?
         N = NOOBLA(K)
         IF( N .GT. 0 ) THEN
C
C           ADRESSE MCN DE LA FIXATION DE LA LIGNE
            MN = LTDELI(LPFIXA,N)
            IF( MN .GT. 0 ) THEN
C
C              LE NUMERO DES 2 SOMMETS DU COTE K
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
C              LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
               NBCOFI = MCN( MN + WBCOFI )
C
               DO 30 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  L = MCN( MN + WUCOFI - 1 + J )
C
C                 LA PENALISATION EST AJOUTEE POUR LE DEPLACEMENT L AU SOMMET K
                  I  = K * 2 - 2 + L
                  IA = I * (I+1) / 2
                  RAIDEU( IA ) = RAIDEU( IA ) + PENALI
C
C                 LA PENALISATION EST AJOUTEE POUR LE DEPLACEMENT L AU SOMMET KK
                  I = KK * 2 - 2 + L
                  I = I * (I+1) / 2
                  RAIDEU( I ) = RAIDEU( I ) + PENALI
 30            CONTINUE
C
            ENDIF
         ENDIF
 35   CONTINUE
C
C
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     -------------------------------------------------------------------
      DO 50 K=1,3
C
C        NUMERO DE POINT DU SOMMET K
         N = NOOBPS(K)
         IF( N .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION?
            MN = LTDEPO(LPFIXA,N)
            IF( MN .GT. 0 ) THEN
C
C              OUI: LE NOMBRE DE COMPOSANTES FIXEES
               NBCOFI = MCN( MN + WBCOFI )
C
               DO 40 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  L = MCN( MN + WUCOFI - 1 + J )
C
C                 LA PENALISATION EST AJOUTEE
                  I = K * 2 - 2 + L
                  I = I * (I+1) / 2
                  RAIDEU( I ) = RAIDEU( I ) + PENALI
C
 40            CONTINUE
C
            ENDIF
         ENDIF
 50   CONTINUE
C
      RETURN
      END
