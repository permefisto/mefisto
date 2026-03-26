      SUBROUTINE TR3P1D( X,      DELTA,  DP,     PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   COND,   CONDUC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CONDUCTIVITE D'UN TETRAEDRE 3P1D
C -----    INTEGRATION AUX SOMMETS DU TETRAEDRE ET DES FACES
C
C ENTREES:
C --------
C X      : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C DELTA  : JACOBIEN DE LA TRANSFORMATION EF REFERENCE -> EF
C DP     : GRADIENT DES POLYNOMES DE BASE AUX SOMMETS DU TETRAEDRE
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
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
C NOOBVO : NUMERO DE L'OBJET VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE
C          DES OBJETS VOLUMES
C COND   : 6 REEL2 POUR LE TENSEUR SYMETRIQUE DE CONDUCTIVITE
C
C SORTIES:
C --------
C CONDUC : MATRICE ELEMENTAIRE SYMETRIQUE DE CONDUCTIVITE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C AJOUTS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  DP(3,4),
     %                  CONDUC(10),
     %                  Rho
      DOUBLE PRECISION  COND(6)
      REAL              X(4,3)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
      INTEGER           NOOBPS(1:4),
     %                  NOOBLA(1:6),
     %                  NOOBSF(1:4)
      INTEGER           NONOFK(3)
      DOUBLE PRECISION  C(6),
     %                  COETEM,
     %                  AUX(12),
     %                  ECHANG,
     %                  PENALI,
     %                  DELTA, DELTAK, DGL(2,3), XYZPI(3)
C
C     ======================
C     CONTRIBUTION DU VOLUME
C     ======================
      MOREE2 = MOTVAR(6)
      CALL AZEROD( 10, CONDUC )

      IF( TESTNL .GT. 0 ) CALL NLDATA0( 4 )
C     RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI ou
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C
C     SI LA CONDUCTIVITE N'EST PAS DECLAREE, SAUT DU CALCUL DE LA CONDUCTIVITE
      IF( LTDEVO(LPCOND,JEU,NOOBVO) .EQ. 0 ) GOTO 35
C
C     CONTRIBUTION DE LA CONDUCTIVITE
C     -------------------------------
      CALL AZEROD( 6, C )
C
      DO 20 K=1,4
C
         IF( TESTNL .GT. 0 ) THEN
C           RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION K=SOMMET K
            TEMPEL=DMCN( (MNTHET-1)/2 + MCN(MNNODL+K-1) )
         ENDIF
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C        D'INTEGRATION K
         XYZPI(1) = X(K,1)
         XYZPI(2) = X(K,2)
         XYZPI(3) = X(K,3)
         CALL RECOND( 4, NOOBVO, 3, XYZPI,
     %                LTDEVO(LPCOND,JEU,NOOBVO), COND )
         DO 10 I=1,6
            C(I) = C(I) + COND(I)
 10      CONTINUE

 20   CONTINUE

      IF( TESTNL .EQ. 9 ) THEN
C        [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
         DO I=1,6
            C(I) = -C(I)
         ENDDO
      ENDIF
C
C     C = CONDUCTIVITE * DELTA * POIDS
      DO 30 I=1,6
         C(I) = C(I) * DELTA / 24D0
 30   CONTINUE
C
C     CONDUC = T(DP) * C * (DP)
      CALL TABA8D( 4, 3, DP, C, CONDUC, AUX )
C
C     CONTRIBUTION DU COEFFICIENT DEVANT LA TEMPERATURE
C     -------------------------------------------------
 35   MNLT = LTDEVO(LPCOET,JEU,NOOBVO)
      IF( MNLT .GT. 0 ) THEN
C
C        ADRESSE DMCN -1 DU VECTEUR GLOBAL TEMPERATURE
         MNT = (MNTHET-1) / MOREE2
         N   = 0
         DO 40 L=1,4
C
            IF( TESTNL .GT. 0 ) THEN
C              NUMERO GLOBAL DU SOMMET L DE L'EF
               NS = MCN(MNNODL-1+L)
C              RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL = DMCN( MNT+NS )
               ONDEPI = 0D0
C
               IF( TESTNL .GE. 6 ) THEN
C                 ONDE COMPLEXE NLSE
C                 PARTIE IMAGINAIRE ACTUELLE AU POINT D'INTEGRATION L de L'ONDE
                  ONDEPI = DMCN( MNT+NBNOEMA+NS )
C
C                 PARTIE REELLE INITIALE AU POINT D'INTEGRATION L de L'ONDE
                  MNT0 = (MNTHET0-1)/MOREE2
                  TEMPEL0 = DMCN( MNT0+NS )
C
C                 PARTIE IMAGINAIRE INITIALE AU POINT D'INTEGRATION L de L'ONDE
                  ONDEPI0 = DMCN( MNT0+NBNOEMA+NS )
               ENDIF
            ENDIF
C
C           RECHERCHE DE LA CAPACITE AU POINT D'INTEGRATION L
            XYZPI(1) = X(L,1)
            XYZPI(2) = X(L,2)
            XYZPI(3) = X(L,3)

            IF( TESTNL .LT. 5 ) THEN
C              COEFFICIENT DEVANT LA TEMPERATURE
               CALL RECOET( 4, NOOBVO, 3, XYZPI, MNLT, COETEM )
            ELSE
C              COEFFICIENT DEVANT L'ONDE NLSE A L'INSTANT TEMPS
               CALL RENLSE( 4, NOOBVO, 3, XYZPI, TEMPS, TEMPEL, ONDEPI,
     %                      MNLT, COETEM )

               IF( TESTNL .EQ. 9 ) THEN
C                 [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
                  COETEM = -COETEM
               ENDIF

               IF( TESTNL .EQ. 6 .OR. TESTNL .EQ. 9 ) THEN
C                 GROSS-PITAEVSKI DEMANDE LA DENSITE DE MASSE Rho
                  MN = LTDEVO(LPMAST,JEU,NOOBVO)
                  IF( MN .GT. 0 ) THEN
                     CALL REMASS( 4, NOOBVO, 3, XYZPI, MN, Rho )
                     COETEM = Rho/PasTemps + COETEM
                  ENDIF
               ENDIF

            ENDIF
C
C           COEF TEMPERATURE  = COETEM * DELTA * POIDS
            N = N + L
            CONDUC(N) = CONDUC(N) + COETEM * DELTA / 24D0
C
 40      CONTINUE
      ENDIF
C
C     ======================
C     CONTRIBUTION DES FACES
C     ======================
      DO 100 K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE UTILISATEUR
            IECHAN = 0
            IF( LTDESU(LPECHA,JEU,NOOB) .GT. 0 ) IECHAN = 1
            IF( LTDESU(LPCONT,JEU,NOOB) .GT. 0  .AND.
     %          PENALI .NE. 0D0            ) IECHAN = 2
            IF( IECHAN .EQ. 0 ) GOTO 100
C
C           UN TABLEAU ECHANGE EXISTE POUR CETTE SURFACE
C           CALCUL DE LA CONTRIBUTION A LA MATRICE DE CONDUCTIVITE
C           ------------------------------------------------------
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C                     CE SONT AUSSI LES POINTS D'INTEGRATION
C
C           RECHERCHE DU JACOBIEN DE G
            N = NONOFK(1)
            DGL(1,1) = X( NONOFK(2), 1 ) - X( N, 1 )
            DGL(2,1) = X( NONOFK(3), 1 ) - X( N, 1 )
            DGL(1,2) = X( NONOFK(2), 2 ) - X( N, 2 )
            DGL(2,2) = X( NONOFK(3), 2 ) - X( N, 2 )
            DGL(1,3) = X( NONOFK(2), 3 ) - X( N, 3 )
            DGL(2,3) = X( NONOFK(3), 3 ) - X( N, 3 )
            CALL JAR2R3( DGL, DELTAK )
C
            DO 80 L=1,3
C
C              LE NUMERO DANS LE TETRAEDRE DU SOMMET L DE LA FACE K
               N  = NONOFK( L )
C
               IF( TESTNL .GT. 0 ) THEN
C                 LA TEMPERATURE AU L DE LA FACE K
                  TEMPEL=DMCN((MNTHET-1)/2+MCN(MNNODL+N-1))
               ENDIF
C
               IF( IECHAN .EQ. 1 ) THEN
C
C                 CALCUL DU COEFFICIENT D'ECHANGE
                  XYZPI(1) = X( N, 1 )
                  XYZPI(2) = X( N, 2 )
                  XYZPI(3) = X( N, 3 )
C                 RECHERCHE DU COEFFICIENT D ECHANGE AU SOMMET L DE LA FACE K
                  CALL REECHA( 3, NOOB, 3, XYZPI,
     %                         LTDESU(LPECHA,JEU,NOOB), ECHANG )
C
               ELSE
C
C                 CALCUL DU COEFFICIENT D'ECHANGE = CONTACT PENALISE
                  ECHANG = PENALI
C
               ENDIF
C
C              SOMMATION AVEC LA MATRICE ELEMENTAIRE
               N = N * ( N + 1 ) / 2
               CONDUC( N ) = CONDUC( N ) + ECHANG * DELTAK / 6D0
 80         CONTINUE
         ENDIF
 100  CONTINUE
C
C     =============================================================
C     CONTRIBUTION DES ARETES A LA CONDITION DE DIRICHLET PENALISEE
C     =============================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 120 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 COEFFICIENT D'ECHANGE = CONTACT PENALISE = PENALI
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 111, 112, 113, 114, 115, 116 ) , K
 111              N = 1
                  L = 2
                  GOTO 118
 112              N = 2
                  L = 3
                  GOTO 118
 113              N = 3
                  L = 1
                  GOTO 118
 114              N = 1
                  L = 4
                  GOTO 118
 115              N = 2
                  L = 4
                  GOTO 118
 116              N = 3
                  L = 4
C
C                 COEFFICIENT DIAGONAL DE LA MATRICE DE CONDUCTIVITE ELEMENTAIRE
 118              N = N * ( N + 1 ) / 2
                  CONDUC(N) = CONDUC(N) + PENALI
C
C                 COEFFICIENT DIAGONAL DE LA MATRICE DE CONDUCTIVITE ELEMENTAIRE
                  N = L * ( L + 1 ) / 2
                  CONDUC(N) = CONDUC(N) + PENALI
               ENDIF
            ENDIF
 120     CONTINUE
C
C     ==============================================================
C     CONTRIBUTION DES SOMMETS A LA CONDITION DE DIRICHLET PENALISEE
C     ==============================================================
         DO 150 K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDEPO(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 COEFFICIENT D'ECHANGE = CONTACT PENALISE = PENALI
C
C                 COEFFICIENT DIAGONAL DE LA MATRICE DE CONDUCTIVITE ELEMENTAIRE
                  N = K * ( K + 1 ) / 2
                  CONDUC(N) = CONDUC(N) + PENALI
C
               ENDIF
            ENDIF
 150     CONTINUE
C
      ENDIF
C
      RETURN
      END
