      SUBROUTINE TR2P1D( X,      PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   COND,   CONDUC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CONDUCTIVITE DU TRIANGLE 2P1D DROIT
C -----    LAGRANGE DE DEGRE 1
C          INTEGRATION NUMERIQUE AUX SOMMETS DU TRIANGLE ET DES ARETES
C
C ENTREES:
C --------
C X      : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
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
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE
C          DES OBJETS SURFACES
C COND   : 6 REEL2 POUR LE TENSEUR SYMETRIQUE DE CONDUCTIVITE
C
C SORTIES:
C --------
C CONDUC : MATRICE ELEMENTAIRE DE CONDUCTIVITE SYMETRIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C AJOUTS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
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
      DOUBLE PRECISION  CONDUC(6)
      DOUBLE PRECISION  ECHANG, PENALI
      REAL              X(3,2)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      INTEGER           NOOBPS(1:3)
      INTEGER           NOOBLA(1:3)
      DOUBLE PRECISION  COND(6), C1, C2, C3, Rho
      DOUBLE PRECISION  DELTA, XYZPI(3)
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32
C
C     ==========================
C     CONTRIBUTION DE LA SURFACE
C     ==========================
      MOREE2 = MOTVAR(6)
      MNT    = 0
      CALL AZEROD( 6, CONDUC )
C
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
C
      DELTA = ABS( X21 * Y31 - X31 * Y21 ) * 6D0

      IF( TESTNL .GT. 0 ) CALL NLDATA0( 3 )
C     RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI ou
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C
C     SI LA CONDUCTIVITE N'EST PAS DECLAREE, SAUT DU CALCUL DE LA CONDUCTIVITE
      IF(  LTDESU(LPCOND,JEU,NOOBSF) .EQ. 0 ) GOTO 18
C
C     CONTRIBUTION DE LA CONDUCTIVITE
C     -------------------------------
      C1 = 0D0
      C2 = 0D0
      C3 = 0D0
C
      DO 10 K=1,3
C
         IF( TESTNL .GT. 0 ) THEN
C           RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION K = SOMMET K
            TEMPEL=DMCN( (MNTHET-1)/MOREE2 + MCN(MNNODL+K-1) )
         ENDIF
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C        D'INTEGRATION K
         XYZPI(1) = X(K,1)
         XYZPI(2) = X(K,2)
         XYZPI(3) = 0D0
C        RECHERCHE DE LA CONDUCTIVITE AU POINT D'INTEGRATION K
         CALL RECOND( 3, NOOBSF, 3, XYZPI,
     %                LTDESU(LPCOND,JEU,NOOBSF), COND )
C
C        LA SOMME DES CONDUCTIVITES AUX 3 SOMMETS
         C1 = C1 + COND(1)
         C2 = C2 + COND(2)
         C3 = C3 + COND(3)
C
 10   CONTINUE

      IF( TESTNL .EQ. 9 ) THEN
C        [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
         C1 = -C1
         C2 = -C2
         C3 = -C3
      ENDIF
C
C     COND = CONDUCTIVITE * POIDS / DELTA
      C1 = C1 / DELTA
      C2 = C2 / DELTA
      C3 = C3 / DELTA
C
      IF( C2 .EQ. 0D0 ) THEN
         CONDUC(1) =   C1 * Y32 * Y32
     %               + C3 * X32 * X32
C
         CONDUC(2) = - C1 * Y32 * Y31
     %               - C3 * X32 * X31
C
         CONDUC(3) =   C1 * Y31 * Y31
     %               + C3 * X31 * X31
C
         CONDUC(4) =   C1 * Y32 * Y21
     %               + C3 * X32 * X21
C
         CONDUC(5) = - C1 * Y31 * Y21
     %               - C3 * X31 * X21
C
         CONDUC(6) =   C1 * Y21 * Y21
     %               + C3 * X21 * X21
C
      ELSE
C
         CONDUC(1) =   C1 * Y32 * Y32
     %               - C2 * X32 * Y32 * 2
     %               + C3 * X32 * X32
C
         CONDUC(2) = - C1 * Y32 * Y31
     %               + C2 * ( X32 * Y31 + Y32 * X31 )
     %               - C3 * X32 * X31
C
         CONDUC(3) =   C1 * Y31 * Y31
     %               - C2 * X31 * Y31 * 2
     %               + C3 * X31 * X31
C
         CONDUC(4) =   C1 * Y32 * Y21
     %               - C2 * ( X32 * Y21 + Y32 * X21 )
     %               + C3 * X32 * X21
C
         CONDUC(5) = - C1 * Y31 * Y21
     %               + C2 * ( X31 * Y21 + Y31 * X21 )
     %               - C3 * X31 * X21
C
         CONDUC(6) =   C1 * Y21 * Y21
     %               - C2 * X21 * Y21 * 2
     %               + C3 * X21 * X21
      ENDIF
C
C     CONTRIBUTION DU COEFFICIENT DEVANT LA TEMPERATURE
C     -------------------------------------------------
 18   MNLT = LTDESU(LPCOET,JEU,NOOBSF)
      IF( MNLT .GT. 0 ) THEN
         DELTA = DELTA / 36D0
         K1    = 0
C        ADRESSE DMCN -1 DU VECTEUR GLOBAL TEMPERATURE
         MNT    = (MNTHET-1) / MOREE2
         DO 20 K=1,3
C
            IF( TESTNL .GT. 0 ) THEN
C              NUMERO GLOBAL DU SOMMET K DE L'EF
               N = MCN(MNNODL-1+K)
C              RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION K
               TEMPEL = DMCN( MNT+N )
               ONDEPI = 0D0
C
               IF( TESTNL .GE. 6 ) THEN
C                 ONDE COMPLEXE NLSE
C                 PARTIE IMAGINAIRE ACTUELLE AU POINT D'INTEGRATION K de L'ONDE
                  ONDEPI = DMCN( MNT+NBNOEMA+N )
C
C                 PARTIE REELLE INITIALE AU POINT D'INTEGRATION K de L'ONDE
                  MNT0 = (MNTHET0-1)/MOREE2
                  TEMPEL0 = DMCN( MNT0+N )
C
C                 PARTIE IMAGINAIRE INITIALE AU POINT D'INTEGRATION K de L'ONDE
                  ONDEPI0 = DMCN( MNT0+NBNOEMA+N )
               ENDIF
            ENDIF
C
C           RECHERCHE DU COEFFICIENT DE LA TEMPERATURE AU SOMMET K
            XYZPI(1) = X(K,1)
            XYZPI(2) = X(K,2)
            XYZPI(3) = 0D0
            IF( TESTNL .LT. 5 ) THEN
C              COEFFICIENT DEVANT LA TEMPERATURE
               CALL RECOET( 3, NOOBSF, 3, XYZPI, MNLT, C1 )
            ELSE
C              COEFFICIENT DEVANT L'ONDE NLSE A L'INSTANT TEMPS
               CALL RENLSE( 3, NOOBSF, 3, XYZPI, TEMPS, TEMPEL, ONDEPI,
     %                      MNLT, C1 )

               IF( TESTNL .EQ. 9 ) THEN
C                 [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
                  C1 = -C1
               ENDIF

               IF( TESTNL .EQ. 6 .OR. TESTNL .EQ. 9 ) THEN
C                 GROSS-PITAEVSKI DEMANDE LA DENSITE DE MASSE Rho
                  MN = LTDESU(LPMAST,JEU,NOOBSF)
                  IF( MN .GT. 0 ) THEN
                     CALL REMASS( 3, NOOBSF, 3, XYZPI, MN, Rho )
                     C1 = Rho/PasTemps + C1
                  ENDIF
               ENDIF
            ENDIF
C
C           COEF TEMPERATURE  = COEF TE  * DELTA * POIDS
            K1 = K1 + K
            CONDUC(K1) = CONDUC(K1) + C1 * DELTA
C
 20      CONTINUE
      ENDIF
C
CCCC
CCCC     CONTRIBUTION EVENTUELLE DE (VITESSE DU FLUIDE * GRADIENT TEMPERATURE)
CCCC     ATTENTION: LA MATRICE N'EST ALORS PLUS SYMETRIQUE!!!...
CCCC     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
CCC      IF( LTDESU(LPVIFL,JEU,NOOBSF) .GT. 0 ) THEN
CCC         CALL REVIFL( 3, NOOBSF, 2, 3, XYZPI,
CCC     %                LTDESU(LPVIFL,JEU,NOOBSF), VITEFL )
CCCC        SUITE A TRAITER ...
CCC      ENDIF
CCC
C
C     =======================
C     CONTRIBUTION DES ARETES
C     =======================
      DO 30 K=1,3
C
C        NO DE LIGNE DE L'ARETE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. SUPPORT D'ECHANGE OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDELI(LPECHA,JEU,NOOB) .GT. 0  ) IECHAN = 1
            IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                 ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 30
C
C           UN TABLEAU ECHANGE OU CONTACT PENALISE EXISTE POUR CETTE LIGNE
C           --------------------------------------------------------------
            IF( K .EQ. 1 ) THEN
               DELTA = X21 * X21 + Y21 * Y21
C              LE NUMERO DU SECOND SOMMET DANS LE TRIANGLE
               K1 = 2
C              LA PLACE DES COEFFICIENTS DIAGONAUX POUR LA CONTRIBUTION DE L'ARE
               M1 = 1
               M2 = 3
            ELSE IF( K .EQ. 2 ) THEN
               DELTA = X32 * X32 + Y32 * Y32
               K1 = 3
C              LA PLACE DES COEFFICIENTS DIAGONAUX POUR LA CONTRIBUTION DE L'ARE
               M1 = 3
               M2 = 6
            ELSE
               DELTA = X31 * X31 + Y31 * Y31
               K1 = 1
C              LA PLACE DES COEFFICIENTS DIAGONAUX POUR LA CONTRIBUTION DE L'ARE
               M1 = 6
               M2 = 1
            ENDIF
            DELTA = SQRT( DELTA ) * 0.5D0
C
C           COEFFICIENT D'ECHANGE AU PREMIER SOMMET DE L'ARETE
C           --------------------------------------------------
            IF( TESTNL .GT. 0 ) THEN
C              LA TEMPERATURE AU SOMMET 1 DE L'ARETE K
               MNT    = (MNTHET-1) / MOREE2
               TEMPEL = DMCN( MNT+MCN(MNNODL+K-1) )
            ENDIF
C
            IF( IECHAN .EQ. 1 ) THEN
C
C              CALCUL DU COEFFICIENT D'ECHANGE
               XYZPI(1) = X(K,1)
               XYZPI(2) = X(K,2)
               XYZPI(3) = 0D0
               CALL REECHA( 2,NOOB, 3, XYZPI,
     %                      LTDELI(LPECHA,JEU,NOOB), ECHANG )
            ELSE
C
C              CALCUL DU COEFFICIENT D'ECHANGE = CONTACT PENALISE
               ECHANG = PENALI
C
            ENDIF
C
C           LA CONTRIBUTION DE L'ARETE A LA MATRICE DE CONDUCTIVITE
            CONDUC(M1) = CONDUC(M1) + ECHANG * DELTA
C
C           COEFFICIENT D'ECHANGE AU SECOND SOMMET DE L'ARETE
C           -------------------------------------------------
            IF( TESTNL .GT. 0 ) THEN
C              LA TEMPERATURE AU SOMMET 2 DE L'ARETE K
               TEMPEL=DMCN( MNT+MCN(MNNODL+K1-1) )
            END IF
C
            IF( IECHAN .EQ. 1 ) THEN
C
C              CALCUL DU COEFFICIENT D'ECHANGE
               XYZPI(1) = X(K1,1)
               XYZPI(2) = X(K1,2)
               XYZPI(3) = 0D0
               CALL REECHA( 2,NOOB, 3, XYZPI,
     %                      LTDELI(LPECHA,JEU,NOOB), ECHANG )
CCC         ELSE
C
C              CALCUL DU COEFFICIENT D'ECHANGE = CONTACT PENALISE
CCC            ECHANG = PENALI
C
            ENDIF
C
C           LA CONTRIBUTION DE L'ARETE A LA MATRICE DE CONDUCTIVITE
            CONDUC(M2) = CONDUC(M2) + ECHANG * DELTA
C
         ENDIF
 30   CONTINUE
C
C     ========================
C     CONTRIBUTION DES SOMMETS
C     ========================
      IF( PENALI .NE. 0D0 ) THEN
         DO 50 K=1,3
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
                  M1 = K * ( K + 1 ) / 2
                  CONDUC(M1) = CONDUC(M1) + PENALI
C
               ENDIF
            ENDIF
 50      CONTINUE
      ENDIF
C
      RETURN
      END
