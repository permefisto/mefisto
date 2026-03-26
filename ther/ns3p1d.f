      SUBROUTINE NS3P1D( DeltaT, X, DELTA, PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     NLSE: CALCUL DU SECOND MEMBRE DES TETRAEDRES 3P1D
C -----           LAGRANGE EN SCHEMA IMPLICITE c'est a dire
C TESTNL=6 REAL PART:
C Rho/dt (Wn+1m-Wn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Vn+1m  +Fr
C TESTNL=6 IMAG PART:
C Rho/dt(-Vn+1m+Vn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Wn+1m+1+Fi
C
C TESTNL=7 REAL PART:
C [M(Rho)] {Wn} + [N(DeltaT beta,V,W)] {Vn+1m} - {DeltaT FOmegaR(tn+1,V,W)}
C TESTNL=7 IMAG PART:
C [M(Rho)] {Vn} - [N(DeltaT beta,V,W)] {Wn+1m} + {DeltaT FOmegaI(tn+1,V,W)}
C
C ENTREES:
C --------
C DeltaT : PAS DE TEMPS DU SCHEMA IMPLICITE EN TEMPS
C X      : COORDONNEES RAYON ET COTE DES 4 POINTS DE L'EF
C          OU X Y DES 4 POINTS DE L'EF
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
C NOOBSF : NUMERO DES SURFACES DES FACES DE L'EF
C NUMISU : NUMERO MINIMAL DES SURFACES UTILISEES
C NUMASU : NUMERO MAXIMAL DES SURFACES UTILISEES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C NOOBVO : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMIVO : NUMERO MINIMAL DES OBJETS SURFACES
C NUMAVO : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C SORTIE :
C --------
C BE     : BE(4,2) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C MODIF  : ALAIN PERRONNET TEXAS A & M University at QATAR  FEVRIER 2011
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  PENALI, BE(4,2)
      REAL              X(4,3)
      INTEGER           NOOBPS(1:4),
     %                  NOOBLA(1:6),
     %                  NOOBSF(1:4)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
      INTEGER           NONOFK(3)
      DOUBLE PRECISION  VN(3), DELTA, DELTAK, XYZ(3),
     %                  DGL(2,3), DGLN, DELTAT, D, NLSECOEF, NLSECOEF0
      DOUBLE PRECISION  MASSE, BETA, FORCE(3), FIXA(2)
      INTEGER           K12(2)
      EQUIVALENCE      (K12(1),K1), (K12(2),K2)
C
C     RECUPERATION DE L'ONDE AUX 4 DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     L'ADRESSE MCN DE VECTEUR EST MNTHET de $MEFISTO/incl/cthet.inc
C     ------------------------------------------------------------------
      CALL NLDATA0( 4 )
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C     EN TOUS LES DL DE L'EF
      MOREE2 =  MOTVAR(6)
      MNON   = (MNTHDL -1) / MOREE2
      MNONn  = (MNTHDLn-1) / MOREE2
      MNON0  = (MNTHDL0-1) / MOREE2
C
C     ==================================
C     CONTRIBUTION DES FORCES VOLUMIQUES
C     ==================================
      DO L = 1, 4
C        COORDONNEES DU POINT D'INTEGRATION L = SOMMET L
         XYZ(1) = X(L,1)
         XYZ(2) = X(L,2)
         XYZ(3) = X(L,3)
C
C        LA VALEUR DE L'ONDE INTERPOLEE AU POINT D'INTEGRATION L
         TEMPEL  = DMCN( MNON     + L )
         ONDEPI  = DMCN( MNON + 4 + L )
         TEMPELn = DMCN( MNONn    + L )
         ONDEPIn = DMCN( MNONn+ 4 + L )
         TEMPEL0 = DMCN( MNON0    + L )
         ONDEPI0 = DMCN( MNON0+ 4 + L )
C
C        DENSITE DE MASSE
C        ----------------
         MASSE = 1D0
         MN = LTDEVO(LPMAST,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            CALL REMASS( 4, NOOBVO, 3, XYZ, MN,  MASSE )
         ENDIF
C
C        BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS ACTUEL
C       -BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS INITIAL
C        ( ex: BETA (V**2+W**2) - BETA (V0**2+W0**2) )
C        -------------------------------------------------------
         BETA = 0D0
         MN   = LTDEVO(LPCOET,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            CALL RENLSE( 4, NOOBVO, 3, XYZ, TEMPS, TEMPEL, ONDEPI, MN,
     %                   NLSECOEF )
            CALL RENLSE( 4, NOOBVO, 3, XYZ, TEMPSINI,TEMPEL0,ONDEPI0,MN,
     %                   NLSECOEF0 )
            BETA = NLSECOEF - NLSECOEF0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDEVO(LPSOUR,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            CALL REFORC( 4, NOOBVO, 2, XYZ(1),XYZ(2),XYZ(3),
     %                   TEMPEL, ONDEPI, 0D0, MN,  FORCE )
         ELSE
            FORCE(1) = 0D0
            FORCE(2) = 0D0
         ENDIF
C
C        CALCUL DU VECTEUR * DELTA / 24D0
         IF( TESTNL .EQ. 6 ) THEN
C
C        REAL PART=Rho/dt (Wn+1m-Wn)-Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2)Vn+1m+F
         BE(L,1)=(MASSE/DeltaT*(ONDEPI-ONDEPIn) -BETA*TEMPEL + FORCE(1))
     %            * DELTA / 24D0
C
C        IMAG PART=Rho/dt(-Vn+1m+Vn)-Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2)Wn+1m+F
         BE(L,1)=(MASSE/DeltaT*(TEMPELn-TEMPEL) -BETA*ONDEPI + FORCE(2))
     %            * DELTA / 24D0
C
         ELSE IF( TESTNL .EQ. 7 ) THEN
C
C        PARTIE REELLE    = DeltaT [N(Beta)] V + [M(Dmasse)] Wn - DeltaT FOmegaR
C        PARTIE IMAGINAIRE= [M(Dmasse)] Vn - DeltaT [N(Beta)] W + DeltaT FOmegaI
C        -----------------------------------------------------------------------
C        LE TERME DeltaT [N(Beta)] V + [M(Dmasse)] W - DeltaT FOmegaR
         BE(L,1)= ( DeltaT * ( BETA*TEMPEL - FORCE(1) ) + MASSE*ONDEPIn)
     %            * DELTA / 24D0
C
C        LE TERME [M(Dmasse)] V - DeltaT [N(Beta)] W + DeltaT FOmegaI
         BE(L,2)= ( MASSE*TEMPELn + DeltaT * ( FORCE(2) - BETA*ONDEPI ))
     %            * DELTA / 24D0
C
         ENDIF
C
      ENDDO
C
C     ===========================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES FACES
C     ===========================================================
      DO 80 K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE UTILISATEUR
            IECHAN = 0
            IF( LTDESU(LPSOUR,JEU,NOOB) .GT. 0 ) IECHAN = 1
            IF( LTDESU(LPCONT,JEU,NOOB) .GT. 0  .AND.
     %          PENALI .NE. 0D0                ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 80
C
C           UN TABLEAU CONTACT EXISTE POUR CETTE FACE K
C           CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C           .....................................................
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
C           LE VECTEUR NORMAL UNITAIRE A LA FACE
            CALL VECNOR( DGL, DGLN, VN )
C
C           CALCUL DU FLUX NORMAL
            DO 60 L=1,3
C
C              LE NUMERO N DANS LE TETRAEDRE DU SOMMET L DE LA FACE K
               N = NONOFK( L )
C
C              LA TEMPERATURE AU SOMMET N
               TEMPEL = DMCN( MNON  +N )
               ONDEPI = DMCN( MNON+4+N )
C
C              LES 3 COORDONNEES DU SOMMET L DE LA FACE K
               XYZ(1) = X(N,1)
               XYZ(2) = X(N,2)
               XYZ(3) = X(N,3)
C
               IF( IECHAN .EQ. 1 ) THEN
C                 FORCE(2) REQUISE
C                 LE VECTEUR NORMAL UNITAIRE EST UTILISE
                  CALL REFORC( 3, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                     VN(1), VN(2), VN(3),
     %                         LTDESU(LPSOUR,JEU,NOOB), FORCE )
C
               ELSE
C                 FIXATION(2) PENALISEE
                  MN =  LTDESU(LPCONT,JEU,NOOB)
                  CALL REFIXA( 3, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1) = 0D0
                  FORCE(2) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA I-EME COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU) = FIXA(I) * PENALI
                  ENDDO
C
               ENDIF
C
               D = DELTAK / 6D0
               IF( TESTNL .EQ. 6 ) THEN
C                 SOMMATION AVEC LE VECTEUR ELEMENTAIRE
                  BE(N,1) = BE(N,1) + D * FORCE(1)
                  BE(N,2) = BE(N,2) + D * FORCE(2)
               ELSE IF( TESTNL .EQ. 7 ) THEN
C                 SOMMATION AVEC LE VECTEUR ELEMENTAIRE
                  D = D * DeltaT
                  BE(N,1) = BE(N,1) - D * FORCE(1)
                  BE(N,2) = BE(N,2) + D * FORCE(2)
               ENDIF
C
 60         CONTINUE
         ENDIF
 80   CONTINUE
C
C     ======================================================
C     CONTRIBUTION DE FORCE OU FIXATION PENALISEE AUX ARETES
C     ======================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 130 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT
C              D'UNE SOURCE OU D'UN CONTACT PENALISE?
               IF(LTDELI(LPSOUR,JEU,NOOB).GT.0 .OR.
     %            LTDELI(LPCONT,JEU,NOOB).GT.0 .AND. PENALI.NE.0D0) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x TEMPERATURE
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 111, 112, 113, 114, 115, 116 ) , K
 111              K1 = 1
                  K2 = 2
                  GOTO 118
 112              K1 = 2
                  K2 = 3
                  GOTO 118
 113              K1 = 3
                  K2 = 1
                  GOTO 118
 114              K1 = 1
                  K2 = 4
                  GOTO 118
 115              K1 = 2
                  K2 = 4
                  GOTO 118
 116              K1 = 3
                  K2 = 4
C
C                 LA LONGUEUR DE L'ARETE K / 2
 118              DELTAK = SQRT( (X(K2,1)-X(K1,1)) ** 2
     %                         + (X(K2,2)-X(K1,2)) ** 2 ) * 0.5D0
C
                  DO M=1,2
C
C                    NUMERO ET COORDONNEES DU SOMMET M DE L'ARETE K
                     KM = K12(M)
                     XYZ(1) = X(KM,1)
                     XYZ(2) = X(KM,2)
                     XYZ(3) = X(KM,3)
C
                     MN = LTDELI(LPSOUR,JEU,NOOB)
                     IF( MN .GT. 0 ) THEN
C
C                       FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
                        MN = LTDELI(LPSOUR,JEU,NOOB)
C                       ATTENTION: ICI OndeR et OndeI SONT PASSES
C                                  A LA PLACE DE VN le VECTEUR NORMAL!
                        TEMPEL = DMCN( MNON  +KM )
                        ONDEPI = DMCN( MNON+4+KM )
                        CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                           TEMPEL, ONDEPI, 0D0,
     %                               MN, FORCE(1) )
                     ENDIF
C
                     MN = LTDELI(LPCONT,JEU,NOOB)
                     IF( PENALI .GT. 0D0 .AND. MN .GT. 0 ) THEN
C                       FIXATION(2) PENALISEE
                        CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                               NBCOFI, FIXA )
C                       NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                        DO I = 1, NBCOFI
C                          LE NUMERO DE LA COMPOSANTE FIXEE
                           NU = MCN( MN + WUCOFI - 1 + I )
                           FORCE(NU) = FIXA(I) * PENALI
                        ENDDO
                     ENDIF
C
C                    LE COEFFICIENT DU SECOND MEMBRE
                     IF( TESTNL .EQ. 6 ) THEN
                        BE(KM,1) = BE(KM,1) + DELTAK * FORCE(1)
                        BE(KM,2) = BE(KM,2) + DELTAK * FORCE(2)
                     ELSE IF( TESTNL .EQ. 7 ) THEN
                        D = DeltaT * DELTAK
                        BE(KM,1) = BE(KM,1) - D * FORCE(1)
                        BE(KM,2) = BE(KM,2) + D * FORCE(2)
                     ENDIF
C
                  ENDDO
C
               ENDIF
            ENDIF
 130     CONTINUE
C
C        =========================================================
C        CONTRIBUTION DES FORCES OU FIXATION PENALISEE AUX SOMMETS
C        =========================================================
         DO K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT
C              D'UNE SOURCE OU D'UN CONTACT PENALISE?
               IF(LTDEPO(LPSOUR,JEU,NOOB).GT.0 .OR.
     %            LTDEPO(LPCONT,JEU,NOOB).GT.0 .AND. PENALI.NE.0D0) THEN
C
C                 OUI: UN TABLEAU FIXATION PENALISE EXISTE POUR CE POINT
C                 FIXATION PENALISE = PENALI x VALEUR
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = X(K,3)
C
                  MN = LTDEPO(LPSOUR,JEU,NOOB)
                  IF( MN .GT. 0 ) THEN
C
C                    FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
C                    ATTENTION: ICI OndeR et OndeI SONT AUSSI PASSES
C                               A LA PLACE DE VN le VECTEUR NORMAL!
C                    RECUPERATION DE LA TEMPERATURE AU SOMMET K
                     TEMPEL = DMCN( MNON  +K )
                     ONDEPI = DMCN( MNON+4+K )
                     CALL REFORC( 1, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                        TEMPEL, ONDEPI, 0D0,
     %                            MN, FORCE(1) )
                  ENDIF
C
                  MN = LTDEPO(LPCONT,JEU,NOOB)
                  IF( PENALI .GT. 0D0 .AND. MN .GT. 0 ) THEN
C                    FIXATION(2) PENALISEE
                     CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                            NBCOFI, FIXA )
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                     DO I = 1, NBCOFI
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        NU = MCN( MN + WUCOFI - 1 + I )
                        FORCE(NU) = FIXA(I) * PENALI
                     ENDDO
                  ENDIF
C
C                 LE COEFFICIENT DU SECOND MEMBRE
                  IF( TESTNL .EQ. 6 ) THEN
                     BE(K,1) = BE(K,1) + FORCE(1)
                     BE(K,2) = BE(K,2) + FORCE(2)
                  ELSE IF( TESTNL .EQ. 7 ) THEN
                     BE(K,1) = BE(K,1) - DeltaT * FORCE(1)
                     BE(K,2) = BE(K,2) + DeltaT * FORCE(2)
                  ENDIF
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
