      SUBROUTINE NS2LAG( DeltaT, D2PI,   NOAXIS, X, PENALI, NBJEUX, JEU,
     %                   NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NBPOLA, NPIA,   POIDSA, POLYA,  DPOLYA,
     %                   NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NBPOLY, NPI,    POLY,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   F1,     F2,     POIDEL,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  NLSE: CALCUL DU SECOND MEMBRE DES EF AXISYMETRIQUES OU 2D
C -----        LAGRANGE ISOPARAMETRIQUES EN SCHEMA IMPLICITE EN TEMPS
C              c'est a dire
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
C D2PI   : 2 FOIS PI
C NOAXIS : 1 SI PROBLEME AXISYMETRIQUE,  0 SINON
C X      : COORDONNEES RAYON ET COTE DES NBPOLY POINTS DE L'EF
C          OU X Y DES NBPOLY POINTS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NBSOMT : NOMBRE DE SOMMETS DE L'EF
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NBPOLA : NOMBRE DE POLYNOMES DE BASE SUR UN COTE DE L'EF
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR UN COTE
C POIDSA : POIDS DES POINTS D INTEGRATION SUR UN COTE
C POLYA  : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION COTE
C DPOLYA : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBCOTE : NOMBRE DES COTES DE L'EF
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES
C
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE
C POLY   : POLY(I, L) = PI(RL, ZL)
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C F1     : RAYON R DES NPI POINTS D'INTEGRATION
C F2     : COTE  Z DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C X      : COORDONNEES RAYON ET COTE DES NBPOLY POINTS DE L'EF
C
C SORTIE :
C --------
C BE     : BE(NBPOLY,2) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET TEXAS A & M University at QATAR      Mars 2011
C AJOUT  :ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Novembre 2013
C23456---------------------------------------------------------------012
      PARAMETER        ( MAXPTI=9, MAXPTA=3 )
C     MAXIMUM POUR LES TYPES D'EF 2D DES POLYNOMES,
C                                    POINTS D'INTEGRATION DE LA SURFACE
C                                    POINTS D'INTEGRATION D'UNE ARETE
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  D2PI, DeltaT
ccc   DOUBLE PRECISION  VN(2)
      DOUBLE PRECISION  POIDSA(NPIA), POLYA(NBPOLA, NPIA),
     %                  DPOLYA(NBPOLA, NPIA),
     %                  POLY(NBPOLY, NPI),
     %                  POIDEL(NPI),
     %                  F1(NPI), F2(NPI), XYZ(3),
     %                  BE(NBPOLY,2),
     %                  PENALI
      REAL              X(NBPOLY, 2)
      INTEGER           NOOBPS( 1:NBSOMT )
      INTEGER           NOOBLA( 1:NBCOTE )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
C
      INTEGER           NOPOAR(3)
      DOUBLE PRECISION  S,  D,  GL(2), DGL(2), PROSCD, DELTA
      DOUBLE PRECISION  MASSE(MAXPTI), BETA(MAXPTI), FORCE(2,MAXPTI),
     %                  V(MAXPTI), W(MAXPTI), Vtn(MAXPTI), Wtn(MAXPTI),
     %                  UK(MAXPTA,2), FIXA(2), NLSECOEF, NLSECOEF0
C
C     RECUPERATION DE L'ONDE AUX NBPOLY DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     L'ADRESSE MCN DE VECTEUR EST MNTHET de $MEFISTO/incl/cthet.inc
C     ------------------------------------------------------------------
      CALL NLDATA0( NBPOLY )
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C     EN TOUS LES DL DE L'EF
      MOREE2 =  MOTVAR(6)
      MNON   = (MNTHDL -1) / MOREE2
C
C     ===================================
C     CONTRIBUTION DES FORCES SURFACIQUES
C     ===================================
      DO L = 1, NPI
C
C        COORDONNEES DU POINT D'INTEGRATION L
         XYZ(1) = F1(L)
         XYZ(2) = F2(L)
         XYZ(3) = 0D0
C
C        LA VALEUR DE L'ONDE INTERPOLEE AU POINT D'INTEGRATION L
         CALL NLDATA1( NBPOLY, POLY(1,L) )
C        PB NON LINEAIRE NLSE:
C        TEMPEL:  PARTIE REELLE     tn+1m AU POINT D'INTEGRATION L de L'ONDE
C        ONDEPI:  PARTIE IMAGINAIRE tn+1m AU POINT D'INTEGRATION L de L'ONDE
C        TEMPELn: PARTIE REELLE     tn    AU POINT D'INTEGRATION L de L'ONDE
C        ONDEPIn: PARTIE IMAGINAIRE tn    AU POINT D'INTEGRATION L de L'ONDE
C        TEMPEL0: PARTIE REELLE     t0    AU POINT D'INTEGRATION L de L'ONDE
C        ONDEPI0: PARTIE IMAGINAIRE t0    AU POINT D'INTEGRATION L de L'ONDE
         V(L)   = TEMPEL
         W(L)   = ONDEPI
         Vtn(L) = TEMPELn
         Wtn(L) = ONDEPIn
C
C        DENSITE DE MASSE
C        ----------------
         MASSE(L) = 1D0
         MN = LTDESU(LPMAST,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REMASS( 3, NOOBSF, 3, XYZ, MN, MASSE(L) )
         ENDIF
C
C        BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS ACTUEL
C       -BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS INITIAL
C        ( ex: BETA (V**2+W**2) - BETA (V0**2+W0**2) )
C        ------------------------------------------------------
         BETA(L) = 0D0
         MN      = LTDESU(LPCOET,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPS, Vtn(L), Wtn(L), MN,
     %                   NLSECOEF )
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPSINI, Vtn(L), Wtn(L),MN,
     %                   NLSECOEF0 )
            BETA(L) = NLSECOEF - NLSECOEF0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN =  LTDESU(LPSOUR,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                 TEMPEL, ONDEPI, 0D0,
     %                   MN, FORCE(1,L) )
         ELSE
            FORCE(1,L) = 0D0
            FORCE(2,L) = 0D0
         ENDIF
C
         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            S = POIDEL(L) * D2PI * F1(L)
         ELSE
C           EF NON AXISYMETRIQUE
            S = POIDEL(L)
         ENDIF
C
C        COEF = COEF * POIDS(l) * DELTA(Bl) (* DeltaT)
         MASSE(L)   = MASSE(L)   * S
         BETA(L)    = BETA(L)    * S
         FORCE(1,L) = FORCE(1,L) * S
         FORCE(2,L) = FORCE(2,L) * S
C
      ENDDO
C
C     CALCUL DU VECTEUR ELEMENTAIRE BE
C     --------------------------------
      DO I = 1, NBPOLY
C
         BE(I,1) = 0D0
         BE(I,2) = 0D0
C
         IF( TESTNL .EQ. 6 ) THEN
            DO L = 1, NPI
C
C              Rho/dt (Wn+1m-Wn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Vn+1m  +F
               BE(I,1) = BE(I,1) + POLY(I,L) *
     %         (MASSE(L)/DeltaT*(W(L)-Wtn(L)) -BETA(L)*V(L) +FORCE(1,L))
C
C              Rho/dt(-Vn+1m+Vn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Wn+1m+1+F
               BE(I,2) = BE(I,2) + POLY(I,L) *
     %         (MASSE(L)/DeltaT*(Vtn(L)-V(L)) -BETA(L)*W(L) +FORCE(2,L))
C
            ENDDO
C
         ELSE IF( TESTNL .EQ. 7 ) THEN
            DO L = 1, NPI
C
C              DeltaT [N(Beta)] V  + [M(Dmasse)] Wn - DeltaT FOmegaR
               BE(I,1) = BE(I,1) + POLY(I,L) *
     %          ( MASSE(L)*Wtn(L) + (BETA(L)*V(L) - FORCE(1,L))*Deltat )
C
C              [M(Dmasse)] Vn -  DeltaT [N(Beta)] W + DeltaT FOmegaI
               BE(I,2) = BE(I,2) + POLY(I,L) *
     %          ( MASSE(L)*Vtn(L) + (FORCE(2,L) - BETA(L)*W(L))*DeltaT )
C
            ENDDO

         ENDIF
C
      ENDDO
C
C     ===========================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES COTES
C     ===========================================================
      DO 50 K = 1, NBCOTE
C
C        LE NUMERO DE LIGNE DU COTE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. SUPPORT DE FLUX OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDELI(LPSOUR,JEU,NOOB) .GT. 0  ) IECHAN = 1
            IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                 ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 50
C
C           UN TABLEAU FORCE ou FIXATION EXISTE POUR CETTE LIGNE
C           CALCUL DE LA CONTRIBUTION DE L'ARETE K A BE
C           ....................................................
C           LE NUMERO DES POINTS DU COTE K
            NOPOAR(1) = K
            IF( K .NE. NBCOTE ) THEN
               NOPOAR(2) = K+1
            ELSE
               NOPOAR(2) = 1
            ENDIF
C           LE NUMERO DU POINT MILIEU
            NOPOAR(3) = K + NBCOTE
C
C           RECUPERATION DE L'ONDE AUX NBPOLA DL DE L'ARETE K
            DO I = 1, NBPOLA
               UK( I, 1 ) = DMCN( MNON+NOPOAR(I) )
               UK( I, 2 ) = DMCN( MNON+NOPOAR(I)+NBPOLY )
            ENDDO
C
            DO L = 1, NPIA
C
C              CALCUL DES COORDONNEES DU POINT D INTEGRATION L
C              ET DU JACOBIEN EN CE POINT
               CALL E22LAG ( NBPOLY, NBPOLA, NOPOAR,
     %                       POLYA(1,L), DPOLYA(1,L),
     %                       X, GL, DGL, DELTA )
C              EN SORTIE GL=LES 2 COORDONNEES DU POINT D'INTEGRATION
C
C              CALCUL DE L'ONDE AU POINT D'INTEGRATION L DE L'ARETE K
               TEMPEL = PROSCD( POLYA(1,L), UK(1,1), NBPOLA )
               ONDEPI = PROSCD( POLYA(1,L), UK(1,2), NBPOLA )
C
C              3 COORDONNEES DU POINT D'INTEGRATION L
               XYZ(1) = GL(1)
               XYZ(2) = GL(2)
               XYZ(3) = 0D0
C
               IF( IECHAN .EQ. 1 ) THEN
C
C                 FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
C                 LE VECTEUR NORMAL UNITAIRE
ccc                  VN(1) =  DGL(2) / DELTA
ccc                  VN(2) = -DGL(1) / DELTA
C                 ATTENTION: ICI OndeR et OndeI SONT PASSES
C                            A LA PLACE DE VN le VECTEUR NORMAL!
                  MN = LTDELI(LPSOUR,JEU,NOOB)
                  CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                     TEMPEL, ONDEPI, 0D0,
     %                         MN, FORCE(1,L) )
C
               ELSE
C
C                 FIXATION(2) PENALISEE POUR CONDITION DE DIRICHLET
                  MN = LTDELI(LPCONT,JEU,NOOB)
                  CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1,L) = 0D0
                  FORCE(2,L) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU,L) = FIXA(I) * PENALI
                  ENDDO
C
               ENDIF
C
               IF( NOAXIS .EQ. 0 ) THEN
C                 EF NON AXISYMETRIQUE
                  D = DELTA * POIDSA(L)
               ELSE
C                 EF AXISYMETRIQUE
                  D = DELTA * POIDSA(L) * D2PI * GL(1)
               ENDIF
C
               IF( TESTNL .EQ. 6 ) THEN
                  DO I=1,NBPOLA
C                    NO DU NOEUD DANS L'EF
                     NU = NOPOAR(I)
                     BE(NU,1) = BE(NU,1) + D * POLYA(I,L) * FORCE(1,L)
                     BE(NU,2) = BE(NU,2) + D * POLYA(I,L) * FORCE(2,L)
                  ENDDO
               ELSE IF( TESTNL .EQ. 7 ) THEN
                  D = D * DeltaT
                  DO I=1,NBPOLA
C                    NO DU NOEUD DANS L'EF
                     NU = NOPOAR(I)
                     BE(NU,1) = BE(NU,1) - D * POLYA(I,L) * FORCE(1,L)
                     BE(NU,2) = BE(NU,2) + D * POLYA(I,L) * FORCE(2,L)
                  ENDDO
               ENDIF
C
            ENDDO
         ENDIF
 50   CONTINUE
C
C     =========================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE AUX SOMMETS
C     =========================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO K = 1, NBSOMT
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
                MN = LTDEPO(LPCONT,JEU,NOOB)
                IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x ONDE
                  TEMPEL = DMCN( MNON          + K )
                  ONDEPI = DMCN( MNON + NBPOLY + K )
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = 0D0
C                 FIXATION(2) PENALISEE
                  CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1,L) = 0D0
                  FORCE(2,L) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU,L) = FIXA(I) * PENALI
                  ENDDO
C
C                 SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
                  IF( NOAXIS .NE. 0 ) THEN
                     FORCE(1,L) = FORCE(1,L) * D2PI * X(K,1)
                     FORCE(2,L) = FORCE(2,L) * D2PI * X(K,1)
                  ENDIF
C
C                 LE COEFFICIENT DU SECOND MEMBRE EST IMPOSE
                  IF( TESTNL .EQ. 6 ) THEN
                     BE(K,1) = BE(K,1) + FORCE(1,L)
                     BE(K,2) = BE(K,2) + FORCE(2,L)
                  ELSE IF( TESTNL .EQ. 7 ) THEN
                     BE(K,1) = BE(K,1) - FORCE(1,L) * DeltaT
                     BE(K,2) = BE(K,2) + FORCE(2,L) * DeltaT
                  ENDIF
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
