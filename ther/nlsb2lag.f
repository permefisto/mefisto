      SUBROUTINE NLSB2LAG( NUELEM, NONOEF, Omega,
     %                     DeltaT, D2PI,   NOAXIS,
     %                     X,      PENALI, NBJEUX, JEU,
     %                     NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                     NBPOLA, NPIA,   POIDSA, POLYA,  DPOLYA,
     %                     NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                     NBPOLY, NPI,    POLY,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     F1,     F2,     POIDEL, DP,
     %                     NBNOMA, Utn,    Utm,
     %                     BE,     IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: TESTNL=6 SCHEMA SEMI-IMPLICITE CALCUL DU SECOND MEMBRE BE
C ----- SUR UN EF LAGRANGE 2D
C     REAL PART:
C      Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn)
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m-Wn+1m)
C     -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
C     +Alfa LAPLACIEN Wn+1m
C
C     IMAG PART:
C     -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn)
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m+Wn+1m)
C     +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
C     -Alfa LAPLACIEN Vn+1m

C ENTREES:
C --------
C NUELEM : NUMERO DE L'EF TRAITE
C NONOEF : NUMERO DES NBPOLY NOEUDS DE L'EF NUELEM
C Omega  : VITESSE ANGULAIRE DE LA ROTATION
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
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C CONDDP : TENSEUR CONDUCTIVITE DE -ALFA LAPLACIEN DE L'EF 2LAG
C          AUX NPI POINTS D'INTEGRATION NUMERIQUE

C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE
C Utn    : U(tn) (NBNOMA,2) U A L'INTANT tn
C Utm    : U(tn+1,m) (NBNOMA,2)  U A L'INTANT tn+1 iteration m
C
C SORTIE :
C --------
C BE     : BE(NBPOLY,2) LE SECOND MEMBRE ELEMENTAIRE
C IERR   : =7 SI PB AXISYMETRIQUE AVEC X=RAYON NEGATIVE OU NUL
C          INCHANGE SI PAS D'ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2014
C23456---------------------------------------------------------------012
      PARAMETER        ( MAXPTI=9, MAXPTA=3 )
C     MAXIMUM POUR LES TYPES D'EF 2D DES POLYNOMES,
C                                    POINTS D'INTEGRATION DE LA SURFACE
C                                    POINTS D'INTEGRATION D'UNE ARETE
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
      include"./incl/gsmenu.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  D2PI, DeltaT, Omega(3), PENALI
ccc   DOUBLE PRECISION  VN(2)
      DOUBLE PRECISION  POIDSA(NPIA), POLYA(NBPOLA, NPIA),
     %                  DPOLYA(NBPOLA, NPIA),
     %                  POLY(NBPOLY, NPI),
     %                  DP(2,NBPOLY,NPI),
     %                  POIDEL(NPI),
     %                  F1(NPI), F2(NPI), XYZ(3),
     %                  BE(NBPOLY,2)
      REAL              X(NBPOLY, 2)
      INTEGER           NONOEF( 1:NBPOLY )
      INTEGER           NOOBPS( 1:NBSOMT )
      INTEGER           NOOBLA( 1:NBCOTE )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )

      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2)
      DOUBLE PRECISION  S, SV, SW, D,  GL(2), DGL(2), PROSCD, DELTA
      DOUBLE PRECISION  Rho(MAXPTI),
     %                  NLSECOEF0(MAXPTI), NLSECOEF(MAXPTI),
     %                  FORCE(2,MAXPTI),
     %                  Vtm(MAXPTI), Wtm(MAXPTI),
     %                  Vtn(MAXPTI), Wtn(MAXPTI),
     %                  F1dYF2dXV(MAXPTI), F1dYF2dXW(MAXPTI),
     %                  VEFtm(8),WEFtm(8), VEFtn(8),WEFtn(8),
     %                  UK(MAXPTA,2), FIXA(2),
     %                  COND(6), CONDP(2,2,MAXPTI)
      INTEGER           NOPOAR(3)

C     INITIALISATION DE LA VITESSE ANGULAIRE Omega/Z SUPPOSEE CONSTANTE
C     -----------------------------------------------------------------
      IF( NUELEM .EQ. 1 ) THEN
         MN = LTDESU(LPVIANT,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
C            RECUPERATION DE LA VITESSE ANGULAIRE Omega SUPPOSEE CONSTANTE
C            VECTEUR(3) DE VITESSE ANGULAIRE AU 1-ER POINT D'INTEGRATION
             CALL REVIAN( 3, NOOBSF, F1(1), F2(1), 0D0,
     %                    LTDESU(LPVIANT,JEU,NOOBSF), Omega )
          ELSE
C            PAS DE ROTATION
             Omega(1) = 0D0
             Omega(2) = 0D0
             Omega(3) = 0D0
          ENDIF
      ENDIF
C
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL tn+1m et tn
C     AUX NBPOLY DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     ------------------------------------------------------------------
      DO I=1,NBPOLY
C        NO GLOBAL DU NOEUD I DE L'EF
         J = NONOEF(I)
C        AU TEMPS tn
         VEFtn(I) = Utn(J,1)
         WEFtn(I) = Utn(J,2)
C        AU TEMPS tn+1m
         VEFtm(I) = Utm(J,1)
         WEFtm(I) = Utm(J,2)
      ENDDO
C
C     ===================================
C     CONTRIBUTION DES FORCES SURFACIQUES
C     ===================================
      DO L = 1, NPI

C        COORDONNEES DU POINT D'INTEGRATION L
         XYZ(1) = F1(L)
         XYZ(2) = F2(L)
         XYZ(3) = 0D0

C        LA VALEUR INTERPOLEE DES ONDES AU POINT D'INTEGRATION bl
C        CALCUL PARTIE REELLE     a tn+1m AU POINT D'INTEGRATION bl de L'ONDE
         Vtm(L) = PROSCD( POLY(1,L), VEFtm, NBPOLY )
C        CALCUL PARTIE IMAGINAIRE a tn+1m AU POINT D'INTEGRATION bl de L'ONDE
         Wtm(L) = PROSCD( POLY(1,L), WEFtm, NBPOLY )
C        CALCUL PARTIE REELLE a tn AU POINT D'INTEGRATION bl de L'ONDE
         Vtn(L) = PROSCD( POLY(1,L), VEFtn, NBPOLY )
C        CALCUL PARTIE REELLE a tn AU POINT D'INTEGRATION bl de L'ONDE
         Wtn(L) = PROSCD( POLY(1,L), WEFtn, NBPOLY )
C
C        DENSITE DE MASSE
C        ----------------
         MN = LTDESU(LPMAST,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            TEMPEL = Vtm(L)
            ONDEPI = Wtm(L)
            CALL REMASS( 3, NOOBSF, 3, XYZ, MN, Rho(L) )
         ELSE
            Rho(L) = 1D0
         ENDIF
C
C        NLSECOEF COEFFICIENT DU TERME NON LINEAIRE AU TEMPS DEMANDE
C        -----------------------------------------------------------
         MN = LTDESU(LPCOET,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            TEMPEL = Vtm(L)
            ONDEPI = Wtm(L)
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPS, Vtm(L), Wtm(L), MN,
     %                   NLSECOEF(L) )

            TEMPEL0 = Vtn(L)
            ONDEPI0 = Wtn(L)
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPSINI, Vtn(L), Wtn(L),MN,
     %                   NLSECOEF0(L) )
         ELSE
            NLSECOEF0(L)= 0D0
            NLSECOEF(L) = 0D0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDESU(LPSOUR,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                 Vtm(L), Wtm(L), 0D0,
     %                   MN, FORCE(1,L) )
         ELSE
            FORCE(1,L) = 0D0
            FORCE(2,L) = 0D0
         ENDIF

C        CALCUL de ( x d/dy - y d/dx ) Vn+1m(bl)
C        CALCUL de ( x d/dy - y d/dx ) Wn+1m(bl)
C        ---------------------------------------
         SV = 0D0
         SW = 0D0
         IF( Omega(1) .NE. 0D0 ) THEN
            DO J=1,NBPOLY
               D  = F1(L) * DP(2,J,L) - F2(L) * DP(1,J,L)
               SV = SV + D * VEFtm(J)
               SW = SW + D * WEFtm(J)
            ENDDO
         ENDIF
C
C        COEFFICIENT ALFA DE - ALFA LAPLACIEN C'EST A DIRE CALCUL DU
C        TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT D'INTEGRATION L
C        si Gross-Pitaevski: COND = 1/2 [Identite]
C        -----------------------------------------------------------
         MN = LTDESU(LPCOND,JEU,NOOBSF)
         TEMPEL = Vtm(L)
         ONDEPI = Wtm(L)
         CALL RECOND( 3, NOOBSF, 3, XYZ, MN, COND )

C        POIDS DE LA FORMULE D'INTEGRATION NUMERIQUE
C        -------------------------------------------
         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            IF( F1(L) .LE. 0D0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'ERREUR: PB AXISYMETRIQUE AVEC EF'
               KERR(2) = 'D''ABSCISSE X=RAYON NEGATIVE OU NULLE'
               CALL LEREUR
               IERR = 7
               RETURN
            ENDIF
            S = POIDEL(L) * D2PI * F1(L)
         ELSE
C           EF NON AXISYMETRIQUE
            S = POIDEL(L)
         ENDIF
C
C        COEF = COEF * POIDS(l) * DELTA(bl)
         Rho(L)       = Rho(L)       * S / DeltaT
         NLSECOEF0(L) = NLSECOEF0(L) * S
         NLSECOEF(L)  = NLSECOEF(L)  * S
         FORCE(1,L)   = FORCE(1,L)   * S
         FORCE(2,L)   = FORCE(2,L)   * S

C        OmegaZ ( x d/dy - y d/dx )Vn+1m(bl)
         F1dYF2dXV(L) = Omega(1) * SV * S
C        OmegaZ ( x d/dy - y d/dx )Wn+1m(bl)
         F1dYF2dXW(L) = Omega(1) * SW * S

C        COND = CONDUCTIVITE * POIDS(l) * DELTA(bl)
         CONDP(1,1,L) = COND(1) * S
         CONDP(2,1,L) = COND(2) * S
         CONDP(1,2,L) = COND(2) * S
         CONDP(2,2,L) = COND(3) * S
C
      ENDDO
C
C     CALCUL DU VECTEUR ELEMENTAIRE BE
C     --------------------------------
      DO I = 1, NBPOLY

C        LE TERME -Alfa LAPLACIEN V ou W INITIALISE BE(I,1:2)
         SV = 0D0
         SW = 0D0
         DO J=1,NBPOLY

            S = 0D0
            DO L = 1, NPI
C              CONDUC = CONDUC + T(DP) * CONDP * (DP) (bl)
               DO K=1,2
                  S = S + DP(K,I,L) * ( CONDP(K,1,L) * DP(1,J,L)
     %                                + CONDP(K,2,L) * DP(2,J,L) )
               ENDDO
            ENDDO

C           -Alfa LAPLACIEN (I,J) Vn+1m (J)
            SV = SV + S * VEFtm(J)
C           +Alfa LAPLACIEN (I,J) Wn+1m (J)
            SW = SW - S * WEFtm(J)

         ENDDO
C        +Alfa LAPLACIEN Wn+1m
         BE(I,1) = SW
C        -Alfa LAPLACIEN Vn+1m
         BE(I,2) = SV

         DO L = 1, NPI

C           REAL PART:
C          +Alfa LAPLACIEN Wn+1m
C           Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m -Wn+1m)
C          -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
            BE(I,1) = BE(I,1) + POLY(I,L)
     %              * ( FORCE(1,L) + FORCE(2,L)
     %                + Rho(L) * ( Vtn(L) + Wtm(L) - Wtn(L) )
     %                + F1dYF2dXV(L) - F1dYF2dXW(L)
     %                - NLSECOEF(L)  * ( Vtm(L) + Wtm(L) )
     %                + NLSECOEF0(L) * Vtm(L)  )

C           IMAG PART:
C          -Alfa LAPLACIEN Vn+1m
C          -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m +Wn+1m)
C          +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
            BE(I,2) = BE(I,2) + POLY(I,L)
     %              * ( FORCE(2,L) - FORCE(1,L)
     %                + Rho(L) * ( -Vtm(L) + Vtn(L) + Wtn(L) )
     %                + F1dYF2dXV(L) + F1dYF2dXW(L)
     %                + NLSECOEF(L)  * ( Vtm(L) - Wtm(L) )
     %                + NLSECOEF0(L) * Wtm(L)  )

         ENDDO

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
C           RECUPERATION DE L'ONDE UK AUX NBPOLA DL DE L'ARETE K
            DO I = 1, NBPOLA
               UK( I, 1 ) = VEFtm( NOPOAR(I) )
               UK( I, 2 ) = WEFtm( NOPOAR(I) )
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
               Vtm(L) = PROSCD( POLYA(1,L), UK(1,1), NBPOLA )
               Wtm(L) = PROSCD( POLYA(1,L), UK(1,2), NBPOLA )
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
     %                                     Vtm(L), Wtm(L), 0D0,
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
                     N = MCN( MN + WUCOFI - 1 + I )
                     FORCE(N,L) = FIXA(I) * PENALI
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
               DO I=1,NBPOLA
C                 NO DU NOEUD DANS L'EF
                  N = NOPOAR(I)
                  BE(N,1) = BE(N,1) + D * POLYA(I,L) * FORCE(1,L)
                  BE(N,2) = BE(N,2) + D * POLYA(I,L) * FORCE(2,L)
               ENDDO
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
                  TEMPEL = VEFtm( K )
                  ONDEPI = WEFtm( K )
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
                     N = MCN( MN + WUCOFI - 1 + I )
                     FORCE(N,L) = FIXA(I) * PENALI
                  ENDDO
C
C                 SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
                  IF( NOAXIS .NE. 0 ) THEN
                     FORCE(1,L) = FORCE(1,L) * D2PI * X(K,1)
                     FORCE(2,L) = FORCE(2,L) * D2PI * X(K,1)
                  ENDIF
C
C                 LE COEFFICIENT DU SECOND MEMBRE EST IMPOSE
                  BE(K,1) = BE(K,1) + FORCE(1,L)
                  BE(K,2) = BE(K,2) + FORCE(2,L)
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
