      SUBROUTINE NS2P1D( DeltaT, D2PI, NOAXIS, X, PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  NLSE: CALCUL DU SECOND MEMBRE DES EF 2P1D AXISYMETRIQUES OU 2D
C -----        LAGRANGE ISOPARAMETRIQUES EN SCHEMA IMPLICITE
C              c'est a dire
C TESTNL=6 REAL PART:
C Rho/dt (Wn+1m-Wn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Vn+1m  +Fr
C TESTNL=6 IMAG PART:
C Rho/dt (-Vn+1m+Vn) -Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2) Wn+1m+1+Fi
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
C X      : COORDONNEES RAYON ET COTE DES 3 SOMMETS DE L'EF
C          OU X Y DES 3 SOMMETS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NOOBPS : NUMERO DE SOMMET DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C SORTIE :
C --------
C BE     : BE(3,2) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR     MARS 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2011
C23456---------------------------------------------------------------012
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
      REAL              X(3, 2)
      INTEGER           NOOBPS( 1:3 )
      INTEGER           NOOBLA( 1:3 )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
C
      DOUBLE PRECISION  XYZ(3), BE(3,2), PENALI, D2PI, VN(2)
      DOUBLE PRECISION  D, S, DELTA, MASSE, BETA, FORCE(2),
     %                  FIXA(2), DeltaT, NLSECOEF, NLSECOEF0
      INTEGER           K12(2)
      EQUIVALENCE      (K12(1),K1), (K12(2),K2)
C
C     RECUPERATION DE L'ONDE AUX 3 DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     L'ADRESSE MCN DE VECTEUR EST MNTHET de $MEFISTO/incl/cthet.inc
C     ------------------------------------------------------------------
      CALL NLDATA0( 3 )
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C     EN TOUS LES DL DE L'EF
      MOREE2 =  MOTVAR(6)
      MNON   = (MNTHDL -1) / MOREE2
      MNONn  = (MNTHDLn-1) / MOREE2
      MNON0  = (MNTHDL0-1) / MOREE2
C
C     ===================================
C     CONTRIBUTION DES FORCES SURFACIQUES
C     ===================================
C     JACOBIEN * POIDS
      DELTA = ABS(  (X(2,1) - X(1,1)) * (X(3,2) - X(1,2))
     %            - (X(3,1) - X(1,1)) * (X(2,2) - X(1,2))  ) / 6D0
C
      DO L = 1, 3
C        COORDONNEES DU POINT D'INTEGRATION L = SOMMET L
         XYZ(1) = X(L,1)
         XYZ(2) = X(L,2)
         XYZ(3) = 0D0
C
C        LA VALEUR DE L'ONDE INTERPOLEE AU POINT D'INTEGRATION L
         TEMPEL  = DMCN( MNON     + L )
         ONDEPI  = DMCN( MNON + 3 + L )
         TEMPELn = DMCN( MNONn    + L )
         ONDEPIn = DMCN( MNONn+ 3 + L )
         TEMPEL0 = DMCN( MNON0    + L )
         ONDEPI0 = DMCN( MNON0+ 3 + L )
C
C        DENSITE DE MASSE
C        ----------------
         MASSE = 1D0
         MN    = LTDESU(LPMAST,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REMASS( 3, NOOBSF, 3, XYZ, MN,  MASSE )
         ENDIF
C
C        BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS ACTUEL
C       -BETA COEFFICIENT DU TERME NON LINEAIRE AU TEMPS INITIAL
C        ( ex: BETA (V**2+W**2) - BETA (V0**2+W0**2) )
C        -------------------------------------------------------
         BETA = 0D0
         MN   = LTDESU(LPCOET,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPS, TEMPEL, ONDEPI, MN,
     %                   NLSECOEF )
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPSINI,TEMPEL0,ONDEPI0,MN,
     %                   NLSECOEF0 )
            BETA = NLSECOEF - NLSECOEF0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDESU(LPSOUR,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                   TEMPEL, ONDEPI, 0D0, MN,  FORCE )
         ELSE
            FORCE(1) = 0D0
            FORCE(2) = 0D0
         ENDIF
C
C        LE COEFFICIENT DELTA * ...
C        --------------------------
         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            S = DELTA * D2PI * X(L,1)
         ELSE
C           EF NON AXISYMETRIQUE
            S = DELTA
         ENDIF
C
C        CALCUL DU VECTEUR * S
         IF( TESTNL .EQ. 6 ) THEN
C
C        REAL PART=Rho/dt (Wn+1m-Wn)-Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2)Vn+1m+F
         BE(L,1)=S*(MASSE/DeltaT*(ONDEPI-ONDEPIn)-BETA*TEMPEL+FORCE(1))
C
C        IMAG PART=Rho/dt(-Vn+1m+Vn)-Beta (Vn+1m**2+Wn+1m**2-V0**2-W0**2)Wn+1m+F
         BE(L,2)=S*(MASSE/DeltaT*(TEMPELn-TEMPEL)-BETA*ONDEPI+FORCE(2))
C
         ELSE IF( TESTNL .EQ. 7 ) THEN
C
C        PARTIE REELLE     = DeltaT [N(Beta)] V  + [M(Dmasse)]      Wn - DeltaT
C        PARTIE IMAGINAIRE = [M(Dmasse)]      Vn - DeltaT [N(Beta)] W  + DeltaT
C        -----------------------------------------------------------------------
C        LE TERME  DeltaT [N(Beta)] V + [M(Dmasse)] Wn - DeltaT FOmegaR
         BE(L,1)=S*( DeltaT*( BETA*TEMPEL - FORCE(1) ) + MASSE*ONDEPIn )
C
C        LE TERME [M(Dmasse)] Vn - DeltaT [N(Beta)] W + DeltaT FOmegaI
         BE(L,2)=S*( MASSE*TEMPELn + DeltaT*( FORCE(2) - BETA*ONDEPI ) )
C
         ENDIF
C
      ENDDO
C
C     =============================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES 3 COTES
C     =============================================================
      MNON = (MNTHDL-1)/MOREE2
      DO 50 K = 1, 3
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
            K1 = K
            IF( K1 .NE. 3 ) THEN
               K2 = K1 + 1
            ELSE
               K2 = 1
            ENDIF
C
C           LE VECTEUR ORTHOGONAL A L'ARETE K
            VN(1) = X(K2,2) - X(K1,2)
            VN(2) = X(K1,1) - X(K2,1)
C
C           LE JACOBIEN
            D = SQRT( VN(1)**2 + VN(2)**2 )
C
C           LE VECTEUR NORMAL UNITAIRE
            VN(1) = VN(1) / D
            VN(2) = VN(2) / D
C
C           LA LONGUEUR DE L'ARETE K / 2
            DELTA = SQRT( (X(K2,1)-X(K,1)) ** 2
     %                  + (X(K2,2)-X(K,2)) ** 2 ) * 0.5D0
C
C           RECUPERATION DE L'ONDE AUX 2 SOMMETS DE L'ARETE K
            DO L = 1, 2
C
C              NO DE 1 A 3 DU SOMMET L DE L'ARETE K
               NSLK = K12(L)
C
C              CALCUL DE L'ONDE AU POINT D'INTEGRATION L DE L'ARETE K
               TEMPEL = DMCN( MNON+NSLK )
               ONDEPI = DMCN( MNON+NSLK+3 )
C
C              3 COORDONNEES DU POINT D'INTEGRATION L
               XYZ(1) = X(NSLK,1)
               XYZ(2) = X(NSLK,2)
               XYZ(3) = 0D0
C
               IF( IECHAN .EQ. 1 ) THEN
C
C                 FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
                  MN = LTDELI(LPSOUR,JEU,NOOB)
C                 ATTENTION: ICI OndeR et OndeI SONT PASSES
C                            A LA PLACE DE VN le VECTEUR NORMAL!
                  CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                     TEMPEL, ONDEPI, 0D0,
     %                         MN, FORCE )
C
               ELSE
C
C                 FIXATION(2) PENALISEE POUR CONDITION NEUMANN OU FOURIER
                  MN = LTDELI(LPCONT,JEU,NOOB)
                  CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1) = 0D0
                  FORCE(2) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU) = FIXA(I) * PENALI
                  ENDDO
C
               ENDIF
C
               IF( NOAXIS .EQ. 0 ) THEN
C                 EF NON AXISYMETRIQUE
                  D = DELTA
               ELSE
C                 EF AXISYMETRIQUE
                  D = DELTA * D2PI * X(NSLK,1)
               ENDIF
C
               IF( TESTNL .EQ. 6 ) THEN
C
                  BE(NSLK,1) = BE(NSLK,1) + D * FORCE(1)
                  BE(NSLK,2) = BE(NSLK,2) + D * FORCE(2)
C
               ELSE IF( TESTNL .EQ. 7 ) THEN
C
                  D = D * DeltaT
                  BE(NSLK,1) = BE(NSLK,1) - D * FORCE(1)
                  BE(NSLK,2) = BE(NSLK,2) + D * FORCE(2)
C
               ENDIF
            ENDDO
         ENDIF
 50   CONTINUE
C
C     =========================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE AUX SOMMETS
C     =========================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO K = 1, 3
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FORCE?
               MNLT = LTDEPO(LPSOUR,JEU,NOOB)
               IF( MNLT .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU SOURCE PENALISE EXISTE POUR CE POINT
C                 FIXATION PENALISE = PENALI x VALEUR
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = 0D0
C
C                 FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
C                 ATTENTION: ICI OndeR et OndeI SONT AUSSI PASSES
C                            A LA PLACE DE VN le VECTEUR NORMAL!
C                 RECUPERATION DE L'ONDE AU SOMMET K
                  TEMPEL = DMCN( MNON  +K )
                  ONDEPI = DMCN( MNON+4+K )
                  CALL REFORC( 1, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                         TEMPEL, ONDEPI, 0D0, MNLT,  FORCE(1) )
C
               ENDIF
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE
C              FIXATION PENALISEE?
               MNLT = LTDEPO(LPCONT,JEU,NOOB)
               IF( MNLT .GT. 0 .AND. PENALI .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x FIXATION DE L'ONDE
                  TEMPEL = DMCN( MNON  +K )
                  ONDEPI = DMCN( MNON+4+K )
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = 0D0
C                 FIXATION(2) PENALISEE
                  CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MNLT,
     %                         NBCOFI, FIXA )
                  FORCE(1) = 0D0
                  FORCE(2) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MNLT + WUCOFI - 1 + I )
                     FORCE(NU) = FIXA(I) * PENALI
                  ENDDO
C
C                 SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
                  IF(NOAXIS .NE. 0) THEN
                     FORCE(1) = FORCE(1) * D2PI * X(K,1)
                     FORCE(2) = FORCE(2) * D2PI * X(K,1)
                  ENDIF
C
C                 LE COEFFICIENT DU SECOND MEMBRE EST IMPOSE
                  IF( TESTNL .EQ. 6 ) THEN
C
                     BE(K,1) = BE(K,1) + FORCE(1)
                     BE(K,2) = BE(K,2) + FORCE(2)
C
                  ELSE IF( TESTNL .EQ. 7 ) THEN
C
                     BE(K,1) = BE(K,1) - FORCE(1) * DeltaT
                     BE(K,2) = BE(K,2) + FORCE(2) * DeltaT
C
                  ENDIF
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
