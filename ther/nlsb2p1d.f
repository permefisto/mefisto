      SUBROUTINE NLSB2P1D( NUELEM, NONOEF, Omega,
     %                     DeltaT, D2PI, NOAXIS, X, PENALI, NBJEUX, JEU,
     %                     NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                     NOOBLA, NUMILI, NUMALI, LTDELI,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     NBNOMA, Utn,    Utm, 
     %                     BE,     IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: TESTNL=6 SCHEMA SEMI-IMPLICITE CALCUL DU SECOND MEMBRE BE
C ----- SUR UN EF TRIANGLE P1 LAGRANGE 2D
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
C
C ENTREES:
C --------
C NUELEM : NUMERO DE L'EF TRAITE
C NONOEF : NUMERO DES 3 NOEUDS SOMMETS DE L'EF NUELEM
C Omega  : VITESSE ANGULAIRE DE LA ROTATION
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

C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE
C Utn    : U(tn) (NBNOMA,2) U A L'INTANT tn
C Utm    : U(tn+1,m) (NBNOMA,2)  U A L'INTANT tn+1 iteration m
C
C SORTIE :
C --------
C BE     : BE(3,2) LE SECOND MEMBRE ELEMENTAIRE
C IERR   : =7 SI PB AXISYMETRIQUE AVEC X=RAYON NEGATIVE OU NUL
C          INCHANGE SI PAS D'ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray     Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
      include"./incl/gsmenu.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      DOUBLE PRECISION  D2PI, DeltaT, Omega(3), PENALI
      REAL              X(3, 2)
      INTEGER           NONOEF( 1:3 )
      INTEGER           NOOBPS( 1:3 )
      INTEGER           NOOBLA( 1:3 )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )

      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2)
      DOUBLE PRECISION  S, D, DELTA, F1dYF2dXV, F1dYF2dXW
      DOUBLE PRECISION  SV, SW, NLSECOEF, NLSECOEF0
      DOUBLE PRECISION  Rho, FORCE(2), FIXA(2), DP(2,3),
     %                  VEFtm(3), WEFtm(3), VEFtn(3), WEFtn(3)
C
      DOUBLE PRECISION  XYZ(3), BE(3,2), VN(2)
      DOUBLE PRECISION  COND(6), CONDUC(6), C1, C2, C3
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32

      INTEGER           K12(2)
      EQUIVALENCE      (K12(1),K1), (K12(2),K2)

C     INITIALISATION DE LA VITESSE ANGULAIRE Omega/Z SUPPOSEE CONSTANTE
C     -----------------------------------------------------------------
      IF( NUELEM .EQ. 1 ) THEN
         MN = LTDESU(LPVIANT,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
C           RECUPERATION DE LA VITESSE ANGULAIRE Omega SUPPOSEE CONSTANTE
C           VECTEUR(3) DE VITESSE ANGULAIRE AU 1-ER SOMMET DU TRIANGLE
            XYZ(1) = X(1,1)
            XYZ(2) = X(1,2)
            CALL REVIAN( 3, NOOBSF, XYZ(1), XYZ(2), 0D0,
     %                   LTDESU(LPVIANT,JEU,NOOBSF), Omega )
C           Si SEUL OMEGA(1) EST DONNE C'EST SUPPOSE ETRE Omega(3)/axeZ
          ELSE
C            PAS DE ROTATION
             Omega(1) = 0D0
             Omega(2) = 0D0
             Omega(3) = 0D0
          ENDIF
      ENDIF
C
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL tn+1m et tn
C     AUX 3 SOMMETS DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     ------------------------------------------------------------------
      DO I = 1, 3
C        NO GLOBAL DU NOEUD I DE L'EF
         J = NONOEF(I)
C        AU TEMPS tn
         VEFtn(I) = Utn(J,1)
         WEFtn(I) = Utn(J,2)
C        AU TEMPS tn+1m
         VEFtm(I) = Utm(J,1)
         WEFtm(I) = Utm(J,2)
      ENDDO

C     ===================================
C     CONTRIBUTION DES FORCES SURFACIQUES
C     ===================================
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)

C     JACOBIEN
      DELTA = ABS( X21 * Y31 - X31 * Y21 )

C     DF-1 DP => DP
      DP(1,1) = -Y32
      DP(1,2) =  Y31
      DP(1,3) = -Y21

      DP(2,1) =  X32
      DP(2,2) = -X31
      DP(2,3) =  X21
C
      DO I = 1, 3

C        COORDONNEES DU POINT D'INTEGRATION I = SOMMET I
         XYZ(1) = X(I,1)
         XYZ(2) = X(I,2)
         XYZ(3) = 0D0

C        DENSITE DE MASSE / DeltaTemps
C        ----------------
         MN = LTDESU(LPMAST,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            TEMPEL = VEFtm(I)
            ONDEPI = WEFtm(I)
            CALL REMASS( 3, NOOBSF, 3, XYZ, MN, Rho )
         ELSE
            Rho = 1D0
         ENDIF
         Rho = Rho / DeltaT
C
C        NLSECOEF COEFFICIENT DU TERME NON LINEAIRE AU TEMPS DEMANDE
C        -----------------------------------------------------------
         MN = LTDESU(LPCOET,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            TEMPEL = VEFtm(I)
            ONDEPI = WEFtm(I)
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPS,   TEMPEL, ONDEPI, MN,
     %                   NLSECOEF )
            TEMPEL0 = VEFtn(I)
            ONDEPI0 = WEFtn(I)
            CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPSINI,TEMPEL0,ONDEPI0,MN,
     %                   NLSECOEF0 )
         ELSE
            NLSECOEF0 = 0D0
            NLSECOEF  = 0D0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDESU(LPSOUR,JEU,NOOBSF)
         IF( MN .GT. 0 ) THEN
            CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                 VEFtm(I), WEFtm(I), 0D0,
     %                   MN, FORCE )
         ELSE
            FORCE(1) = 0D0
            FORCE(2) = 0D0
         ENDIF

C        CALCUL de ( x d/dy - y d/dx ) Vn+1m(bl)
C        CALCUL de ( x d/dy - y d/dx ) Wn+1m(bl)
C        ---------------------------------------
         F1dYF2dXV = 0D0
         F1dYF2dXW = 0D0
         DO J = 1, 3
            D  = X(I,1) * DP(2,J) - X(I,2) * DP(1,J)
            F1dYF2dXV = F1dYF2dXV + D * VEFtm(J)
            F1dYF2dXW = F1dYF2dXW + D * WEFtm(J)
         ENDDO

C        OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )(bl)
         D = Omega(1) / DELTA
         F1dYF2dXV = F1dYF2dXV * D
C        OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )(bl)
         F1dYF2dXW = F1dYF2dXW * D

         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            IF( X(I,1) .LE. 0D0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'ERREUR: PB AXISYMETRIQUE AVEC EF'
               KERR(2) = 'D''ABSCISSE X=RAYON NEGATIVE OU NULLE'
               CALL LEREUR
               IERR = 7
               RETURN
            ENDIF
            S = DELTA * D2PI * X(I,1) / 6D0
         ELSE
C           EF NON AXISYMETRIQUE
            S = DELTA / 6D0
         ENDIF

C        COEFFICIENT ALFA DE - ALFA LAPLACIEN C'EST A DIRE CALCUL DU
C        TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT D'INTEGRATION L
C        si Gross-Pitaevski: COND = 1/2 [Identite]
C        -----------------------------------------------------------
         C1 = 0D0
         C2 = 0D0
         C3 = 0D0
C
         DO J=1,3
C
C           RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C           D'INTEGRATION J
            XYZ(1) = X(J,1)
            XYZ(2) = X(J,2)
            XYZ(3) = 0D0

C           RECHERCHE DE LA CONDUCTIVITE AU POINT D'INTEGRATION J
            TEMPEL = VEFtm(J)
            ONDEPI = WEFtm(J)
            CALL RECOND( 3, NOOBSF, 3, XYZ,
     %                   LTDESU(LPCOND,JEU,NOOBSF), COND )
C
C           LA SOMME DES CONDUCTIVITES AUX 3 SOMMETS
            C1 = C1 + COND(1)
            C2 = C2 + COND(2)
            C3 = C3 + COND(3)
C
         ENDDO
C
C        COND = CONDUCTIVITE * POIDS / DELTA
         D  = DELTA * 6D0
         C1 = C1 / D
         C2 = C2 / D
         C3 = C3 / D
C
         IF( C2 .EQ. 0D0 ) THEN
            CONDUC(1) =   C1 * Y32 * Y32
     %                  + C3 * X32 * X32
C
            CONDUC(2) = - C1 * Y32 * Y31
     %                  - C3 * X32 * X31
C
            CONDUC(3) =   C1 * Y31 * Y31
     %                  + C3 * X31 * X31
C
            CONDUC(4) =   C1 * Y32 * Y21
     %                  + C3 * X32 * X21
C
            CONDUC(5) = - C1 * Y31 * Y21
     %                  - C3 * X31 * X21
C
            CONDUC(6) =   C1 * Y21 * Y21
     %                  + C3 * X21 * X21
C
         ELSE
C
            CONDUC(1) =   C1 * Y32 * Y32
     %                  - C2 * X32 * Y32 * 2
     %                  + C3 * X32 * X32
C
            CONDUC(2) = - C1 * Y32 * Y31
     %                  + C2 * ( X32 * Y31 + Y32 * X31 )
     %                  - C3 * X32 * X31
C
            CONDUC(3) =   C1 * Y31 * Y31
     %                  - C2 * X31 * Y31 * 2
     %                  + C3 * X31 * X31
C
            CONDUC(4) =   C1 * Y32 * Y21
     %                  - C2 * ( X32 * Y21 + Y32 * X21 )
     %                  + C3 * X32 * X21
C
            CONDUC(5) = - C1 * Y31 * Y21
     %                  + C2 * ( X31 * Y21 + Y31 * X21 )
     %                  - C3 * X31 * X21
C
            CONDUC(6) =   C1 * Y21 * Y21
     %                  - C2 * X21 * Y21 * 2
     %                  + C3 * X21 * X21
         ENDIF

C        LE TERME -Alfa LAPLACIEN V ou W INITIALISE BE(I,1:2)
         GOTO( 11, 12, 13 ), I
 11      SV = CONDUC(1)*VEFtm(1) +CONDUC(2)*VEFtm(2) +CONDUC(4)*VEFtm(3)
         SW = CONDUC(1)*WEFtm(1) +CONDUC(2)*WEFtm(2) +CONDUC(4)*WEFtm(3)
         GOTO 20
 12      SV = CONDUC(2)*VEFtm(1) +CONDUC(3)*VEFtm(2) +CONDUC(5)*VEFtm(3)
         SW = CONDUC(2)*WEFtm(1) +CONDUC(3)*WEFtm(2) +CONDUC(5)*WEFtm(3)
         GOTO 20
 13      SV = CONDUC(4)*VEFtm(1) +CONDUC(5)*VEFtm(2) +CONDUC(6)*VEFtm(3)
         SW = CONDUC(4)*WEFtm(1) +CONDUC(5)*WEFtm(2) +CONDUC(6)*WEFtm(3)

C        REAL PART:
C          +Alfa LAPLACIEN Wn+1m
C           Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m -Wn+1m)
C          -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
 20      BE(I,1) = -SW + S 
     %              * ( FORCE(1) + FORCE(2)
     %                + Rho * ( VEFtn(I) + WEFtm(I) - WEFtn(I) )
     %                + F1dYF2dXV - F1dYF2dXW
     %                - NLSECOEF  * ( VEFtm(I) + WEFtm(I) )
     %                + NLSECOEF0 * VEFtm(I)  )

C        IMAG PART:
C          -Alfa LAPLACIEN Vn+1m
C          -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m +Wn+1m)
C          +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
         BE(I,2) = SV + S
     %              * ( FORCE(2) - FORCE(1)
     %                + Rho * ( -VEFtm(I) + VEFtn(I) + WEFtn(I) )
     %                + F1dYF2dXV + F1dYF2dXW
     %                + NLSECOEF  * ( VEFtm(I) - WEFtm(I) )
     %                + NLSECOEF0 * WEFtm(I)  )

      ENDDO
C
C     =============================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES 3 COTES
C     =============================================================
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
               TEMPEL = VEFtm( NSLK )
               ONDEPI = WEFtm( NSLK )
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
               BE(NSLK,1) = BE(NSLK,1) + D * FORCE(1)
               BE(NSLK,2) = BE(NSLK,2) + D * FORCE(2)

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
                  TEMPEL = VEFtm( K )
                  ONDEPI = WEFtm( K )
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
                  TEMPEL = VEFtm( K )
                  ONDEPI = WEFtm( K )
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
                  BE(K,1) = BE(K,1) + FORCE(1)
                  BE(K,2) = BE(K,2) + FORCE(2)
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
