      SUBROUTINE GP3BE2P1( DeltaT, Rho, Omega,  Alfa, Beta, Force,
     %                     NONOEF, X,   NBNOMA, Utn,  Utm, 
     %                     BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GPE: TESTNL=10 ASSEMBLAGE DANS BG DES 4 COEFFICIENTS DES
C ----- 2 POLYNOMES DU 3-EME DEGRE D'UN EF TRIANGLE P1 LAGRANGE 2D

C -  Beta Vn+1m+1**3 
C -  Beta Wn+1m * Vn+1m+1**2
C + {  Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Vn+1m
C + { -Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Wn+1m
C + Rho/dt (Wn-Vn) - Fr -Fi = 0

C - Beta Wn+1m+1**3 
C + Beta Vn+1m * Wn+1m+1**2
C + {  Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Wn+1m
C - { -Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Vn+1m
C - Rho/dt (Wn+Vn) + Fr -Fi = 0

C ENTREES:
C --------
C DeltaT : PAS DE TEMPS DU SCHEMA IMPLICITE EN TEMPS
C Rho    : DENSITE DE MASSE DE L'EF
C Omega  : VITESSE ANGULAIRE DE LA ROTATION
C Alfa   : COEFFICIENT DU LAPLACIEN
C Beta   : COEFFICIENT DE | V**2 + W**2 |
C Force  : 2 COMPOSANTES DE LA FORCE SUPPOSEE CONSTANTE
C NONOEF : NUMERO DES 3 NOEUDS SOMMETS DE L'EF
C X      : COORDONNEES RAYON ET COTE DES 3 SOMMETS DE L'EF
C          OU X Y DES 3 SOMMETS DE L'EF
C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE
C Utn    : U(tn) (NBNOMA,2) U A L'INTANT tn
C Utm    : U(tn+1,m) (NBNOMA,2)  U A L'INTANT tn+1 iteration m
C
C SORTIE :
C --------
C BG     : BG(NBNOMA,0:3,2) LES 4 COEFFICIENTS DES 2 POLYNOMES DE DEGRE 3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Juin 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      DOUBLE PRECISION  DeltaT, Rho, Omega(3), Alfa, Beta, Force(2)
      REAL              X(3,2)
      INTEGER           NONOEF(3)

      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2), BG(NBNOMA,0:3,2)
      DOUBLE PRECISION  DP(2,3), VEFtm(3), WEFtm(3), VEFtn(3), WEFtn(3)
C
      DOUBLE PRECISION  XYZ(3), RhoDel, V, DPV1,DPV2, DPW1,DPW2
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32,
     %                  DELTA, DELTAS6, D6Beta
C
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL tn+1m et tn
C     AUX 3 SOMMETS DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     ------------------------------------------------------------------
      DO I = 1, 3
C        NO GLOBAL DU NOEUD I DE L'EF
         K = NONOEF(I)
C        AU TEMPS tn
         VEFtn(I) = Utn(K,1)
         WEFtn(I) = Utn(K,2)
C        AU TEMPS tn+1m
         VEFtm(I) = Utm(K,1)
         WEFtm(I) = Utm(K,2)
      ENDDO

C     ===================================
C     CONTRIBUTION DES FORCES SURFACIQUES
C     ===================================
      XYZ(3) = 0D0
      RhoDel = Rho / DeltaT

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
      DP(1,1) = -Y32 / DELTA
      DP(1,2) =  Y31 / DELTA
      DP(1,3) = -Y21 / DELTA

      DP(2,1) =  X32 / DELTA
      DP(2,2) = -X31 / DELTA
      DP(2,3) =  X21 / DELTA

C     COEFFICIENTS CALCULES UNE SEULE FOIS
      DELTAS6 = DELTA / 6D0
      D6Beta  = DELTAS6 * BETA

      DO I = 1, 3

C        NO GLOBAL DU NOEUD I DE L'EF
         K = NONOEF(I)

C        COORDONNEES DU POINT D'INTEGRATION I = SOMMET I
         XYZ(1) = X(I,1)
         XYZ(2) = X(I,2)

C        LE POTENTIEL QUARTIC V = r2/2 +r4/4
         V = XYZ(1)**2 + XYZ(2)**2
         V = V * ( 0.5D0 + V * 0.25D0 )

         DPV1 = DP(1,1)*VEFtm(1) +DP(1,2)*VEFtm(2) +DP(1,3)*VEFtm(3)
         DPV2 = DP(2,1)*VEFtm(1) +DP(2,2)*VEFtm(2) +DP(2,3)*VEFtm(3)

         DPW1 = DP(1,1)*WEFtm(1) +DP(1,2)*WEFtm(2) +DP(1,3)*WEFtm(3)
         DPW2 = DP(2,1)*WEFtm(1) +DP(2,2)*WEFtm(2) +DP(2,3)*WEFtm(3)

C -  Beta Vn+1m+1**3 
C -  Beta Wn+1m * Vn+1m+1**2
C + {  Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Vn+1m
C + { -Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Wn+1m
C + Rho/dt (Wn-Vn) - Fr -Fi = 0
C        COEFFICIENT CONSTANT DU POLYNOME EN V
         BG(K,0,1) = BG(K,0,1) + DELTAS6 *
     %             ( ( -RhoDel - V - Beta * WEFtm(I)**2 ) * WEFtm(I)
     %              +( Alfa * DP(1,I) - Omega(1) * xyz(2) ) * DPW1
     %              +( Alfa * DP(2,I) + Omega(1) * xyz(1) ) * DPW2
     %              + RhoDel * ( WEFtn(I) - VEFtn(I) )
     %              - Force(1) - Force(2)
     %              +( Alfa * DP(1,I) + Omega(1) * xyz(2) ) * 
     %                  (DPV1-DP(1,I)*VEFtm(I))
     %              +( Alfa * DP(2,I) - Omega(1) * xyz(1) ) *
     %                  (DPV2-DP(2,I)*VEFtm(I)) )

C        COEFFICIENT DU MONOME V
         BG(K,1,1) = BG(K,1,1) + DELTAS6 *
     %             ( RhoDel - V - Beta * WEFtm(I)**2
     %             +( Alfa * DP(1,I) + Omega(1) * xyz(2) ) * DP(1,I)
     %             +( Alfa * DP(2,I) - Omega(1) * xyz(1) ) * DP(2,I) )

C        COEFFICIENT DU MONOME V**2
         BG(K,2,1) = BG(K,2,1) - D6Beta * WEFtm(I)

C        COEFFICIENT DU MONOME V**3
         BG(K,3,1) = BG(K,3,1) - D6Beta

C - Beta Wn+1m+1**3 
C + Beta Vn+1m * Wn+1m+1**2
C + {  Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Wn+1m
C - { -Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Vn+1m
C - Rho/dt (Wn+Vn) + Fr -Fi = 0
C        COEFFICIENT CONSTANT DU POLYNOME EN W
         BG(K,0,2) = BG(K,0,2) + DELTAS6 *
     %             ( ( RhoDel + V + Beta * VEFtm(I)**2 ) * VEFtm(I)
     %              +(-Alfa * DP(1,I) + Omega(1) * xyz(2) ) * DPV1
     %              +(-Alfa * DP(2,I) - Omega(1) * xyz(1) ) * DPV2
     %              - RhoDel * ( VEFtn(I) + WEFtn(I) )
     %              + Force(1) - Force(2)
     %              +( Alfa * DP(1,I) + Omega(1) * xyz(2) ) * 
     %                  (DPW1-DP(1,I)*WEFtm(I))
     %              +( Alfa * DP(2,I) - Omega(1) * xyz(1) ) *
     %                  (DPW2-DP(2,I)*WEFtm(I)) )

C        COEFFICIENT DU MONOME W
         BG(K,1,2) = BG(K,1,2) + DELTAS6 *
     %             ( RhoDel - V - Beta * VEFtm(I)**2
     %             +( Alfa * DP(1,I) + Omega(1) * xyz(2) ) * DP(1,I)
     %             +( Alfa * DP(2,I) - Omega(1) * xyz(1) ) * DP(2,I) )

C        COEFFICIENT DU MONOME W**2
         BG(K,2,2) = BG(K,2,2) + D6Beta * VEFtm(I)

C        COEFFICIENT DU MONOME W**3
         BG(K,3,2) = BG(K,3,2) - D6Beta

      ENDDO

      RETURN
      END
