      SUBROUTINE E32LAG( NBPOLY, NPI, POIDS, POLY, DPOLY, DDPOLY, X,
     %                   F1, F2, POIDEL, DP, DFM1,
     %                   DDP, DDFM21, DDFM22 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F1,F2  = COORDONNEES DES POINTS D INTEGRATION
C                         NUMERIQUES DE L ELEMENT COURANT
C                DFM1   = (DF) -1 AUX POINTS D INTEGRATION
C                DP     = DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C                POIDEL = PRODUIT DU POIDS PAR DELTA PAR R
C                         AUX POINTS D'INTEGRATION
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES DE BASE
C NPI    : NOMBRE DE POINTS D'INTEGRATION SUR L ELEMENT
C POIDS  : VALEUR DES NPI POIDS
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C DDPOLY : DDPOLY(I,J,L)=DERIVEE DDPJ/DDXI(POINT L) SELON L'ORDRE
C          DDXX DDXY DDYY I=1,3 DE PJ AU POINT L
C X      : COORDONNEES DES POINTS DE L'EF
C
C PARAMETRES RESULTATS:
C ---------------------
C F1,F2  : COORDONNEES XX ET YY DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C DFM1   : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D INTEGRATION
C DDP    : DERIVEES SECONDES DES POLYNOMES DE BASE AUX NPI POINTS
C          D'INTEGRATION DE L'ELEMENT FINI COURANT
C DDFM21 : MATRICE BLOC 21 DE LA MATRICE INVERSE POUR PASSER
C          DES DERIVEES SECONDES SUR EF REFERENCE A CELLES SUR EF COURANT
C DDFM22 : MATRICE BLOC DIAGONALE 22 DE LA MATRICE INVERSE POUR PASSER
C          DES DERIVEES SECONDES SUR EF REFERENCE A CELLES SUR EF COURANT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1995
C23456---------------------------------------------------------------012
      DOUBLE PRECISION POIDS(NPI),POLY(NBPOLY,NPI),
     %                 DPOLY(2,NBPOLY,NPI),F1(NPI),F2(NPI),
     %                 POIDEL(NPI),DP(2,NBPOLY,NPI),DFM1(4,NPI),
     %                 DDPOLY(3,NBPOLY,NPI),DDP(3,NBPOLY,NPI),
     %                 DDFM21(3,2,NPI),DDFM22(3,3,NPI)
      DOUBLE PRECISION XX,YY,DF1,DF2,DF3,DF4,D,DABS,
     %                 DD21(3,2),DXDY(3,3),DAUX(3,2)
      REAL             X(NBPOLY,2)
C
C     F1,DFM1,DELTA
C     =============
      DO 2 L=1,NPI
         XX = 0.D0
         YY = 0.D0
         DF1= 0.D0
         DF2= 0.D0
         DF3= 0.D0
         DF4= 0.D0
         DD21(1,1) = 0D0
         DD21(2,1) = 0D0
         DD21(3,1) = 0D0
         DD21(1,2) = 0D0
         DD21(2,2) = 0D0
         DD21(3,2) = 0D0
         DO 3 I=1,NBPOLY
C
            XX  = XX + POLY(I,L) * X(I,1)
            YY  = YY + POLY(I,L) * X(I,2)
C
C           DF1 = DF1/DX
            DF1 = DF1 + DPOLY(1,I,L) * X(I,1)
C           DF2 = DF1/DY
            DF2 = DF2 + DPOLY(2,I,L) * X(I,1)
C           DF3 = DF2/DX
            DF3 = DF3 + DPOLY(1,I,L) * X(I,2)
C           DF4 = DF2/DY
            DF4 = DF4 + DPOLY(2,I,L) * X(I,2)
C
C           LES DERIVEES SECONDES DE -F AU POINT L
            DD21(1,1) = DD21(1,1) - DDPOLY(1,I,L) * X(I,1)
            DD21(2,1) = DD21(2,1) - DDPOLY(2,I,L) * X(I,1)
            DD21(3,1) = DD21(3,1) - DDPOLY(3,I,L) * X(I,1)
C
            DD21(1,2) = DD21(1,2) - DDPOLY(1,I,L) * X(I,2)
            DD21(2,2) = DD21(2,2) - DDPOLY(2,I,L) * X(I,2)
            DD21(3,2) = DD21(3,2) - DDPOLY(3,I,L) * X(I,2)
    3    CONTINUE
C
         F1(L) = XX
         F2(L) = YY
C
         D     = DF1 * DF4 - DF2 * DF3
         DFM1(1,L) = DF4 / D
         DFM1(2,L) =-DF2 / D
         DFM1(3,L) =-DF3 / D
         DFM1(4,L) = DF1 / D
C
         POIDEL(L) = DABS(D) * POIDS(L)
C
C        LES CARRES DES DERIVEES PREMIERES
C        DXDY(1,1) = DF1/DX * DF1/DX
         DXDY(1,1) = DF1 * DF1
C        DXDY(2,1) = DF1/DX * DF1/DY
         DXDY(2,1) = DF1 * DF2
C        DXDY(3,1) = DF1/DY * DF1/DY
         DXDY(3,1) = DF2 * DF2
C
C        DXDY(1,2) = DF1/DX * DF2/DX + DF1/DX * DF2/DX
         DXDY(1,2) = 2 * DF1 * DF3
C        DXDY(2,2) = DF1/DX * DF2/DY + DF2/DX * DF1/DY
         DXDY(2,2) = DF1 * DF4 + DF2 * DF3
C        DXDY(3,2) = DF1/DY * DF2/DY + DF1/DY * DF2/DY
         DXDY(3,2) = 2 * DF2 * DF4
C
C        DXDY(1,3) = DF2/DX * DF2/DX
         DXDY(1,3) = DF3 * DF3
C        DXDY(2,3) = DF2/DX * DF2/DY
         DXDY(2,3) = DF3 * DF4
C        DXDY(3,3) = DF2/DY * DF2/DY
         DXDY(3,3) = DF4 * DF4
C
C        DP = DF-1 * DP SUR L'ELEMENT FINI UNITE
C        =======================================
         DO 4 I=1,NBPOLY
            DP(1,I,L) = DFM1(1,L) * DPOLY(1,I,L) +
     &                  DFM1(3,L) * DPOLY(2,I,L)
            DP(2,I,L) = DFM1(2,L) * DPOLY(1,I,L) +
     &                  DFM1(4,L) * DPOLY(2,I,L)
    4    CONTINUE
C
C        INVERSION DE LA MATRICE DDF
C        ===========================
C        DDF-1 AU POINT L => D22**-1
         CALL M33INV( DXDY, D, DDFM22(1,1,L) )
CCC         CALL DFM1R3( DXDY, DDFM22(1,1,L), IERR )
C
C        D22**-1 * -D21
         CALL AB0D( 3,3,2, DDFM22(1,1,L), DD21, DAUX )
C
C        D22**-1 * -D21 * D11**-1
         CALL AB0D( 3,2,2, DAUX, DFM1(1,L), DDFM21(1,1,L) )
C
C        DDP = DDF-1 * (DP, DDP) SUR L'ELEMENT DE REFERENCE
C        ==================================================
         DO 6 I=1, NBPOLY
            DDP(1,I,L) = DDFM21(1,1,L) *  DPOLY(1,I,L) +
     &                   DDFM21(1,2,L) *  DPOLY(2,I,L) +
     &                   DDFM22(1,1,L) * DDPOLY(1,I,L) +
     &                   DDFM22(1,2,L) * DDPOLY(2,I,L) +
     &                   DDFM22(1,3,L) * DDPOLY(3,I,L)
            DDP(2,I,L) = DDFM21(2,1,L) *  DPOLY(1,I,L) +
     &                   DDFM21(2,2,L) *  DPOLY(2,I,L) +
     &                   DDFM22(2,1,L) * DDPOLY(1,I,L) +
     &                   DDFM22(2,2,L) * DDPOLY(2,I,L) +
     &                   DDFM22(2,3,L) * DDPOLY(3,I,L)
            DDP(3,I,L) = DDFM21(3,1,L) *  DPOLY(1,I,L) +
     &                   DDFM21(3,2,L) *  DPOLY(2,I,L) +
     &                   DDFM22(3,1,L) * DDPOLY(1,I,L) +
     &                   DDFM22(3,2,L) * DDPOLY(2,I,L) +
     &                   DDFM22(3,3,L) * DDPOLY(3,I,L)
    6    CONTINUE
    2 CONTINUE
      END
