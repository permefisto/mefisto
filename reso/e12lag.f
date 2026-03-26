      SUBROUTINE E12LAG( NBPOLY, NPI, POIDS, POLY, DPOLY, X,
     %                   F1, F2, POIDEL, DP, DFM1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F1,F2  = COORDONNEES DES POINTS D INTEGRATION
C                         NUMERIQUES DE L ELEMENT COURANT
C                DFM1   = (DF) -1 AUX POINTS D INTEGRATION
C                DP     = DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C                POIDEL = PRODUIT DU POIDS PAR DELTA
C                         AUX POINTS D'INTEGRATION
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES DE BASE
C NPI    : NOMBRE DE POINTS D'INTEGRATION SUR L ELEMENT
C POIDS  : VALEUR DES NPI POIDS
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C X      : COORDONNEES DES POINTS DE L'EF
C
C PARAMETRES RESULTATS:
C ---------------------
C F1,F2  : COORDONNEES XX ET YY DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C DFM1   : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D INTEGRATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      DOUBLE PRECISION POIDS(NPI), POLY(NBPOLY,NPI),
     %                 DPOLY(2,NBPOLY,NPI), F1(NPI), F2(NPI),
     %                 POIDEL(NPI), DP(2,NBPOLY,NPI), DFM1(4,NPI)
      DOUBLE PRECISION XX, YY, DF1, DF2, DF3, DF4, D
      INTRINSIC        ABS
      REAL             X(NBPOLY,2)
C
C     F1,F2,DFM1,DELTA
C     ================
      DO L=1,NPI
         XX = 0.D0
         YY = 0.D0
         DF1= 0.D0
         DF2= 0.D0
         DF3= 0.D0
         DF4= 0.D0
         DO I=1,NBPOLY
            XX  = XX + POLY(I,L) * X(I,1)
            YY  = YY + POLY(I,L) * X(I,2)
            DF1 = DF1 + DPOLY(1,I,L) * X(I,1)
            DF2 = DF2 + DPOLY(2,I,L) * X(I,1)
            DF3 = DF3 + DPOLY(1,I,L) * X(I,2)
            DF4 = DF4 + DPOLY(2,I,L) * X(I,2)
         ENDDO
         F1(L) = XX
         F2(L) = YY
         D     = DF1 * DF4 - DF2 * DF3
         DFM1(1,L) = DF4 / D
         DFM1(2,L) =-DF2 / D
         DFM1(3,L) =-DF3 / D
         DFM1(4,L) = DF1 / D
         POIDEL(L) = ABS(D) * POIDS(L)
C
C        DP = DF-1 * DP SUR L'ELEMENT FINI UNITE
C        =======================================
         DO I=1,NBPOLY
            DP(1,I,L) = DFM1(1,L) * DPOLY(1,I,L) +
     %                  DFM1(3,L) * DPOLY(2,I,L)
            DP(2,I,L) = DFM1(2,L) * DPOLY(1,I,L) +
     %                  DFM1(4,L) * DPOLY(2,I,L)
         ENDDO
      ENDDO
C
      RETURN
      END
