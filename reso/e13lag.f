      SUBROUTINE E13LAG( NBPOLY, NPI,    POIDS, POLY, DPOLY, X,
     %                   F,      POIDEL, DP,    DFM1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F        = COORDONNEES DES POINTS D INTEGRATION
C                           NUMERIQUES DE L EF COURANT
C                DFM1     = (DF) -1 AUX POINTS D INTEGRATION
C                DP       = DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
C                POIDEL   = PRODUIT DU POIDS PAR DELTA
C                           AUX POINTS D'INTEGRATION
C
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POLYNOMES DE BASE DE L'EF VOLUMIQUE
C NPI    : NOMBRE DE POINTS D'INTEGRATION SUR L EF VOLUMIQUE
C POIDS  : VALEUR DES NPI POIDS
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C X      : COORDONNEES DES POINTS DE L'EF
C
C SORTIES:
C --------
C F      : COORDONNEES XX YY ZZ DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L'EF COURANT
C DFM1   : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D'INTEGRATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  POIDS(NPI),
     %                  POLY(NBPOLY,NPI),
     %                  DPOLY(3,NBPOLY,NPI),
     %                  F(NPI,3),
     %                  POIDEL(NPI),
     %                  DP(3,NBPOLY,NPI),
     %                  DFM1(3,3,NPI)
      REAL              X(NBPOLY,3)
      DOUBLE PRECISION  DF(3,3), D, DABS
C
cccC *********************************************************************
cccC SENSIBILITE DE L'INTEGRATION du PB de l'HYDROGEN  Int( 1/r ) sur boule
cccC     PROJECTION DES NOEUDS FRONTIERE SUR LA SPHERE
ccc      DOUBLE PRECISION RADIUS, R  (Article Extrapolation.pdf ...)
ccc      RADIUS = 60D0
ccc      DO 5 I=1,NBPOLY
ccc         R = SQRT( X(I,1)**2 + X(I,2)**2 +  X(I,3)**2 )
ccc         IF( R .GT. RADIUS*0.98D0 ) THEN
cccC           PROJECTION SUR LA SPHERE CENTREE DE RAYON RADIUS
ccc            X(I,1) = X(I,1) * RADIUS / R
ccc            X(I,2) = X(I,2) * RADIUS / R
ccc            X(I,3) = X(I,3) * RADIUS / R
ccc         ENDIF
ccc 5    CONTINUE
cccC **********************************************************************
C
      DO 100 L=1,NPI
C
C        F1,DFM1,DELTA  AU POINT D'INTEGRATION L
C        =======================================
         DO 20 K=1,3
C           LA COORDONNEE K DE F AU POINT D'INTEGRATION L
            F(L,K) = 0.D0
C           LA COLONNE K DE DF AU POINT D'INTEGRATION L
            DF(1,K) = 0.D0
            DF(2,K) = 0.D0
            DF(3,K) = 0.D0
            DO 10 I=1,NBPOLY
               D       = X(I,K)
               F(L,K)  = F(L,K)  + POLY(I,L)    * D
               DF(1,K) = DF(1,K) + DPOLY(1,I,L) * D
               DF(2,K) = DF(2,K) + DPOLY(2,I,L) * D
               DF(3,K) = DF(3,K) + DPOLY(3,I,L) * D
 10         CONTINUE
 20      CONTINUE
C
C        [DF]-1 AU POINT D'INTEGRATION L
         CALL M33INV( DF, D, DFM1(1,1,L) )
C
C        POIDS * | DELTA | AU POINT D'INTEGRATION L
         POIDEL(L) = DABS(D) * POIDS(L)
C
C        DP(XL) = DF-1(XCL) * DPC(XCL) AU POINT D'INTEGRATION L
C        ======================================================
         DO 60 I=1,NBPOLY
            DO 50 K=1,3
               DP(K,I,L) = 0.D0
               DO 40 J=1,3
                  DP(K,I,L) = DP(K,I,L) + DFM1(K,J,L) * DPOLY(J,I,L)
 40            CONTINUE
 50         CONTINUE
 60      CONTINUE
 100  CONTINUE
C
      RETURN
      END
