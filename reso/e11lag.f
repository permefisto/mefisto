      SUBROUTINE E11LAG( NBPOLY, NPI, POIDS, POLY, DPOLY, X,
     %                   F1,     POIDEL, DP, DF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F1     = COORDONNEES DES POINTS D INTEGRATION
C                         NUMERIQUES DE L ELEMENT FINI COURANT
C                DFM1   = (DF) -1 AUX POINTS D INTEGRATION
C                DP     = DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C                POIDEL = PRODUIT DU POIDS PAR DELTA PAR R
C                         AUX POINTS D'INTEGRATION
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES DE BASE
C NPI    : NOMBRE DE POINTS D'INTEGRATION SUR L ELEMENT FINI
C POIDS  : VALEUR DES NPI POIDS DES POINTS D'INTEGRATION
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C X      : COORDONNEE DES NBPOLY POINTS DE L'EF
C
C PARAMETRES RESULTATS:
C ---------------------
C F1     : COORDONNEE XX DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C          AUX POINTS D'INTEGRATION NUMERIQUE INTERNES
C DF     : DETERMINANT DE LA MATRICE JACOBIENNE DE F: e ref -> e
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     JUIN 2009
C23456---------------------------------------------------------------012
      DOUBLE PRECISION POIDS(NPI), POLY(NBPOLY,NPI),
     %                 DPOLY(NBPOLY,NPI), F1(NPI),
     %                 POIDEL(NPI), DP(NBPOLY,NPI)
      DOUBLE PRECISION XX, DF
      REAL             X(NBPOLY)
C
C     DF = he
C     =======
      DF = X(2) - X(1)
C
      DO 8 L=1,NPI
C
C        F1(bl)
C        ======
         XX = 0.D0
         DO 2 I=1,NBPOLY
            XX = XX + POLY(I,L) * X(I)
 2       CONTINUE
         F1(L) = XX
C
C        POIDS(l) * DF
C        =============
         POIDEL(L) = ABS(DF) * POIDS(L)
C
C        DP = DF-1 * DP SUR L'ELEMENT FINI UNITE
C        =======================================
         DO 4 I=1,NBPOLY
            DP(I,L) = DPOLY(I,L) / DF
 4       CONTINUE
C
 8    CONTINUE
C
      RETURN
      END
