      SUBROUTINE E16Q1C( NBPOLY, NPI,    POIDS, POLY, DPOLY, X,
     %                   F,      POIDEL, DP,    DFM1, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER POUR LE 6CUBE
C -----    F        = COORDONNEES DES POINTS D INTEGRATION
C                     NUMERIQUES DE L EF COURANT
C          DFM1     = (DF) -1 AUX POINTS D INTEGRATION
C          DP       = DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
C          POIDEL   = PRODUIT DU POIDS PAR DELTA
C                     AUX POINTS D'INTEGRATION
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POLYNOMES DE BASE DE L'EF 6Q1C
C NPI    : NOMBRE DE POINTS D'INTEGRATION SUR L EF VOLUMIQUE
C POIDS  : VALEUR DES NPI POIDS
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C X      : COORDONNEES DES POINTS DE L'EF
C
C SORTIES:
C --------
C F      : COORDONNEES XYZUVW  DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L'EF COURANT
C DFM1   : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D'INTEGRATION
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              X(NBPOLY,6)
      DOUBLE PRECISION  POIDS(         NPI),
     %                  POLY(   NBPOLY,NPI),
     %                  DPOLY(6,NBPOLY,NPI),
     %                  F(NPI,6),
     %                  POIDEL(NPI),
     %                  DP(6,NBPOLY,NPI),
     %                  DFM1(6,6,NPI)
      DOUBLE PRECISION  DF(6,6), DELTA, DABS, D
ccc
ccc      double precision  dfdfm1(6,6)   verification de l'inversion
C
C      LES 6 XYZUVW DES NBPOLY NOEUDS=SOMMETS
CCC         DO 31 J=1,NBPOLY
CCC            WRITE(IMPRIM,10031) (J,K,X(J,K),K=1,6)
CCC 31      CONTINUE
CCC10031    FORMAT(6('  XYZUVW(',I2,',',I1,')=',D15.8))
C
      DO 100 L=1,NPI
C
C        F1,DFM1,DELTA  AU POINT D'INTEGRATION L
C        =======================================
         DO 30 K=1,6
C           LA COORDONNEE K DE F AU POINT D'INTEGRATION L
            F(L,K) = 0.D0
C           LA COLONNE K DE DF AU POINT D'INTEGRATION L
            DO 5 J=1,6
               DF(J,K) = 0.D0
 5          CONTINUE
            DO 20 I=1,NBPOLY
               D       = X(I,K)
               F(L,K)  = F(L,K) + POLY(I,L) * D
               DO 10 J=1,6
                  DF(J,K) = DF(J,K) + DPOLY(J,I,L) * D
 10            CONTINUE
 20         CONTINUE
 30      CONTINUE
C
CCC         WRITE(IMPRIM,10029) (I,L,POLY(I,L),I=1,NBPOLY)
CCC10029    FORMAT(6('  POLY(',I2,',',I2,')=',D15.8))
C
CCC         WRITE(IMPRIM,10032) (L,K,F(L,K),K=1,6)
CCC10032    FORMAT(6('  F(',I2,',',I1,')=',D15.8))
C
CCC         DO 33 J=1,6
CCC            WRITE(IMPRIM,10033) (J,K,DF(J,K),K=1,6)
CCC 33      CONTINUE
CCC10033    FORMAT(6('  DF(',I1,',',I1,')=',D15.8))
CCC         DO 34 K=1,NBPOLY
CCC            WRITE(IMPRIM,10034) (J,K,L,DPOLY(J,K,L),J=1,6)
CCC 34      CONTINUE
CCC10034    FORMAT(6('  DPOLY(',I1,',',I2,',',I2,')=',D15.8))
C
C        [DF]-1 AU POINT D'INTEGRATION L
         CALL M66INV( DF,  DFM1(1,1,L), DELTA, IERR )
ccc
ccc         L'inverse DFM1 est obtenue avec une bonne precision 10**-30
ccc         do 22 i=1,6
ccc            do 23 j=1,6
ccc               dfdfm1(i,j) = 0D0
ccc               do 24 k=1,6
ccc                  dfdfm1(i,j) = dfdfm1(i,j) + df(i,k) * dfm1(k,j,l)
ccc 24            continue
ccc 23         continue
ccc            print 10000,(i,j,dfdfm1(i,j),j=1,6)
ccc10000       format(6('  Id(',I1,',',I1,')=',D15.7))
ccc 22      continue
c
         IF( IERR .NE. 0 ) RETURN
C
C        | DELTA | * POIDS  AU POINT D'INTEGRATION L
         POIDEL(L) = DABS(DELTA) * POIDS(L)
C
C        DP(XL) = DF-1(XCL) * DPC(XCL) AU POINT D'INTEGRATION L
C        ======================================================
         DO 60 I=1,NBPOLY
            DO 50 K=1,6
               DP(K,I,L) = 0.D0
               DO 40 J=1,6
                  DP(K,I,L) = DP(K,I,L) + DFM1(K,J,L) * DPOLY(J,I,L)
 40            CONTINUE
 50         CONTINUE
 60      CONTINUE
 100  CONTINUE
C
      RETURN
      END
