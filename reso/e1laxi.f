      SUBROUTINE E1LAXI( D2PI,NBPOLY,NPI,POIDS,POLY,DPOLY,X,
     %                   F1,F2,POIDEL,DP,DFM1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F1,F2= COORDONNEES DES POINTS D INTEGRATION
C                       NUMERIQUES DE L ELEMENT COURANT
C                DFM1 = (DF) -1 AUX POINTS D INTEGRATION
C                DP   = DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C                POIDEL : PRODUIT DU POIDS PAR DELTA PAR R
C                         AUX POINTS D'INTEGRATION
C
C PARAMETRES D ENTREE :
C ---------------------
C D2PI   : 2 FOIS PI
C NBPOLY : NBRE DE POLYNOMES DE BASE
C NPI    : NBRE DE POINTS D'INTEGRATION SUR L ELEMENT FINI
C POIDS  : VALEUR DES NPI POIDS
C POLY   : POLY(I,J)=VALEUR DE PI(POINT J D INTEGRATION)
C DPOLY  : DPOLY(I,J,L)=DERIVEE DPJ/DXI(POINT L)
C X      : RAYON ET COTE DES NBPOLY POINTS DE L ELEMENT
C
C PARAMETRES RESULTATS:
C ---------------------
C F1,F2  : COORDONNEES R ET Z DES NPI POINTS D INTEGRATION DE L ELEMENT
C POIDEL : DELTA * POIDS(NPI) * D2PI * R(NPI) DES NPI POINTS D INTEGRATI
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C DFM1   : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D INTEGRATION
C
C NBPOLY : NBRE DE POLYNOMES DE BASE SI PAS D'ERREUR
C          0 SI EF AVEC X<=0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE PARIS   JUILLET 1989
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION POIDS(NPI),POLY(NBPOLY,NPI),
     %                 DPOLY(2,NBPOLY,NPI),D2PI,F1(NPI),F2(NPI),
     %                 POIDEL(NPI),DP(2,NBPOLY,NPI),DFM1(4,NPI)
      DOUBLE PRECISION RR,ZZ,DF1,DF2,DF3,DF4,D,DABS,R,Z
      REAL             X(NBPOLY,2)
C
C     F1,DFM1,DELTA
C     =============
      DO 2 L=1,NPI
         RR = 0.D0
         ZZ = 0.D0
         DF1= 0.D0
         DF2= 0.D0
         DF3= 0.D0
         DF4= 0.D0
         DO 3 I=1,NBPOLY
            R = X(I,1)
            Z = X(I,2)
            RR  = RR + POLY(I,L) * R
            ZZ  = ZZ + POLY(I,L) * Z
            DF1 = DF1 + DPOLY(1,I,L) * R
            DF2 = DF2 + DPOLY(2,I,L) * R
            DF3 = DF3 + DPOLY(1,I,L) * Z
            DF4 = DF4 + DPOLY(2,I,L) * Z
    3    CONTINUE
C
         IF( RR .LE. 0D0 ) THEN
C           EF AXISYMETRIQUE AVEC X <= 0   INTERDIT !
            NBLGRC(NRERR) = 2
            KERR(1) = 'ERREUR: PB AXISYMETRIQUE AVEC EF '
            KERR(2) = 'D''ABSCISSE X=R NEGATIVE OU NULLE'
            CALL LEREUR
C           MISE A ZERO POUR TEST EN SORTIE
            NBPOLY = 0
            RETURN
         ENDIF
C
         F1(L) = RR
         F2(L) = ZZ
         D     = DF1 * DF4 - DF2 * DF3
         DFM1(1,L) = DF4 / D
         DFM1(2,L) =-DF2 / D
         DFM1(3,L) =-DF3 / D
         DFM1(4,L) = DF1 / D
         POIDEL(L) = DABS(D) * POIDS(L) * D2PI * RR
C
C        DP = DF-1 * DP SUR L'ELEMENT FINI UNITE
C        =======================================
         DO 4 I=1,NBPOLY
            DP(1,I,L) = DFM1(1,L) * DPOLY(1,I,L) +
     &                  DFM1(3,L) * DPOLY(2,I,L)
            DP(2,I,L) = DFM1(2,L) * DPOLY(1,I,L) +
     &                  DFM1(4,L) * DPOLY(2,I,L)
    4    CONTINUE
    2 CONTINUE
      END
