       SUBROUTINE E1LTDI (NBPOLY,NPI,POIDS,POLY,DPOLY,X,IP,
     %                    F1,F2,F3,POIDEL,DP,DFM1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER
C -----          F1,F2,F3= COORDONNEES DES POINTS D'INTEGRATION NUMERIQUES DE L'
C                DFM1    = (DF) -1 AUX POINTS D'INTEGRATION
C                DP      = DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
C                POIDEL  = PRODUIT DU POIDS PAR DELTA AUX POINTS D'INTEGRATION
C
C PARAMETRES D'ENTREE :
C ---------------------
C NBPOLY  : NBRE DE POLYNOMES DE BASE
C NPI     : NBRE DE POINTS D'INTEGRATION SUR L ELEMENT
C POIDS   : VALEUR DES NPI POIDS
C POLY    : POLY(I,J) =VALEUR DE PI (POINT J D'INTEGRATION)
C DPOLY   : DPOLY(I,J,L0 =DERIVEE DPJ/DXI(POINT L)
C X       : COORDONNEES DES NBPOLY POINTS DE L'ELEMENT
C
C PARAMETRES RESULTATS :
C-----------------------
C IP       : IP(J) = POSITION DU J-EME D.L. COMPOSANTE PAR COMPOSANTE DANS L'ORD
C F1,F2,F3 : COORDONNEES X,Y,Z PES NPI POINTS D'INTEGRATION DE L'ELEMENT
C POIDEL   : DELTA * POIDS(NPI) DES NPI POINTS D'INTEGRATION
C DP       : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C DFM1     : MATRICE JACOBIENNE INVERSE AUX NPI POINTS D'INTRGRATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ZEBIC,DUTERTE
C...............................................................................
        DOUBLE PRECISION POIDS (NPI) ,POLY (NBPOLY,NPI),
     %                   DPOLY (3,NBPOLY,NPI),
     %                   F1(NPI), F2(NPI), F3(NPI),
     %                   POIDEL(NPI), DP(3,NBPOLY,NPI),
     %                   DFM1(3,3,NPI),DF(3,3)
C
        DOUBLE PRECISION XX,YY,ZZ,D,DABS
        REAL             X(NBPOLY,3)
        INTEGER          IP(*)
C                              EN FAIT IP(3*NBPOLY)  !
C       IP(I) = POSITION DU I-EME D.L. PAR COMPOSANTE DANS L'ORDRE PAR NOEUD PAR
C       ========================================================================
C
C CALCUL DE IP
C
        J=1
        DO 1 I=1,NBPOLY
           IP(I) = J
           IP(I+NBPOLY) = J+1
           IP(I+2*NBPOLY) = J+2
           J = J+3
 1      CONTINUE
C
        DO 2 L=1,NPI
           XX = 0.D0
           YY = 0.D0
           ZZ = 0.D0
           CALL AZER0D(6,DF)
           DO 3 J=1,NBPOLY
              XX = XX + POLY(J,L) * X(J,1)
              YY = YY + POLY(J,L) * X(J,2)
              ZZ = ZZ + POLY(J,L) * X(J,3)
              DO 10 I=1,3
                 DF(I,1) = DF(I,1) + DPOLY(I,J,L) * X(J,1)
                 DF(I,2) = DF(I,2) + DPOLY(I,J,L) * X(J,2)
                 DF(I,3) = DF(I,3) + DPOLY(I,J,L) * X(J,3)
 10           CONTINUE
 3         CONTINUE
           F1(L) = XX
           F2(L) = YY
           F3(L) = ZZ
           D =   DF(1,1) * (DF(2,2)*DF(3,3) - DF(2,3)*DF(3,2))
     %         - DF(2,1) * (DF(1,2)*DF(3,3) - DF(1,3)*DF(3,2))
     %         + DF(3,1) * (DF(1,2)*DF(2,3) - DF(2,2)*DF(1,3))
C
           DFM1(1,1,L) = (1/D) * DF(1,1) * (DF(2,2)*DF(3,3)
     %                     - DF(2,3)*DF(3,2))
           DFM1(1,2,L) = -(1/D) * DF(2,1) * (DF(1,2)*DF(3,3)
     %                      - DF(1,3)*DF(3,2))
           DFM1(1,3,L) = (1/D) * DF(3,1) * (DF(1,2)*DF(2,3)
     %                       - DF(1,3)*DF(2,2))
C
           DFM1(2,1,L) = -(1/D) * DF(1,2) * (DF(2,1)*DF(3,3)
     %                   - DF(3,1)*DF(2,3))
           DFM1(2,2,L) = (1/D) * DF(2,2) * (DF(1,1)*DF(3,3)
     %                    - DF(1,3)*DF(3,1))
           DFM1(2,3,L) = -(1/D) * DF(3,2) * (DF(1,1)*DF(3,2)
     %                     - DF(2,2)*DF(3,1))
C
           DFM1(3,1,L) = (1/D) * DF(1,3) * (DF(2,1)*DF(3,2)
     %                   - DF(2,2)*DF(3,1))
           DFM1(3,2,L) = -(1/D) * DF(2,3) * (DF(1,1)*DF(3,2)
     %                    - DF(1,2)*DF(3,1))
           DFM1(3,3,L) = (1/D) * DF(3,3) * (DF(1,1)*DF(2,2)
     %                     - DF(1,2)*DF(2,1))
C
           POIDEL(L) = DABS(D) * POIDS(L)
C
C          CALCUL DE ( DFM1 * DP )
           DO 11 J=1,NBPOLY
              DO 12 I=1,3
                 DP(I,J,L)=0.D0
 12           CONTINUE
 11        CONTINUE
           DO 4 J=1,NBPOLY
              DO 13 I=1,3
                 DO 14 K=1,3
                    DP(I,J,L) = DP(I,J,L) + (DFM1(K,I,L)*DPOLY(K,J,L))
 14              CONTINUE
 13           CONTINUE
 4         CONTINUE
 2      CONTINUE
        END
