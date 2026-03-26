      SUBROUTINE ERLAXI( D2PI,   PENALI, NBSOMT, NBCOTE,
     &                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     &                   NOOBLA, NUMILI, NUMALI, LTDELI,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   NBPOLA, NBPOLY, NPI,    POLY,
     &                   IP,     F1,     F2,     POIDEL, DP,
     &                   G1,     G2,     G3,
     &                   ELAS,   AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE RAIDEUR DES ELEMENTS AXISYMETRIQUES
C -----
C
C ENTREES:
C --------
C D2PI   : 2 FOIS PI
C PENALI : 1/EPSILON DE LA PENALISATION DU DEPLACEMENT IMPOSE
C NBSOMT : NOMBRE DE SOMMETS DE L ELEMENT FINI
C NBCOTE : NOMBRE DES COTES  DE L ELEMENT FINI
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL  DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL  DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES POINTS DE L'OBJET
C
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES LIGNES DE L'OBJET
C NUMALI : NUMERO MAXIMAL DES LIGNES DE L'OBJET
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES LIGNES DE L'OBJET
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES SURFACES DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C
C NBPOLA : NOMBRE DE POLYNOMES DE BASE SUR UN COTE DE L'ELEMENT FINI
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR UN COTE
C POIDSA : POIDS DES POINTS D INTEGRATION SUR UN COTE
C POLYA  : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION COTE
C
C NBPOLY : NOMBRE DE POLYNOMES DE L'ELEMENT FINI COMPLET
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE
C POIDS  : LES NPI POIDS DE LA FORMULE D INTEGRATION
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I,  L)= P(I) (XL)
C
C IP     : IP(J) = POSITION DI J-EME D.L. COMPOSANTE PAR COMPOSANTE
C          DANS LA NUMEROTATION PAR NOEUDS
C F1     : COORDONNEES XX DES NPI POINTS D INTEGRATION DE L ELEMENT
C F2     : COORDONNEES YY DES NPI POINTS D INTEGRATION DE L ELEMENT
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C G1,G2,G3: TABLEAUX AUXILIAIRES
C
C SORTIES:
C --------
C ELAS   : TENSEUR SYMETRIQUE DE L ELASTICITE AXISYMETRIQUE
C AE     : LA MATRICE DE RAIDEUR SYMETRIQUE PLEINE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE PARIS    JUILLET 1989
C ......................................................................
      include"./incl/a___fixation.inc"
      include"./incl/donela.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      DOUBLE PRECISION D2PI
      INTEGER          LTDEPO( 1:MXDOEL, NUMIPO:NUMAPO )
      INTEGER          LTDELI( 1:MXDOEL, NUMILI:NUMALI )
      INTEGER          LTDESU( 1:MXDOEL, NUMISU:NUMASU )
      INTEGER          NOOBPS( 1:NBSOMT )
      INTEGER          NOOBLA( 1:NBCOTE )
      DOUBLE PRECISION PENALI,
     %                 POLY(NBPOLY,NPI),
     %                 F1(NPI),F2(NPI),POIDEL(NPI),
     %                 DP(2,NBPOLY,NPI),
     %                 G1(NBPOLY,NBPOLY),
     %                 G2(NBPOLY),
     %                 G3(NBPOLY,2),
     %                 AE(*)
      DOUBLE PRECISION AUX,A11(1),A12(2),A13(2),A22(2,2),
     %                 A23(2,2),A33(3),ELAS(10),COEFDE(2,9)
      INTEGER          IP(*),NOPOAR(3)
C
C     INITIALISATION DE AE
C     --------------------
      CALL AZEROD( NBPOLY*(2*NBPOLY+1) , AE )
C
C     CONSTRUCTION DES MATRICES   (A11 , A12 , A13)
C     -------------------------   (      A22 , A23)  =  (TD * E * D)
C                                 (            A33)
      DO 10 L=1,NPI
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE L'ELASTICITE AU POINT
C        D'INTEGRATION L
         CALL REELAS( 3,NOOBSF,23,F1(L),F2(L),0D0,LTDESU(LPYOUN,NOOBSF),
     %                ELAS)
C
         AUX      = POIDEL(L) / F1(L)
         A11(1)   = AUX * ELAS(6) / F1(L)
         A12(1)   = AUX * ELAS(5)
         A12(2)   = AUX * ELAS(9)
         A13(1)   = AUX * ELAS(9)
         A13(2)   = AUX * ELAS(4)
         AUX      = POIDEL(L)
         A22(1,1) = AUX * ELAS(3)
         A22(2,1) = AUX*ELAS(8)
         A22(1,2) = A22(2,1)
         A22(2,2) = AUX * ELAS(10)
         A23(1,1) = A22(2,1)
         A23(2,1) = A22(2,2)
         A23(1,2) = AUX * ELAS(2)
         A23(2,2) = AUX * ELAS(7)
         A33(1)   = A22(2,2)
         A33(2)   = A23(2,2)
         A33(3)   = AUX * ELAS(1)
C
C        CALCUL DE G1
C        ------------
         CALL AB0D(1,1,NBPOLY,A11,POLY(1,L),G2)
         CALL AB1D(1,2,NBPOLY,A12,DP(1,1,L),G2)
         CALL TAB0D(NBPOLY,1,NBPOLY,POLY(1,L),G2,G1)
C
         CALL TAB0D(2,1,NBPOLY,A12,POLY(1,L),G3)
         CALL AB1D(2,2,NBPOLY,A22,DP(1,1,L),G3)
         CALL TAB1D(NBPOLY,2,NBPOLY,DP(1,1,L),G3,G1)
C        ON PLONGE G1 DANS AE
         CALL PLONAD(-1,NBPOLY,NBPOLY,1,1,IP,G1,G1,AE)
C
C        CALCUL DE G2
C        ------------
         CALL TAB0D(NBPOLY,1,2,POLY(1,L),A13,G3)
         CALL TAB1D(NBPOLY,2,2,DP(1,1,L),A23,G3)
         CALL AB0D(NBPOLY,2,NBPOLY,G3,DP(1,1,L),G1)
C        ON PLONGE G2 DANS AE
         CALL PLONAD(-1,NBPOLY,NBPOLY,1,2,IP,G1,G1,AE)
C
C        CALCUL DE G3
C        ------------
         CALL TABA8D(NBPOLY,2,DP(1,1,L),A33,G1,G3)
C        ON PLONGE G3 DANS AE
         CALL PLONAD(1,NBPOLY,NBPOLY,2,2,IP,G1,G1,AE)
C
 10   CONTINUE
C
C     CONTRIBUTION DU COEFFICIENT DU DEPLACEMENT SUR LA SURFACE
C     SOUS LA FORME DE COEFFICIENTS DIAGONAUX DU TYPE DE LA MASSE
C     ===========================================================
C     LA SURFACE SUPPORTE T ELLE UN COEFFICIENT DU DEPLACEMENT?
C     ADRESSE MCN DU COEFFICIENT DU DEPLACEMENT
      MN = LTDESU(LPCOED,NOOBSF)
      IF( MN .GT. 0 ) THEN
         CALL EN2LAG( D2PI,   1,
     &                NOOBSF, NUMISU, NUMASU, LTDESU,
     &                NBPOLY, NPI,    POLY,
     &                F1, F2, POIDEL, COEFDE,
     &                AE )
      ENDIF
C
C     CONTRIBUTION DE LA PENALISATION DES DEPLACEMENTS IMPOSES
C     ========================================================
      IF( PENALI .EQ. 0D0 ) RETURN
C
C     CONTRIBUTION DES ARETES
C     -----------------------
      DO 400 K=1,NBCOTE
C
         IF( NOOBLA(K) .GT. 0 ) THEN
C
            MN = LTDELI(LPFIXA,NOOBLA(K))
            IF( MN .GT. 0 ) THEN
C
C              LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT SUR LE COTE K
               NBCOFI = MCN( MN + WBCOFI )
C
C              LE NUMERO DES POINTS DU COTE K
               NOPOAR(1) = K
               IF( K .NE. NBCOTE ) THEN
                   NOPOAR(2) = K+1
               ELSE
                   NOPOAR(2) = 1
               ENDIF
C              LE NUMERO DU POINT MILIEU
               NOPOAR(3) = K + NBCOTE
C
C              CALCUL DES FIXATIONS PENALISEES SUR LE COTE K
               DO 360 I=1,NBPOLA
C
C                 LE NUMERO DU NOEUD I DE L'ARETE K
                  II = NOPOAR(I)
C
                  DO 350 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     L = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL FIXE
                     L = II * 2 - 2 + L
C
C                    LE NUMERO DU COEFFICIENT DIAGONAL L
                     N = L * ( L + 1 ) / 2
                     AE(N) = AE(N) + PENALI
C
 350              CONTINUE
 360           CONTINUE
            ENDIF
         ENDIF
C
 400  CONTINUE
C
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     -------------------------------------------------------------------
      DO 450 K=1,NBSOMT
C
C        NUMERO DE POINT DU SOMMET K
         N = NOOBPS(K)
         IF( N .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION?
            MN = LTDEPO(LPFIXA,N)
            IF( MN .GT. 0 ) THEN
C
C              OUI: LE NOMBRE DE COMPOSANTES FIXEES
               NBCOFI = MCN( MN + WBCOFI )
C
               DO 440 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  L = MCN( MN + WUCOFI - 1 + J )
C
C                 LES SOMMETS SONT NUMEROTES AVANT LES MILIEUX DES COTES
C                 LE NUMERO DU DL FIXE
                  L = K * 2 - 2 + L
C
C                 LE NUMERO DU COEFFICIENT DIAGONAL L
                  M = L * ( L + 1 ) / 2
                  AE(M) = AE(M) + PENALI
C
 440           CONTINUE
C
            ENDIF
         ENDIF
 450  CONTINUE
C
      RETURN
      END
