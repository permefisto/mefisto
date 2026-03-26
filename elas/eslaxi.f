      SUBROUTINE ESLAXI( D2PI,   X,      PENALI, NBSOMT, NBCOTE, NDSM,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NBPOLA, NPIA,   POIDSA, POLYA,  DPOLYA,
     %                   NBPOLY, NPI,    POLY,
     %                   NUELEM, NBELEM, NUNDEL,
     %                   MNTEMP, NTDLT,  TEMPER,
     %                   ELAS,   FOMEGA, FGAMMA, G1,G2,  F1,F2,
     %                   POIDEL, DP,     IP,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES SECONDS MEMBRES DES ELEMENTS AXISYMETRIQUES
C -----    DE DEGRE 1 OU 2 SUR UN TRIANGLE OU UN QUADRANGLE
C
C
C ENTREES:
C --------
C D2PI   : 2 FOIS PI
C X      : 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C PENALI : 1/EPSILON DE LA PENALISATION DU DEPLACEMENT IMPOSE
C          LA VALEUR DOIT ETRE SUFFISAMMENT GRANDE POUR ECRASER CELLES
C          DU VECTEUR ELEMENTAIRE ET GLOBAL
C
C NBSOMT : NOMBRE DE SOMMETS DE L ELEMENT FINI
C NBCOTE : NOMBRE DES COTES  DE L ELEMENT FINI
C NDSM   : NOMBRE DE SECONDS MEMBRES
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL  DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL  DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES POINTS DE L'OBJET
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES LIGNES DE L'OBJET
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES SURFACES DE L'OBJET
C
C NBPOLA : NOMBRE DE POLYNOMES DE BASE SUR UN COTE DE L'ELEMENT FINI
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR UN COTE
C POIDSA : POIDS DES POINTS D INTEGRATION SUR UN COTE
C POLYA  : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION COTE
C DPOLYA : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBPOLY : NOMBRE DE POLYNOMES DE L'ELEMENT FINI COMPLET
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE
C POIDS  : LES NPI POIDS DE LA FORMULE D INTEGRATION
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I,  L)= P(I) (XL)
C DPOLY  : IDEM POUR LES DERIVEES DE CES POLYNOMES
C          DPOLY(I,J,L)=DP(J)/DX(I) (XL)
C
C NUELEM : NUMERO DE L'ELEMENT FINI
C NBELEM : NOMBRE D'ELEMENTS FINIS
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C MNTEMP : >0 CALCUL DEMANDE DES CONTRAINTES THERMIQUES,=<0 SINON
C NTDLT  : NOMBRE DE NOEUDS THERMIQUES OU DL THERMIQUES DU MAILLAGE
C TEMPER : TEMPERATURES AUX NOEUDS DU MAILLAGE TEMPER(NTDLT,NDSM)
C ELAS   : TENSEUR DE L ELASTICITE SYMETRIQUE CALCULE DANS CE SP
C FOMEGA : EFFORTS DE VOLUME CONSTANT PAR SOUS-DOMAINE (NDSM,2)
C          TABLEAU REMPLI DANS CE SOUS PROGRAMME
C FGAMMA  : PRESSION OU CONTRAINTE INITIALE DECLAREE FGAMMA(NDSM,3)
C          TABLEAU REMPLI DANS CE SP
C
C G1,G2  : TABLEAUX AUXILIAIRES
C F1     : RAYON R DES NPI POINTS D'INTEGRATION
C F2     : COTE  Z DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI)  DES NPI POINTS D INTEGRATION
C DP     : DP(2,NBPOLY,NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS DE BASE ISOPARAMETRIQUES
C IP     : IP(J)=POSITION DU J-EME D.L. COMPOSANTE PAR COMPOSANTE DANS L
C          NUMEROTATION NOEUD PAR NOEUD
C
C SORTIES:
C --------
C BE     : BE(NDSM,2*NBPOLY) LES NDSM SECONDS MEMBRES ELEMENTAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1999
C23456---------------------------------------------------------------012
      include"./incl/donela.inc"
      include"./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      REAL              X(NBPOLY,2)
      DOUBLE PRECISION  D2PI, PENALI
      DOUBLE PRECISION  POIDSA(NPIA),
     %                  POLYA(NBPOLA,NPIA),
     %                  DPOLYA(NBPOLA,NPIA),
     %                  POLY(NBPOLY,NPI),
     %                  POIDEL(NPI),
     %                  DP(2,NBPOLY,NPI),
     %                  F1(NPI),
     %                  F2(NPI),
     %                  FOMEGA(NDSM,2),
     %                  FGAMMA(NDSM,3),
     %                  TEMPER(NTDLT,NDSM),
     %                  ELAS(10),
     %                  A(5),
     %                  G1(NBPOLY,2),
     %                  G2(NBPOLY,2,NBPOLY),
     %                  BE(NDSM,*)
C
      INTEGER           NOOBPS(NBSOMT), NOOBLA(NBCOTE), NOOBSF
      INTEGER           LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
C
      INTEGER           IP(*),NUNDEL(1:NBELEM,1:NBPOLY)
      INTEGER           NOPOAR(3)
      DOUBLE PRECISION  D,DD,GL(2),DGL(2),DELTA,XD,YD
      DOUBLE PRECISION  TEMPIN,DILATA
C
C     INITIALISATION DE BE
      CALL AZEROD(NDSM*2*NBPOLY,BE)
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES
C     ------------------------------------
      IF( LTDESU(LPFORC,NOOBSF) .GT. 0 ) THEN
         DO 15 L=1,NPI
C           LA VALEUR DES EFFORTS SURFACIQUES EN CE POINT
            CALL REFORC(3,NOOBSF,2, F1(L),F2(L),0D0, 0D0,0D0,0D0,
     %                  LTDESU(LPFORC,NOOBSF),FOMEGA)
            DO 10 I=1,NBPOLY
               DO 9 J=1,2
                  I1 = IP(I + NBPOLY * (J-1))
                  DO 8 N=1,NDSM
                     BE(N,I1) = BE(N,I1) +
     %                          POIDEL(L)*POLY(I,L)*FOMEGA(N,J)
    8             CONTINUE
    9          CONTINUE
   10       CONTINUE
   15    CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES EFFORTS SUR LES COTES
C     ----------------------------------------
      DO 50 K=1,NBCOTE
         IF( NOOBLA(K) .GT. 0 ) THEN
            IF( LTDELI(LPFORC,NOOBLA(K)) .GT. 0 ) THEN
C              LE NUMERO DES POINTS DU COTE
               NOPOAR(1) = K
               IF( K .NE. NBCOTE ) THEN
                  NOPOAR(2) = K+1
               ELSE
                  NOPOAR(2) = 1
               ENDIF
C              LE NUMERO DU POINT MILIEU
               NOPOAR(3) = K + NBCOTE
               DO 40 L=1,NPIA
C                 CALCUL DES COORDONNEES DU POINT D INTEGRATION
C                 ET DU JACOBIEN
                  CALL E22LAG( NBPOLY,NBPOLA,NOPOAR,
     %                        POLYA(1,L),DPOLYA(1,L),
     %                        X,GL,DGL,DELTA)
C                 EN SORTIE GL=LES 2 COORDONNEES DU POINT D'INTEGRATION
C
C                 LE VECTEUR NORMAL UNITAIRE
                  A(1) =  DGL(2) / DELTA
                  A(2) = -DGL(1) / DELTA
C
C                 CALCUL DES FORCES EN CE POINT DU COTE K
                  CALL REFORC( 2,NOOBLA(K),2,
     %                         GL(1),GL(2),0D0, A(1),A(2),0D0,
     %                         LTDELI(LPFORC,NOOBLA(K)), FGAMMA)
C
                  D = DELTA * POIDSA(L) * D2PI * GL(1)
                  DO 30 I=1,NBPOLA
                     II = NOPOAR(I)
                     I1 = IP(II)
                     II = IP(II+NBPOLY)
                     DD = D * POLYA(I,L)
                     DO 20 N=1,NDSM
                        BE(N,I1) = BE(N,I1) + DD * FGAMMA(N,1)
                        BE(N,II) = BE(N,II) + DD * FGAMMA(N,2)
 20                  CONTINUE
 30               CONTINUE
 40            CONTINUE
            ENDIF
         ENDIF
 50   CONTINUE
C
C     CONTRIBUTION DES CONTRAINTES INITIALES
C     --------------------------------------
      IF( LTDESU(LPCOIN,NOOBSF) .GT. 0 ) THEN
         DO 85 L=1,NPI
C           LA VALEUR DES CONTRAINTES INITIALES EN CE POINT
            CALL RECOIN( 3,NOOBSF,4,F1(L),F2(L),0.D0,
     %                   LTDESU(LPCOIN,NOOBSF), A )
C           A =  SZZ SRR STETA SRZ
C
C           T(D) * A
            A(5) = A(1)
            A(1) = A(3) / F1(L)
            A(3) = A(4)
C
C           CALCUL DU TABLEAU TP * U/R = T(P) A(1)
            CALL TAB0D(NBPOLY,1,1,POLY(1,L),A(1),G1(1,1))
C
C           CALCUL DU TABLEAU T(DP)*T(S11 S12)
C           CALCUL DE G1(NBPOLY,1) = T(DP) T(A(2),A(3))
            CALL TAB1D(NBPOLY,2,1,DP(1,1,L),A(2),G1(1,1))
C
C           CALCUL DU TABLEAU T(DP)*T(S12 S22)
C           CALCUL DE G1(NBPOLY,2) = T(DP) T(A(4),A(5))
            CALL TAB0D(NBPOLY,2,1,DP(1,1,L),A(4),G1(1,2))
C
C           AJOUT AUX SECONDS MEMBRES ELEMENTAIRES
            DO 80 I=1,NBPOLY
               DO 75 J=1,2
                  I1 = IP(I + NBPOLY * (J-1))
C                 LA COMPOSANTE DES CONTRAINTES INITIALES
                  D  = G1(I,J)
                  DO 70 N=1,NDSM
                     BE(N,I1) = BE(N,I1) - D
   70             CONTINUE
   75          CONTINUE
   80       CONTINUE
   85    CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES CONTRAINTES THERMIQUES
C     ----------------------------------------
      IF( MNTEMP .LE. 0 ) RETURN
      IF( LTDESU(LPDILA,NOOBSF) .LE. 0 ) RETURN
C
C     LE VECTEUR DE DILATATION THERMIQUE VAUT ICI (DILATA,DILATA,DILATA,0.)
      CALL AZEROD(NBPOLY*2*NBPOLY,G2)
      DO 150 L=1,NPI
C
C        LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
         CALL REELAS( 3,NOOBSF,23,F1(L),F2(L),0D0,
     %                LTDESU(LPYOUN,NOOBSF), ELAS )
C
C        LE COEFFICIENT DE DILATATION THERMIQUE
         CALL REDILA( 3,NOOBSF,F1(L),F2(L),0D0,LTDESU(LPDILA,NOOBSF),
     %                DILATA,TEMPIN )
C
C        CALCUL DE A = T(D) (E) (DILATA) (PT BL)
         D = POIDEL(L) * DILATA
         A(1) = D * ( ELAS(4) + ELAS(5) + ELAS(6) ) / F1(L)
         A(2) = D * ( ELAS(2) + ELAS(3) + ELAS(5) )
         A(3) = D * ( ELAS(7) + ELAS(8) + ELAS(9) )
         A(4) = D * ( ELAS(1) + ELAS(2) + ELAS(4) )
C
C        CALCUL DE G1(NBPOLY,1) = T(P) (A(1))
         CALL TAB0D(NBPOLY,1,1,POLY(1,L),A(1),G1(1,1) )
C
C        CALCUL DE G1(NBPOLY,1) = T(P) (A(1)) + T(DP) (A(2),A(3))
         CALL TAB1D(NBPOLY,2,1,DP(1,1,L),A(2),G1(1,1) )
C
C        CALCUL DE G1(NBPOLY,2) = T(DP) (A(3),A(4))
         CALL TAB0D(NBPOLY,2,1,DP(1,1,L),A(3),G1(1,2) )
C
C        CALCUL DE G2 = (G1) (P)
         CALL AB1D(NBPOLY*2,1,NBPOLY,G1,POLY(1,L),G2  )
 150  CONTINUE
C
C     LE SECOND MEMBRE DU AUX CONTRAINTES THERMIQUES
      DO 200 I=1,NBPOLY
         DO 180 J=1,2
            I1 = IP( I + NBPOLY * ( J - 1 ) )
            DO 170 N=1,NDSM
               D = 0D0
               DO 160 L=1,NBPOLY
                  D = D + G2(I,J,L)*(TEMPER(NUNDEL(NUELEM,L),N)-TEMPIN)
 160           CONTINUE
               BE(N,I1) = BE(N,I1) + D
 170        CONTINUE
 180     CONTINUE
 200  CONTINUE
C
      IF( PENALI .EQ. 0D0 ) RETURN
C
C     CONTRIBUTION DES ARETES A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     ------------------------------------------------------------------
      DO 400 K=1,NBCOTE
C
          IF( NOOBLA(K) .GT. 0 ) THEN
C
              MN = LTDELI(LPFIXA,NOOBLA(K))
              IF( MN .GT. 0 ) THEN
C
C                LE NUMERO DES POINTS DU COTE K
                 NOPOAR(1) = K
                 IF( K .NE. NBCOTE ) THEN
                    NOPOAR(2) = K+1
                 ELSE
                    NOPOAR(2) = 1
                 ENDIF
C                LE NUMERO DU POINT MILIEU
                 NOPOAR(3) = K + NBCOTE
C
C                CALCUL DES FIXATIONS PENALISEES SUR LE COTE K
                 DO 360 I=1,NBPOLA
C
C                   LE NUMERO DU NOEUD I DE L'ARETE K
                    II = NOPOAR(I)
                    XD = X(II,1)
                    YD = X(II,2)
                    CALL REFIXA( 2, NOOBLA(K),  XD, YD, 0D0, MN,
     &                           NBCOFI, FGAMMA )
C                   NBCOFI LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT SUR LE
C
                    DO 350 J=1,NBCOFI
C
C                      LE NUMERO DE LA COMPOSANTE FIXEE
                       L = MCN( MN + WUCOFI - 1 + J )
C
C                      LE NUMERO DU DL FIXE
                       L = II * 2 - 2 + L
C
                       DO 340 N=1,NDSM
                          BE(N,L) = BE(N,L) + PENALI * FGAMMA(1,J)
 340                   CONTINUE
C
 350                CONTINUE
 360             CONTINUE
             ENDIF
          ENDIF
C
 400  CONTINUE
C
C
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     -------------------------------------------------------------------
      DO 450 K=1,NBSOMT
C
C        NUMERO DE POINT DU SOMMET K
         NL = NOOBPS(K)
         IF( NL .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION?
            MN = LTDEPO( LPFIXA, NL )
            IF( MN .GT. 0 ) THEN
C
C              CALCUL DE LA FIXATION AU SOMMET K
               XD = X(K,1)
               YD = X(K,2)
               CALL REFIXA( 1, NL, XD, YD, 0.D0, MN, NBCOFI, FGAMMA )
C
               DO 440 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  L = MCN( MN + WUCOFI - 1 + J )
C
C                 LA PENALISATION EST AJOUTEE
                  L = K * 2 - 2 + L
C
                  DO 435 N=1,NDSM
                     BE(N,I) = BE(N,I) + PENALI * FGAMMA(1,J)
 435              CONTINUE
C
 440           CONTINUE
C
            ENDIF
         ENDIF
 450  CONTINUE
C
      RETURN
      END
