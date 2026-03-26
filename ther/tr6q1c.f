      SUBROUTINE TR6Q1C( NBJEUX, JEU,    NBPOLY, NPI,    POLY,
     &                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     &                   F,      POIDEL, DP,
     &                   COND,   COEFTE, CONDUC)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CONDUCTIVITE D'UN EF 6D
C -----    6CUBE  6Q1C  LAGRANGE DE DEGRE 1 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POLYNOMES DE L'EF VOLUME
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DE L'EF VOLUME
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I,L)=P(I) (XL)
C
C NOOBVO : NUMERO DE L'OBJET VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE
C          DES OBJETS VOLUMES
C
C F      : COORDONNEES XX YY ZZ DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION SUR
C          L'EF COURANT
C
C SORTIES:
C --------
C COND   : TENSEUR SYMETRIQUE  DE CONDUCTIVITE
C COEFTE : COEFFICIENT DE LA TEMPERATURE AUX POINTS D'INTEGRATION
C CONDUC : MATRICE ELEMENTAIRE DE CONDUCTIVITE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      DOUBLE PRECISION  POLY(NBPOLY,NPI),
     %                  F(NPI,6),
     %                  POIDEL(NPI)
      DOUBLE PRECISION  DP(6,NBPOLY,NPI),
     %                  CONDUC(NBPOLY*(NBPOLY+1)/2)
      DOUBLE PRECISION  COND(21), CONDP(6,6), COEFTE(NPI)
      DOUBLE PRECISION  PROSCD, S, XYZPI(6)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
C
C     INITIALISATION A ZERO DE LA MATRICE DE CONDUCTIVITE ELEMENTAIRE
C     ---------------------------------------------------------------
      CALL AZEROD( NBPOLY*(NBPOLY+1)/2, CONDUC )
C
C     ======================
C     CONTRIBUTION DU VOLUME
C     ======================
      IF( TESTNL .GE. 1 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN = (MNTHET-1)/2
         DO 2 I=1,NBPOLY
            DMCN( (MNTHDL-1)/2+I )=DMCN( MN+MCN(MNNODL+I-1) )
 2       CONTINUE
      ENDIF
C
C     SI LA CONDUCTIVITE N'EST PAS DECLAREE, SAUT DU CALCUL DE LA CONDUCTIVITE
      IF( LTDEVO(LPCOND,JEU,NOOBVO) .EQ. 0 ) GOTO 35
C
      DO 30 L=1,NPI
C
         IF( TESTNL .GE. 1 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
         ENDIF
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C        D'INTEGRATION L
         DO 4 I=1,6
            XYZPI(I) = F(L,I)
 4       CONTINUE
         CALL RECOND( 4, NOOBVO, 6, XYZPI, LTDEVO(LPCOND,JEU,NOOBVO),
     %                COND )
C
C        COND = CONDUCTIVITE * (POIDS * DELTA)
         N = 0
         DO 8 I=1,6
            DO 6 J=1,I
               N = N + 1
               CONDP(I,J) = COND(N) * POIDEL(L)
               CONDP(J,I) = COND(N) * POIDEL(L)
 6          CONTINUE
 8       CONTINUE
C
C        CONDUC = CONDUC + T(DP) * CONDP * (DP)
         M = 0
         DO 25 I=1,NBPOLY
            DO 20 J=1,I
               S = 0D0
               DO 10 K=1,6
                  S = S + DP(K,I,L) * ( CONDP(K,1) * DP(1,J,L)
     %                                + CONDP(K,2) * DP(2,J,L)
     %                                + CONDP(K,3) * DP(3,J,L)
     %                                + CONDP(K,4) * DP(4,J,L)
     %                                + CONDP(K,5) * DP(5,J,L)
     %                                + CONDP(K,6) * DP(6,J,L) )
 10            CONTINUE
               M = M + 1
               CONDUC(M) = CONDUC(M) + S
 20         CONTINUE
 25      CONTINUE
C
 30   CONTINUE
C
C     CONTRIBUTION DU COEFFICIENT DE LA TEMPERATURE
 35   IF( LTDEVO(LPCOET,JEU,NOOBVO) .GT. 0 ) THEN
C
         DO 60 L=1,NPI
C
            IF( TESTNL .GE. 1 ) THEN
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL=PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
            ENDIF
C
C           RECHERCHE DU COEFFICIENT DE LA TEMPERATURE AU POINT D'INTEGRATION L
            DO 40 I=1,6
               XYZPI(I) = F(L,I)
 40         CONTINUE
            CALL RECOET( 4, NOOBVO, 6, XYZPI,
     %                   LTDEVO(LPCOET,JEU,NOOBVO), COEFTE(L) )
C
C           COEF TEMPERATURE = COEFTE * POIDS * DELTA
            COEFTE(L) = COEFTE(L) * POIDEL(L)
C
 60      CONTINUE
C
C        SOMMATION AVEC LA MATRICE ELEMENTAIRE
         M = 0
         DO 90 I=1,NBPOLY
            DO 80 J=1,I
               S = 0.D0
               DO 70 L=1,NPI
                  S = S + COEFTE(L) * POLY(I,L) * POLY(J,L)
 70            CONTINUE
               M = M + 1
               CONDUC(M) = CONDUC(M) + S
 80         CONTINUE
 90      CONTINUE
      ENDIF
C
      RETURN
      END
