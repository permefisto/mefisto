      SUBROUTINE TM6Q1C( NBJEUX, JEU,    NBPOLY, NPI,    POLY,
     &                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     &                   F,      POIDEL,
     &                   CAPA,   NCODSM, CAPAEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE D'UN EF 6D
C -----    6CUBE  6Q1C  LAGRANGE DE DEGRE 1 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POLYNOMES DE L'ELEMENT FINI
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DANS LE VOLUME
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I, L)= P(I) (XL)
C
C NOOBVO : NUMERO DU VOLUME DE L'EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES UTILISES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES UTILISES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES OBJETS VOLUMES
C F      : COORDONNEES XX YY ZZ DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D'INTEGRATION DE L'EF
C
C SORTIES:
C --------
C CAPA   : CAPACITE EN CHAQUE POINT D'INTEGRATION DE L'EF
C NCODSM : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE DE CAPACITE
C          1 CAR SYMETRIQUE PLEINE
C CAPAEF : MATRICE DE CAPACITE STOCKEE   1
C                                        2 3
C                                        4 5 6
C                                        7 ...
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  CAPA(NPI),
     %                  CAPAEF(1),
     %                  S,
     %                  XYZPI(6),
     %                  PROSCD
      DOUBLE PRECISION  POLY(NBPOLY, NPI),
     %                  F(NPI,6),
     %                  POIDEL(NPI)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
C
C     CONTRIBUTION DU VOLUME A LA MATRICE DE CAPACITE ELEMENTAIRE
C     ===========================================================
      IF( TESTNL .GE. 1 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN = (MNTHET-1)/2
         DO 2 I=1,NBPOLY
            DMCN((MNTHDL-1)/2+I)=DMCN(MN+MCN(MNNODL+(I-1)))
 2       CONTINUE
      ENDIF
C
      DO 10 L=1,NPI
C
         IF( TESTNL .GE. 1 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
         ENDIF
C
C        RECHERCHE DE LA CAPACITE AU POINT D'INTEGRATION L
         DO 5 I=1,6
            XYZPI(I) = F(L,I)
 5       CONTINUE
         CALL RECAPA( 4, NOOBVO, 6, XYZPI,
     %                LTDEVO(LPMAST,JEU,NOOBVO),
     %                LTDEVO(LPCHMA,JEU,NOOBVO),
     %                CAPA(L) )
C
C        CAPACITE = CAPACITE * POIDS * DELTA
         CAPA(L) = CAPA(L) * POIDEL(L)
C
 10   CONTINUE
C
C     LE CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE
C     ------------------------------------------------
      N = 0
      DO 40 I=1,NBPOLY
         DO 30 J=1,I
            S = 0.D0
            DO 20 L=1,NPI
               S = S + CAPA(L) * POLY(I,L) * POLY(J,L)
 20         CONTINUE
            N = N + 1
            CAPAEF(N) = S
 30      CONTINUE
 40   CONTINUE
C
C     MATRICE ELEMENTAIRE DE CAPACITE SYMETRIQUE PLEINE
      NCODSM = 1
      RETURN
      END
