      SUBROUTINE TM3LAG( NBJEUX, JEU,    NBPOLY, NPI,    POLY,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   F,      POIDEL,
     %                   CAPA,   NCODSM, CAPAEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE DES EF 3D
C -----    LAGRANGE DE DEGRE 1 OU 2 ISOPARAMETRIQUES
C
C ENTREES:
C --------
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NBPOLY : NOMBRE DE POLYNOMES DE L'ELEMENT FINI VOLUMIQUE
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DANS LE VOLUME
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C           POLY(I, L)= P(I) (XL)
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
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
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
     %                  S
      DOUBLE PRECISION  POLY(NBPOLY, NPI),
     %                  F(NPI,3),
     %                  POIDEL(NPI),
     %                  XYZPI(3),
     %                  PROSCD
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
C
C     CONTRIBUTION DU VOLUME A LA MATRICE DE CAPACITE ELEMENTAIRE
C     ===========================================================
      IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN = (MNTHET-1)/2
         DO I=1,NBPOLY
            DMCN((MNTHDL-1)/2+I)=DMCN(MN+MCN(MNNODL+(I-1)))
         ENDDO
      ENDIF
C
      DO L=1,NPI
C
         IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
         ELSE
            TEMPEL = 0D0
         ENDIF
         ONDEPI = 0D0
C
C        RECHERCHE DE LA CAPACITE AU POINT D'INTEGRATION L
         XYZPI(1) = F(L,1)
         XYZPI(2) = F(L,2)
         XYZPI(3) = F(L,3)
         CALL RECAPA( 4, NOOBVO, 3, XYZPI,
     %                LTDEVO(LPMAST,JEU,NOOBVO),
     %                LTDEVO(LPCHMA,JEU,NOOBVO),
     %                CAPA(L) )
C
C        CAPACITE = CAPACITE * POIDS * DELTA
         CAPA(L) = CAPA(L) * POIDEL(L)
C
      ENDDO
C
C     LE CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE
C     ------------------------------------------------
      N = 0
      DO 40 J=1,NBPOLY
         DO 30 I=1,J
            S = 0.D0
            DO 20 L=1,NPI
               S = S + CAPA(L) * POLY(J,L) * POLY(I,L)
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
