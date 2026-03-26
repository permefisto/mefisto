      SUBROUTINE TM2LAG( D2PI,   NOAXIS, NBJEUX, JEU,
     %                   NBPOLY, NPI,    POLY,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   F1, F2, POIDEL,
     %                   CAPA,   NCODSM, CAPAEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE  Rho Cp
C -----    ( ou encore NOMMEE CAPACITE THERMIQUE )
C          DES EF 2D de LAGRANGE DE DEGRE 1 OU 2 ISOPARAMETRIQUES
C
C ENTREES:
C --------
C D2PI   : 2 FOIS PI
C NOAXIS : 1 SI PROBLEME AXISYMETRIQUE, 0 SINON
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NBPOLY : NBRE DE POLYNOMES DE L'ELEMENT SURFACE
C NPI    : NBRE DE POINTS D INTEGRATION NUMERIQUE SUR LA SURFACE
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C           POLY(I,L)= P(I) (XL)
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES OBJETS SURFACES
C F1     : COORDONNEES XX DES NPI POINTS D INTEGRATION DE L'EF
C F2     : COORDONNEES YY DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C
C SORTIES:
C --------
C CAPA   : CAPACITE EN CHAQUE POINT D'INTEGRATION
C NCODSM : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE DE CAPACITE
C          1 CAR SYMETRIQUE PLEINE
C CAPAEF : MATRICE DE CAPACITE STOCKEE   1
C                                        2 3
C                                        4 5 6
C                                        7 ...
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1993
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
      DOUBLE PRECISION  D2PI, DELTA, CAPA(NPI), CAPAEF(1:*), S
      DOUBLE PRECISION  POLY(NBPOLY, NPI),
     %                  F1(NPI), F2(NPI), POIDEL(NPI), PROSCD, XYZPI(3)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
C
C     CONTRIBUTION DE LA SURFACE
C     ==========================
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
C        RECHERCHE DE LA CAPACITE THERMIQUE MASSIQUE AU POINT D'INTEGRATION L
         XYZPI(1) = F1(L)
         XYZPI(2) = F2(L)
         XYZPI(3) = 0D0
         CALL RECAPA( 3, NOOBSF, 3, XYZPI,
     %                LTDESU(LPMAST,JEU,NOOBSF),
     %                LTDESU(LPCHMA,JEU,NOOBSF),
     %                CAPA(L) )
C
         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            DELTA = POIDEL(L) * D2PI * F1(L)
         ELSE
C           EF NON AXISYMETRIQUE
            DELTA = POIDEL(L)
         ENDIF
C
C        CAPACITE = CAPACITE * POIDS * DELTA
         CAPA(L) = CAPA(L) * DELTA
C
      ENDDO
C
C     LE CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE
C     ------------------------------------------------
      N = 0
      DO J=1,NBPOLY
         DO I=1,J
            S = 0.D0
            DO L=1,NPI
               S = S +  CAPA(L) * POLY(J,L) * POLY(I,L)
            ENDDO
            N = N + 1
            CAPAEF(N) = S
         ENDDO
      ENDDO
C
C     MATRICE ELEMENTAIRE DE CAPACITE SYMETRIQUE PLEINE
      NCODSM = 1
      RETURN
      END
