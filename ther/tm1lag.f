      SUBROUTINE TM1LAG( NBJEUX, JEU,    NBPOLY, NPI,    POLY,
     &                   NOOBLA, NUMILI, NUMALI, LTDELI,
     &                   F1,     POIDEL,
     &                   CAPA,   NCODSM, CAPAEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE DES EF
C -----    SEGMENT 1D LAGRANGE DE DEGRE 1 OU 2
C
C ENTREES:
C --------
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C NBPOLY : NBRE DE POLYNOMES DE L'ELEMENT LIGNE
C NPI    : NBRE DE POINTS D INTEGRATION NUMERIQUE SUR LA LIGNE
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C           POLY(I,L)= P(I) (XL)
C
C NOOBLA : NUMERO DE L'OBJET LIGNE DE CET ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES OBJETS LIGNES
C F1     : COORDONNEES XX DES NPI POINTS D INTEGRATION DE L'EF
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
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     JUIN 2009
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
      DOUBLE PRECISION  CAPA(NPI), CAPAEF(1:*), S
      DOUBLE PRECISION  POLY(NBPOLY, NPI),
     %                  F1(NPI), POIDEL(NPI), PROSCD, XYZPI(3)
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
C
C     CONTRIBUTION DE LA LIGNE
C     ========================
      IF( TESTNL .GE. 1 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN = (MNTHET-1)/2
         DO 5 I=1,NBPOLY
            DMCN((MNTHDL-1)/2+I) = DMCN( MN + MCN(MNNODL+I-1) )
 5       CONTINUE
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
         XYZPI(1) = F1(L)
         XYZPI(2) = 0D0
         XYZPI(3) = 0D0
         CALL RECAPA( 2, NOOBLA, 3, XYZPI,
     %                LTDELI(LPMAST,JEU,NOOBLA),
     %                LTDELI(LPCHMA,JEU,NOOBLA),
     %                CAPA(L) )
C
C        CAPACITE = CAPACITE * POIDS * DELTA
         CAPA(L) = CAPA(L) *  POIDEL(L)
C
 10   CONTINUE
C
C     LE CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE
C     ------------------------------------------------
      N = 0
      DO 40 J=1,NBPOLY
         DO 30 I=1,J
            S = 0.D0
            DO 20 L=1,NPI
               S = S +  CAPA(L) * POLY(J,L) * POLY(I,L)
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
