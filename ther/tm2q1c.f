      SUBROUTINE TM2Q1C( X,      NBJEUX, JEU,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   NCODSM, CAPAEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DIAGONALE DE CAPACITE CALORIFIQUE
C -----    POUR UN QUADRANGLE 2D DE TYPE QUAD 2Q1C ET UNE INTEGRATION
C          NUMERIQUE AUX 4 SOMMETS DU QUADRANGLE
C
C ENTREES:
C --------
C X      : LES 2 COORDONNEES DES 4 SOMMETS DE L'ELEMENT FINI
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C NOOBSF : NUMERO DE SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES SURFACES DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES SURFACES DE L'OBJET
C
C SORTIE :
C --------
C NCODSM : CODE DE STOCKAGE DE LA MATRICE DE CAPACITE ELEMENTAIRE
C          0 CAR DIAGONALE
C CAPAEF : MATRICE DE CAPACITE STOCKEE SOUS FORME DIAGONALE
C          (1,1) (2,2) (3,3) (4,4)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1998
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), DMCN(1))
C
      DOUBLE PRECISION  CAPAEF(4)
      REAL              X(4,2)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
C
      DOUBLE PRECISION  S(2,4)
      DATA              S/ 0D0,0D0, 1D0,0D0, 1D0,1D0, 0D0,1D0 /
C
      DOUBLE PRECISION  DELTA, XYZPI(3), X1234, Y1234
C
C     CONTRIBUTION DE LA SURFACE
C     ==========================
      X1234 = X(1,1) - X(2,1) + X(3,1) - X(4,1)
      Y1234 = X(1,2) - X(2,2) + X(3,2) - X(4,2)
C
      DO 10 L=1,4
C
         IF( TESTNL .GE. 1 ) THEN
C           RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=DMCN((MNTHET-1)/2+MCN(MNNODL+L-1))
         ENDIF
C
C        RECHERCHE DE LA CAPACITE AU POINT D'INTEGRATION L
         XYZPI(1) = X(L,1)
         XYZPI(2) = X(L,2)
         XYZPI(3) = 0D0
         CALL RECAPA( 3, NOOBSF, 3, XYZPI,
     %                LTDESU(LPMAST,JEU,NOOBSF),
     %                LTDESU(LPCHMA,JEU,NOOBSF),
     %                CAPAEF(L) )
C
C        LE JACOBIEN AU SOMMET L
         DELTA = ( X(2,1) - X(1,1) + S(2,L) * X1234 ) *
     %           ( X(4,2) - X(1,2) + S(1,L) * Y1234 ) -
     %           ( X(2,2) - X(1,2) + S(2,L) * Y1234 ) *
     %           ( X(4,1) - X(1,1) + S(1,L) * X1234 )
C
C        CAPACITE  = CAPACITE  * ABS(DELTA) * POIDS
         CAPAEF(L) = CAPAEF(L) * ABS(DELTA) * 0.25D0
C
 10   CONTINUE
C
C     MATRICE DE CAPACITE ELEMENTAIRE DIAGONALE
      NCODSM = 0
      RETURN
      END
