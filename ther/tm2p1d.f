      SUBROUTINE TM2P1D( X, NBJEUX, JEU, NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   NCODSM, CAPAEF)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CAPACITE CALORIFIQUE
C -----    POUR UN TRIANGLE 2D DE TYPE TRIA 2P1D
C
C ENTREES:
C --------
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C X      : LES 2 COORDONNEES DES 3 POINTS DE L'ELEMENT FINI
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES OBJETS SURFACES
C
C SORTIE :
C --------
C NCODSM : CODE DE STOCKAGE DE LA MATRICE DE CAPACITE ELEMENTAIRE
C          0 CAR DIAGONALE
C CAPAEF : MATRICE DE CAPACITE STOCKEE SOUS FORME DIAGONALE
C          (1,1) (2,2) (3,3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
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
      DOUBLE PRECISION  CAPAEF(3)
      REAL              X(3,2)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
C
      DOUBLE PRECISION  DELTA, XYZPI(3)
C
C     CONTRIBUTION DE LA SURFACE
C     ==========================
C     JACOBIEN * POIDS
      DELTA = ABS(  (X(2,1) - X(1,1)) * (X(3,2) - X(1,2))
     %            - (X(3,1) - X(1,1)) * (X(2,2) - X(1,2))  ) / 6D0
C
      ONDEPI = 0D0
      DO L=1,3
C
         IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C           RECHERCHE DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL = DMCN((MNTHET-1)/2+MCN(MNNODL+L-1))
         ELSE
            TEMPEL = 0D0
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
C        CAPACITE  = CAPACITE  * JACOBIEN * POIDS
         CAPAEF(L) = CAPAEF(L) * DELTA
C
      ENDDO
C
C     MATRICE DE CAPACITE ELEMENTAIRE DIAGONALE
      NCODSM = 0
      RETURN
      END
