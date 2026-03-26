      SUBROUTINE EM2Q1C( X,      NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NCODSM, MASSEF)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE MASSE (DIAGONALE)
C -----    POUR UN QUADRANGLE 2D DE TYPE QUAD 2Q1C
C
C ENTREES:
C --------
C X      : LES 2 COORDONNEES DES 4 SOMMETS DU QUADRANGLE
C NOOBSF : NUMERO DE LA SURFACE DE CE QUADRANGLE
C NUMISU : NUMERO MINIMAL DES SURFACES DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES MASSE
C          DES SURFACES DE L'OBJET
C
C SORTIES:
C --------
C NCODSM : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE DE MASSE
C          0 DIAGONALE
C MASSEF : MATRICE DE MASSE DIAGONALE (1,1) (2,2) ... (8,8)
C          RANGEE PAR NOEUDS ( U1S1,U2S1, U1S2,U2S2, ... , U1S4,U2S4 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1998
C23456---------------------------------------------------------------012
      include"./incl/donela.inc"
      DOUBLE PRECISION  MASSEF(8), RHO
      REAL              X(4, 2)
      INTEGER           LTDESU( 1:MXDOEL, NUMISU:NUMASU )
      DOUBLE PRECISION  DELTA, S1234(2)
      DOUBLE PRECISION  S(2,4), XYZPI(3)
      DATA              S/ 0D0,0D0, 1D0,0D0, 1D0,1D0, 0D0,1D0 /
C
C     CALCUL DES COEFFICIENTS DIAGONAUX DE LA MASSE ELEMENTAIRE
C     =========================================================
      S1234(1) = X(1,1) - X(2,1) + X(3,1) - X(4,1)
      S1234(2) = X(1,2) - X(2,2) + X(3,2) - X(4,2)
C
      NC = 1
      DO 10 L=1,4
C
C        RECHERCHE DE LA DENSITE SURFACIQUE DE MASSE
C        AU POINT D'INTEGRATION L
         XYZPI(1) = X(L,1)
         XYZPI(2) = X(L,2)
         XYZPI(3) = 0D0
        CALL REMASS( 3, NOOBSF, 3, XYZPI,
     %                LTDESU(LPMASS,NOOBSF), RHO )
C
C        LE JACOBIEN AU SOMMET L
         DELTA = ( X(2,1) - X(1,1) + S(2,L) * S1234(1) ) *
     %           ( X(4,2) - X(1,2) + S(1,L) * S1234(2) ) -
     %           ( X(2,2) - X(1,2) + S(2,L) * S1234(2) ) *
     %           ( X(4,1) - X(1,1) + S(1,L) * S1234(1) )
C
C        MASSE = MASSE * POIDS * DELTA
         MASSEF(NC  ) = RHO * DELTA * 0.25D0
         MASSEF(NC+1) = MASSEF(NC)
         NC = NC + 2
C
 10   CONTINUE
C
C     CODE MATRICE ELEMENTAIRE DE MASSE DIAGONALE
      NCODSM = 0
      RETURN
      END
