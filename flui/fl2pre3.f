      SUBROUTINE FL2PRE3( NBPOE,  X,      NONOEF,
     %                    NTDLVP, NBVECT, VXYZPR, PREMIN, NDDLNV,
     %                    VOLTET, INTPRE,
     %                    NOFOPR, NBNOPR, PREEXA, INTPER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE DES NBVECT (PRESSION-PRESSION MIN)**2
C -----           DE L'INTEGRALE DU VOLUME DU TETRAEDRE
C          SUR UN TETRAEDRE DE PRESSION POLYNOME DE LAGRANGE DE DEGRE 1
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C X      : 3 COORDONNEES DES NBPOE NOEUDS DU TETRAEDRE
C NONOEF : NUMEROS DES NBPOE NOEUDS DU TETRAEDRE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYZPR : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C PREMIN : PRESSION MINIMALE DE CHACUN DES NBVECT VECTEURS SOLUTION
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C NOFOPR : NUMERO DE LA FONCTION UTILISATEUR PRESSION_EXACTE(t,x,y,z)
C          0 SI ELLE N'EXISTE PAS => INTPER N'EST PAS CALCULE
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C PREEXA : PREEXA(NBNOPR,NBVECT) PRESSION EXACTE EN CHAQUE NOEUD TEMPS
C
C SORTIES:
C --------
C VOLTET : VOLUME DU TETRAEDRE
C                                                         (Pe1)
C INTPRE : Som  Jacobien (Pe1,Pe2,Pe3,Pe4) [int Pi Pj DX] (Pe2)
C          e de E                                         (Pe3)
C                                                         (Pe4)
C                                                            (P1E-P1C)
C INTPER : Som Jacobien (P1E-P1C,...,P4E-P4C) [int Pi Pj DX] (  ...  )
C          e de E                                            (P4E-P4C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC &SAINT PIERRE DU PERRAY JANVIER 2012
C23456---------------------------------------------------------------012
      INTEGER           NONOEF(NBPOE), NDDLNV(0:*)
      REAL              X(NBPOE,3)
      DOUBLE PRECISION  VXYZPR(NTDLVP,NBVECT), PREEXA(NBNOPR,NBVECT),
     %                  INTPRE(NBVECT), INTPER(NBVECT), PREMIN(NBVECT),
     %                  VOLTET, INTPEF, INTPERR,
     %                  PresExI, PresExJ, PREI, PREJ, PREMI
      DOUBLE PRECISION  X12, Y12, Z12, X13, Y13, Z13, X14, Y14, Z14
      DOUBLE PRECISION  D, E, F
      INTRINSIC         SQRT
C
      X12 = X(1,1) - X(2,1)
      Y12 = X(1,2) - X(2,2)
      Z12 = X(1,3) - X(2,3)
C
      X13 = X(1,1) - X(3,1)
      Y13 = X(1,2) - X(3,2)
      Z13 = X(1,3) - X(3,3)
C
      X14 = X(1,1) - X(4,1)
      Y14 = X(1,2) - X(4,2)
      Z14 = X(1,3) - X(4,3)
C
      D = Y14 * Z13 - Y13 * Z14
      E = Y12 * Z14 - Y14 * Z12
      F = Y13 * Z12 - Y12 * Z13
C
C     6 x (VOLUME DU TETRAEDRE) = JACOBIEN
      VOLTET = ABS( X12 * D + X13 * E + X14 * F )
C
C                                               1/60   1/120  1/120 1/120
C     COEFFICIENTS DE LA MATRICE Int Pi Pj dX = 1/120  1/60   1/120 1/120
C                                               1/120  1/120  1/60  1/120
C                                               1/120  1/120  1/120 1/60
C     => INTEGRALE EXACTE POUR Pi COORDONNEES BARYCENTRIQUES
C
      DO NV=1, NBVECT
C
C        PRESSION MINIMALE DU VECTEUR NV
         PREMI   = PREMIN(NV)
C
         INTPEF  = 0D0
         INTPERR = 0D0
C
         DO J=1,4
C
C           NUMERO DU SOMMET J DU TETRAEDRE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE PRESSION AU NOEUD NS
            NDLJ = NDDLNV( NSJ )
C
C           PRESSION AU SOMMET J
            PREJ = VXYZPR(NDLJ,NV)
C
C           LA PRESSION_EXACTE(t,x,y,z) AU NOEUD J
            IF( NOFOPR .GT. 0 ) PresExJ = PREEXA( NSJ, NV )
C
            DO I=1,4
C
C              NUMERO DU SOMMET I DU TETRAEDRE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE PRESSION AU NOEUD NS
               NDLI = NDDLNV( NSI )
C
C              PRESSION AU SOMMET I
               PREI = VXYZPR(NDLI,NV)
C
C              NORME L1 DE PRESSION**2 INTEGREE EXACTEMENT
               IF( I .EQ. J ) THEN
C
C                 CARRE NORME L2 DE LA PRESSION SUR L'EF INTEGREE EXACTEMENT
                  INTPEF = INTPEF + ( (PREJ-PREMI)**2 ) * 2D0
C
                  IF( NOFOPR .GT. 0 ) THEN
C                    PRESSION_EXACTE(temps,x,y,z) AU NOEUD J DEJA CALCULE
C                    ERREUR Pression = PressionExacte - PressionCalculee
                     INTPERR = INTPERR + ( (PresExJ-PREJ)**2 ) * 2D0
                  ENDIF
C
               ELSE
C
C                 CARRE NORME L2 DE LA PRESSION SUR L'EF INTEGREE EXACTEMENT
                  INTPEF = INTPEF + (PREI-PREMI) * (PREJ-PREMI)
C
                  IF( NOFOPR .GT. 0 ) THEN
C                    LA PRESSION_EXACTE(t,x,y,z) AU NOEUD J
                     PresExI = PREEXA( NSI, NV )
C                    ERREUR Pression = PressionExacte - PressionCalculee
                     INTPERR = INTPERR + (PresExI-PREI) * (PresExJ-PREJ)
                  ENDIF
C
               ENDIF
C
            ENDDO
C
         ENDDO
C
C        NORME L1 DE PRESSION**2 INTEGREE EXACTEMENT => NORME L2 AU CARRE
         INTPRE(NV) = INTPRE(NV) + INTPEF  * VOLTET / 120D0
         INTPER(NV) = INTPER(NV) + INTPERR * VOLTET / 120D0
C
      ENDDO
C
C     VOLUME DU TETRAEDRE
      VOLTET = VOLTET / 6D0
C
      RETURN
      END
