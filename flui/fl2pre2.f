      SUBROUTINE FL2PRE2( NBPOE,  X,      NONOEF,
     %                    NTDLVP, NBVECT, VXYPR,  PREMIN, NDDLNV,
     %                    SURFTR, INTPRE,
     %                    NOFOPR, NBNOPR, PREEXA, INTPER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE DES NBVECT (PRESSION-PRESSION MIN)**2
C -----           DE L'INTEGRALE DE LA SURFACE DU TRIANGLE
C          SUR UN TRIANGLE DE PRESSION (POLYNOME DE LAGRANGE DE DEGRE 1)
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C X      : 3 COORDONNEES DES NBPOE NOEUDS DU TRIANGLE
C NONOEF : NUMEROS DES NBPOE NOEUDS DU TRIANGLE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYPR  : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C PREMIN : PRESSION MINIMALE DE CHACUN DES NBVECT VECTEURS SOLUTION
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE   0:NBNOVI
C NOFOPR : NUMERO DE LA FONCTION UTILISATEUR PRESSION_EXACTE(t,x,y,z)
C          0 SI ELLE N'EXISTE PAS => INTPER N'EST PAS CALCULE
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C PREEXA : PREEXA(NBNOPR,NBVECT) PRESSION EXACTE EN CHAQUE NOEUD TEMPS
C
C SORTIES:
C --------
C SURFTR : SURFACE DU TRIANGLE
C                                                      (Pe1)
C INTPRE : Som   Jacobien (Pe1,Pe2,Pe3) [int Pi Pj DX] (Pe2)
C          e de E                                      (Pe3)
C                                                                (P1E-P1C)
C INTPER : Som Jacobien (P1E-P1C,P2E-P2C,P3E-P3C) [int Pi Pj DX] (P2E-P2C)
C          e de E                                                (P3E-P3C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC &SAINT PIERRE DU PERRAY JANVIER 2012
C23456---------------------------------------------------------------012
      INTEGER           NONOEF(NBPOE), NDDLNV(0:*)
      REAL              X(NBPOE,2)
      DOUBLE PRECISION  VXYPR(NTDLVP,NBVECT), PREEXA(NBNOPR,NBVECT),
     %                  PREMIN(NBVECT), INTPRE(NBVECT),  INTPER(NBVECT),
     %                  SURFTR,  INTPEF, INTPERR,
     %                  PresExI, PresExJ, PREMI, PR
C
C     2 x (SURFACE DU TRIANGLE) = JACOBIEN
      SURFTR = ABS( (X(2,1) - X(1,1)) * (X(3,2) - X(1,2))
     %            - (X(3,1) - X(1,1)) * (X(2,2) - X(1,2)) )
C
C                                               1/12  1/24  1/24
C     COEFFICIENTS DE LA MATRICE Int Pi Pj dX = 1/24  1/12  1/24
C                                               1/24  1/24  1/12
C
C     avec INTEGRALE EXACTE POUR Pi COORDONNEES BARYCENTRIQUES
C
      DO NV=1, NBVECT
C
C        PRESSION MINIMALE DU VECTEUR NV
         PREMI   = PREMIN(NV)
C
         INTPEF  = 0D0
         INTPERR = 0D0
C
C        CALCUL DE L'INTEGRALE SUR L'EF
         DO J=1,3
C
C           NUMERO DU SOMMET J DU TRIANGLE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE PRESSION AU NOEUD NSJ
            NDLJ = NDDLNV( NSJ )
C
C           LA PRESSION_EXACTE(t,x,y,z) AU NOEUD J
            IF( NOFOPR .GT. 0 ) PresExJ = PREEXA( NSJ, NV )
C
            DO I=1,3
C
C              NUMERO DU SOMMET I DU TRIANGLE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE PRESSION AU NOEUD NSI
               NDLI = NDDLNV( NSI )
C
               IF( I .EQ. J ) THEN
C
C                 CARRE DE LA NORME L2 DE LA PRESSION-PREMIN
C                 SUR L'EF INTEGREE EXACTEMENT
                  PR = VXYPR(NDLI,NV) - PREMI
                  INTPEF = INTPEF + (PR**2) * 2D0
C
                  IF( NOFOPR .GT. 0 ) THEN
C                    PRESSION_EXACTE(temps,x,y,z) AU NOEUD J DEJA CALCULE
C                    ERREUR Pression = PressionExacte - PressionCalculee
                     INTPERR=INTPERR + ((PresExJ-VXYPR(NDLJ,NV))**2)*2D0
                  ENDIF
C
               ELSE
C
C                 CARRE DE LA NORME L2 DE LA PRESSION-PREMIN
C                 SUR L'EF INTEGREE EXACTEMENT
                  INTPEF = INTPEF + (VXYPR(NDLI,NV)-PREMI)
     %                            * (VXYPR(NDLJ,NV)-PREMI)
C
                  IF( NOFOPR .GT. 0 ) THEN
C                    LA PRESSION_EXACTE(t,x,y,z) AU NOEUD J
                     PresExI = PREEXA( NSI, NV )
C                    ERREUR Pression = PressionExacte - PressionCalculee
                     INTPERR = INTPERR + (PresExI-VXYPR(NDLI,NV))
     %                                 * (PresExJ-VXYPR(NDLJ,NV))
                  ENDIF
C
               ENDIF
C
            ENDDO
C
         ENDDO
C
C        CARRE DE LA NORME L2 DE LA PRESSION-PREMIN
C        INTEGREE EXACTEMENT SUR L'EF
         INTPRE(NV) = INTPRE(NV) + INTPEF  * SURFTR / 24D0
         INTPER(NV) = INTPER(NV) + INTPERR * SURFTR / 24D0
C
      ENDDO
C
C     SURFACE DU TRIANGLE
      SURFTR = SURFTR / 2D0
C
      RETURN
      END
