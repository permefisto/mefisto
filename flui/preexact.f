      SUBROUTINE PREEXACT( NOFOPR, NBNOPR, NBVECT, TIMES,  XYZN,
     %                     PREEXA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE TABLEAU DE LA PRESSION EXACTE(NOEUDS,TEMPS)
C -----
C ENTREES:
C --------
C NOFOPR : NUMERO DE LA FONCTION UTILISATEUR PRESSION_EXACTE(t,x,y,z)
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C NBVECT : NOMBRE TOTAL DE VECTEURS PRESSION PRESSION
C TIMES  : NBVECT TEMPS DU CALCUL DE LA PRESSION
C XYZN   : 3 COORDONNEES DES NBNOPR NOEUDS DU MAILLAGE
C
C SORTIE :
C --------
C PREEXA : PRESSION_EXACTE(NOEUD,TEMPS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  FEVRIER 2012
C234567...............................................................12
      DOUBLE PRECISION  PREEXA(NBNOPR,NBVECT), DPARAF(4)
      REAL              TIMES(NBVECT), XYZN(3,NBNOPR)
C
      DO K=1,NBVECT
C
C        LE TEMPS K
         TEMPS = TIMES(K)
C
         DO N=1,NBNOPR
C
C           PRESSION_EXACTE(t,x,y,z,nocomp) ou
C           EXACT_VELOCITY(t,x,y,z,nocomp) AU NOEUD ET AU TEMPS K
            DPARAF(1) = TEMPS
            DPARAF(2) = XYZN(1,N)
            DPARAF(3) = XYZN(2,N)
            DPARAF(4) = XYZN(3,N)
C
C           LA PRESSION EXACTE AU NOEUD I AU TEMPS K
            CALL FONVAL( NOFOPR, 4, DPARAF, NCODEV, PREEXA(N,K) )
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
