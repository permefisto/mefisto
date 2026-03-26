      SUBROUTINE TRDISY( N, A, B, D, U, V, X )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RESOUDRE UN SYSTEME LINEAIRE A N INCONNUES AVEC UNE
C -----    MATRICE TRIDIAGONALE SYMETRIQUE
C
C ENTREES:
C --------
C N      : NOMBRE DE LIGNES DE LA MATRICE ET DU SECOND MEMBRE
C A      : SOUS DIAGONALE DE LA MATRICE QUI EST SYMETRIQUE
C B      : DIAGONALE DE LA MATRICE
C D      : SECOND-MEMBRE
C
C AUXILIAIRE:
C -----------
C U,V    : 2 VECTEUR DE N COMPOSANTES REELLES
C
C SORTIES:
C --------
C X     : LES N COMPOSANTES DE LA SOLUTION DU SYSTEME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1995
C2345X7..............................................................012
      REAL   A(N), B(N), D(N), U(N), V(N), X(N)
C
      N1 = N - 1
C
C     RESOLUTION PAR SUBSTITUTION
C     XI-1 = UI-1 - VI-1 * XI PORTE DANS LA LIGNE I ...
C     DANS LA DERNIERE LIGNE XN=... ET REMONTEE
C
C     CALCUL DES UI VI POUR I=1,...,N-1
C     =================================================
      V(1) = A(2) / B(1)
      U(1) = D(1) / B(1)
      DO 10 K=2,N1
         V(K) = A(K+1) / ( B(K) - A(K) * V(K-1) )
         U(K) = ( D(K) - A(K) * U(K-1) ) * V(K) / A(K+1)
 10   CONTINUE
C
C     REMONTEE : XI-1 = UI-1 - VI-1 * XI
C     ==================================
      X(N) = ( D(N) - A(N) * U(N1) ) / ( B(N) - A(N) * V(N1) )
      DO 20 K=N1,1,-1
         X(K) = U(K) - V(K) * X(K+1)
 20   CONTINUE
      END
