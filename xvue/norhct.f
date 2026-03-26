        SUBROUTINE NORHCT( U, V, X, Y, Z,  NORMAL, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR AU POINT (U,V) DES 3 COMPOSANTES DE LA NORMALE
C ----- ORTHOGONALE AUX 2 TANGENTES EN U ET V DE L'INTERPOLATION
C       DU TRIANGLE RECTANGLE UNITE HSIEH-CLOUGH-TOCHER REDUIT
C       (REDUIT => DERIVEE NORMALE P1 SUR LES 3 COTES)
C       ET POUR LES 9 DEGRES DE LIBERTE STOCKES DANS X Y Z
C       CF ARTICLE DE M. BERNARDOU et K.HASSAN (1981)
C              ^
C              | V
C              |
C             S3=(0,1)
C              X
C              | \ \
C              |  \   ^
C              -   \    \
C              |    \     \
C              |     \      \
C              | E2   \       \
C              |       \    E1  \
C              ^     /   \        ^
C              |   /         \      \
C              | /      E3        \   \
C              X----->------------<-----X ----> U
C             S1=(0,0)                S2=(1,0)
C
C ENTREES:
C --------
C U, V   : VALEURS DES 2 PARAMETRES SUR LE TRIANGLE RECTANGLE UNITE
C          ( 0=<U,V<=1  ET V<=1-U )
C X, Y, Z: VALEURS DES 9 DEGRES DE LIBERTE RANGEES SELON K POUR X, Y, Z
C K=1  F(S1),        K=2  F(S2),        K=3  F(S3),
C K=4 DF(S1)(S2-S1), K=5 DF(S1)(S3-S1),
C K=6 DF(S2)(S3-S2), K=7 DF(S2)(S1-S2),
C K=8 DF(S3)(S1-S3), K=9 DF(S3)(S2-S3)
C
C SORTIES:
C --------
C NORMAL : 3 COMPOSANTES DE LA NORMALE UNITAIRE AU POINT (U,V)
C IERR   : 1 SI LES 2 TANGENTES SONT COLINAIRES
C          0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1998
C2345X7..............................................................012
      REAL    U, V, X(9), Y(9), Z(9), TGU(3), TGV(3), NORMAL(3)
C
C     LES 3 COMPOSANTES DE LA TG EN U
      CALL TGUHCT( U, V, X, Y, Z,  TGU(1), TGU(2), TGU(3) )
C
C     LES 3 COMPOSANTES DE LA TG EN V
      CALL TGVHCT( U, V, X, Y, Z,  TGV(1), TGV(2), TGV(3) )
C
C     LES 3 COMPOSANTES DE LA NORMALE UNITAIRE
      CALL PROVER( TGU, TGV, NORMAL )
C
      S = SQRT( NORMAL(1)**2 + NORMAL(2)**2 + NORMAL(3)**2 )
      IF( S .EQ. 0 ) THEN
C        ERREUR: 2 TANGENTES COLINEAIRES
         IERR = 1
      ELSE
         IERR = 0
         NORMAL(1) = NORMAL(1) / S
         NORMAL(2) = NORMAL(2) / S
         NORMAL(3) = NORMAL(3) / S
      ENDIF
      END
