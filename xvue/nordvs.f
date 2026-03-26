        SUBROUTINE NORDVS( U, V, X, Y, Z,  NORMAL, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE VECTEUR NORMAL UNITAIRE AU POINT (U,V) DU CARRE UNITE
C ----- DE L'INTERPOLATION C1 DE Fraeijs de VEUBEKE SANDER REDUITE
C       (DERIVEE NORMALE P1 SUR LES 4 COTES)
C
C              ^ V
C              |
C              |
C              |
C             S4=(0,1)               S3=(1,1)
C              X----->------------<---X
C              | \     p1+p2+p3     / |
C              |   \              /   |
C              -     \          /     -
C              |       \      /       |
C              |          \ /         |
C              | p1+p3   /   \  p1+p2 |
C              |       /       \      |
C              ^     /           \    ^
C              |   /               \  |
C              | /        p1         \|
C              X----->----------<-----X  ------------> U
C             S1=(0,0)               S2=(1,0)
C
C ENTREES:
C --------
C U,V    : VALEURS DES 2 PARAMETRES SUR LE CARRE UNITE ( 0<U,V<1 )
C X, Y, Z: VALEURS DES 12 DEGRES DE LIBERTE RANGEES SELON K POUR X Y Z
C
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C
C SORTIES:
C --------
C NORMAL : 3 COMPOSANTES DE LA NORMALE UNITAIRE AU POINT (U,V)
C IERR   : 1 SI LES 2 TANGENTES SONT COLINAIRES
C          0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1998
C2345X7..............................................................012
      REAL    U, V, X(12), Y(12), Z(12), TGU(3), TGV(3), NORMAL(3)
C
C     LES 3 COMPOSANTES DE LA TG EN U
      CALL TGUDVS( U, V, X, Y, Z,  TGU(1), TGU(2), TGU(3) )
C
C     LES 3 COMPOSANTES DE LA TG EN V
      CALL TGVDVS( U, V, X, Y, Z,  TGV(1), TGV(2), TGV(3) )
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
