        SUBROUTINE DDUDVS( U, V, DUDVSR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR AU POINT (U,V) DU CARRE UNITE DE LA DERIVEE
C ----- PAR RAPPORT AU PARAMETRE U DES 12 FONCTIONS DE BASE DE
C       L'INTERPOLATION C1 DE Fraeijs de VEUBEKE SANDER REDUITE
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
C
C SORTIES:
C --------
C DUDVSR : VALEUR AU POINT (U,V) DES 12 FONCTIONS DE BASE SELON L'ORDRE K
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL  U, V, DUDVSR(12)
C
      UU = U * U
      UV = U * V
      VV = V * V
C
      DUDVSR( 1) =
     &    -6 * U
     &    +3 * VV
     &    +6 * UU
      DUDVSR( 2) =
     &     6 * U
     &    -3 * VV
     &    -6 * UU
      DUDVSR( 3) =
     &     3 * VV
      DUDVSR( 4) =
     &    -3 * VV
      DUDVSR( 5) =
     &       1
     &    -  4 * U
     &    +0.5 * VV
     &    +  3 * UU
      DUDVSR( 6) =
     &    -        V
     &    +        VV
      DUDVSR( 7) =
     &             V
     &    -        VV
      DUDVSR( 8) =
     &       2 * U
     &    -0.5 * VV
     &    -  3 * UU
      DUDVSR( 9) =
     &     0.5 * VV
      DUDVSR(10) =
     &             VV
      DUDVSR(11) =
     &    -        VV
      DUDVSR(12) =
     &    -0.5 * VV
C
       IF( V .GT. 1.-U ) THEN
C         CALCUL DE P2
      DUDVSR( 1) = DUDVSR( 1)
     &    -3
     &    +6 * V
     &    +6 * U
     &    -6 * UV
     &    -3 * VV
     &    -3 * UU
      DUDVSR( 2) = DUDVSR( 2)
     &    +3
     &    -6 * V
     &    -6 * U
     &    +6 * UV
     &    +3 * VV
     &    +3 * UU
      DUDVSR( 3) = DUDVSR( 3)
     &    -3
     &    +6 * V
     &    +6 * U
     &    -6 * UV
     &    -3 * VV
     &    -3 * UU
      DUDVSR( 4) = DUDVSR( 4)
     &    +3
     &    -6 * V
     &    -6 * U
     &    +6 * UV
     &    +3 * VV
     &    +3 * UU
      DUDVSR( 5) = DUDVSR( 5)
     &    -1.5
     &    +  2 * V
     &    +  3 * U
     &    -  2 * UV
     &    -0.5 * VV
     &    -1.5 * UU
      DUDVSR( 6) = DUDVSR( 6)
     &    -1
     &    +2   * V
     &    + U
     &    - UV
     &    - VV
      DUDVSR( 7) = DUDVSR( 7)
     &    +1
     &    -2   * V
     &    -        U
     &    +        UV
     &    +        VV
      DUDVSR( 8) = DUDVSR( 8)
     &    +0.5
     &    -        V
     &    -  2 * U
     &    +  2 * UV
     &    +0.5 * VV
     &    +1.5 * UU
      DUDVSR( 9) = DUDVSR( 9)
     &    -0.5
     &    +        V
     &    +  2 * U
     &    -  2 * UV
     &    -0.5 * VV
     &    -1.5 * UU
      DUDVSR(10) = DUDVSR(10)
     &    +        V
     &    -        UV
     &    -        VV
      DUDVSR(11) = DUDVSR(11)
     &    -        V
     &    +        UV
     &    +        VV
      DUDVSR(12) = DUDVSR(12)
     &    +1.5
     &    -  2 * V
     &    -  3 * U
     &    +  2 * UV
     &    +0.5 * VV
     &    +1.5 * UU
       ENDIF
C
       IF( V .GT. U ) THEN
      DUDVSR( 1) = DUDVSR( 1)
     &    +6 * UV
     &    -3 * VV
     &    -3 * UU
      DUDVSR( 2) = DUDVSR( 2)
     &    -6 * UV
     &    +3 * VV
     &    +3 * UU
      DUDVSR( 3) = DUDVSR( 3)
     &    +6 * UV
     &    -3 * VV
     &    -3 * UU
      DUDVSR( 4) = DUDVSR( 4)
     &    -6 * UV
     &    +3 * VV
     &    +3 * UU
      DUDVSR( 5) = DUDVSR( 5)
     &    -        V
     &    +        U
     &    +  2 * UV
     &    -0.5 * VV
     &    -1.5 * UU
      DUDVSR( 6) = DUDVSR( 6)
     &    +        V
     &    -        U
     &    +        UV
     &    -        VV
      DUDVSR( 7) = DUDVSR( 7)
     &    -        V
     &    +        U
     &    -        UV
     &    +        VV
      DUDVSR( 8) = DUDVSR( 8)
     &    -  2 * UV
     &    +0.5 * VV
     &    +1.5 * UU
      DUDVSR( 9) = DUDVSR( 9)
     &    +  2 * UV
     &    -0.5 * VV
     &    -1.5 * UU
      DUDVSR(10) = DUDVSR(10)
     &    +        UV
     &    -        VV
      DUDVSR(11) = DUDVSR(11)
     &    -        UV
     &    +        VV
      DUDVSR(12) = DUDVSR(12)
     &    +        V
     &    -        U
     &    -  2 * UV
     &    +0.5 * VV
     &    +1.5 * UU
      ENDIF
      END
