        SUBROUTINE DDVDVS( U, V, DVDVSR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR AU POINT (U,V) DU CARRE UNITE DE LA DERIVEE
C ----- PAR RAPPORT AU PARAMETRE V DES 12 FONCTIONS DE BASE DE
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
C DVDVSR : VALEUR AU POINT (U,V) DES 12 FONCTIONS DE BASE SELON L'ORDRE K
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL  U, V, DVDVSR(12)
C
      UU = U * U
      UV = U * V
      VV = V * V
C
      DVDVSR( 1) =
     &    -6 * V
     &    +6 * UV
     &    +3 * VV
      DVDVSR( 2) =
     &    -6 * UV
     &    +3 * VV
      DVDVSR( 3) =
     &     6 * UV
     &    -3 * VV
      DVDVSR( 4) =
     &     6 * V
     &    -6 * UV
     &    -3 * VV
      DVDVSR( 5) =
     &    -      V
     &    +      UV
      DVDVSR( 6) =
     &             1
     &    -        U
     &    -3   * V
     &    +2   * UV
     &    +1.5 * VV
      DVDVSR( 7) =
     &             U
     &    -        V
     &    -2   * UV
     &    +1.5 * VV
      DVDVSR( 8) =
     &    -      UV
      DVDVSR( 9) =
     &           UV
      DVDVSR(10) =
     &       2 * UV
     &    -1.5 * VV
      DVDVSR(11) =
     &       2 * V
     &    -  2 * UV
     &    -1.5 * VV
      DVDVSR(12) =
     &           V
     &    -      UV
C
       IF( V .GT. 1.-U ) THEN
C         CALCUL DE P2
      DVDVSR( 1) = DVDVSR( 1)
     &    -3
     &    +6 * U
     &    +6 * V
     &    -3 * UU
     &    -6 * UV
     &    -3 * VV
      DVDVSR( 2) = DVDVSR( 2)
     &    +3
     &    -6 * U
     &    -6 * V
     &    +3 * UU
     &    +6 * UV
     &    +3 * VV
      DVDVSR( 3) = DVDVSR( 3)
     &    -3
     &    +6 * U
     &    +6 * V
     &    -3 * UU
     &    -6 * UV
     &    -3 * VV
      DVDVSR( 4) = DVDVSR( 4)
     &    +3
     &    -6 * U
     &    -6 * V
     &    +3 * UU
     &    +6 * UV
     &    +3 * VV
      DVDVSR( 5) = DVDVSR( 5)
     &    -1
     &    +2 * U
     &    +      V
     &    -      UU
     &    -      UV
      DVDVSR( 6) = DVDVSR( 6)
     &    -1.5
     &    +2   * U
     &    +3   * V
     &    -0.5 * UU
     &    -  2 * UV
     &    -1.5 * VV
      DVDVSR( 7) = DVDVSR( 7)
     &    +1.5
     &    -2   * U
     &    -3   * V
     &    +0.5 * UU
     &    +  2 * UV
     &    +1.5 * VV
      DVDVSR( 8) = DVDVSR( 8)
     &    -        U
     &    +        UU
     &    +        UV
      DVDVSR( 9) = DVDVSR( 9)
     &    +        U
     &    -        UU
     &    -        UV
      DVDVSR(10) = DVDVSR(10)
     &    -0.5
     &    +        U
     &    +  2 * V
     &    -0.5 * UU
     &    -  2 * UV
     &    -1.5 * VV
      DVDVSR(11) = DVDVSR(11)
     &    +0.5
     &    -        U
     &    -  2 * V
     &    +0.5 * UU
     &    +  2 * UV
     &    +1.5 * VV
      DVDVSR(12) = DVDVSR(12)
     &    +1
     &    -2   * U
     &    -        V
     &    +        UU
     &    +        UV
       ENDIF
C
       IF( V .GT. U ) THEN
      DVDVSR( 1) = DVDVSR( 1)
     &    +3 * UU
     &    -6 * UV
     &    +3 * VV
      DVDVSR( 2) = DVDVSR( 2)
     &    -3 * UU
     &    +6 * UV
     &    -3 * VV
      DVDVSR( 3) = DVDVSR( 3)
     &    +3 * UU
     &    -6 * UV
     &    +3 * VV
      DVDVSR( 4) = DVDVSR( 4)
     &    -3 * UU
     &    +6 * UV
     &    -3 * VV
      DVDVSR( 5) = DVDVSR( 5)
     &    -        U
     &    +        V
     &    +        UU
     &    -        UV
      DVDVSR( 6) = DVDVSR( 6)
     &    +        U
     &    -        V
     &    +0.5 * UU
     &    -  2 * UV
     &    +1.5 * VV
      DVDVSR( 7) = DVDVSR( 7)
     &    -        U
     &    +        V
     &    -0.5 * UU
     &    +  2 * UV
     &    -1.5 * VV
      DVDVSR( 8) = DVDVSR( 8)
     &    -        UU
     &    +        UV
      DVDVSR( 9) = DVDVSR( 9)
     &    +        UU
     &    -        UV
      DVDVSR(10) = DVDVSR(10)
     &    +0.5 * UU
     &    -  2 * UV
     &    +1.5 * VV
      DVDVSR(11) = DVDVSR(11)
     &    -0.5 * UU
     &    +  2 * UV
     &    -1.5 * VV
      DVDVSR(12) = DVDVSR(12)
     &    +        U
     &    -        V
     &    -        UU
     &    +        UV
      ENDIF
      END
