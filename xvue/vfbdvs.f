        SUBROUTINE VFBDVS( U, V, FBDVSR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR AU POINT (U,V) DU CARRE UNITE DES
C ----- 12 FONCTIONS DE BASE DE L'INTERPOLATION C1 DE
C       Fraeijs de VEUBEKE SANDER REDUITE
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
C FdVeSa : FdVeSa(10,3,12) TABLEAU DES COEFFICIENTS DES 12 FONCTIONS DE BASE
C          FdVeSa(I,J,K) I =NUMERO DU COEFFICIENT DU POLYNOME (1 a 10)
C                        J =NUMERO DU POLYNOME pj
C                        K =NUMERO DU DEGRE DE LIBERTE
C
C AVEC L'ORDRE SUIVANT POUR CHAQUE INDICE :
C P(U,V)  = SOMME pkl U**k V**l  = SOMME pi U**k V**l
C I =1 p00, I =2 p10, I =3 p01,  I =4 p11, I =5 p20, I =6 p02,
C I =7 p21, I =8 p12, I =9 p30, I =10 p03
C
C J =1 2 3 COMME INDIQUE SUR LA FIGURE CI DESSUS
C J =1  => P1 AGIT SUR LES 4 TRIANGLES
C J =2  => P2 AGIT SUR LE TRIANGLE S2S3S4
C J =3  => P3 AGIT SUR LE TRIANGLE S1S3S4
C
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C
C ENTREES:
C --------
C U,V    : VALEURS DES 2 PARAMETRES SUR LE CARRE UNITE ( 0<U,V<1 )
C
C SORTIES:
C --------
C FBDVSR : VALEUR AU POINT (U,V) DES 12 FONCTIONS DE BASE SELON L'ORDRE K
C K=1 F(S1),         K=2  F(S2),         K=3  F(S3),         K=4  F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C2345X7..............................................................012
      REAL  U, V, FBDVSR(12)
C
      UV = U * V
      UU = U * U
      VV = V * V
C
      UUU = U  * UU
      VVV = V  * VV
      UVV = U  * VV
      UUV = UU * V
C
      FBDVSR( 1) =
     &     1
     &    -3 * UU
     &    -3 * VV
     &    +3 * UVV
     &    +2 * UUU
     &    +    VVV
      FBDVSR( 2) =
     &     3 * UU
     &    -3 * UVV
     &    -2 * UUU
     &    +    VVV
      FBDVSR( 3) =
     &     3 * UVV
     &    - VVV
      FBDVSR( 4) =
     &     3 * VV
     &    -3 * UVV
     &    -    VVV
      FBDVSR( 5) =
     &           U
     &    -2   * UU
     &    -0.5 * VV
     &    +0.5 * UVV
     &    +      UUU
      FBDVSR( 6) =
     &           V
     &    -      UV
     &    -1.5 * VV
     &    +      UVV
     &    +0.5 * VVV
      FBDVSR( 7) =
     &           UV
     &    -0.5 * VV
     &    -      UVV
     &    +0.5 * VVV
      FBDVSR( 8) =
     &           UU
     &    -0.5 * UVV
     &    -      UUU
      FBDVSR( 9) =
     &     0.5 * UVV
      FBDVSR(10) =
     &           UVV
     &    -0.5 * VVV
      FBDVSR(11) =
     &           VV
     &    -      UVV
     &    -0.5 * VVV
      FBDVSR(12) =
     &     0.5 * VV
     &    -0.5 * UVV
C
       IF( V .GT. 1.-U ) THEN
C         CALCUL DE P2
      FBDVSR( 1) = FBDVSR( 1)
     &    +1
     &    -3 * U
     &    -3 * V
     &    +6 * UV
     &    +3 * UU
     &    +3 * VV
     &    -3 * UUV
     &    -3 * UVV
     &    -    UUU
     &    -    VVV
      FBDVSR( 2) = FBDVSR( 2)
     &    -1
     &    +3 * U
     &    +3 * V
     &    -6 * UV
     &    -3 * UU
     &    -3 * VV
     &    +3 * UUV
     &    +3 * UVV
     &    +    UUU
     &    +    VVV
      FBDVSR( 3) = FBDVSR( 3)
     &    +1
     &    -3 * U
     &    -3 * V
     &    +6 * UV
     &    +3 * UU
     &    +3 * VV
     &    -3 * UUV
     &    -3 * UVV
     &    -    UUU
     &    -    VVV
      FBDVSR( 4) = FBDVSR( 4)
     &    -1
     &    +3 * U
     &    +3 * V
     &    -6 * UV
     &    -3 * UU
     &    -3 * VV
     &    +3 * UUV
     &    +3 * UVV
     &    +    UUU
     &    +    VVV
      FBDVSR( 5) = FBDVSR( 5)
     &    +0.5
     &    -1.5 * U
     &    -      V
     &    +2   * UV
     &    +1.5 * UU
     &    +0.5 * VV
     &    -      UUV
     &    -0.5 * UVV
     &    -0.5 * UUU
      FBDVSR( 6) = FBDVSR( 6)
     &    +0.5
     &    -      U
     &    -1.5 * V
     &    +2   * UV
     &    +0.5 * UU
     &    +1.5 * VV
     &    -0.5 * UUV
     &    -      UVV
     &    -0.5 * VVV
      FBDVSR( 7) = FBDVSR( 7)
     &    -0.5
     &    +      U
     &    +1.5 * V
     &    -2   * UV
     &    -0.5 * UU
     &    -1.5 * VV
     &    +0.5 * UUV
     &    +      UVV
     &    +0.5 * VVV
      FBDVSR( 8) = FBDVSR( 8)
     &    +0.5 * U
     &    -      UV
     &    -      UU
     &    +      UUV
     &    +0.5 * UVV
     &    +0.5 * UUU
      FBDVSR( 9) = FBDVSR( 9)
     &    -0.5 * U
     &    +      UV
     &    +      UU
     &    -      UUV
     &    -0.5 * UVV
     &    -0.5 * UUU
      FBDVSR(10) = FBDVSR(10)
     &    -0.5 * V
     &    +      UV
     &    +      VV
     &    -0.5 * UUV
     &    -      UVV
     &    -0.5 * VVV
      FBDVSR(11) = FBDVSR(11)
     &    +0.5 * V
     &    -      UV
     &    -      VV
     &    +0.5 * UUV
     &    +      UVV
     &    +0.5 * VVV
      FBDVSR(12) = FBDVSR(12)
     &    -0.5
     &    +1.5 * U
     &    +      V
     &    -2   * UV
     &    -1.5 * UU
     &    -0.5 * VV
     &    +      UUV
     &    +0.5 * UVV
     &    +0.5 * UUU
       ENDIF
C
       IF( V .GT. U ) THEN
      FBDVSR( 1) = FBDVSR( 1)
     &    +3 * UUV
     &    -3 * UVV
     &    -    UUU
     &    +    VVV
      FBDVSR( 2) = FBDVSR( 2)
     &    -3 * UUV
     &    +3 * UVV
     &    +    UUU
     &    -    VVV
      FBDVSR( 3) = FBDVSR( 3)
     &    +3 * UUV
     &    -3 * UVV
     &    -    UUU
     &    +    VVV
      FBDVSR( 4) = FBDVSR( 4)
     &    -3 * UUV
     &    +3 * UVV
     &    +    UUU
     &    -    VVV
      FBDVSR( 5) = FBDVSR( 5)
     &    -      UV
     &    +0.5 * UU
     &    +0.5 * VV
     &    +      UUV
     &    -0.5 * UVV
     &    -0.5 * UUU
      FBDVSR( 6) = FBDVSR( 6)
     &    +      UV
     &    -0.5 * UU
     &    -0.5 * VV
     &    +0.5 * UUV
     &    -      UVV
     &    +0.5 * VVV
      FBDVSR( 7) = FBDVSR( 7)
     &    -        UV
     &    +0.5 * UU
     &    +0.5 * VV
     &    -0.5 * UUV
     &    +      UVV
     &    -0.5 * VVV
      FBDVSR( 8) = FBDVSR( 8)
     &    -      UUV
     &    +0.5 * UVV
     &    +0.5 * UUU
      FBDVSR( 9) = FBDVSR( 9)
     &    +      UUV
     &    -0.5 * UVV
     &    -0.5 * UUU
      FBDVSR(10) = FBDVSR(10)
     &    +0.5 * UUV
     &    -      UVV
     &    +0.5 * VVV
      FBDVSR(11) = FBDVSR(11)
     &    -0.5 * UUV
     &    +      UVV
     &    -0.5 * VVV
      FBDVSR(12) = FBDVSR(12)
     &    +        UV
     &    -0.5 * UU
     &    -0.5 * VV
     &    -      UUV
     &    +0.5 * UVV
     &    +0.5 * UUU
      ENDIF
      END
