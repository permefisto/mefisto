        SUBROUTINE DDUHCT( U, V, DDUFB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR AU POINT (U,V) DE LA DERIVEE PAR RAPPORT
C ----- AU PARAMETRE U DES 9 FONTIONS DE BASE DE L'EF
C       DU TRIANGLE RECTANGLE UNITE HSIEH-CLOUGH-TOCHER REDUIT
C       (REDUIT => DERIVEE NORMALE P1 SUR LES 3 COTES)
C
C       CF ARTICLE DE M. BERNARDOU et K.HASSAN (1981)
C
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
C ENTREE :
C --------
C U,V    : VALEURS DES 2 PARAMETRES SUR LE TRIANGLE RECTANGLE UNITE
C          ( 0=<U,V<=1  ET V<=1-U )
C
C SORTIE :
C --------
C DDUFB  : LA DERIVEE PAR RAPPORT A U DES 9 FONCTIONS DE BASE
C          DU HCT REDUIT RANGEES SELON K
C          K=1  F(S1),        K=2  F(S2),        K=3  F(S3),
C          K=4 DF(S1)(S2-S1), K=5 DF(S1)(S3-S1),
C          K=6 DF(S2)(S3-S2), K=7 DF(S2)(S1-S2),
C          K=8 DF(S3)(S1-S3), K=9 DF(S3)(S2-S3)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1996
C2345X7..............................................................012
      REAL    U, V, DDUFB(9)
      REAL    L1, L2, L3
C
      IF( V*2 .GE. 1-U .AND. V .GE. 1-2*U ) THEN
C
C     LE POINT (U,V) EST DANS LE SOUS TRIANGLE E1
C     ===========================================
C     LES COORDONNEES BARYCENTRIQUES DU POINT (U,V) DU TRIANGLE RECTANGLE UNITE
      L1 = 1.0 - U - V
      L2 = U
      L3 = V
C
      CB1  = -3 * L1 * L1
      CB2  =  3 * L2 * L2
CCC   CB3  =  0
      CB4  = -2 * L1 * L3
      CB5  =(-2 * L2 + L1 ) * L1
      CB6  =( 2 * L1 - L2 ) * L2
      CB7  =  2 * L2 * L3
      CB8  =      L3 * L3
      CB9  = -    L3 * L3
      CB10 =(     L1 - L2 ) * L3
C
      DDUFB(1) = -CB1 + 6 * (CB4 + CB5)
      DDUFB(2) =  CB1 + CB2 - 1.5 * (CB4+CB5) + 3 * (CB6+CB7+CB10)
      DDUFB(3) =  CB1 - 1.5 * (CB4+CB5) + 3 * (CB8+CB9+CB10)
      DDUFB(5) =  2 * CB4 + 0.5 * (CB5-CB1)
      DDUFB(4) =  0.5  * (CB4-CB1) + 2 * CB5
      DDUFB(7) =  0.5  * (CB1-CB4) - CB5 + CB6 + CB10
      DDUFB(6) =  0.25 * (CB5-CB4) + CB7 + 0.5 * CB10
      DDUFB(9) =  0.25 * (CB4-CB5) + CB8 + 0.5 * CB10
      DDUFB(8) =  0.5  * (CB1-CB5) - CB4 + CB9 + CB10
C
      RETURN
      ENDIF
C
      IF( V .GE. U .AND. V .LE. 1-2*U ) THEN
C
C     LE POINT (U,V) EST DANS LE SOUS TRIANGLE E2
C     ===========================================
C     LES COORDONNEES BARYCENTRIQUES DU POINT (U,V) DU TRIANGLE RECTANGLE UNITE
C     PERMUTEES CIRCULAIREMENT
      L1 = U
      L2 = V
      L3 = 1.0 - U - V
C
      CB1  =  3 * L1 * L1
CCC   CB2  =  0
      CB3  = -3 * L3 * L3
      CB4  =( 2 * L3 - L1 ) * L1
      CB5  =  2 * L1 * L2
      CB6  =      L2 * L2
      CB7  = -    L2 * L2
      CB8  = -2 * L3 * L2
      CB9  =(-2 * L1 + L3 ) * L3
      CB10 =( L3 - L1 ) * L2
C
      DDUFB(2) = 0.5 * CB1 + 3 * CB4 + 4.5 * CB5
      DDUFB(3) =-0.5 * CB1 + 1.5 * CB5 + 3 * (CB6+CB7)
      DDUFB(1) = CB1 + CB3 + 3 * (CB8+CB9-CB5) + 6 * CB10
      DDUFB(7) = 0.5  * (CB4+CB5)
      DDUFB(6) =-0.25 *  CB1 + 0.5 * CB4 + 1.25 * CB5
      DDUFB(9) = 0.25 * (CB1-CB5) - 0.5 * CB4 + CB6 + CB10
      DDUFB(8) = 0.5  * (CB4-CB1) + CB5 + CB7 - CB10
      DDUFB(5) = 0.5  * (CB1-CB4) - CB5 + CB8 + 2 * CB10
      DDUFB(4) = 0.5  * (CB4-CB5) + CB9 + CB10
C
      RETURN
      ENDIF
C
C     LE POINT (U,V) EST DANS LE SOUS TRIANGLE E3
C     ===========================================
C     LES COORDONNEES BARYCENTRIQUES DU POINT (U,V) DU TRIANGLE RECTANGLE UNITE
C     PERMUTEES CIRCULAIREMENT
      L1 = V
      L2 = 1.0 - U - V
      L3 = U
C
CCC   CB1  = 0
      CB2  = -3 * L2 * L2
      CB3  =  3 * L3 * L3
      CB4  =      L1 * L1
      CB5  = -    L1 * L1
      CB6  = -2 * L2 * L1
      CB7  =(-2 * L3 + L2 ) * L2
      CB8  =( 2 * L2 - L3 ) * L3
      CB9  =  2 * L3 * L1
      CB10 =( L2 - L3 ) * L1
C
      DDUFB(3) =  4.5 * CB4 + 3 * CB5
      DDUFB(1) =  CB2 + 3 * (CB6+CB7-CB4) + 6 * CB10
      DDUFB(2) =  CB3 + 1.5 * CB4 + 3 * (CB8+CB9)
      DDUFB(9) =  1.25 * CB4 + 0.5 * CB5
      DDUFB(8) =  0.5  * (CB4+CB5)
      DDUFB(5) =  0.5  * (CB5-CB4) + CB6 + CB10
      DDUFB(4) = -0.5  * CB5 - CB4 + CB7 + 2 * CB10
      DDUFB(7) =  CB4 + 0.5 * CB5 + CB8 - CB10
      DDUFB(6) = -0.25 * CB4 - 0.5 * CB5 + CB9 + CB10
C
      RETURN
      END
