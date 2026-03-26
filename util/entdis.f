      SUBROUTINE ENTDIS( BASE , IXYZ , L )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DU CARRE DE LA DISTANCE A L'ORIGINE DU POINT IXYZ
C -----  CALCUL EXACT EN ENTIERS MULTI-MOTS
C
C ENTREES :
C ---------
C BASE    : BASE DU CALCUL ENTIER MULTI-MOTS
C IXYZ    : 3 COORDONNEES ENTIERES DU POINT  ( TOUTES >=0 )
C
C SORTIES :
C ---------
C L       : DISTANCE CARRE A L'ORIGINE
C           ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NL
C           L = SOMME ( L(I) * BASE ** I )   ET  L(-1) = SIGNE( L )
C                I=0                             L(-2) = NL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER  BASE,IXYZ(3),L(-2:*)
      INTEGER  A(-2:0),B(-2:2)
C
      A(-2) = 0
      A(-1) = 1
      A( 0) = IXYZ(1)
C
C     L MAX = 3 * ( BASE - 1 ) ** 2 > BASE ** 2 => B(-2:2)
C
C     L = IXYZ(1) * IXYZ(1)
      CALL ENTMUL( BASE , A , A , L )
C
C     L = L + IXYZ(2) * IXYZ(2) + IXYZ(3) * IXYZ(3)
      DO 10 I=2,3
C
C        B = IXYZ(I) * IXYZ(I)
         A(0) = IXYZ( I )
         CALL ENTMUL( BASE , A , A , B )
C
C        L = L + IXYZ(I) * IXYZ(I)
         CALL ENTSOM( BASE , L , B , L )
 10   CONTINUE
C
C     LA SOMME DES CARRES EST POSITIVE
      L(-1) = 1
      END
