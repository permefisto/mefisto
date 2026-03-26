      SUBROUTINE POIDTEP5( NBPOID, POTEP5 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES POIDS DE LA FORMULE DE QUADRATURE A 15 POINTS
C -----  SUR LE TETRAEDRE DE REFERENCE EXACTE POUR P5
C
C SORTIES:
C --------
C NBPOID : NOMBRE DE POINTS D'INTEGRATION=15
C POTEP5 : VALEURS DES POIDS DE LA FORMULE DE QUADRATURE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    Juin 2007
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  POTEP5(NBPOID), S, RAC
      INTRINSIC         SQRT
C
      NBPOID = 15
C
      POTEP5(1) = 8D0/405D0
C
      RAC = SQRT( 15D0 )
      S   = (2665D0 + 14D0 * RAC) / 226800D0
      J   = 0
C
 5    DO 10 I=1,4
         K         = J + 2
         POTEP5(K) = S
         K         = J + 3
         POTEP5(K) = S
         K         = J + 4
         POTEP5(K) = S
         K         = J + 5
         POTEP5(K) = S
 10   CONTINUE
C
      IF( J .EQ. 4 ) GOTO 20
      S  = (2665D0 - 14D0 * RAC) / 226800D0
      J  = 4
      GOTO 5
C
 20   S  = 5D0 / 567D0
      DO 30 I=10,15
         POTEP5(I) = S
 30   CONTINUE
C
      RETURN
      END
