      SUBROUTINE PER2AR( I , J , NOARET , NLARET )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PERMUTATION DES 2 ARETES I ET J DANS LES TABLEAUX
C -----    NOARET NLARET
C
C ENTREES:
C --------
C I J    : NUMERO DES ARETES A PERMUTER
C
C ENTREES ET SORTIES:
C -------------------
C NOARET :  LE NUMERO POOL DES ARETES DU CUBE
C NLARET :  LE NUMERO DE GROUPE DES ARETES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  MARS 1987
C.....................................................................
      INTEGER  NOARET(1:*),NLARET(1:*)
C
C     PERMUTATION DES ARETES I ET J
      NN        = NOARET(I)
      NOARET(I) = NOARET(J)
      NOARET(J) = NN
      NN        = NLARET(I)
      NLARET(I) = NLARET(J)
      NLARET(J) = NN
      END
