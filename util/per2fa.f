      SUBROUTINE PER2FA( I , J , NOFACE , NSFACE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PERMUTATION DES 2 FACES I ET J DANS LES TABLEAUX
C -----    NOFACE NSFACE
C
C ENTREES:
C --------
C I J    : NUMERO DES FACES A PERMUTER
C
C ENTREES ET SORTIES:
C -------------------
C NOFACE :  LE NUMERO POOL DES FACES DU CUBE
C NSFACE :  LE NUMERO DE GROUPE DES FACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  MARS 1987
C.....................................................................
      INTEGER  NOFACE(0:*),NSFACE(1:*)
C
C     PERMUTATION DES FACES I ET J
      NN        = NOFACE(I)
      NOFACE(I) = NOFACE(J)
      NOFACE(J) = NN
      NN        = NSFACE(I)
      NSFACE(I) = NSFACE(J)
      NSFACE(J) = NN
      END
