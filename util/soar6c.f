      SUBROUTINE SOAR6C( NBARET, NOSOAR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FOURNIR LE NUMERO DES SOMMETS DES ARETES D'UN 3-CUBE D'UN 6-CUBE
C -----
C
C SORTIES:
C --------
C NBARET : 12=NOMBRE D'ARETES du 3-CUBE D'UN 6-CUBE
C NOSOAR : NOSOAR(I,J) NO DU I-EME SOMMET DE L'ARETE J D'UN 6-CUBE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  TAMU & UPMC LJLL             SEPTEMBRE 2005
C2345X7..............................................................012
      INTEGER  NOSOAR(2,12)
C
C     LES ARETES
      NBARET = 12
C
      NOSOAR(1,1) = 1
      NOSOAR(2,1) = 2
      NOSOAR(1,2) = 2
      NOSOAR(2,2) = 4
      NOSOAR(1,3) = 4
      NOSOAR(2,3) = 3
      NOSOAR(1,4) = 3
      NOSOAR(2,4) = 1
C
      NOSOAR(1,5) = 1
      NOSOAR(2,5) = 5
      NOSOAR(1,6) = 2
      NOSOAR(2,6) = 6
      NOSOAR(1,7) = 4
      NOSOAR(2,7) = 8
      NOSOAR(1,8) = 3
      NOSOAR(2,8) = 7
C
      NOSOAR(1, 9) = 5
      NOSOAR(2, 9) = 6
      NOSOAR(1,10) = 6
      NOSOAR(2,10) = 8
      NOSOAR(1,11) = 8
      NOSOAR(2,11) = 7
      NOSOAR(1,12) = 7
      NOSOAR(2,12) = 5
C
      RETURN
      END
