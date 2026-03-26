      SUBROUTINE ENTIDE( NXYZ , IXYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRANSFORMATION DES COORDONNEES ENTIERES EN COORDONNEES ENTIERES
C -----  VARIANTE DE MISE INITIALE EN COORDONNEES ENTIERES
C
C ENTREES :
C ---------
C NXYZ    : 3 COORDONNEES ENTIERES
C
C SORTIES :
C ---------
C IXYZ    : 3 COORDONNEES ENTIERES CHACUNE SUR UN MOT ENTIER MACHINE
C
C POUR PLUS DE PRECISION VOIR LE SOUS-PROGRAMME ENTDET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  JUILLET 88
C23456...............................................................012
      INTEGER    NXYZ(3),IXYZ(3)
C
      IXYZ(1) = NXYZ(1)
      IXYZ(2) = NXYZ(2)
      IXYZ(3) = NXYZ(3)
      END
