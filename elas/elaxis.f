      SUBROUTINE ELAXIS ( NOTELM, SURFAM, ELAS )
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT  :   CALCUL DE LA MATRICE D ELASTICITE ( 4*4 ) POUR UN
C  ------   PROBLEME AXISYMETRIQUE
C
C  ENTREES:
C  --------
C  NOTELM : OPTION DE TRAITEMENT DU MILIEU
C  SURFAM : SI NOTELM(1)>0  SURFAM(1)=YOUNG SURFAM(2)=POISSON
C           SINON          SURAFM(1..10) = ELAS TENSEUR D ELASTICITE
C
C  SORTIE :
C  --------
C  ELAS   : LES 10 COEFFICIENTS DU TENSEUR SYMETRIQUE D ELASTICITE
C           ISOTROPE OU ANISOTROPE
C ......................................................................
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE  UPMC     PARIS   AVRIL 84
C ......................................................................
      DOUBLE PRECISION SURFAM(10),ELAS(10),D
      INTEGER          NOTELM(1)
C
C     MATERIAU HOMOGENE ISOTROPE
C     --------------------------
      IF ( NOTELM(1) .LE. 0 ) GO TO 1
      D        = SURFAM(1) / (1.D0 + SURFAM(2))
      ELAS(10) = D * 0.5D0
      D        = D / (1.D0 - 2.D0 * SURFAM(2))
      ELAS( 1) = D * ( 1.D0 - SURFAM(2))
      ELAS( 3) = ELAS(1)
      ELAS( 6) = ELAS(1)
      ELAS( 2) = D * SURFAM(2)
      ELAS( 4) = ELAS(2)
      ELAS( 5) = ELAS(2)
      ELAS( 7) = 0.D0
      ELAS( 8) = 0.D0
      ELAS( 9) = 0.D0
      GOTO 1000
C
C     MATERIAU ANISOTROPE
C     -------------------
    1 DO 2 I=1,10
         ELAS(I) = SURFAM(I)
    2 CONTINUE
C
 1000 RETURN
      END
