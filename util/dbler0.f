      DOUBLE PRECISION FUNCTION DBLER0( R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  METTRE A ZERO LES CHIFFRES 7,8, ... DE LA MANTISSE DE R
C -----      R =  RR * 10**NBCHIR   AVEC 1 <= RR < 10
C
C ENTREES:
C --------
C R      : LE NOMBRE REEL
C
C SORTIES:
C --------
C DBLER0 : LE DOUBLE PRECISION ARRONDI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1992
C2345X7..............................................................012
      DOUBLE PRECISION D
C
C     L'EXPOSANT DE 10 DU NOMBRE REEL R TEL QUE
C     R = RR * 10**NBCHIR   AVEC 1 <= RR < 10
      N = NBCHIR( R )
C
C     PASSAGE AU DOUBLE PRECISION
      D = DBLE( R )
      D = D * ( 10D0 ** (-N+6) )
C
C     ARRONDI DANS UN ENTIER
      I = NINT( D )
C
C     RETOUR AU BON EXPOSANT
      DBLER0 = I * ( 10D0 ** (N-6) )
      END
