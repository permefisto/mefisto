      INTEGER FUNCTION NBPRMT( NOCOPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE NOMBRE DE PARAMETRES DU CODE OPERATION NOCOPE
C -----
C
C ENTREE :
C --------
C NOCOPE : LE CODE OPERATION ( CF $MEFISTO/incl/lu.inc )
C
C SORTIE :
C --------
C NBPRMT : NOMBRE DE PARAMETRES DE CETTE OPERATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      IF( NOCOPE .LT. 100 ) THEN
         NBPRMT = 0
      ELSE IF( NOCOPE  .LT. 200 ) THEN
         NBPRMT = 1
      ELSE IF( NOCOPE .LT. 1000 ) THEN
         NBPRMT = 2
      ELSE
         NBPRMT = NOCOPE / 1000
      ENDIF
      END
