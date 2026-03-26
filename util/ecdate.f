      SUBROUTINE ECDATE( NTABL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DONNER LA DATE ACTUELLE EN POSITION 1 ET 2 DU TABLEAU NTABL
C -----
C
C SORTIE:
C -------
C NTABL : NTABL(1:2) DOUBLE PRECISION STOCKE DANS UN ENTIER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS          MARS 1996
C.....................................................................
      INTEGER           NTABL(1:2)
C
C     LE REEL DOUBLE PRECISION EST EN EQUIVALENCE SUR 2 ENTIERS
      INTEGER           NDATE(1:2)
      DOUBLE PRECISION  DATE
      EQUIVALENCE      (DATE,NDATE(1))
C
C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
      CALL  SECONDES1970( DATE )
C
C     CETTE DATE EST STOCKEE DANS 2 ENTIERS SANS CONVERSION
      NTABL(1) = NDATE(1)
      NTABL(2) = NDATE(2)
C
CCC      PRINT *, 'ECDATE: DATE=',DATE,' NDATE=',NDATE(1),' ',NDATE(2)
      RETURN
      END
