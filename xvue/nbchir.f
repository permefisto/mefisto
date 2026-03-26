      INTEGER FUNCTION NBCHIR( R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER L'EXPOSANT DE 10 D'UN NOMBRE REEL TEL QUE
C -----      R =  RR * 10**NBCHIR   AVEC 1 <= RR < 10
C
C ENTREES:
C --------
C R      : LE NOMBRE REEL
C
C SORTIES:
C --------
C NBCHIR : L'EXPOSANT DE 10 DU NOMBRE REEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1992
C2345X7..............................................................012
C
C     PROTECTION DU CAS ZERO ET DE R
      NBCHIR = 0
      IF( R .EQ. 0 ) RETURN
C
      RR = ABS( R )
C
 10   IF( RR .GE. 10.0 ) THEN
C
C        REEL PLUS GRAND QUE 10
         NBCHIR = NBCHIR + 1
         RR     = RR / 10
         GOTO 10
C
      ELSE IF( RR .LT. 1.0 ) THEN
C
C        REEL PLUS PETIT QUE 0.1
         NBCHIR = NBCHIR - 1
         RR     = RR * 10
         GOTO 10
      ENDIF
C
C     REEL COMPRIS ENTRE 1 ET 10
      END
