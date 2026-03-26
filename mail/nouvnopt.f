      INTEGER FUNCTION NOUVNOPT( NP, NONOPT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NUMERO DE POINT A TRAVERS LE TABLEAU NONOPT
C -----    INDIQUANT PAR VALEUR NEGATIVE LE BON NUMERO DE POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Decembre 2011
C2345X7..............................................................012
      INTEGER NONOPT(1:*)
C
      NOUVNOPT = NP
 10   NOUVNOPT = NONOPT( ABS(NOUVNOPT) )
      IF( NOUVNOPT .LT. 0 ) GOTO 10
C
      RETURN
      END
