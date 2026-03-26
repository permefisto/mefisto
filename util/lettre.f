      LOGICAL FUNCTION LETTRE( CAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER VRAI SI LE 1-ER CARACTERE DE CAR EST UNE LETTRE
C -----           FAUX SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      CHARACTER*1 CAR
C
      LETTRE = ( CAR .GE. 'a' .AND. CAR .LE. 'z' ) .OR.
     %         ( CAR .GE. 'A' .AND. CAR .LE. 'Z' )
      END
