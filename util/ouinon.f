      SUBROUTINE OUINON( KMOT , NON )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE MOT KMOT ET RETOURNER 0 SI NON , 1 SI OUI
C -----
C SORTIES :
C ---------
C KMOT    : LA CHAINE DE CARACTERES NON BLANC LUE
C NON     : 1 SI 1-ER CARACTERE NON BLANC EST '1' , 'O' , 'Y'
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      CHARACTER*(*) KMOT
      CHARACTER*1   CAR
C
C     LECTURE DU MOT
      NCVALS = 0
      CALL LIRCAR(NCVALS, KMOT )
      IF( NCVALS .EQ. -1 ) RETURN
      CAR = KMOT(1:1)
C
C     TEST SUR LE PREMIER CARACTERE
      IF( CAR .EQ. '1' .OR. CAR .EQ. 'O' .OR. CAR .EQ. 'Y' ) THEN
         NON = 1
      ELSE
         NON = 0
      ENDIF
      END
