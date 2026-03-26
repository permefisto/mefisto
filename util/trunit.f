      SUBROUTINE TRUNIT( NF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER UN NUMERO D'UNITE LOGIQUE LIBRE
C -----

C SORTIE :
C --------
C NF     : NUMERO DE L'UNITE LOGIQUE LIBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      LOGICAL LOPEN

      NF = 20
 1    NF = NF + 1
      INQUIRE( UNIT=NF , OPENED=LOPEN )
      IF( LOPEN ) GOTO 1

C     L'UNITE NF EST LIBRE
      RETURN
      END
