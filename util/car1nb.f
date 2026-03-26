      SUBROUTINE CAR1NB( NL , NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NL,NC) DANS KTD DU 1-ER CARACTERE NON BLANC
C ----- APRES LE CARACTERE COURANT (NL,NC)
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C
C SORTIES :
C ---------
C NL , NC : POSITION DANS KTD DU PREMIER CARACTERE NON BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     LE CARACTERE APRES
 10   CALL CARAPR( NL , NC )
C
C     LE CARACTERE EST-IL BLANC ?
      IF( KTD(NL)(NC:NC) .EQ. ' ' ) GOTO 10
      RETURN
      END
