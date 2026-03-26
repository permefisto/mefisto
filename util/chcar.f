      SUBROUTINE CHCAR( CAR,  NL, NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NL,NC) DU PROCHAIN CARACTERE NON BLANC
C ----- QUI DOIT ETRE UN CAR
C
C ENTREES :
C ---------
C CAR     : LE CARACTERE A RETROUVER
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C SORTIES :
C ---------
C NL ,NC  : POSITION DE CAR     NL=0 SI LE 1-ER CARCTERE N'EST PAS CAR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      CHARACTER*1  CAR
C
C     POSITION DU 1-ER CARACTERE NON BLANC
      CALL CAR1NB( NL , NC )
      IF( KTD(NL)(NC:NC) .NE. CAR ) THEN
C        CE N'EST PAS CAR
         NL = 0
      ENDIF
C
      RETURN
      END
