      SUBROUTINE CARAPR( NL , NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA POSITION DANS KTD DU CARACTERE QUI SUIT
C ----- (NL,NC) EN ENTREE
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE ACTUEL
C
C SORTIES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI LE SUIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      IF( NC .GE. NCKTD ) THEN
         NC = 1
         NL = NL + 1
      ELSE
         NC = NC + 1
      ENDIF
      END
