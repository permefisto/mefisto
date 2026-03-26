      SUBROUTINE CARAVA( NL , NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA POSITION DANS KTD DU CARACTERE QUI PRECEDE
C ----- (NL,NC) EN ENTREE
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE ACTUEL
C
C SORTIES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      IF( NC .LE. 1 ) THEN
         NC = NCKTD
         NL = NL - 1
      ELSE
         NC = NC - 1
      ENDIF
      END
