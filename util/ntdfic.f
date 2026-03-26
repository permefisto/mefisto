      SUBROUTINE NTDFIC( NMTD, NMFIC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER UN NOM DE TD EN NOM DE FICHIER CORRESPONDANT
C ----- ~>POINT>>DEFINITION   DEVIENT   ~/td/d/a_point__definition
C
C ENTREES :
C ---------
C NMTD    : NOM DU TABLEAU DESCRIPTEUR ISSU DE DICOTD
C
C SORTIES :
C ---------
C NMFIC   : NOM DU FICHIER DANS LE REPERTOIRE TD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/homdir.inc"
C.......................................................................
      CHARACTER*(*)  NMTD,NMFIC
      CHARACTER*160  KNOM
C
      KNOM = NMTD
C
C     CHANGEMENT DES > EN /
 20   L = INDEX( KNOM , '>' )
      IF( L .GT. 0 ) THEN
         KNOM(L:L) = '_'
         GOTO 20
      ENDIF
C
C     LE REPERTOIRE DE DEPART
      NMFIC = HOMDIR // '/td/d/a' // KNOM(2:160)
C
      RETURN
      END
