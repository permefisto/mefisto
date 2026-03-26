      SUBROUTINE NMDTMS( KNOMTS , NMTOBJ , NMOBJT , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOM DU TYPE ET LE NOM DE PLSVO A PARTIR
C ----- DU NOM D'UN TMS
C       EXEMPLE : TS : ~>POINT>P1>DEFINITION
C                 => NMTOBJ = 'POINT'
C                    NMOBJT = 'P1'
C
C ENTREE :
C --------
C KNOMTS : NOM DU TABLEAU MS
C
C SORTIES :
C ---------
C NMTOBJ : NOM DU TYPE  POINTS ou LIGNES ou ...
C NMOBJT : NOM DU PLSVO... OU '    ' SINON
C IERR   : =0 SI PAS D'ERREUR
C          =1 NOM DE TMS INCORRECT SANS >
C          =2 NOM DE TMS SANS NOM DE PLSVO    ( EXEMPLE: '~>POINT>' )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1991
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     KNOMTS,NMTOBJ,NMOBJT
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     KNOMTS DEBUTE T IL PAR ~   ?
      N1 = 1
      IF( KNOMTS( 1: 1) .EQ. '~' ) N1 = 2
      IF( KNOMTS(N1:N1) .EQ. '>' ) N1 = N1 + 1
C
C     RECHERCHE DE >
      L  = LEN( KNOMTS )
      N2 = INDEX( KNOMTS(N1:L) , '>' )
      IF( N2 .LE. 0 ) THEN
         IERR = 1
         RETURN
      ENDIF
      N2 = N2 + N1 - 1
C
C     LE NOM DU TYPE
      NMTOBJ = KNOMTS(N1:N2-1)
C
C     RECHERCHE DE > SUIVANT
      N1 = N2 + 1
      IF( N1 .GT. L ) THEN
         IERR = 2
         RETURN
      ENDIF
      N2 = INDEX( KNOMTS(N1:L) , '>' )
      IF( N2 .GT. 0 ) THEN
C
C        LE DERNIER CARACTERE DU PLSVO PRECEDE >
         N2 = N2 + N1 - 2
C
      ELSE
C
C        EXISTE T IL UN BLANC FINAL?
         N2 = INDEX( KNOMTS(N1:L) , ' ' )
         IF( N2 .GT. 0 ) THEN
C           OUI
            N2 = N2 + N1 - 2
         ELSE
C           NON : DERNIER CARACTERE DU TMS
            N2  = L
         ENDIF
      ENDIF
      NMOBJT = KNOMTS(N1:N2)
      IERR = 0
      END
