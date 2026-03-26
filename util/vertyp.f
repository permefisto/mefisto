      SUBROUTINE VERTYP( NUTYOB , NUOBJT , NUALGO , LADEFI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   VERIFIER SI LE TYPE DE PLSVO EST CORRECT
C -----   AFFICHAGE D'UN DIAGNOSTIC SINON
C
C ENTREES :
C ---------
C NUTYOB : NUMERO DU TYPE DE PLSVO   1:POINT, 2:LIGNE, ...
C NUOBJT : NUMERO DU PLSVO DANS SON LEXIQUE
C NUALGO : NUMERO DE L'ALGORITHME DE DEFINITION DU MAILLAGE
C LADEFI : NUMERO LU DANS LA DEFINITION DU PLSVO
C
C SORTIES :
C ---------
C IERR   : CODE D'ERREUR   0 SI PAS d'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      CHARACTER*24    KNOM
      CHARACTER*10    NMTYOB,KTYOBJ
C
      IF( LADEFI .NE. NUALGO ) THEN
         NBLGRC(NRERR) = 2
         KTYOBJ  = NMTYOB( NUTYOB )
         L       = NUDCNB( KTYOBJ )
         CALL NMOBNU( KTYOBJ , NUOBJT , KNOM )
         KERR(1) = KTYOBJ(1:L) // ' ' // KNOM
         KERR(2) = 'TYPE DEFINITION INCORRECT'
         CALL LEREUR
         IERR = 1
      ELSE
         IERR = 0
      ENDIF
      END
