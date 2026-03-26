      SUBROUTINE FICAKC( NFFICA , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : METTRE DANS LE BUFFER KFICA LE NOM KNOM
C -----
C
C ENTREES :
C --------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C KNOM    : NOM A PORTER DANS LE BUFFER KFICA
C           EN FAIT SEULEMENT LES CARACTERES JUSQU'AU PREMIER BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     KNOM
C
C     LE NOMBRE DE CARACTERES A METTRE DANS KFICA
      NBC = INDEX( KNOM , ' ' )
      IF( NBC .LE. 0 ) THEN
         NBC = LEN( KNOM )
      ENDIF
C
C     LE BUFFER SERAIT-IL INSUFFISANT?
      IF( LCFICA + NBC .GT. NCFICA ) THEN
C        OUI . IL EST PORTE SUR LE FICHIER
         CALL FICAWR( NFFICA )
      ENDIF
C
C     REMPLISSAGE DANS LE BUFFER KFICA
      KFICA(LCFICA+1:LCFICA+NBC) = ' ' // KNOM(1:NBC-1)
      LCFICA = LCFICA + NBC
      END
