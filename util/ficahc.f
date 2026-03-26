      SUBROUTINE FICAHC( NFFICA , CAR4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE DANS LE BUFFER KFICA LES 4 CARACTERES CAR4
C -----
C ENTREES :
C --------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C CAR4    : 4 CARACTERES A METTRE DANS LE BUFFER KFICA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   MARS 1989
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*4       CAR4
C
C     LE BUFFER SERAIT-IL INSUFFISANT?
      IF( LCFICA + 5 .GT. NCFICA ) THEN
C        OUI . IL EST PORTE SUR LE FICHIER
         CALL FICAWR( NFFICA )
      ENDIF
C
C     REMPLISSAGE DANS LE BUFFER KFICA
      KFICA( LCFICA+1 : LCFICA+5 ) = ' ' // CAR4
C
C     MISE A JOUR DU POINTEUR SUR LE DERNIER CARACTERE
C     DANS LE BUFFER KFICA
      LCFICA = LCFICA + 5
C
      RETURN
      END
