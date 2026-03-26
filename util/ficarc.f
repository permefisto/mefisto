      SUBROUTINE FICARC( NFFICA , REEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE DANS LE BUFFER KFICA LE REEL REEL
C -----
C
C ENTREES :
C --------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C REEL    : VALEUR REELLE A METTRE DANS LE BUFFER KFICA
C       EN FAIT SEULEMENT LES CARACTERES NON BLANC PRECEDE D'UN BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              REEL
      CHARACTER*15      KREEL
C
C     LA CONVERSION DU REEL EN CARACTERES
      WRITE( KREEL , '(E15.7)' ) REEL
C
C     ELIMINATION DES BLANCS A GAUCHE
      DO 20 I1=1,15
         IF( KREEL(I1:I1) .NE. ' ' ) GOTO 30
 20   CONTINUE
C
C     ELIMINATION DES BLANCS A DROITE
 30   DO 40 I2=15,I1+1,-1
         IF( KREEL(I2:I2) .NE. ' ' ) GOTO 50
 40   CONTINUE
C
C     LE BUFFER SERAIT-IL INSUFFISANT?
 50   NBC = I2 - I1 + 2
      IF( LCFICA + NBC .GT. NCFICA ) THEN
C        OUI . IL EST PORTE SUR LE FICHIER
         CALL FICAWR( NFFICA )
      ENDIF
C
C     REMPLISSAGE DANS LE BUFFER KFICA
      KFICA( LCFICA+1 : LCFICA+NBC ) = ' ' // KREEL( I1 : I2 )
C
C     MISE A JOUR DU POINTEUR SUR LE DERNIER CARACTERE
C     DANS LE BUFFER KFICA
      LCFICA = LCFICA + NBC
      END
