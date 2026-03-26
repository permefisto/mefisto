      SUBROUTINE FICADC( NFFICA , REEL2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE DANS LE BUFFER KFICA LE REEL DOUBLE PRECISION REEL2
C -----
C
C ENTREES :
C --------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C REEL2   : VALEUR REEL2LE A METTRE DANS LE BUFFER KFICA
C       EN FAIT SEULEMENT LES CARACTERES NON BLANC PRECEDE D'UN BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  REEL2
      CHARACTER*25      KREEL2
C
C     LA CONVERSION DU REEL2 EN CARACTERES
      WRITE( KREEL2 , '(D25.17)' ) REEL2
C
C     ELIMINATION DES BLANCS A GAUCHE
      DO 20 I1=1,25
         IF( KREEL2(I1:I1) .NE. ' ' ) GOTO 30
 20   CONTINUE
C
C     ELIMINATION DES BLANCS A DROITE
 30   DO 40 I2=25,I1+1,-1
         IF( KREEL2(I2:I2) .NE. ' ' ) GOTO 50
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
      KFICA( LCFICA+1 : LCFICA+NBC ) = ' ' // KREEL2( I1 : I2 )
C
C     MISE A JOUR DU POINTEUR SUR LE DERNIER CARACTERE
C     DANS LE BUFFER KFICA
      LCFICA = LCFICA + NBC
      END
