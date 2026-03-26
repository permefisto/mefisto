      SUBROUTINE FICAXC( NFFICA , XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE DANS LE BUFFER KFICA LES 3 COORDONNEES XYZ
C -----
C
C ENTREES :
C --------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER DES CARACTERES
C XYZ     : VALEUR DES 3 COORDONNEES A METTRE DANS LE BUFFER KFICA
C       EN FAIT SEULEMENT LES CARACTERES NON BLANC PRECEDE D'UN BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              XYZ(3)
      CHARACTER*15      KXYZ
C
      DO 100 N=1,3
C
C        LA CONVERSION DU XYZ EN CARACTERES
         WRITE( KXYZ , '(E15.7)' ) XYZ(N)
C
C        ELIMINATION DES BLANCS A GAUCHE
         DO 20 I1=1,15
            IF( KXYZ(I1:I1) .NE. ' ' ) GOTO 30
 20      CONTINUE
C
C        ELIMINATION DES BLANCS A DROITE
 30      DO 40 I2=15,I1+1,-1
            IF( KXYZ(I2:I2) .NE. ' ' ) GOTO 50
 40      CONTINUE
C
C        LE BUFFER SERAIT-IL INSUFFISANT?
 50      NBC = I2 - I1 + 2
         IF( LCFICA + NBC .GT. NCFICA ) THEN
C           OUI . IL EST PORTE SUR LE FICHIER
            CALL FICAWR( NFFICA )
         ENDIF
C
C        REMPLISSAGE DANS LE BUFFER KFICA
         KFICA( LCFICA+1 : LCFICA+NBC ) = ' ' // KXYZ( I1 : I2 )
C
C        MISE A JOUR DU POINTEUR SUR LE DERNIER CARACTERE
C        DANS LE BUFFER KFICA
         LCFICA = LCFICA + NBC
 100  CONTINUE
      END
