      SUBROUTINE RNTAIL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RENOMMER LA FONCTION TAILLE_IDEALE(x,y,z)
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1997
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*24   KNOM
C
C     LA FONCTION TAILLE_IDEALE EXISTE T ELLE?
      LANG = 0
      CALL LXNMNO( NTFONC, 'TAILLE_IDEALE', NOFOTI, MNFONC )
      IF( NOFOTI .LE. 0 ) THEN
C        LA FONCTION EDGE_LENGTH EXISTE T ELLE?
         LANG = 1
         CALL LXNMNO( NTFONC, 'EDGE_LENGTH', NOFOTI, MNFONC )
      ENDIF
C
      IF( NOFOTI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='LA FONCTION TAILLE-IDEALE(X,Y,Z) N''EXISTE PAS'
         ELSE
            KERR(1)='The FUNCTION EDGE_LENGTH(X,Y,Z) DOES NOT EXIST'
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     LA FONCTION TAILLE_IDEALE ou EDGE_LENGTH existe
      CALL INVITE( 82 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C
C     LA FONCTION TAILLE_IDEALE EXISTE T ELLE?
      IF( LANG .EQ. 0 ) THEN
C        OUI: LE NOUVEAU NOM EST IMPOSE
         CALL LXNMNM( NTFONC, 'TAILLE_IDEALE', KNOM )
      ELSE
C        OUI: LE NOUVEAU NOM EST IMPOSE
         CALL LXNMNM( NTFONC, 'EDGE_LENGTH', KNOM )
      ENDIF
C
 9000 RETURN
      END
