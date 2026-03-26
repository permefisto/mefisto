      SUBROUTINE DFTRAN
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR LES TRANSFORMATIONS C-A-D
C------ CREER ET LIRE LE TABLEAU ~>TRANSFO>...>DEFINITION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*24      KNOMTR
      CHARACTER*160     KNMTD,KNMTS
C
C     NOM_DE_LA_TRANSFORMATION  LE CARACTERE @ POUR FINIR
 100  CALL INVITE( 63 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , KNOMTR )
      IF( NCVALS .EQ. -1 ) GOTO 9000
      ILEXIS = 1
C
C     RECHERCHE DU NOM DE LA TRANSFORMATION DANS LE LEXIQUE
C     DES TRANSFORMATIONS
 150  CALL LXLXOU( NTTRAN , KNOMTR , NTLXTR , MNLXTR )
C
C     S'IL N'EXISTE PAS IL EST CREE
      IF( NTLXTR .LE. 0 ) THEN
         ILEXIS = 0
         CALL LXLXDC( NTTRAN , KNOMTR , 24 , 8 )
         GOTO 150
      ENDIF
C
C     MODIFICATION OU CREATION DE LA DEFINITION DE LA TRANSFORMATION
C     ==============================================================
      KNMTD = '~>TRANSFO>>DEFINITION'
      I     = INDEX( KNOMTR , ' ' )
      IF( I .LE. 1 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOM INCORRECT DE TRANSFORMATION: '//KNOMTR
         CALL LEREUR
         GOTO 100
      ENDIF
      KNMTS = '~>TRANSFO>' // KNOMTR(1:I-1) // '>DEFINITION'
      CALL MOTSTD( KNMTD(1:NUDCNB(KNMTD)), KNMTS(1:NUDCNB(KNMTS)),
     %             NRETOU)
      IF( NRETOU .NE. 0 ) THEN
C        DESTRUCTION SI LA TRANSFORMATION N'EXISTAIT PAS AVANT
         IF( ILEXIS .EQ. 0 ) CALL LXLXDS( NTTRAN , KNOMTR )
      ENDIF
      GOTO 100
C
C     FIN DE DEFINITION DES TRANSFORMATIONS
 9000 RETURN
      END
