      SUBROUTINE TUER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SUIVI DES TMS OBJETS FONCTIONS DE MEFISTO
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C2345X7..............................................................012
      include"./incl/pp.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*160     KNOM
C
C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
 10   CALL LIMTCL( 'tuer' , NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
      GOTO( 100 , 200 , 300 , 400 )  , NMTCL
C
C     DESTRUCTION D'UN TMS
C     ====================
 100  CALL INVITE( 56 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C     RECHERCHE DU NUMERO DE TMS DE CE TMS
      CALL NUMTMS( KNOM , NUTMS )
C     DESTRUCTION DU TMS
      IF( NUTMS .GT. 0 ) THEN
         CALL TAMSDS( NUTMS )
      ENDIF
      GOTO 10
C
C     DESTRUCTION D'UN OBJET
C     ======================
C     LECTURE DE SON TYPE
 200  CALL LIMTCL( 'typ_objt' , NUTYOB )
      IF( NUTYOB .LE. 0 ) GOTO 10
C     LE LEXIQUE DE CE TYPE D'OBJET
      NTOBJ = NTMN( NUTYOB )
C
      GOTO( 21, 22, 23, 24, 25, 26, 27 ),NUTYOB
 21   CALL INVITE( 51 )
      GOTO 29
 22   CALL INVITE( 40 )
      GOTO 29
 23   CALL INVITE( 42 )
      GOTO 29
 24   CALL INVITE( 60 )
      GOTO 29
 25   CALL INVITE( 45 )
      GOTO 29
 26   CALL INVITE( 38 )
      GOTO 29
 27   CALL INVITE( 44 )
C
 29   NCVALS = 0
      CALL LIRCAR( NCVALS , KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C     RECHERCHE DU NOM DANS LE LEXIQUE
      CALL LXNMNO( NTOBJ , KNOM , NUMOB , MNOBJ )
C
C     L'OBJET A T IL ETE RETROUVE ?
      IF( NUMOB .LE. 0 ) THEN
C        NON.L'OBJET N'EXISTE PAS
         NUMOB = NUDCNB( KNOM )
         KERR(1) = KNOM(1:NUMOB)  // ' NON RETROUVE'
         NBLGRC(NRERR) = 1
         CALL LERESU
         CALL LXIM0( MNOBJ )
      ELSE
C        OUI.L'OBJET EST DETRUIT
         CALL LXLXDS( NTOBJ , KNOM )
      ENDIF
      GOTO 10
C
C     DESTRUCTION DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C     ===============================================
C     RECHERCHE DU NOM TAILLE_IDEALE DANS LE LEXIQUE
 300  CALL LXNMNO( NTFONC , 'TAILLE_IDEALE' , NUMOB , MNOBJ )
C
C     LE LEXIQUE DE CETTE FONCTION A T IL ETE RETROUVE?
      IF( NUMOB .LE. 0 ) THEN
C        NON.LA FONCTION N'EXISTE PAS
         KERR(1) = 'TAILLE_IDEALE(X,Y,Z) NON RETROUVEE'
         NBLGRC(NRERR) = 1
         CALL LERESU
         CALL LXIM0( MNOBJ )
      ELSE
C        OUI.LA FONCTION EST DETRUITE
         CALL LXLXDS( NTFONC , 'TAILLE_IDEALE' )
      ENDIF
      GOTO 10
C
C     DESTRUCTION DE LA FONCTION EDGE_LENGTH(X,Y,Z)
C     =============================================
C     RECHERCHE DU NOM EDGE_LENGTH DANS LE LEXIQUE
 400  CALL LXNMNO( NTFONC , 'EDGE_LENGTH' , NUMOB , MNOBJ )
C
C     LE LEXIQUE DE CETTE FONCTION A T IL ETE RETROUVE?
      IF( NUMOB .LE. 0 ) THEN
C        NON.LA FONCTION N'EXISTE PAS
         KERR(1) = 'EDGE_LENGTH(X,Y,Z) UNFOUND'
         NBLGRC(NRERR) = 1
         CALL LERESU
         CALL LXIM0( MNOBJ )
      ELSE
C        OUI.LA FONCTION EST DETRUITE
         CALL LXLXDS( NTFONC , 'EDGE_LENGTH' )
      ENDIF
      GOTO 10
C
 9000 RETURN
      END
