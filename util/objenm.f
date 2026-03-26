      SUBROUTINE OBJENM( NMTOBJ, NMOBJT, NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE LE NOM D'UN OBJET
C -----    RETOURNER SON NOM APRES VERIFICATION DE SON EXISTENCE
C          DANS LE LEXIQUE DE SON TYPE D'OBJETS
C
C ENTREES:
C --------
C NMTOBJ : CHAINE DE CARACTERES DEFINISSANT LE TYPE DE L'OBJET
C          'POINT' OU 'LIGNE' OU 'SURFACE' OU 'VOLUME' OU 'OBJET'
C
C SORTIES:
C --------
C NMOBJT : NOM DE L'OBJET
C NUOBJT : > 0 NUMERO DE L'OBJET DANS SON LEXIQUE
C          = 0 SI NOM NON RETROUVE
C          =-1 SI LE NOM EST '@' POUR RETOUR SANS CALCUL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE   PARIS  DECEMBRE 1985
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      CHARACTER*(*)     NMOBJT
      CHARACTER*(*)     NMTOBJ
C
C     LE NUMERO DU TABLEAU LEXIQUE DANS LE COMMON / NTMNLT /
      NUOBJT = 0
      NTYOB  = NTYOBJ( NMTOBJ )
      IF( NTYOB .LE. 0 ) RETURN
      NTOBJT = NTMN( NTYOB )
C
C     LECTURE DU NOM DE L'OBJET
 10   GOTO( 21, 22, 23, 24, 25, 26, 27 ),NTYOB
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
      CALL LIRCAR( NCVALS , NMOBJT )
      IF( NCVALS .EQ. -1 ) THEN
         NUOBJT = -1
         RETURN
      ENDIF
C
C     RECHERCHE DU NOM DANS LE LEXIQUE
      CALL LXNMNO( NTOBJT, NMOBJT, NUOBJT, MNOBJT )
C
C     L'OBJET A T IL ETE RETROUVE ?
      IF( NUOBJT .LE. 0 ) THEN
C
C        NON. L'OBJET N'EXISTE PAS
C        -------------------------
         N = NUDCNB( NMOBJT )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19000) NMTOBJ, NMOBJT(1:N)
         ELSE
            WRITE(IMPRIM,29000) NMTOBJ, NMOBJT(1:N)
         ENDIF
19000    FORMAT(1X,A,1X,A,' NON RETROUVE')
29000    FORMAT(1X,A,1X,A,' NOT RECOVERED')
C
C        CONSTRUCTION DU TEXTE DU RECTANGLE DE L'ERREUR
         NBLGRC(NRERR) = 2
         KERR(1) =  NMTOBJ // ' : ' // NMOBJT(1:N)
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = '... NOM NON RETROUVE PARMI ...'
         ELSE
            KERR(2) = '... NAME NOT RECOVERED AMONG ...'
         ENDIF
C
C        AFFICHAGE DES NOMS DU LEXIQUE AVEC APPEL DU SP LEREUR
         CALL LXIM0( MNOBJT )
C
         IF( INTERA .GE. 3 ) GOTO 10
         NUOBJT = 0
C
      ENDIF
C
      RETURN
      END
