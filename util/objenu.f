      SUBROUTINE OBJENU( NMTOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE LE NOM D'UN OBJET
C -----    RETOURNER SON NUMERO DANS LE LEXIQUE OBJET
C
C ENTREES:
C --------
C NMTOBJ : CHAINE DE CARACTERES DEFINISSANT LE TYPE DE L'OBJET
C          'POINT' OU 'LIGNE' OU 'SURFACE' OU 'VOLUME' OU 'OBJET'
C
C SORTIE :
C --------
C NUOBJT : NUMERO DANS LE LEXIQUE DE L'OBJET DE NOM LU
C          =>  0 SI TRAITEMENT BATCH ET NOM NON RETROUVE
C          => -1 SI LE NOM EST '<' POUR RETOUR SANS CALCUL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE   PARIS  DECEMBRE 1985
C..............................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / NOMNOM / KNOM
      CHARACTER*24      KNOM
      CHARACTER*(*)     NMTOBJ
      CHARACTER*8       MMTOBJ
C
C     LE NUMERO DU TABLEAU LEXIQUE DANS LE COMMON / NTMNLT /
      N = NTYOBJ( NMTOBJ )
      NTOBJT = NTMN( N )
C
C     LECTURE DU NOM DE L'OBJET
 10   GOTO( 21, 22, 23, 24, 25 ),N
 21   CALL INVITE( 51 )
      GOTO 29
 22   CALL INVITE( 40 )
      GOTO 29
 23   CALL INVITE( 42 )
      GOTO 29
 24   CALL INVITE( 60 )
      GOTO 29
 25   CALL INVITE( 45 )
 29   NCVALS = 0
      CALL LIRCAR( NCVALS , KNOM )
      IF( NCVALS .EQ. -1 ) THEN
C         LE NOM DE L'OBJET EST '@' => RETOUR
          NUOBJT = -1
          RETURN
      ENDIF
C
C     RECHERCHE DU NOM DANS LE LEXIQUE
      CALL LXNMNO( NTOBJT , KNOM , NUOBJT , MNOBJT )
C
C     L'OBJET A T IL ETE RETROUVE ?
      IF( NUOBJT .LE. 0 ) THEN
C
C        NON.L'OBJET N'EXISTE PAS
C        -------------------------
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  NMTOBJ//' INCONNU:'//KNOM
            WRITE(IMPRIM,19000) NMTOBJ,KNOM
19000       FORMAT(1X,A7,1X,A,' NON RETROUVE PARMI ...')
         ELSE
            MMTOBJ = NMTOBJ
            KERR(1) =  MMTOBJ // ' UNKNOWN:' // KNOM
            IF( INDEX( MMTOBJ, 'LIGNE' ) .GT. 0 ) MMTOBJ = 'LINE'
            IF( INDEX( MMTOBJ, 'OBJET' ) .GT. 0 ) MMTOBJ = 'OBJECT'
            WRITE(IMPRIM,29000) MMTOBJ, KNOM
29000       FORMAT(1X,A7,1X,A,' NOT RECOVERED AMONG ...')
         ENDIF
         CALL LEREUR
         CALL LXIM0( MNOBJT )
         IF( INTERA .GE. 3 ) GOTO 10
         NUOBJT = 0
      ENDIF
C
      RETURN
      END
