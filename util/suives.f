      SUBROUTINE SUIVES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MODIFIER LES AFFICHAGES OU MODE D'ENTREE DES DONNEES
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/lu.inc"
      include"./incl/pilect.inc"
      include"./incl/trvari.inc"
      CHARACTER*(NBCALI) NOM
      COMMON / MSIMTA /  NOIMPR
      COMMON / BAVARD /  IBAVAR
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NUNITE(29)
      LOGICAL            LOGIQ
C
C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
 10   CALL LIMTCL( 'entrsort' , NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
      GOTO( 100 , 200 , 300 , 400 , 500 , 600 ) , NMTCL
C
C     MODE D''ENTREE DES DONNEES
C     ==========================
 100  CALL LIMTCL( 'mode_es' , NMTCL1 )
C     0: 'Sans INTERACTION (BATCH aveugle)'
C     1: 'InterActif avec DESSINS sans CLAVIER sans SOURIS'
C     3: 'InterActif avec DESSINS avec CLAVIER avec SOURIS' ) ;
      IF( NMTCL1 .LT. 0 ) GOTO 10

      IF( NMTCL1 .EQ. 0 ) THEN
C        0: 'Sans INTERACTION (BATCH aveugle)'
         INTERA = 0
         INTERB(LHLECT) = INTERA
         GOTO 10
      ENDIF

      GOTO( 110, 10, 110, 10, 10 ) , NMTCL1
 110  INTERA = MIN(3,NMTCL1)
C     RE-OUVERTURE DU LOGICIEL LU ET xvue SELON LA NOUVELLE INTER ACTIVITE
      INTERB(LHLECT) = INTERA
CCCC     FERMETURE DES TRACES ET SAISIES DE VOIR
CCC      CALL XVFERMER
CCCC     REOUVERTURE DU LOGICIEL xvue
CCC      CALL XVINIT
CCC      CALL LUOU
      GOTO 10
C
C     BAVARDAGE
C     =========
 200  CALL INVITE( 35 )
      NCVALS = 4
      N      = IBAVAR
      CALL LIRENT( NCVALS , N )
      IF( N .LE. 0 ) GOTO 10
      IBAVAR = MAX( 0 , MIN( N , 10 ) )
      GOTO 10
C
C     ADRESSAGE
C     =========
 300  IF( NOIMPR .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: AFFICHAGE DES ADRESSAGES'
         CALL LERESU
         NOIMPR = 1
      ELSE
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: FIN D''AFFICHAGE DES ADRESSAGES'
         CALL LERESU
         NOIMPR = 0
      ENDIF
      GOTO 10
C
C     ECHO
C     ====
 400  IF( LUIMPR .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: DEBUT D''ECHO DES DONNEES'
         CALL LERESU
         LUIMPR = 1
      ELSE
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: FIN D''ECHO DES DONNEES'
         CALL LERESU
         LUIMPR = 0
      ENDIF
      GOTO 10
C
C     REDIRECTION DES AFFICHAGES
C     ==========================
 500  CALL LIMTCL( 'affiche' , NMTCL1 )
      IF( NMTCL1 .LE. 0 ) GOTO 10
      GOTO( 510 , 520 ) , NMTCL1
C
C     ECRAN
 510  IMPRIM = IINFO( 'IMPRIMANTE INITIALE' )
      GOTO 10
C
C     FICHIER DE NOM ?
 520  CALL INVITE( 49 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , NOM )
      IF( NCVALS .LE. 0 ) GOTO 10
C     RECHERCHE D'UN NUMERO D'UNITE LIBRE
      CALL TRUNIT( NF )
C     L'IMPRIMANTE DEVIENT LE FICHIER NF
      CALL MINUSC( NOM )
      OPEN(FILE=NOM,UNIT=NF,ACCESS='SEQUENTIAL',
     %     FORM='FORMATTED',STATUS='UNKNOWN',
     %     IOSTAT= N )
      IF( N .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) =  'LU: NOM INCORRECT DE FICHIER'
         KERR(2) = NOM
         CALL LEREUR
         GOTO 520
      ENDIF
C     LE NOUVEAU FICHIER SUPPORT DES AFFICHAGES
      IMPRIM = NF
      GOTO 10
C
C     LECTURE
C     =======
 600  CALL LIMTCL( 'lecteur' , NMTCL1 )
      IF( NMTCL1 .LE. 0 ) GOTO 10
      GOTO( 610 , 620 , 630 ) , NMTCL1
C
C     1: 'CLAVIER classique'
 610  LECTEU = IINFO( 'LECTEUR INITIAL' )
      LHLECT = 1
      INTERA = INTERB( LHLECT )
      GOTO 10
C
C     2: 'FICHIER a son DEBUT (avec rembobinage)'
C     FICHIER DE NOM ?
 620  CALL INVITE( 49 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , NOM )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL MINUSC( NOM )
      INQUIRE( FILE=NOM , OPENED=LOGIQ , NUMBER=NF , IOSTAT=N )
      IF( LOGIQ ) GOTO 622
C     RECHERCHE D'UN NUMERO D'UNITE LIBRE
      CALL TRUNIT( NF )
      OPEN(FILE=NOM,UNIT=NF,ACCESS='SEQUENTIAL',
     %     FORM='FORMATTED',STATUS='OLD',
     %     IOSTAT= N )
 622  IF( N .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) =  'LU: NOM INCONNU DE FICHIER'
         KERR(2) = NOM
         CALL LEREUR
         GOTO 620
      ENDIF
      IF (LHLECT.GE.MXLECT) THEN
          NBLGRC(NRERR) = 1
          KERR(1) ='PILE SATUREE DES LECTEURS. AUGMENTER MXLECT'
          CALL LEREUR
          GOTO 10
      ENDIF
      REWIND NF
      LECTEU = NF
      LHLECT = LHLECT + 1
      LPLECT(LHLECT) = LECTEU
C     PAS DE VISUALISATION GRAPHIQUE LORS DE LECTURE SUR FICHIER
      INTERA = 1
      INTERB(LHLECT) = 1
      GOTO 9000
C
C     3: 'FICHIER dans l''ETAT de son dernier abandon' ) ;
C     FICHIER DE NOM ?
 630  CALL INVITE( 49 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , NOM )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL MINUSC( NOM )
      INQUIRE( FILE=NOM , OPENED=LOGIQ , NUMBER=NF , IOSTAT=N )
      IF( N .NE. 0 .OR. NF .LE. 0 .OR. (.NOT.LOGIQ) ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NOM
         KERR(2) = 'LU: FICHIER NON OUVERT'
         CALL LEREUR
         GOTO 600
      ENDIF
      IF (LHLECT.GE.MXLECT) THEN
         NBLGRC(NRERR) = 1
         KERR(1) ='PILE SATUREE DES LECTEURS. AUGMENTER MXLECT'
         CALL LEREUR
         GOTO 10
      ENDIF
      LECTEU = NF
      LHLECT = LHLECT + 1
      LPLECT(LHLECT) = LECTEU
C     PAS DE VISUALISATION GRAPHIQUE LORS DE LECTURE SUR FICHIER
      INTERA = 1
      INTERB(LHLECT) = 1
      GOTO 9000
C
 9000 RETURN
      END
