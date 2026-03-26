      SUBROUTINE DFLIGN
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR LES LIGNES C-A-D
C------ CREER ET LIRE LE TABLEAU ~>LIGNE>...>DEFINITION
C       GENERER       LE TABLEAU ~>LIGNE>...>XYZSOMMET
C       GENERER       LE TABLEAU ~>LIGNE>...>NSEF
C       TRACER LES LIGNES SI LA CONSOLE EST INTERACTIVE GRAPHIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1988
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/xyzext.inc"
      include"./incl/pilect.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
C
C     DECLARATION DU SUPER-TABLEAU NUMERIQUE MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
C     LES VARIABLES LOCALES
      CHARACTER*24  KNOMLI
      CHARACTER*80  KNMTD,KNMTS
C
C     SAISIE DU NOM_DE_LA_LIGNE OU  LE CARACTERE @ POUR FINIR
C     =======================================================
C     AFFICHAGE DE 'NOM DE LA LIGNE' COMME INVITE
 100  CALL INVITE( 41 )
C     LECTURE D'UNE CHAINE DE CARACTERES SELON LE LANGAGE LU
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMLI )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C
C     LE RECTANGLE HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
C     AFFICHAGE DANS LE RECTANGLE HISTORIQUE DU NOM DE LA LIGNE
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'LIGNE: ' // KNOMLI
      ELSE
         KHIST(1) = 'LINE: ' // KNOMLI
      ENDIF
      CALL LHISTO
      ILEXIS = 1
C
C     RECHERCHE DU NOM DE LA LIGNE DANS LE LEXIQUE DES LIGNES
 150  IF( NTLIGN .LE. 0 ) GOTO 9000
      CALL LXLXOU( NTLIGN, KNOMLI, NTLXLI, MNLXLI )
C
      IF( NTLXLI .LE. 0 ) THEN
C        LA LIGNE N'EXISTE PAS. ELLE EST CREEE
         ILEXIS = 0
         CALL LXLXDC( NTLIGN, KNOMLI, 24, 8 )
         GOTO 150
      ENDIF
C
C     MODIFICATION OU CREATION DU TMS DE DEFINITION DE LA LIGNE
      KNMTD = '~>LIGNE>>DEFINITION'
      I     = INDEX( KNOMLI, ' ' )
      IF( I .EQ. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT DE LIGNE: '//KNOMLI
         ELSE
            KERR(1) = 'INCORRECT NAME of LINE: '//KNOMLI
         ENDIF
         CALL LEREUR
         GOTO 100
      ELSE IF( I .EQ. 0 ) THEN
         I = 25
      ENDIF
C
C     ========================================================
C     SAISIE SELON LE TD DES DONNEES DE DEFINITION DE LA LIGNE
C     ========================================================
      KNMTS = '~>LIGNE>' // KNOMLI(1:I-1) // '>DEFINITION'
      CALL MOTSTD( KNMTD, KNMTS, NRETOU )
C
      IF( NRETOU .NE. 0 ) THEN
C        PB => DESTRUCTION SI LA LIGNE N'EXISTAIT PAS AVANT
         IF( ILEXIS .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'DESTRUCTION TMS:'//KNMTS
            CALL LXLXDS( NTLIGN, KNOMLI )
         ENDIF
         GOTO 100
      ENDIF
C
C     REOUVERTURE DU TMS '~>LIGNE>NOM_DE_LA_LIGNE>DEFINITION'
      CALL LXTSOU( NTLXLI, 'DEFINITION', NTDFLI, MNDFLI )
      IF( NTDFLI .LE. 0 ) THEN
C        TABLEAU 'DEFINITION' INEXISTANT
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIGNE NON DEFINIE: '//KNOMLI
         ELSE
            KERR(1) = 'LINE NOT DEFINED: '//KNOMLI
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
C
C     MISE A JOUR DANS LE CAS D'UNE DONNEE DIRECTE DES TABLEAUX
C     'XYZSOMMET' ET 'NSEF'
      NUTYLI = MCN( MNDFLI + WUTYLI )
      IF( NUTYLI .EQ. 10 ) THEN
         CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
C        MISE A JOUR DE SA DATE DE CREATION
         CALL ECDATE( MCN( MNSOLI ) )
         CALL LXTSOU( NTLXLI, 'NSEF', NTNSEF, MNNSEF )
C        MISE A JOUR DE SA DATE DE CREATION
         CALL ECDATE( MCN( MNNSEF) )
      ENDIF
C
C     LE CALCUL DES VARIABLES DU TABLEAU 'XYZSOMMET' DE LA LIGNE
C                             DU TABLEAU 'NSEF'      DE LA LIGNE
      IF( INTERA .LE. 0 ) THEN
C        SI CONSOLE NON INTERACTIVE GRAPHIQUE => CALCUL DIRECT
         CALL MAILEX( 'LIGNE', KNOMLI,
     %                 NTNSEF, MNNSEF, NTSOLI, MNSOLI, IERR )
         IF( IERR .NE. 0 ) THEN
C           DESTRUCTION DU LEXIQUE DE LA LIGNE NON MAILLEE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ATTENTION LA LIGNE ' // KNOMLI
               KERR(2) = 'EST DETRUITE SUITE A PB LORS DU MAILLAGE'
            ELSE
               KERR(1) = 'ATTENTION: the LINE ' // KNOMLI
               KERR(2) = 'is DELETED from the PROBLEM SEEN BEFORE'
            ENDIF
            CALL LEREUR
            CALL LXLXDS( NTLIGN, KNOMLI )
         ENDIF
      ELSE
C        SI CONSOLE INTERACTIVE GRAPHIQUE
C        AFFICHAGE DU TABLEAU DEFINITION
         CALL AFTSTD( KNMTS )
C        => CALCUL INDIRECT LORS DE LA DEMANDE DE TRACE DE LA LIGNE
         NBC = INDEX( KNOMLI, ' ' ) - 1
         IF( NBC .LE. 0 ) NBC = LEN( KNOMLI )
         CALL NUOBNM( 'LIGNE', KNOMLI, NULIGN )
         CALL T1MOBJ( 'LIGNE', KNOMLI(1:NBC), NULIGN )
C
         IF( INIEXT .EQ. 2 ) THEN
C           TRACE D'UNE LIGNE QUI A MODIFIE MIN MAX COOEXT
            CALL EFFACEMEMPX
            CALL ITEMS0
            CALL VISEE0
            CALL TRAXES
            CALL TRTOBJ( 'POINT', 0 )
            CALL TRTOBJ( 'LIGNE', 0 )
            CALL TRFINS( 'LIGNE' )
         ENDIF
      ENDIF
      GOTO 100
C
C     FIN DE DEFINITION DES LIGNES
C     ============================
C     LE RECTANGLE HISTORIQUE EST EFFACE
 9000 CALL RECTEF( NRHIST )
      END
