      SUBROUTINE DFPOIN
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR LES POINTS C-A-D
C------ CREER ET LIRE LE TABLEAU ~>POINT>...>DEFINITION
C       GENERER       LE TABLEAU ~>POINT>...>XYZSOMMET
C       TRACER LES LIGNES SI LA CONSOLE EST INTERACTIVE GRAPHIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C....................................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/nisafr.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      CHARACTER*24      KNOMPO
      CHARACTER*80      KNMTD,KNMTS
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE MODE DE SAISIE DES POINTS
      CALL LIMTCL( 'saisi_pt', I )
      IF( I .EQ. -1 ) GOTO 9000
      IF( INTERA .LT. 3 ) THEN
         IF( I .NE. 1 ) THEN
CCC            NBLGRC(NRERR) = 3
CCC            KERR(1) = 'SAISIE AVEC SOURIS INTERDITE'
CCC            KERR(2) = 'EN MODE SANS SOURIS'
CCC            KERR(3) = 'FRAPPE DES VALEURS FORCEE'
CCC            CALL LEREUR
C           PASSAGE FORCE PAR LE TYPE 1
            I  = 1
         ENDIF
      ELSE
C        INTERACTIVITE AVEC SOURIS
         IF( I .GT. 1 ) THEN
C           LA FRAPPE SUR LE FICHIER FRAPPE EST SUPPRIMEE
            NISAFR = 0
         ENDIF
      ENDIF
C
      GOTO( 100, 2000, 3000 ), I
C
C     MODE FRAPPE DES VALEURS
C     =======================
C     NOM_DU_POINT  LE CARACTERE @ POUR FINIR
 100  CALL INVITE( 55 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMPO )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C
C     TRACE DU NOM DU POINT
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'POINT: ' // KNOMPO
      CALL LHISTO
      ILEXIS = 1
C
C     RECHERCHE DU NOM DU POINT DANS LE LEXIQUE DES POINTS
 150  IF( NTPOIN .LE. 0 ) GOTO 9000
      CALL LXLXOU( NTPOIN, KNOMPO, NTLXPO, MNLXPO )
C
C     S'IL N'EXISTE PAS IL EST CREE
      IF( NTLXPO .LE. 0 ) THEN
         ILEXIS = 0
         CALL LXLXDC( NTPOIN, KNOMPO, 24, 8 )
         GOTO 150
      ENDIF
C
C     MODIFICATION OU CREATION DE LA DEFINITION DU POINT
C     ==================================================
      KNMTD = '~>POINT>>DEFINITION'
      I     = INDEX( KNOMPO, ' ' )
      IF( I .EQ. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT DE POINT: ' // KNOMPO
         ELSE
            KERR(1) = 'INCORRECT NAME of POINT: ' // KNOMPO
         ENDIF
         CALL LEREUR
         GOTO 100
      ELSE IF( I .EQ. 0 ) THEN
         I = 25
      ENDIF
      KNMTS = '~>POINT>' // KNOMPO(1:I-1) // '>DEFINITION'
      CALL MOTSTD( KNMTD, KNMTS, NRETOU )
      IF( NRETOU .NE. 0 ) THEN
C        DESTRUCTION SI LE POINT N'EXISTAIT PAS AVANT
         IF( ILEXIS .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'DESTRUCTION TMS:'//KNMTS
            CALL LXLXDS( NTPOIN, KNOMPO )
         ENDIF
         GOTO 100
      ENDIF
C
C     CREATION DU POINT SELON LE TABLEAU DEFINITION
C     OUVERTURE DU TABLEAU  'DEFINITION'
      CALL LXTSOU( NTLXPO, 'DEFINITION', NTDFPO, MNDFPO )
      IF( NTDFPO .LE. 0 ) THEN
C        TABLEAU 'DEFINITION' INEXISTANT
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINT NON DEFINI: '//KNOMPO
         ELSE
            KERR(1) = 'POINT NOT DEFINED: '//KNOMPO
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
C
C     LA GENERATION DU TABLEAU 'XYZSOMMET' DU POINT
C     ---------------------------------------------
      IF( INTERA .LE. 0 ) THEN
C        SI CONSOLE NON INTERACTIVE GRAPHIQUE
C        ALORS GENERATION DU TABLEAU  'XYZSOMMET'
         CALL MAILEX( 'POINT', KNOMPO,
     %                 NT, MN, NTSOPO, MNSOPO, IERR )
         IF( IERR .NE. 0 ) THEN
C           DESTRUCTION DU LEXIQUE DU POINT NON MAILLE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ATTENTION LE POINT ' // KNOMPO
               KERR(2) = 'EST DETRUIT SUITE A PB LORS DU MAILLAGE'
            ELSE
               KERR(1) = 'ATTENTION: the POINT ' // KNOMPO
               KERR(2) = 'is DELETED from the PROBLEM SEEN BEFORE'
            ENDIF
            CALL LEREUR
            CALL LXLXDS( NTPOIN, KNOMPO )
         ENDIF
      ELSE
C        SI CONSOLE INTERACTIVE GRAPHIQUE ALORS
C        AFFICHAGE DU TABLEAU DEFINITION
         CALL AFTSTD( KNMTS )
C        GENERATION DU TABLEAU 'XYZSOMMET'
C        LORS DU TRACE AUTOMATIQUE DE TOUS LES POINTS
         CALL ITEMS0
         CALL EFFACEMEMPX
         CALL TRTOBJ( 'POINT', 0 )
         CALL TRAXES
         CALL TRFINS( 'POINT' )
      ENDIF
      GOTO 100
C
C     MODE DE SAISIE DE POINTS 2D
C     ===========================
 2000 CALL SAPT2D
      GOTO 9000
C
C     MODE DE SAISIE DE POINTS 3D
C     ===========================
 3000 CALL SAPT3D
C
C     FIN DE DEFINITION DES POINTS ET TRACE DE TOUS LES POINTS
C     ========================================================
 9000 IF( INTERA .GE. 1 ) THEN
         CALL EFFACEMEMPX
         CALL ITEMS0
         CALL VISEE0
         CALL TRTOBJ( 'POINT', 0 )
         CALL TRAXES
         CALL TRFINS( 'POINT' )
      ENDIF
C
C     LA SAUVEGARDE SUR LE FICHIER FRAPPE EST RESTAUREE
      NISAFR = 1
      CALL RECTEF( NRHIST )
      RETURN
      END
