      SUBROUTINE DFOBJE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  DEFINIR LES OBJETS c-a-d
C------  CREER ET LIRE LE TABLEAU ~>OBJET>...>DEFINITION
C        TRACER LES OBJETS SI LA CONSOLE EST INTERACTIVE GRAPHIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1989
C....................................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      CHARACTER*10      NMTYOB, KNOMTY
      CHARACTER*24      KNOMOB, KNOM
      CHARACTER*80      KNMTD,  KNMTS
      INTEGER           SDOUNO

C     NOM_DE_L'OBJET  LE CARACTERE @ POUR FINIR
 100  CALL INVITE( 48 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMOB )
      IF( NCVALS .EQ. -1 ) GOTO 9000

C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
      ILEXIS = 1

C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      IF( NTOBJE .LE. 0 ) GOTO 9000
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )

      IF( NTLXOB .LE. 0 ) THEN
C        L'OBJET KNOMOB N'EXISTE PAS: IL EST CREE
         ILEXIS = 0
         CALL LXLXDC( NTOBJE, KNOMOB, 24, 8 )
         CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
         GOTO 200
      ELSE
C        L'OBJET KNOMOB EXISTE:
         IF( LHLECT .GT. 1 ) GOTO 200

C        L'OBJET EST TRACE PUIS MODIFIE
 5       IF( INTERA .GE. 1 ) THEN
C           TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
            CALL EFFACEMEMPX
            CALL VISEE0
            CALL T1OBJE( KNOMOB )
         ENDIF

C        OUVERTURE DU TABLEAU 'DEFINITION' DE L'OBJET
         CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
         IF( NTDFOB .LE. 0 ) THEN
C           TABLEAU 'DEFINITION' INEXISTANT
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DEFINITION INCONNUE de l''OBJET: ' // KNOMOB
            ELSE
               KERR(1) = 'UNKNOWN DEFINITION of OBJECT: ' // KNOMOB
            ENDIF
            CALL LEREUR
            GOTO 200
         ENDIF
         SDOUNO = MCN( MNDFOB+WDOUNO )
         IF( SDOUNO .NE. 0 ) GOTO 200
         NBDOBJ = MCN(MNDFOB+WBDOBJ)
         IF( NBDOBJ .LE. 0 ) GOTO 200

         CALL LIMTCL( 'modifobj', NMTCL0 )
         IF( NMTCL0 .LT. 0 ) RETURN
         IF( NMTCL0 .EQ. 0 ) GOTO 200
         GOTO( 10, 20, 30 ), NMTCL0

C        SUPPRIMER cet OBJET
C        -------------------
 10      CALL LXLXDS( NTOBJE, KNOMOB )
         IF( INTERA .GE. 1 ) CALL EFFACEMEMPX
         RETURN

C        AJOUTER un PLSV a L'OBJET
C        -------------------------
 20      CALL LIRTPLSV( NOTYPE, KNOM, NUPLSV )
         IF( NOTYPE .LT. 0 ) RETURN

C        RECHERCHE du PLSV DANS L'OBJET
         MNTYOB = MNDFOB + WTYOBJ
         DO 25 I=1,NBDOBJ
C           LE NUMERO DU TYPE DE PLSV
            IF( MCN( MNTYOB ) .EQ. NOTYPE ) THEN
C              LE NUMERO DU PLSV DANS SON LEXIQUE
               IF( MCN( MNTYOB + 1 ) .EQ. NUPLSV ) THEN
C                 LE I-EME PLSV A AJOUTER EXISTE DEJA DANS L'OBJET
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LE PLSV ' // KNOM
                     KERR(2) = 'EXISTE DEJA DANS OBJET' // KNOMOB
                  ELSE
                     KERR(1) = 'The PLSV ' // KNOM
                     KERR(2) = 'ALREADY EXISTS in OBJECT' // KNOMOB
                  ENDIF
                  CALL LEREUR
                  GOTO 20
               ENDIF
            ENDIF
            MNTYOB = MNTYOB + 2
 25      CONTINUE

C        PLSV NON RETROUVE DANS L'OBJET => IL EST AJOUTE EN AUGMENTANT DEFINITIO
         L = MOTVAR( NUMTYP('TYPEOBJET') )
         CALL TAMSAU( NTDFOB, WTYOBJ + NBDOBJ*L + L )
         CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
         MNTYOB = MNDFOB + WTYOBJ + NBDOBJ * L
C        AJOUT DU PLSV DANS LA DEFINITION
         MCN( MNTYOB     ) = NOTYPE
         MCN( MNTYOB + 1 ) = NUPLSV
         NBDOBJ = NBDOBJ + 1
         MCN(MNDFOB+WBDOBJ) = NBDOBJ
C        MISE EN ORDRE V, S, L, P
         MNTYOB = MNDFOB + WTYOBJ
         CALL VSLPOB( NBDOBJ, MNTYOB )
C        LA DATE DU TMS DEFINITION DE L'OBJET
         CALL ECDATE( MCN(MNDFOB) )
         GOTO 5

C        RETIRER un PLSV de l'OBJET
C        --------------------------
 30      CALL LIRTPLSV( NOTYPE, KNOM, NUPLSV )
         IF( NOTYPE .LT. 0 ) RETURN

C        RECHERCHE du PLSV DANS L'OBJET
         MNTYOB = MNDFOB + WTYOBJ
         DO 35 I=1,NBDOBJ
C           LE NUMERO DU TYPE DE PLSV
            IF( MCN( MNTYOB ) .EQ. NOTYPE ) THEN
C              LE NUMERO DU PLSV DANS SON LEXIQUE
               IF( MCN( MNTYOB + 1 ) .EQ. NUPLSV ) THEN
C                 LE I-EME PLSV EST A SUPPRIMER
C                 DECALAGE VERS LE DEBUT DU TABLEAU
                  MN = MNTYOB
                  DO 32 J=I+1,NBDOBJ
                     MCN( MN     ) = MCN( MN + 2 )
                     MCN( MN + 1 ) = MCN( MN + 3 )
                     MN = MN + 2
 32               CONTINUE
C                 UN PLSV DE MOINS DANS L'OBJET
                  MCN(MNDFOB+WBDOBJ) = NBDOBJ - 1
C                 LA DATE DU TMS DEFINITION DE L'OBJET
                  CALL ECDATE( MCN(MNDFOB) )
                  GOTO 5
               ENDIF
            ENDIF
            MNTYOB = MNTYOB + 2
 35      CONTINUE
C        PLSV NON RETROUVE DANS L'OBJET
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PLSV: ' // KNOMOB // ' NON RETROUVE DANS OBJET'
         ELSE
            KERR(1) = 'PLSV: ' // KNOMOB // ' UNKNOWN in OBJECT'
         ENDIF
         CALL LEREUR
         GOTO 5

      ENDIF

C     MODIFICATION OU CREATION DE LA DEFINITION DE L'OBJET
C     ====================================================
 200  KNMTD = '~>OBJET>>DEFINITION'
      I = INDEX( KNOMOB, ' ' )
      IF( I .EQ. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT D''OBJET: '//KNOMOB
         ELSE
            KERR(1) = 'INCORRECT NAME of OBJECT: '//KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 100
      ELSE IF( I .EQ. 0 ) THEN
         I = 25
      ENDIF
      KNMTS = '~>OBJET>' // KNOMOB(1:I-1) // '>DEFINITION'
      CALL MOTSTD( KNMTD, KNMTS, NRETOU )
      IF( NRETOU .NE. 0 ) THEN
C        DESTRUCTION SI L'OBJET N'EXISTAIT PAS AVANT
         IF( ILEXIS .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'DESTRUCTION TMS:' // KNMTS
            CALL LXLXDS( NTOBJE, KNOMOB )
         ENDIF
         GOTO 100
      ENDIF

C     CREATION DE L'OBJET SELON LE TABLEAU DEFINITION
C     OUVERTURE DU TABLEAU  'DEFINITION'
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
      IF( NTDFOB .LE. 0 ) THEN
C        TABLEAU 'DEFINITION' INEXISTANT
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET NON DEFINI: ' // KNOMOB
         ELSE
            KERR(1) = 'UNDEFINED OBJECT: ' // KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF

C     TRAITEMENT CLASSIQUE(-1) SANS RENUMEROTATION DES NOEUDS
C     TRAITEMENT CLASSIQUE( 0) AVEC RENUMEROTATION DES NOEUDS
C     TRAITEMENT PAR SOUS-DOMAINES(1)
C     TRAITEMENT PAR JOINTS(2) ?
      IF( MCN( MNDFOB+WDOUNO ) .GT. 0 ) THEN
C        TRAITEMENT CLASSIQUE FORCE
         MCN( MNDFOB+WDOUNO ) = 0
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SOUS-DOMAINE NON OPERATIONNEL'
            KERR(2) = 'TRAITEMENT CLASSIQUE FORCE'
            KERR(3) = 'AVEC RENUMEROTATION DES NOEUDS'
         ELSE
            KERR(1) = 'SUB-DOMAIN NOT OPERATIONAL'
            KERR(2) = 'CLASSICAL TRAITMENT FORCED'
            KERR(3) = 'WITH NODES RENUMBERING'
         ENDIF
         CALL LERESU
      ENDIF
      SDOUNO = MCN( MNDFOB+WDOUNO )

C     PROTECTION POUR LES OBJETS AVEC VOLUMES
C     SOUS-DOMAINES ET JOINTS INTERDITS
      MN = MNDFOB + WTYOBJ
      DO 220 I=1,MCN(MNDFOB+WBDOBJ)
C        LE NO DU TYPE ET LE NUMERO DE L'OBJET
         NTY = MCN( MN     )
         IF( NTY .EQ. 4 ) THEN
C           EXISTENCE D'UN VOLUME => RESOLUTION CLASSIQUE DU SYSTEME
            IF( SDOUNO .GT. 0 ) GOTO 218
CCC         ELSE IF( NTY .EQ. 5 ) THEN
CCCC           COHERENCE DE LA VARIABLE SDOUNO AVEC LES AUTRES OBJETS
CCC            NUO = MCN( MN + 1 )
CCC            CALL LXNLOU( NTOBJE, NUO, NTOB, MNOB )
CCC            IF( NTOB .GT. 0 ) THEN
CCC               CALL LXTSOU( NTOB, 'DEFINITION', NT, MN )
CCC               IF( NT .GT. 0 ) THEN
CCC                  IF( MCN( MN + WDOUNO ) .NE. SDOUNO ) GOTO 218
CCCC                 LE SOUS-OBJET N'EST PAS IMPOSE TRAITEMENT CLASSIQUE
CCCC                 CAR CELA NE SEMBLE PAS NECESSAIRE
CCC               ENDIF
CCC            ENDIF
         ENDIF
         GOTO 219

C        INCOHERENCE : TRAITEMENT CLASSIQUE IMPOSE
 218     SDOUNO = 0
         MCN( MNDFOB+WDOUNO ) = 0
         NBLGRC(NRERR) = 4
         KERR(1) = 'ATTENTION: OBJET AVEC UN VOLUME'
         KERR(2) = 'SOUS-DOMAINE NON OPERATIONNEL'
         KERR(3) = 'TRAITEMENT CLASSIQUE FORCE'
         KERR(4) = 'AVEC RENUMEROTATION DES NOEUDS'
         CALL LERESU

C        LE PLSVO SUIVANT
 219     MN  = MN + 2
 220  CONTINUE
C
C     MISE EN ORDRE V, S, L, P
      MNOBPR = MNDFOB + WTYOBJ
      CALL VSLPOB( MCN(MNDFOB+WBDOBJ), MNOBPR )

C     AFFICHAGE DU TABLEAU DEFINITION
      CALL AFTSTD( KNMTS )

      IF( INTERA .GE. 1 ) THEN
C        SI CONSOLE INTERACTIVE GRAPHIQUE ALORS
C        TRACE AUTOMATIQUE DE L'OBJET
         CALL VISEE0
C        LE TRACE DES OBJETS DE L'OBJET
         MN = MNDFOB + WTYOBJ
         DO 300 I=1,MCN(MNDFOB+WBDOBJ)
C           LE NO DU TYPE ET LE NUMERO DE L'OBJET
            NTY = MCN( MN     )
            NUO = MCN( MN + 1 )
C           LE NOM DU TYPE
            KNOMTY = NMTYOB( NTY )
C           LE NOM KNOM DE L'OBJET
            CALL NMOBNU( KNOMTY, NUO,  KNOM )
            CALL T1MOBJ( KNOMTY, KNOM, NUO )
            MN = MN + 2
 300     CONTINUE
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'FIN de DEFINITION de l''OBJET: ', KNOMOB
      ELSE
         WRITE(IMPRIM,*) 'END of DEFINITION of OBJECT: ', KNOMOB
      ENDIF
      WRITE(IMPRIM,19000)
19000 FORMAT(100('=')/)
      GOTO 100

C     FIN DE DEFINITION DES OBJETS
 9000 CALL RECTEF( NRHIST )

      RETURN
      END
