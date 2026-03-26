      SUBROUTINE MOTYVM( NOIDEN, NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MODIFIER OU CREER PAR MENU LA VALEUR DE
C -----     L'IDENTIFICATEUR NOIDEN DU TABLEAU IDENT
C
C ENTREES :
C ---------
C NOIDEN  : NUMERO DE L'IDENTIFICATEUR
C
C SORTIES :
C ---------
C NRETOU  : 0 SI LA VALEUR EST CREEE
C           1 SI LA VALEUR EST ABANDONNEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        mars 1990
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/typobj.inc"
      include"./incl/nbcamo.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/langue.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*10      NMTYOB, KNOMOB
      CHARACTER*80      KNOM
      CHARACTER*84      KLGFRA
      CHARACTER*4       KAUX
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C     CE QUI SUIT POUR RECUPERER LES 2 EVENTUELS MOTS D'UN REEL2
      REAL              R(2)
      DOUBLE PRECISION  D
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
      EQUIVALENCE     (D,R(1))
C
C     INITIALISATIONS
      KNOM   = ' '
      NRETOU = 0
      INITIK = 0
C
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
C     ====================================================
      NBV    = 1
C     POSITION DANS LE TAS
      N = IDENT(2,NOIDEN)
      IF( N .GT. 0 ) THEN
C        C'EST UN TABLEAU
C        LE NOMBRE DE SES INDICES
         M = LETAS( N )
         DO 10 I=1,M
            NBV = NBV * ( LETAS(N+2) - LETAS(N+1) + 1 )
            N   = N + 2
 10      CONTINUE
      ENDIF
C
C     LE TABLEAU MC EST IL ASSEZ GRAND ?
C
C     LE NOMBRE DE MOTS DE LA VARIABLE ET DU TABLEAU DECLARE
      L2MOTS = NBVATC(LHTMS)
C     LE TYPE DE L'IDENTIFICATEUR
      NOTYPE = IDENT(1,NOIDEN)
C     LE NOMBRE DE MOTS D'UNE VARIABLE DE L'IDENTIFICATEUR
      MOTS   = MOTVAR( NOTYPE )
C
C     LE NOMBRE DE MOTS ACTUELS DU TABLEAU
      L1MOTS = LDTS(LHTMS) + NBV * MOTS
C
      IF( L1MOTS .GT. L2MOTS ) THEN
C        LE TABLEAU MC EST TROP PETIT.IL EST AUGMENTE
         L = L1MOTS + 32
         CALL TNMCAU( 'MOTS', NBVATC(LHTMS),  L ,
     %                NBVATC(LHTMS), MCTAMS(LHTMS) )
         NBVATC(LHTMS) = L
      ENDIF
C
C     DEMANDE SELON LE TYPE
C     =====================
      IF( NOTYPE .LE. 0 .OR. NOTYPE .GT. NBTYPV ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(2)(1:3),'(I3)') NOTYPE
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MOTYVM:'//KERR(2)(1:3)//' TYPE INCORRECT'
         ELSE
            KERR(1) = 'MOTYVM:'//KERR(2)(1:3)//' INCORRECT TYPE'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     LA POSITION DU 1-ER CARACTERE DE L'IDENTIFICATEUR
      NL = IDENT(4,NOIDEN)
      NC = IDENT(5,NOIDEN)
C
C     LA LIGNE DE FRAPPE CONTIENT EN COMMENTAIRE LE NOM DE L'IDENTIFICATEUR
C     SUIVI DE L'EVENTUEL TEXTE COMMENTAIRE
      IF( LHLECT .EQ. 1 ) THEN
         L      = NUDCNB( KTD(NL)(NC:NCKTD) )
         KLGFRA = '{ ' //  KTD(NL)(NC:NC+L-1) // ' }'
         INITIK = 0
      ENDIF
C
      GOTO( 9000, 22  , 9000, 100, 300, 400, 9000, 9000, 9000,
     %      9000, 200 , 500 , 600 ), NOTYPE
C          LOGIQUE   => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL      => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2 => 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPE_OBJET=>13
C
C     CARACTERE4
C     ----------
 22   CALL CHMOT( 'caractere', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCHAI( NLD,NCD, NLF,NCF )
      CALL AFLIGN
      DO 40 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 30
C        LECTURE DES 4 CARACTERES
         NCVALS = 0
         CALL LIRCAR( NCVALS, KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
         MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = ICHARX( KNOM(1:NBCAMO) )
 30      LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 40   CONTINUE
      GOTO 8000
C
C   ENTIER {(<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })}
C   -------------------------------------------------------------------
C     RECHERCHE DE entier
 100  CALL CHMOT( 'entier', NL,NC, NLD,NCD, NL1,NC1 )
      NL = NL1
      NC = NC1
C     RECHERCHE ;
      CALL CHMOT( ';', NL,NC, NLD,NCD, NLPV,NCPV )
C     RECHERCHE ( AVANT ;
      CALL CHMOSC( '(', NLPV,NCPV, NL,NC, NLD,NCD, NLF,NCF )
C
      IF( NLD .GT. 0 ) THEN
C
C        CAS: (<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })
         NL = NLD
         NC = NCD - 1
C        RECHERCHE DE ( ET )
         CALL CHPAOF( NL,NC, NLD,NCD, NLF,NCF )
C        LE CARACTERE (
         NL = NLD
         NC = NCD
C        LA LIGNE DE MENU AVANT LA 1-ERE OPTION
         NLMEN0 = 0
C
C        PARCOURS DES DIFFERENTES OPTIONS DEBUTANT
C        TOUTES PAR INTERVAL_ENT :
C        PASSAGE FORCE A LA LIGNE
 110     NBLGRC(NRMENU) = NBLGRC(NRMENU) + 1
         MDLGRC(NRMENU) = 0
C
         CALL CHINTR( NL,NC,NL1,NC1,NL2,NC2,N1,N2 )
C        NL1,NC1 : POSITION DU PREMIER CARACTERE DE L'INTERVALLE
C                  NL1 = 0 SI INTERVALLE_E INCORRECT
C        NL2,NC2 : POSITION DANS KTD DU DERNIER CARACTERE DE INTERVAL_ENT
C        N1 ,N2  : VALEURS INITIALE ET FINALE DE L'INTERVALLE
C                  N1=N2 SI INTERVALLE_E = ID_ENTIER
         IF( NL1 .LE. 0 ) THEN
C           INTERVALLE INCORRECT
            NBLGRC(NRERR) = 3
            KERR(2) = KTD(NL)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='MOTYVM: INTERVALLE INCORRECT '
               KERR(3) ='TD INCORRECT . A REDEFINIR'
            ELSE
               KERR(1) ='MOTYVM: INCORRECT INTERVAL'
               KERR(3) ='INCORRECT TD. DEFINE AGAIN'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
         NL = NL2
         NC = NC2
C        RECHERCHE DU :
         CALL CHCAR( ':', NL, NC )
C        RECHERCHE DE LA CHAINE '...'
         CALL CHCHAI( NL,NC,NL1,NC1,NL2,NC2 )
         NL = NL2
         NC = NC2
C        RECHERCHE DE , OU )
         CALL CAR1NB( NL, NC )
         IF( KTD(NL)(NC:NC) .NE. ',' .AND.
     %       KTD(NL)(NC:NC) .NE. ')' ) THEN
C          ( INTERVAL_ENT : CHAINE {,INTERVAL_ENT : CHAINE} ) INCORRECT
             NBLGRC(NRERR) = 3
             KERR(1) ='MOTYVM:(... , ... ) INCORRECT '
             KERR(2) = KTD(NL)
             IF( LANGAG .EQ. 0 ) THEN
                KERR(3) ='TD INCORRECT . A REDEFINIR'
             ELSE
                KERR(3) ='INCORRECT TD. DEFINE AGAIN'
             ENDIF
             CALL LEREUR
             GOTO 9900
         ENDIF
C
C        MISE DANS LE MENU DE CETTE OPTION
         CALL AFCHAI( NL1,NC1, NL,NC-1 )
         CALL AFLIGN
C        LE NUMERO ASSOCIE A CETTE OPTION
         IF( N1 .EQ. N2 ) THEN
            DO 140 I=NLMEN0+1,NBLGRC(NRMENU)
               NAMENU( I ) = N1
 140        CONTINUE
         ELSE
C           UN INTERVALLE DE PLUSIEURS VALEURS
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'MOTYVM: INTERVALLE AU LIEU D''UN ENTIER'
               KERR(2) = 'TD INCORRECT . A REDEFINIR'
            ELSE
               KERR(1) = 'MOTYVM: INTERVAL in PLACE of an INTEGER'
               KERR(2) = 'INCORRECT TD. DEFINE AGAIN'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C        LA DERNIERE LIGNE DE CETTE OPTION
         NLMEN0 = NBLGRC(NRMENU)
         IF( KTD(NL)(NC:NC) .NE. ')' ) GOTO 110
C
C        TOUTES LES OPTIONS DU MENU ONT ETE VUES
C        LECTURE DES VALEURS
         LIMENU = 1
      ELSE
C        LECTURE SIMPLE D'UN ENTIER
         LIMENU = 0
      ENDIF
C
      DO 190 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 180
C        LECTURE DE L'ENTIER
 172     NCVALS = 0
         CALL LIRENT( NCVALS, NVALEN )
         IF( NCVALS .EQ. -1 ) GOTO 9900
C
         IF( LIMENU .NE. 0 ) THEN
C           VERIFICATION DE LA VALIDITE DE L'ENTIER DANS LE MENU
            DO 175 KK=1,NBLGRC(NRMENU)
               J = KK
               IF( NVALEN .EQ. NAMENU(J) ) THEN
C                 VALEUR LUE CORRECTE
C                 SAUVEGARDE DE LA SIGNIFICATION DE L'ENTIER FRAPPE
C                 SUR LE FICHIER
                  IF( LHLECT .EQ. 1 ) THEN
C                    RECHERCHE DE L'EVENTUELLE SECONDE LIGNE DU MENU
C                    CORRESPONDANT A CETTE VALEUR
                     IF( NVALEN .EQ. NAMENU(J+1) ) J=J+1
                     L      = NUDCNB( KMENU(J) )
                     KLGFRA = '{ ' // KMENU(J)(1:L) // ' }'
                     CALL SANSDBL( KLGFRA, L )
                     WRITE(NFFRAP,*) KLGFRA(1:L)
                     INITIK = 1
                  ENDIF
                  GOTO 179
               ENDIF
 175        CONTINUE
C           LA VALEUR EST INCORRECTE
            NBLGRC(NRERR) = 2
            WRITE( KERR(MXLGER)(1:12), '(I12)' ) NVALEN
            IF( INTERA .EQ. 0 ) THEN
               KERR(1) = 'ENTIER LU INCORRECT' //  KERR(MXLGER)(1:12)
               KERR(2) = 'REDONNER CET ENTIER'
            ELSE
               KERR(1) = 'INCORRECT READ INTEGER' //  KERR(MXLGER)(1:12)
               KERR(2) = 'GIVE AGAIN THIS INTEGER'
            ENDIF
            CALL LEREUR
C           REPRISE DE LA LECTURE DE L'ENTIER
            GOTO 172
         ENDIF
C
C        LA VALEUR EST COPIEE DANS LE TABLEAU MS COPIE EN MC
 179     MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
 180     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 190  CONTINUE
C     LE RESULTAT EST MIS SOUS FORME DE CARACTERES
      KNOM = ' '
      WRITE( KNOM(1:12), '(I12)' ) NVALEN
      GOTO 8000
C
C     ^NOM_LEXIQUE
C     ------------
 200  CALL CHMOT( '^', NL,NC, NLD,NCD, NLF,NCF )
      NL = NLF
      NC = NCF
      CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
C
C     LE NUMERO DU TAMS DU LEXIQUE PERE DU TS
C     FUSION DU NOM DU TD ET TS
      CALL FUSNOM( KTD(NLD)(NCD:NCF), NOMTS(LHTMS),  KNOM, NCODER )
      IF( NCODER .EQ. 0 ) THEN
C        NOMS NON FUSIONNABLES
         GOTO 9900
      ENDIF
C
      IF( NCODER .EQ. 1 ) THEN
C        OUVERTURE DU LEXIQUE DU TD
         CALL TNOUVR( KNOM, KAUX, NTLX, MNLX )
      ELSE
C        RECHERCHE ET OUVERTURE DU LEXIQUE FILS
         CALL LXNTPN( KNOM, NTLX, MNLX, N1, N2 )
      ENDIF
      IF( NTLX .LE. 0 ) THEN
C        LE LEXIQUE N'EXISTE PAS
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='MOTYVM: LEXIQUE INCONNU'
         ELSE
            KERR(1) ='MOTYVM: UNKNOWN LEXICON'
         ENDIF
         KERR(2) = KNOM
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     LECTURE DES NOMS DANS LE LEXIQUE  AVEC ^ LE CARACTERE PRECEDANT
      CALL AFCHAI( NLD,NCD-1, NLF,NCF )
      CALL AFLIGN
      DO 240 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 235
C        LECTURE DU NOM DANS LE LEXIQUE NTLX
         NCVALS = 0
         CALL LIRLEX( NTLX, NCVALS, KNOM, NVALEN )
         IF( NCVALS .LE. 0 ) GOTO 9900
C
C        Y A T IL RECURSIVITE ?
         IF( NETOBR .NE. 0 ) THEN
            IF( NTLOBR .EQ. NTLX .AND. NUMOBR .EQ. NVALEN ) THEN
C               OUI : LE LEXIQUE S'APPELLE LUI MEME
                NBLGRC(NRERR) = 2
                IF( LANGAG .EQ. 0 ) THEN
                   KERR(1) = 'RECURSIVITE INTERDITE'
                   KERR(2) = 'LEXIQUE SE REFERANT A LUI MEME'
                ELSE
                   KERR(1) = 'FORBIDDEN RECURSIVITY'
                   KERR(2) = 'LEXICON REFERS ITSELF'
                ENDIF
                CALL LEREUR
C               SIMULATION D'ABANDON UTILISATEUR POUR CONSERVER
C               LE TMS INITIAL EVENTUEL
                GOTO 9900
            ENDIF
         ENDIF
C
C        SAUVEGARDE DU NOM DU LEXIQUE SUR LE FICHIER FRAPPE
         IF( LHLECT .EQ. 1 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
             KLGFRA = '{ dans le LEXIQUE ' // KTD(NLD)(NCD:NCF) // ' }'
            ELSE
             KLGFRA = '{ in LEXICON ' // KTD(NLD)(NCD:NCF) // ' }'
            ENDIF
            CALL SANSDBL( KLGFRA, L )
            WRITE(NFFRAP,*) KLGFRA(1:L)
            INITIK = 1
         ENDIF
C
C        LA VALEUR EST COPIEE DANS LE TABLEAU MS COPIE EN MC
         MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
 235     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 240  CONTINUE
      GOTO 8000
C
C     REEL
C     ----
 300  CALL CHMOT( 'reel', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCHAI( NLD,NCD, NLF,NCF )
      CALL AFLIGN
      DO 330 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 320
C        LECTURE DU REEL SIMPLE
         NCVALS = 0
         CALL LIRRSP( NCVALS, RMCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) )
         IF( NCVALS .EQ. -1 ) GOTO 9900
 320     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 330  CONTINUE
C     LE RESULTAT EST MIS SOUS FORME DE CARACTERES
      KNOM = ' '
      WRITE( KNOM(1:25),'(G25.17)') RMCN(MCTAMS(LHTMS)+LDTS(LHTMS)-MOTS)
      GOTO 8000
C
C     REEL2
C     -----
 400  CALL CHMOT( 'reel2', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCHAI( NLD,NCD, NLF,NCF )
      CALL AFLIGN
      DO 430 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 420
C        LECTURE DU REEL DOUBLE
         NCVALS = 0
         CALL LIRRDP( NCVALS, RMCN( MCTAMS(LHTMS) + LDTS(LHTMS) )  )
         IF( NCVALS .EQ. -1 ) GOTO 9900
 420     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 430  CONTINUE
C     LE RESULTAT EST MIS SOUS FORME DE CARACTERES
      R(1) = RMCN( MCTAMS(LHTMS)+LDTS(LHTMS)-MOTS )
      IF( MOTS .GT. 1 ) THEN
         R(2) = RMCN( MCTAMS(LHTMS)+LDTS(LHTMS)-MOTS + 1 )
      ENDIF
      KNOM = ' '
      WRITE( KNOM(1:25) ,'(G25.17)') D
      GOTO 8000
C
C     XYZ
C     ---
 500  CALL CHMOT( 'xyz', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCHAI( NLD,NCD, NLF,NCF )
      CALL AFLIGN
      DO 530 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 525
C        LECTURE DE XYZ
         NCVALS = 0
         CALL LIRXYZ( NCVALS, RMCN( MCTAMS(LHTMS)+LDTS(LHTMS) ) )
         IF( NCVALS .EQ. -1 ) GOTO 9900
 525     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 530  CONTINUE
C     LE RESULTAT EST MIS SOUS FORME DE CARACTERES
      KNOM = ' '
      WRITE( KNOM(1:25),'(G25.17)') RMCN(MCTAMS(LHTMS)+LDTS(LHTMS)-MOTS)
      GOTO 8000
C
C     TYPEOBJET
C     ---------
 600  CALL CHMOT( 'typeobjet', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCHAI( NLD,NCD, NLF,NCF )
      CALL AFLIGN
      DO 670 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 650
C        LE NUMERO DU TYPE DE L'OBJET
C        LECTURE DU TYPE D'OBJET
 610     NCVALS = 0
         CALL INVITE( 91 )
         CALL LIRCAR( NCVALS, KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
         N1 = NTYOBJ( KNOM )
         IF( N1 .EQ. 0 ) THEN
C           ERREUR A REDONNER
            GOTO 610
         ENDIF
C        OUVERTURE DU LEXIQUE DANS LE LEXIQUE ADAM
         KNOMOB = NMTYOB( N1 )
         CALL LXLXOU( NTADAM, KNOMOB, NTLX, MNLX )
C        LECTURE DU NOM DE L'OBJET DANS SON LEXIQUE
         GOTO(611, 612, 613, 614, 615 ),N1
 611     CALL INVITE( 51 )
         GOTO 630
 612     CALL INVITE( 40 )
         GOTO 630
 613     CALL INVITE( 42 )
         GOTO 630
 614     CALL INVITE( 60 )
         GOTO 630
 615     CALL INVITE( 45 )
C
C        LECTURE DU NOM
 630     NCVALS = 0
         CALL LIRCAR( NCVALS, KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
C        OUVERTURE DANS LE LEXIQUE
         CALL LXNMNO( NTLX, KNOM, NVALEN, N2 )
         IF( NVALEN .LE. 0 ) THEN
C           NOM INCORRECT NON RETROUVE DANS LE LEXIQUE
            IF (LHLECT.GT. 1 .AND. LHLECT.LE.MXLECT ) THEN
C              LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C              AFFICHAGE DU RETOUR AU CLAVIER
               CALL QUIFIC
            ELSE
               NBLGRC(NRERR) = 0
            ENDIF
            N2 = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KERR(NBLGRC(NRERR)+1) = 'NOM INCONNU: ' // KNOM(1:N2)
               KERR(NBLGRC(NRERR)+2) = 'A choisir PARMI:'
            ELSE
               KERR(NBLGRC(NRERR)+1) = 'UNKNOWN NAME: ' // KNOM(1:N2)
               KERR(NBLGRC(NRERR)+2) = 'To be chosen AMONG:'
            ENDIF
            NBLGRC(NRERR) = NBLGRC(NRERR) + 2
            CALL LXIM0( MNLX )
            GOTO 630
         ENDIF
C        LES VALEURS SONT COPIEES DANS LE TABLEAU MS COPIE EN MC
         N2 = MCTAMS(LHTMS) + LDTS(LHTMS)
         MCN( N2     ) = N1
         MCN( N2 + 1 ) = NVALEN
 650     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 670  CONTINUE
C
C     SORTIE AVEC EFFACEMENT DU MENU
 8000 CALL RECTEF( NRMENU )
      NUIDEN = 0
      GOTO 9995
C
C     ERREUR TYPE INCORRECT ICI
 9000 NBLGRC(NRERR) = 1
      WRITE(KERR(2)(1:3),'(I3)') NOTYPE
      KERR(1) =  'MOTYVM: TYPE'//KERR(2)(1:3)//' INCORRECT'
      CALL LEREUR
C
C     ABANDON DE LA LECTURE OU ERREUR
 9900 NRETOU = 1
C     LE MENU EST EFFACE
      CALL RECTEF( NRMENU )
C     L'INVITE EST EFFACEE
      CALL RECTEF( NRINVI )
C
C     SAUVEGARDE SUR LE FICHIER D'UN COMMENTAIRE ECLAIRANT LA DONNEE
 9995 NBLGRC(NRMENU) = 0
      IF( LHLECT .EQ. 1 ) THEN
         IF( NRETOU.EQ.0 .AND. INITIK.EQ.0 .AND. NOMUET.NE.0 ) THEN
            CALL SANSDBL( KLGFRA, L )
            WRITE(NFFRAP,*) KLGFRA(1:L)
         ENDIF
      ENDIF

      RETURN
      END
