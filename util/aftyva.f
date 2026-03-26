      SUBROUTINE AFTYVA( NOIDEN, NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER NBV VARIABLES DEFINIES PAR L'IDENTIFICATEUR NOIDEN
C -----
C
C ENTREES :
C ---------
C NOIDEN  : NUMERO DE L'IDENTIFICATEUR A AFFICHER
C
C SORTIES :
C ---------
C NRETOU  : 0 SI PAS D'ERREUR RENCONTREE >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      PARAMETER  (MXIND=5)
      include"./incl/nbcamo.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
      include"./incl/langue.inc"
C.......................................................................
      CHARACTER*80      KNOM
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*10      NMTYOB,KNM
      CHARACTER*4       KAUX
      CHARACTER*66      KFORMA
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NBIND, MIMXIN(1:2,1:MXIND)
C
C     A PRIORI PAS D'ERREUR
      NRETOU = 0
C
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
      NBV    = 1
      NBIND  = 0
C
C     POSITION DANS LE TAS
      N = IDENT(2,NOIDEN)
      IF( N .GT. 0 ) THEN
C
C        C'EST UN TABLEAU
C        LE NOMBRE DE SES INDICES
         NBIND  = LETAS( N )
C        NBV LE NOMBRE DE VARIABLES DU TABLEAU
         DO 10 I=1,NBIND
            IF( LETAS(N+1) .GT. LETAS(N+2) ) THEN
C              TABLEAU VIDE
               RETURN
            ENDIF
C           LA VALEUR MINIMALE DE L'INDICE I
            MIMXIN(1,I) = LETAS(N+1)
C           LA VALEUR MAXIMALE DE L'INDICE I
            MIMXIN(2,I) = LETAS(N+2)
            NBV = NBV * ( LETAS(N+2) - LETAS(N+1) + 1 )
            N   = N + 2
 10      CONTINUE
      ENDIF
C
C     LE TYPE DE LA VARIABLE A AFFICHER
      NOTYPE = IDENT(1,NOIDEN)
C
C     LE TYPE DE LA VARIABLE A AFFICHER
      IF( NOTYPE .LE. 0 .OR. NOTYPE .GT. NBTYPV ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'AFTYVA: TYPE de VARIABLE INCORRECT'
         CALL LEREUR
         NRETOU = 1
         GOTO 9999
      ENDIF
C
C     LA POSITION DU 1-ER CARACTERE DE L'IDENTIFICATEUR
      NL = IDENT(4,NOIDEN)
      NC = IDENT(5,NOIDEN)
C     NOMBRE DE MOTS DE LA VARIABLE
      MOT = MOTVAR( NOTYPE )
C
      GOTO( 9000, 20  , 9000, 100, 300, 400, 9000, 9000, 9000,
     %      9000, 200 , 500 , 600 ), NOTYPE
C          LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPEOBJET=>13
C          TMS      => 21
C
C     CARACTERE4
C     ==========
 20   CALL CHMOT( 'caractere', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCAR( '=' )
      IF( NOMUET .NE. 0 ) THEN
         IF( NBV .LE. 1 ) THEN
            CALL AFCARX( MCN( MNTAMS(LHTMS)+LDTS(LHTMS) ) )
         ELSE
C           AFFICHAGE SUR UNE LIGNE DU MAXIMUM DE CAR4
            CALL AFLIGN
            KFORMA(1:9) = '(10(''('',I'
            NBC = NBCHIF( NBV )
            WRITE(KFORMA(10:10),'(I1)') NBC
            KFORMA(11:16) = ','':'',A'
            WRITE(KFORMA(17:17),'(I1)') NBCAMO
            KFORMA(18:25) = ','')'',))'
            M = MNTAMS(LHTMS) + LDTS(LHTMS) - 1
            WRITE(IMPRIM,KFORMA(1:25)) (I,MCN(M+I),I=1,NBV)
         ENDIF
      ENDIF
      LDTS(LHTMS) = LDTS(LHTMS) + MOT * NBV
      GOTO 9999
C
C   ENTIER {(<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })}
C   ===================================================================
C     RECHERCHE DE entier
 100  CALL CHMOT( 'entier', NL,NC, NLD,NCD, NL1,NC1 )
      NL = NL1
      NC = NC1
C     RECHERCHE ;
      CALL CHMOT( ';', NL,NC, NLD,NCD, NLPV,NCPV )
C     RECHERCHE ( APRES entier ET AVANT ;
      CALL CHMOSC( '(', NLPV,NCPV, NL,NC, NLD,NCD, NLF,NCF )
      IF( NLD .GT. 0 ) THEN
C        IL EXISTE UNE (
         NL = NLD
         NC = NCD - 1
C        RECHERCHE ( ET )
         CALL CHPAOF( NL,NC, NLD,NCD, NLF,NCF )
C
         DO 170 I=1,NBV
            IF( NOMUET .EQ. 0 ) GOTO 150
C           LA VALEUR ENTIERE A AFFICHER
            NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
C
C           VERIFICATION QUE NVALEN EST UNE VALEUR PERMISE
C           LE CARACTERE (
            NL = NLD
            NC = NCD
 110        CALL CHINTR( NL,NC,NL1,NC1,NL2,NC2,N1,N2 )
            IF( NL1 .LE. 0 ) THEN
C              VALEUR NON COMPRISE DANS L'INTERVALLE DE DEFINITION
               NBLGRC(NRERR) = 3+NL-NLD
               WRITE(KERR(2)(1:12),'(I12)') NVALEN
               KERR(1) ='AFTYVA: ENTIER ' // KERR(2)(1:12)
               KERR(2) ='DIFFERENT DES VALEURS LICITES'
               DO 111 J=NLD,NL
                  KERR(3+J-NLD) = KTD(J)
 111           CONTINUE
               CALL LEREUR
               NRETOU = 1
               GOTO 9999
            ENDIF
C
            IF( NVALEN .LT. N1 .OR. NVALEN .GT. N2 ) THEN
C              VALEUR NON DANS L'INTERVALLE. ON PASSE A LA SUIVANTE
               NL = NL2
               NC = NC2
C              RECHERCHE DU :
               CALL CHCAR( ':', NL, NC )
C              RECHERCHE DE LA CHAINE '...'
               CALL CHCHAI( NL,NC,NL1,NC1,NL2,NC2 )
               NL = NL2
               NC = NC2
C              RECHERCHE DE , OU )
               CALL CAR1NB( NL, NC )
               IF( KTD(NL)(NC:NC) .NE. ',' .AND.
     %             KTD(NL)(NC:NC) .NE. ')'      ) THEN
C                 ( INTERVAL_ENT : CHAINE {,INTERVAL_ENT : CHAINE} )
C                 INCORRECT
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'AFTYVA: INCORRECT ( ... , ... )'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  NRETOU = 2
                  GOTO 9999
               ENDIF
               GOTO 110
            ENDIF
C
C           LA VALEUR EST DANS L'INTERVALLE ACTUEL
            CALL AFCAR( '=' )
            CALL AFENTI( NVALEN )
C
C           LA CHAINE EXPLICATIVE DE CETTE VALEUR
            NLL = NL2
            NCC = NC2
C           NLL NCC CARACTERE DE FIN DE <INTERVAL_ENT>
C           RECHERCHE DE :
            CALL CHCAR( ':', NLL, NCC )
            IF( NLL .GT. 0 ) THEN
C              RECHERCHE DE LA CHAINE
               CALL CHCHAI( NLL,NCC, NL1,NC1, NL2,NC2 )
               IF( NL1 .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'AFTYVA: CHAINE INEXISTANTE'
                  ELSE
                     KERR(1) = 'AFTYVA: UNKNOWN STRING of CHARACTERS'
                  ENDIF
                  KERR(2) = KTD(NLL)
                  CALL LEREUR
                  NRETOU = 1
                  GOTO 9999
               ENDIF
C              AFFICHAGE DE LA CHAINE
               CALL AFCHAI( NL1,NC1, NL2,NC2 )
            ENDIF
            CALL AFLIGN
C
 150        LDTS(LHTMS) = LDTS(LHTMS) + MOT
 170     CONTINUE
C
      ELSE
C        IL N'EXISTE PAS DE ( )
         CALL AFCAR( '=' )
         CALL AFTAEN(NBIND, MIMXIN, MOT, MCN(MNTAMS(LHTMS)+LDTS(LHTMS)))
      ENDIF
      GOTO 9999
C
C     ^NOM_LEXIQUE
C     ============
 200  CALL CHMOT( '^', NL,NC, NLD,NCD, NLF,NCF )
      NL = NLF
      NC = NCF
      CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
C
C     LE NUMERO DU TAMS DU LEXIQUE PERE DU TS
C     FUSION DU NOM DU TD ET TS
      CALL FUSNOM( KTD(NLD)(NCD:NCF), NOMTS(LHTMS), KNOM, NCODER )
      IF( NCODER .EQ. 0 ) THEN
C        NOMS NON FUSIONNABLES
         NRETOU = 1
         GOTO 9999
      ENDIF
      LKNOM = INDEX( KNOM, '"' )
      IF( LKNOM .LE. 0 ) THEN
         LKNOM = NUDCNB( KNOM )
      ELSE
         LKNOM = LKNOM - 1
      ENDIF
      IF( NCODER .EQ. 1 ) THEN
C        OUVERTURE DU LEXIQUE DU TD
         CALL TNOUVR( KNOM(1:LKNOM), KAUX, NTLX, MNLX )
      ELSE
C        RECHERCHE ET OUVERTURE DU LEXIQUE FILS
         CALL LXNTPN( KNOM(1:LKNOM), NTLX, MNLX, N1, N2 )
      ENDIF
      IF( NTLX .LE. 0 ) THEN
C        LE LEXIQUE N'EXISTE PAS
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='AFTYVA: LEXIQUE INCONNU'
         ELSE
            KERR(1) ='AFTYVA: UNKNOWN LEXICON'
         ENDIF
         KERR(2) = KNOM
         CALL LEREUR
         NRETOU = 1
         GOTO 9999
      ENDIF
C
C     AFFICHAGE DES NOMS DANS LE LEXIQUE
      CALL AFCAR( '=' )
      IF( NBV .GT. 0 ) CALL AFLIGN
      DO 240 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 235
C        ( NO VARIABLE :
         CALL AFCAR( '(' )
         CALL AFENTI( I )
         CALL AFCAR( ':' )
C        OUVERTURE DANS LE LEXIQUE
         NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
C        L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
         MN = MNLX + MCN( MNLX ) * NVALEN
C        LE NOM CONVERTI D'ENTIERS EN CARACTERES
C        KNOM OFFRE UNE PROTECTION DE DEBORDEMENT
         CALL ENTNOM( MCN(MNLX+2), MCN(MN), KNOM )
         J = INDEX( KNOM, ' ' ) + 1
         KNOM(J:J) = ')'
         CALL AFCAR( KNOM(1:J) )
         CALL AFLIGN
 235     LDTS(LHTMS) = LDTS(LHTMS) + MOT
 240  CONTINUE
      GOTO 9999
C
C     REEL
C     ====
 300  CALL CHMOT( 'reel', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCAR( '=' )
      CALL AFTAR1( NBIND, MIMXIN, MOT, RMCN(MNTAMS(LHTMS)+LDTS(LHTMS)) )
      GOTO 9999
C
C     REEL2
C     =====
 400  CALL CHMOT( 'reel2', NL,NC, NLD,NCD, NLF,NCF )
      CALL AFCAR( '=' )
      CALL AFTAR2( NBIND, MIMXIN, MOT, RMCN(MNTAMS(LHTMS)+LDTS(LHTMS)) )
      GOTO 9999
C
C     XYZ
C     ===
 500  CALL CHMOT( 'xyz', NL,NC, NLD,NCD, NLF,NCF )
C     LA LIGNE ACTUELLE EST AFFICHEE
      CALL AFLIGN
C     AFFICHAGE SUR UNE LIGNE DES 3 COORDONNEES DE CHAQUE XYZ
      KFORMA(1:11) = '( '' XYZ '',I'
      NBC = NBCHIF( NBV )
      WRITE(KFORMA(12:12),'(I1)') NBC
      KFORMA(13:26) = ','' :'',3G14.7) '
      DO 530 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 525
C        AFFICHAGE SUR UNE LIGNE DES 3 COORDONNEES
         NVALEN  = MNTAMS(LHTMS) + LDTS(LHTMS)
         WRITE(IMPRIM,KFORMA(1:25)) I,(RMCN(NVALEN+L),L=0,2)
 525     LDTS(LHTMS) = LDTS(LHTMS) + MOT
 530  CONTINUE
C     UNE LIGNE BLANCHE
      WRITE(IMPRIM,*)
      GOTO 9999
C
C     TYPEOBJET
C     =========
 600  CALL CHMOT( 'typeobjet', NL,NC, NLD,NCD, NLF,NCF )
C     AFFICHAGE DES NOMS DANS LE LEXIQUE
      CALL AFCAR( '=' )
      IF( NBV .GT. 1 ) CALL AFLIGN
      DO 640 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 635
C        ( NO VARIABLE :
         CALL AFCAR( '(' )
         CALL AFENTI( I )
         CALL AFCAR( ':' )
C        OUVERTURE DANS LE LEXIQUE DE L'OBJET
         N2  = MNTAMS(LHTMS) + LDTS(LHTMS)
         N1  = MCN( N2 )
         KNM = NMTYOB( N1 )
         J = INDEX( KNM, 'LIGNE' )
         IF( J .GT. 0 ) KNM(J:J+4)= 'LINE '
         J = INDEX( KNM, 'OBJET' )
         IF( J .GT. 0 ) KNM(J:J+5)= 'OBJECT'
         J = INDEX( KNM, ' ' )
         CALL AFCAR( KNM(1:J) )
         CALL LXNLOU( NTADAM, N1, NTLX, MNLX )
         NVALEN = MCN( N2 + 1 )
C        L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
         MN = MNLX + MCN( MNLX ) * NVALEN
C        LE NOM CONVERTI D'ENTIERS EN CARACTERES
C        KNOM OFFRE UNE PROTECTION DE DEBORDEMENT
         CALL ENTNOM( MCN(MNLX+2), MCN(MN), KNOM )
         J = INDEX( KNOM, ' ' ) + 1
         KNOM(J:J) = ')'
         CALL AFCAR( KNOM(1:J) )
         CALL AFLIGN
 635     LDTS(LHTMS) = LDTS(LHTMS) + MOT
 640  CONTINUE
      GOTO 9999
C
C     ERREUR TYPE INCORRECT ICI
 9000 NBLGRC(NRERR) = 1
      WRITE(KERR(2)(1:4),'(I4)') NOTYPE
      KERR(1) = 'AFTYVA:TYPE' // KERR(2)(1:4) // ' INCORRECT'
      CALL LEREUR
      NRETOU = 1
C
 9999 RETURN
      END
