      SUBROUTINE NCTYVA( NFFICA , NOIDEN , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : METTRE SUR LE FICHIER NFFICA LES NBV  VARIABLES
C ----- DEFINIES PAR L'IDENTIFICATEUR NOIDEN
C
C ENTREES :
C ---------
C NFFICA  : NUMERO D'UNITE LOGIQUE DU FICHIER CARACTERE
C NOIDEN  : NUMERO DE L'IDENTIFICATEUR A AFFICHER
C
C SORTIES :
C ---------
C NRETOU  : 0 SI PAS D'ERREUR RENCONTREE >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*80      KNOM
      CHARACTER*4       CHARX,CAR4
      CHARACTER*10      NMTYOB
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
C
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
      NBV    = 1
      NRETOU = 0
C
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
C     LE TYPE DE LA VARIABLE A AFFICHER
      NOTYPE = IDENT(1,NOIDEN)
      IF( NOTYPE .LE. 0 .OR. NOTYPE .GT. NBTYPV ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
         KERR(1) = 'NCTYVA:TYPE INCORRECT '//KERR(MXLGER)(1:4)
         CALL LEREUR
         NRETOU = 1
         RETURN
      ENDIF
C
C     LA POSITION DU 1-ER CARACTERE DE L'IDENTIFICATEUR
      NL = IDENT(4,NOIDEN)
      NC = IDENT(5,NOIDEN)
C     NOMBRE DE MOTS DE LA VARIABLE
      MOT = MOTVAR( NOTYPE )
C
      GOTO( 9000 ,  50 , 9000 , 100 , 300 , 400 , 9000 , 9000 , 9000,
     %      9000 , 200 , 500  , 600 ) , NOTYPE
C          LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPEOBJET=>13 TMS      => 21
C
C     4 CARACTERES
C     ============
 50   DO 55 I=1,NBV
         CAR4 = CHARX( MCN( MNTAMS(LHTMS) + LDTS(LHTMS) ) )
         CALL FICAHC( NFFICA , CAR4 )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 55   CONTINUE
      RETURN
C
C   ENTIER {(<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })}
C   ===================================================================
C     RECHERCHE DE entier
 100  CALL CHMOT( 'entier' , NL,NC, NLD,NCD, NL1,NC1 )
      NL = NL1
      NC = NC1
C     RECHERCHE ;
      CALL CHMOT( ';' , NL,NC, NLD,NCD, NLPV,NCPV )
C     RECHERCHE ( AVANT ;
      CALL CHMOSC( '(' , NLPV,NCPV, NL,NC, NLD,NCD, NLF,NCF )
      IF( NLD .GT. 0 ) THEN
C        IL EXISTE UNE (
         NL = NLD
         NC = NCD - 1
         CALL CHPAOF( NL,NC, NLD,NCD, NLF,NCF )
C
         DO 170 I=1,NBV
C           LA VALEUR A TRANSFERER
            NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
C
C           VERIFICATION QUE NVALEN EST UNE VALEUR PERMISE
C           LE CARACTERE (
            NL = NLD
            NC = NCD
C
 110        CALL CHINTR( NL,NC,NL1,NC1,NL2,NC2,N1,N2 )
            IF( NL1 .LE. 0 ) THEN
C              VALEUR NON COMPRISE DANS L'INTERVALLE DE DEFINITION
               NBLGRC(NRERR) = 3+NL-NLD
               WRITE(KERR(2)(1:12),'(I12)') NVALEN
               KERR(1) ='NCTYVA: ENTIER ' // KERR(2)(1:12)
               KERR(2) ='DIFFERENT DES VALEURS LICITES'
               DO 111 J=NLD,NL
                  KERR(3+J-NLD) = KTD(J)
 111           CONTINUE
               CALL LEREUR
               NRETOU = 1
               RETURN
            ENDIF
C
            IF( NVALEN .LT. N1 .OR. NVALEN .GT. N2 ) THEN
C              VALEUR NON DANS L'INTERVALLE. ON PASSE A LA SUIVANTE
               NL = NL2
               NC = NC2
C              RECHERCHE DU :
               CALL CHCAR( ':' , NL , NC )
C              RECHERCHE DE LA CHAINE '...'
               CALL CHCHAI( NL,NC,NL1,NC1,NL2,NC2 )
               NL = NL2
               NC = NC2
C              RECHERCHE DE , OU )
               CALL CAR1NB( NL , NC )
               IF( KTD(NL)(NC:NC) .NE. ',' ) THEN
C                 ( INTERVAL_ENT : CHAINE {,INTERVAL_ENT : CHAINE} )
C                 INCORRECT
                  NBLGRC(NRERR) = 21
                  KERR(1) = 'NCTYVA:(... , ... ) INCORRECT '
                  KERR(2) =  KTD(NL)
                  CALL LEREUR
                  NRETOU = 2
                  RETURN
               ENDIF
               GOTO 110
            ENDIF
C
C           LA VALEUR EST DANS L'INTERVALLE ACTUEL
C
C           LA CHAINE EXPLICATIVE DE CETTE VALEUR
            NLL = NL2
            NCC = NC2
C           NLL NCC CARACTERE DE FIN DE <INTERVAL_ENT>
C           RECHERCHE DE :
            CALL CHCAR( ':' , NLL , NCC )
C           RECHERCHE DE LA CHAINE
            CALL CHCHAI( NLL,NCC, NL1,NC1, NL2,NC2 )
            IF( NL1 .LE. 0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'NCTYVA: CHAINE INEXISTANTE'
               KERR(2) = KTD(NLL)
               CALL LEREUR
               NRETOU = 1
               RETURN
            ENDIF
C
C           ENTIER PORTE DANS LE BUFFER KFICA
            CALL FICAEC( NFFICA , NVALEN )
            LDTS(LHTMS) = LDTS(LHTMS) + MOT
 170     CONTINUE
C
      ELSE
C        IL N'EXISTE PAS DE ( )
         DO 190 I=1,NBV
            NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
            CALL FICAEC( NFFICA , NVALEN )
            LDTS(LHTMS) = LDTS(LHTMS) + MOT
 190     CONTINUE
      ENDIF
      RETURN
C
C     ^NOM_LEXIQUE
C     ============
 200  CALL CHMOT( '^' , NL,NC, NLD,NCD, NLF,NCF )
      NL = NLF
      NC = NCF
      CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
C
C     LE NUMERO DU TAMS DU LEXIQUE PERE DU TS
C     FUSION DU NOM DU TD ET TS
      CALL FUSNOM( KTD(NLD)(NCD:NCF) , NOMTS(LHTMS) , KNOM , NCODER )
      IF( NCODER .EQ. 0 ) THEN
C        NOMS NON FUSIONNABLES
         NRETOU = 1
         RETURN
      ENDIF
C
      IF( NCODER .EQ. 1 ) THEN
C        OUVERTURE DU LEXIQUE DU TD
         CALL TNOUVR( KNOM , CAR4 , NTLX , MNLX )
      ELSE
C        RECHERCHE ET OUVERTURE DU LEXIQUE FILS
         CALL LXNTPN( KNOM , NTLX , MNLX , N1 , N2 )
      ENDIF
      IF( MNLX .LE. 0 ) THEN
C        LE LEXIQUE N'EXISTE PAS
         NBLGRC(NRERR) = 2
         KERR(1) = 'NCTYVA: LEXIQUE INCONNU'
         KERR(2) = KNOM
         CALL LEREUR
         NRETOU = 1
         RETURN
      ENDIF
C
C     CONVERSION DES NOMS DANS LE LEXIQUE
      DO 240 I=1,NBV
C        OUVERTURE DANS LE LEXIQUE
         NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
         IF( NVALEN .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) =  KNOM
            KERR(2) ='NON DANS LE LEXIQUE PERE'
            CALL LEREUR
            CALL LXIM0( MNLX )
            NRETOU = 1
         ENDIF
C
C        L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
         MN = MNLX + MCN( MNLX ) * NVALEN
C        LE NOM CONVERTI D'ENTIERS EN CARACTERES
C        KNOM OFFRE UNE PROTECTION DE DEBORDEMENT
         CALL ENTNOM( MCN(MNLX+2) , MCN(MN) , KNOM )
C
C        ATTENTION: LE NOM EST ECRIT , PAS LE NUMERO
C                   EN EFFET LE NUMERO DANS LE LEXIQUE EST LOCAL
C                   LE NOM LUI RESTE GLOBAL
C
         CALL FICAKC( NFFICA , KNOM )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 240  CONTINUE
      RETURN
C
C     REEL
C     ====
 300  DO 330 I=1,NBV
         R = RMCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
         CALL FICARC( NFFICA , R )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 330  CONTINUE
      RETURN
C
C     REEL2
C     =====
 400  DO 430 I=1,NBV
         CALL FICADC( NFFICA , RMCN( MNTAMS(LHTMS) + LDTS(LHTMS) )   )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 430  CONTINUE
      RETURN
C
C     XYZ
C     ===
 500  DO 530 I=1,NBV
         NVALEN  = MNTAMS(LHTMS) + LDTS(LHTMS)
         CALL FICAXC( NFFICA , RMCN( NVALEN ) )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 530  CONTINUE
      RETURN
C
C     TYPEOBJET
C     =========
C     CONVERSION DES NOMS DANS LE LEXIQUE
 600  DO 640 I=1,NBV
C        OUVERTURE DANS LE LEXIQUE
         NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
         IF( NVALEN .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) =  KNOM
            KERR(2) = 'NON DANS LE LEXIQUE ADAM'
            CALL LEREUR
            NRETOU = 1
         ENDIF
C        LE TYPE DE L'OBJET
         CALL FICACK( NFFICA , NMTYOB(NVALEN) )
C        OUVERTURE DE SON LEXIQUE
         CALL LXNLOU( NTADAM , NVALEN , NTLX , MNLX )
C
C        LE NOM DE L'OBJET
         NVALEN = MCN( MNTAMS(LHTMS) + LDTS(LHTMS) + 1 )
         IF( NVALEN .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOM
            KERR(2) ='NON DANS LE LEXIQUE ADAM'
            CALL LEREUR
            CALL LXIM0( MNLX )
            NRETOU = 1
         ENDIF
C
C        L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
         MN = MNLX + MCN( MNLX ) * NVALEN
C        LE NOM CONVERTI D'ENTIERS EN CARACTERES
C        KNOM OFFRE UNE PROTECTION DE DEBORDEMENT
         CALL ENTNOM( MCN(MNLX+2) , MCN(MN) , KNOM )
C
C        ATTENTION: LE NOM EST ECRIT , PAS LE NUMERO
C                   EN EFFET LE NUMERO DANS LE LEXIQUE EST LOCAL
C                   LE NOM LUI RESTE GLOBAL
C
         CALL FICAKC( NFFICA , KNOM )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 640  CONTINUE
      RETURN
C
C     ERREUR TYPE INCORRECT ICI
 9000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
      KERR(1) = 'NCTYVA:TYPE INCORRECT '//KERR(MXLGER)(1:4)
      CALL LEREUR
      NRETOU = 1
      END
