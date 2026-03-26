      SUBROUTINE CNTYVA( NFFICA , NOIDEN , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CONVERTIR LES CARACTERES DE L'IDENTIFICATEUR NOIDEN LU
C ----- DANS LE BUFFER KFICA EN VALEUR NUMERIQUE
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
      include"./incl/langue.inc"
      include"./incl/td.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*80      KNOM
      CHARACTER*4       CAR4
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
         KERR(1) = 'CNTYVA:'//KERR(MXLGER)(1:4)//' TYPE INCORRECT'
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
         CALL FICACH( NFFICA , CAR4 , NRETOU )
         MCN( MNTAMS(LHTMS) + LDTS(LHTMS) ) = ICHARX( CAR4 )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 55   CONTINUE
C
C   ENTIER {(<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })}
C   ===================================================================
C     RECHERCHE ;
 100  DO 190 I=1,NBV
         CALL FICACE( NFFICA , MCN( MNTAMS(LHTMS)+LDTS(LHTMS) ) ,
     %                NRETOU )
         IF( NRETOU .NE. 0 ) RETURN
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 190  CONTINUE
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
         CALL LXNTPO( KNOM , NTLX , MNLX , N1 , N2 )
      ENDIF
      IF( NTLX .LE. 0 ) THEN
C        LE LEXIQUE N'EXISTE PAS
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CNTYVA: LEXIQUE INCONNU'
         ELSE
            KERR(1) = 'CNTYVA: UNKNOWN LEXICON'
         ENDIF
         KERR(2) = KNOM
         CALL LEREUR
         NRETOU = 1
         RETURN
      ENDIF
C
C     CONVERSION DES NOMS DANS LE LEXIQUE
      DO 240 I=1,NBV
C
C        ATTENTION: LE NOM EST ECRIT , PAS LE NUMERO
C                   EN EFFET LE NUMERO DANS LE LEXIQUE EST LOCAL
C                   LE NOM LUI RESTE GLOBAL
C        LE NOM DANS LE LEXIQUE
         CALL FICACK( NFFICA , KNOM , NRETOU )
C        RECHERCHE DE SON NUMERO DANS LE LEXIQUE
 210     CALL LXNMNO( NTLX , KNOM , NVALEN , MNLX )
         IF( NVALEN .LE. 0 ) THEN
C           IL EST DECLARE LEXIQUE DE 8 NOMS DE 24 CARACTERES
            CALL LXLXDC( NTLX , KNOM , 24 , 8 )
            GOTO 210
         ENDIF
C        LE NUMERO DANS LE LEXIQUE
         MCN( MNTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 240  CONTINUE
      RETURN
C
C     REEL
C     ====
 300  DO 330 I=1,NBV
         CALL FICACR( NFFICA , RMCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
     %              , NRETOU )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 330  CONTINUE
      RETURN
C
C     REEL2
C     =====
 400  DO 430 I=1,NBV
         CALL FICACD( NFFICA , RMCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
     %              , NRETOU )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 430  CONTINUE
      RETURN
C
C     XYZ
C     ===
 500  DO 530 I=1,NBV
         CALL FICACX( NFFICA , RMCN( MNTAMS(LHTMS) + LDTS(LHTMS) )
     %              , NRETOU )
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 530  CONTINUE
      RETURN
C
C     TYPEOBJET
C     =========
C     CONVERSION DES NOMS DANS LE LEXIQUE EN ENTIERS
 600  DO 640 I=1,NBV
C
C        ATTENTION: LE NOM EST LU , PAS LE NUMERO
C                   EN EFFET LE NUMERO DANS LE LEXIQUE EST LOCAL
C                   LE NOM LUI RESTE GLOBAL
C
C        LE NOM DU TYPE DE L'OBJET
         CALL FICACK( NFFICA , KNOM , NRETOU )
C        RECHERCHE DE SON NUMERO DANS LE LEXIQUE ADAM
 610     CALL LXNMNO( NTADAM , KNOM , NVALEN , N1 )
         IF( NVALEN .LE. 0 ) THEN
C           IL EST DECLARE LEXIQUE DE 8 NOMS DE 24 CARACTERES
            CALL LXLXDC( NTADAM , KNOM , 24 , 8 )
            GOTO 610
         ENDIF
C        OUVERTURE DU LEXIQUE DU TYPE DE L'OBJET
         CALL LXNLOU( NTADAM , NVALEN , NTLX , MNLX )
C        LE NUMERO DANS LE LEXIQUE
         MCN( MNTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
C        LE NOM DE L'OBJET
         CALL FICACK( NFFICA , KNOM , NRETOU )
 620     CALL LXNMNO( NTLX   , KNOM , NVALEN , MNLX )
         IF( NVALEN .LE. 0 ) THEN
C           IL EST DECLARE LEXIQUE DE 8 NOMS DE 24 CARACTERES
            CALL LXLXDC( NTLX , KNOM , 24 , 8 )
            GOTO 620
         ENDIF
C        LE NUMERO DANS LE LEXIQUE
         MCN( MNTAMS(LHTMS) + LDTS(LHTMS) + 1 ) = NVALEN
C
         LDTS(LHTMS) = LDTS(LHTMS) + MOT
 640  CONTINUE
      RETURN
C
C     ERREUR TYPE INCORRECT ICI
 9000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
      KERR(1) = 'CNTYVA:TYPE'//KERR(MXLGER)(1:4)//' INCORRECT'
      CALL LEREUR
      NRETOU = 1
      END
