      SUBROUTINE MOTYVA( NOIDEN , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MODIFIER OU CREER LA VALEUR DE L'IDENTIFICATEUR NOIDEN
C -----     DU TABLEAU IDENT
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
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/typobj.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/msvaau.inc"
      include"./incl/pilect.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*10      NMTYOB
      CHARACTER*80      KNOM,KNOMI
      CHARACTER*84      KLGFRA
      CHARACTER*4       KAUX
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
C
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
      INITIK = 0
      NRETOU = 0
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
         CALL TNMCAU( 'MOTS' , NBVATC(LHTMS) ,  L ,
     %                NBVATC(LHTMS) , MCTAMS(LHTMS) )
         NBVATC(LHTMS) = L
      ENDIF
C
C     DEMANDE SELON LE TYPE
      IF( NOTYPE .LE. 0 .OR. NOTYPE .GT. NBTYPV ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MOTYVA: '//KERR(MXLGER)(1:4)
     %              //' TYPE INCORRECT'
         ELSE
            KERR(1) = 'MOTYVA:'//KERR(MXLGER)(1:4)
     %              //' INCORRECT TYPE'
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
      GOTO( 9000 , 22   , 9000 , 100 , 300 , 400 , 9000 , 9000 , 9000,
     %      9000 , 200  , 500  , 600 ) , NOTYPE
C          LOGIQUE   => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL      => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2 => 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPE_OBJET=>13
C
C     CARACTERE4
 22   CALL CHMOT( 'caractere' , NL,NC, NLD,NCD, NLF,NCF )
      DO 40 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 30
         NCVALS = 0
         CALL LIRCAR( NCVALS , KNOMI )
         IF( NCVALS .EQ. -1 ) GOTO 9900
         MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = ICHARX( KNOMI )
 30      LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 40   CONTINUE
      GOTO 9995
C
C   ENTIER {(<INTERVAL_ENT> : <CHAINE> {, <INTERVAL_ENT> : <CHAINE> })}
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
         CALL AFCHAI( NLD,NCD, NLF,NCF )
         CALL AFLIGN
C
         DO 170 I=1,NBV
            IF( NOMUET .EQ. 0 ) GOTO 150
            NCVALS = 0
            CALL LIRENT(NCVALS, NVALEN )
            IF( NCVALS .EQ. -1 ) GOTO 9900
C
C           VERIFICATION QUE NVALEN EST UNE VALEUR PERMISE
C           LE CARACTERE ( POUR DEMARRER LA RECHERCHE
            NL = NLD
            NC = NCD
C
 110        CALL CHINTR( NL,NC,NL1,NC1,NL2,NC2,N1,N2 )
            IF( NL1 .LE. 0 ) THEN
C              VALEUR NON COMPRISE DANS L'INTERVALLE DE DEFINITION
               NBLGRC(NRERR) = 3+NL-NLD
               WRITE(KERR(2)(1:12),'(I12)') NVALEN
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='MOTYVA: ENTIER ' // KERR(2)(1:12)
                  KERR(2) ='DIFFERENT DES VALEURS LICITES'
               ELSE
                  KERR(1) ='MOTYVA: INTEGER ' // KERR(2)(1:12)
                  KERR(2) ='VALUE INCORRECT'
               ENDIF
               DO 111 J=NLD,NL
                  KERR(3+J-NLD) = KTD(J)
 111           CONTINUE
               CALL LEREUR
               GOTO 9900
            ENDIF
C
            IF( NVALEN .LT. N1 .OR. NVALEN .GT. N2 ) THEN
C              VALEUR INCORRECTE.ON PASSE A LA SUIVANTE
               NL = NL2
               NC = NC2
C              RECHERCHE DU :
               CALL CHCAR( ':' , NL , NC )
C              RECHERCHE DE LA CHAINE '...'
               CALL CHCHAI( NL,NC,NL1,NC1,NL2,NC2 )
               NL = NL2
               NC = NC2
C              RECHERCHE DE ,
               CALL CAR1NB( NL , NC )
               IF( KTD(NL)(NC:NC) .NE. ',' ) THEN
C                 LA FIN DES OPTIONS EST ATTEINTE SANS VALEUR CORRECTE
C                 ( INTERVAL_ENT : CHAINE {,INTERVAL_ENT : CHAINE} )
C                 INCORRECT OU VALEUR LUE INCORRECTE
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:12),'(I12)') NVALEN
                  KERR(1) = KIDENT(NOIDEN)//'='//KERR(MXLGER)(1:12)
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(2) = 'VALEUR INCORRECTE => A REDONNER'
                  ELSE
                     KERR(2) = 'INCORRECT VALUE => GIVE AGAIN'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
               GOTO 110
            ENDIF
C
C           LA VALEUR EST CORRECTE
            MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
C           SAUVEGARDE DE LA SIGNIFICATION DE L'ENTIER FRAPPE SUR LE FICHIER
            IF( LHLECT .EQ. 1 ) THEN
C              RECHERCHE DU :
               CALL CHCAR( ':' , NL2 , NC2 )
               NLL = NL2
               NCC = NC2
C              RECHERCHE DE LA CHAINE '...'
               CALL CHCHAI( NLL,NCC,NL1,NC1,NL2,NC2 )
               L      = NUDCNB( KTD(NL2) )
               KLGFRA = '{ ' // KTD(NL2)(NC2:L) // ' }'
               CALL SANSDBL( KLGFRA, L )
               WRITE(NFFRAP,*) KLGFRA(1:L)
               INITIK = 1
            ENDIF
C
 150        LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 170     CONTINUE
C
      ELSE
C        IL N'EXISTE PAS DE ( )
         DO 190 I=1,NBV
            IF( NOMUET .EQ. 0 ) GOTO 180
            NCVALS = 0
            CALL LIRENT(NCVALS, NVALEN )
            IF( NCVALS .EQ. -1 ) GOTO 9900
C           LA VALEUR EST COPIEE DANS LE TABLEAU MS COPIE EN MC
            MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
 180        LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 190     CONTINUE
      ENDIF
      GOTO 9995
C
C     ^NOM_LEXIQUE
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
         GOTO 9900
      ENDIF
C
      IF( NCODER .EQ. 1 ) THEN
C        OUVERTURE DU LEXIQUE DU TD
         CALL TNOUVR( KNOM , KAUX , NTLX , MNLX )
      ELSE
C        RECHERCHE ET OUVERTURE DU LEXIQUE FILS
         CALL LXNTPN( KNOM , NTLX , MNLX , N1 , N2 )
      ENDIF
      IF( NTLX .LE. 0 ) THEN
C        LE LEXIQUE N'EXISTE PAS
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='MOTYVA: LEXIQUE INCONNU'
         ELSE
            KERR(1) ='MOTYVA: UNKNOWN LEXIQUE'
         ENDIF
         KERR(2) = KNOM
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     LECTURE DES NOMS DANS LE LEXIQUE
      DO 240 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 235
         IF( LANGAG .EQ. 0 ) THEN
            KNOMI = 'NOM DANS LE LEXIQUE '
     %            // KTD(NLD)(NCD:NCF)//' =?'
         ELSE
            KNOMI = 'NAME in the LEXIQUE '
     %            // KTD(NLD)(NCD:NCF)//' =?'
         ENDIF
         CALL INVITD( KNOMI )
C
 232     NCVALS = 0
         CALL LIRCAR( NCVALS , KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
C
C        OUVERTURE DE KNOM DANS LE LEXIQUE
         CALL LXNMNO( NTLX , KNOM , NVALEN , MNLXK )
         IF( NVALEN .LE. 0 ) THEN

C           NOM INCORRECT NON RETROUVE DANS LE LEXIQUE
            NBLGRC( NRERR ) = 0
            N2 = INDEX( KNOM, ' ' )
            IF( N2 .LE. 0 ) N2=1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(NBLGRC(NRERR)+1) = 'NOM INCONNU: ' // KNOM(1:N2)
            KERR(NBLGRC(NRERR)+2) = 'A choisir PARMI:'
            ELSE
            KERR(NBLGRC(NRERR)+1) = 'UNKNOWN NAME: '// KNOM(1:N2)
            KERR(NBLGRC(NRERR)+2) = 'MUST be CHOSEN AMONG:'
            ENDIF
            NBLGRC(NRERR) = NBLGRC(NRERR) + 2
            CALL LXIM0( MNLXK )

            IF (LHLECT.GT. 1 .AND. LHLECT.LE.MXLECT ) THEN
C              LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C              AFFICHAGE DU RETOUR AU CLAVIER
               CALL QUIFIC
            ENDIF

            GOTO 232
         ENDIF
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
                   KERR(2) = 'LEXIQUE REFERS TO ITSELF'
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
             KLGFRA = '{ DANS LE LEXIQUE ' // KTD(NLD)(NCD:NCF) // ' }'
             ELSE
             KLGFRA = '{ IN THE LEXICON ' // KTD(NLD)(NCD:NCF) // ' }'
             ENDIF
             CALL SANSDBL( KLGFRA, L )
             WRITE(NFFRAP,*) KLGFRA(1:L)
             INITIK = 1
         ENDIF
C
C        LA VALEUR EST COPIEE DANS LE TABLEAU MS COPIE EN MC
         MCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = NVALEN
C
 235     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 240  CONTINUE
      GOTO 9995
C
C     REEL
 300  CALL CHMOT( 'reel' , NL,NC, NLD,NCD, NLF,NCF )
      DO 330 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 320
         NCVALS = 0
         CALL LIRRSP( NCVALS , R )
         IF( NCVALS .EQ. -1 ) GOTO 9900
         RMCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) = R
C
 320     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 330  CONTINUE
      GOTO 9995
C
C     REEL2
 400  CALL CHMOT( 'reel2' , NL,NC, NLD,NCD, NLF,NCF )
      DO 430 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 420
         NCVALS = 0
         CALL LIRRDP( NCVALS, RMCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) )
         IF( NCVALS .EQ. -1 ) GOTO 9900
C
 420     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 430  CONTINUE
      GOTO 9995
C
C     XYZ
 500  CALL CHMOT( 'xyz' , NL,NC, NLD,NCD, NLF,NCF )
      DO 530 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 525
         NCVALS = 0
         CALL LIRXYZ( NCVALS ,  RMCN( MCTAMS(LHTMS) + LDTS(LHTMS) ) )
         IF( NCVALS .EQ. -1 ) GOTO 9900
 525     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 530  CONTINUE
      GOTO 9995
C
C     TYPEOBJET
 600  CALL CHMOT( 'typeobjet' , NL,NC, NLD,NCD, NLF,NCF )
      DO 670 I=1,NBV
         IF( NOMUET .EQ. 0 ) GOTO 650
C        LE NUMERO DU TYPE DE L'OBJET
 610     CALL INVITE( 90 )
         NCVALS = 0
         CALL LIRCAR( NCVALS, KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
         N1 = NTYOBJ( KNOM )
         IF( N1 .EQ. 0 ) THEN
C           ERREUR A REDONNER
            GOTO 610
         ENDIF
C
C        OUVERTURE DU LEXIQUE DANS LE LEXIQUE ADAM
         CALL LXLXOU( NTADAM , NMTYOB( N1 ) , NTLX , MNLX )
C        LECTURE DU NOM DE L'OBJET DANS SON LEXIQUE
         IF( LHLECT .EQ. 1 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               KNOMI = 'NOM DE '//NMTYOB(N1)//'=?'
            ELSE
               KNOMI = NMTYOB(N1)//'''s NAME=?'
            ENDIF
            CALL INVITD( KNOMI )
         ENDIF
 632     NCVALS = 0
         CALL LIRCAR( NCVALS , KNOM )
         IF( NCVALS .EQ. -1 ) GOTO 9900
C
C        OUVERTURE DANS LE LEXIQUE
         CALL LXNMNO( NTLX , KNOM , NVALEN , N2 )
         IF( NVALEN .LE. 0 ) THEN
            L = INDEX( KNOM, ' ' )
            IF( LANGAG .EQ. 0 ) THEN
            KERR(NBLGRC(NRERR)+1) = 'NOM INCONNU: ' // KNOM(1:L)
            KERR(NBLGRC(NRERR)+2) = 'A CHOISIR PARMI:'
            ELSE
            KERR(NBLGRC(NRERR)+1) = 'UNKNOWN NAME: '// KNOM(1:L)
            KERR(NBLGRC(NRERR)+2) = 'MUST BE CHOSEN AMONG:'
            ENDIF
            NBLGRC(NRERR) = NBLGRC(NRERR) + 2
            CALL LXIM0( N2 )
            GOTO 632
         ENDIF
C
C        LES VALEURS SONT COPIEES DANS LE TABLEAU MS COPIE EN MC
         N2 = MCTAMS(LHTMS) + LDTS(LHTMS)
         MCN( N2     ) = N1
         MCN( N2 + 1 ) = NVALEN
C
 650     LDTS(LHTMS) = LDTS(LHTMS) + MOTS
 670  CONTINUE
      GOTO 9995
C
C     ERREUR TYPE INCORRECT ICI
 9000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'MOTYVA:TYPE'//KERR(MXLGER)(1:4)//' INCORRECT'
      ELSE
         KERR(1) = 'MOTYVA: INCORRECT TYPE '//KERR(MXLGER)(1:4)
      ENDIF
      CALL LEREUR
C
C     ABANDON DE LA LECTURE
 9900 NRETOU = 1
C
C     SAUVEGARDE SUR LE FICHIER D'UN COMMENTAIRE ECLAIRANT LA DONNEE
 9995 IF( LHLECT .EQ. 1 ) THEN
         IF( NRETOU.EQ.0 .AND. INITIK.EQ.0 .AND. NOMUET.NE.0 ) THEN
            CALL SANSDBL( KLGFRA, L )
            WRITE(NFFRAP,*) KLGFRA(1:L)
         ENDIF
      ENDIF
C
      RETURN
      END
