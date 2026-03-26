      SUBROUTINE ALEXPD( NLD   , NCD   , NLMAX,  NCMAX,
     %                   NCODEV, DBLVAL, NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : ANALYSE LEXICALE DE LA DONNEE DANS KLG
C ----- TRADUCTION DES CARACTERES LUS EN UNE SUITE CHAINEE D'OPERATIONS
C
C       CF ~LU/GRAMMAIRE DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NLD,NCD  : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C MODIFIES :
C ----------
C NLMAX,NCMAX : POSITION DANS KLG DU DERNIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NLD,NCD: POSITION DANS KLG DU DERNIER CARACTERE
C          EFFECTIVEMENT TRAITE
C NCODEV : 0 DBLVAL N'EST PAS INITIALISEE
C          1 DBLVAL EST INITIALISEE
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C NRETOU :  0 SI PAS D'ERREUR
C          >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON /UNITES/LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NOAFTS,NUNIT(26)
      DOUBLE PRECISION  DBLVAL
      CHARACTER*1       CAR
C
C     LA HAUTEUR DE LA PILE DES EXPRESSIONS
      LHEXPD = 1
      NBEXPD = 1
      NOPER  = 0
      NRETOU = 0
      NCODEV = 0
C     LE NOMBRE DE ( NON FERMEES DE L'EXPRESSION NBEXPR
      NBPAOU( LHEXPD ) = 0
C     LA POSITION DANS NCOPER DE LA DERNIERE ( OU ) RENCONTREE
      NDERP2( LHEXPD ) = 0
C     L'ETAT ACTUEL AVANT LA PROCHAINE OPERATION
      NETAT ( LHEXPD ) = 0
C     LA PREMIERE OPERATION DE L'EXPRESSION
      LSEXPD(1,NBEXPD) = NOPER + 1
      LSEXPD(2,NBEXPD) = 0
      LNEXPD(LHEXPD)   = 1
C
      NL = NLD
      NC = NCD
      IF( KLG(NL)(NC:NC) .EQ. ' ' ) CALL CARPNB( NL , NC )
C
C     AFFICHER NE PEUT ETRE QU'AU DEBUT DE L'EXPRESSION ARITHMETIQUE
      IF( KLG(NL)(NC:NC+7) .EQ. 'AFFICHER' .OR.
     %    KLG(NL)(NC:NC+6) .EQ. 'DISPLAY' ) THEN
C
C        AFFICHER ou DISPLAY EXP_ARITH {, EXP_ARITH } ;
         IF( LANGAG .EQ. 0 ) THEN
            KLG(NL)(NC:NC+8) = 'AFFICHER('
         ELSE
            KLG(NL)(NC:NC+7) = 'DISPLAY('
         ENDIF
C        UNE ) EST INTERCALEE JUSTE AVANT NLMAX,NCMAX
C        DECALAGE DE 1 CARACTERE VERS LA GAUCHE
         IF( KLG(NLMAX)(NCMAX:NCMAX) .EQ. ' ' ) THEN
            NLF = NLMAX
            NCF = NCMAX
         ELSE
C           DECALAGE POUR AVOIR UN BLANC EN CETTE POSITION NLMAX,NCMAX
            NLF = NLMAX
            NCF = NCMAX
            CALL DEPKLG( NLF , NCF , NLMAX , NCMAX , 1 )
         ENDIF
         KLG(NLF)(NCF:NCF) = ')'
C        PSEUDO FONCTION UTILISATEUR DE NUMERO 0
         NOFONC = 0
         NETAT(LHEXPD) = 6
         NOPER = NOPER + 1
         NCOPER(1,NOPER) = NETAT(LHEXPD)
         NCOPER(2,NOPER) = 0
         NCOPER(3,NOPER) = NOFONC
         NCOPER(4,NOPER) = 0
C        PASSAGE EN FIN D'IDENTIFICATEUR DE FONCTION
         NC = NC + 7
CCCC
CCCC        TEMOIN D'INSTRUCTION AFFICHER NOM de VARIABLE ou EXPRESSION
CCC         NOAFTS = 1
         GOTO 100
      ENDIF
C
C     NLD,NCD 1-ER CARACTERE DE L'EXPRESSION
      NL = NLD
      NC = NCD
C
C     ANALYSE DE L'OPERATION SUIVANTE
C     ===============================
C     SI CARACTERE BLANC PASSAGE AU PROCHAIN CARACTERE NON BLANC
 100  IF( KLG(NL)(NC:NC) .EQ. ' ' ) CALL CARPNB( NL , NC )
      IF( NL .GT. NLMAX ) GOTO 8000
      IF( NL .EQ. NLMAX .AND. NC .GT. NCMAX ) GOTO 8000
C
C     LE CARACTERE A TRAITER
 105  CAR = KLG(NL)(NC:NC)
C
      IF( CAR .EQ. ';' ) GOTO 8000
      IF( CAR .EQ. '(' ) THEN
C
C        ATTENTION : 2 TYPES DE (     ( ) ET (  ,  ,  )
C                    CELLES DES EXPRESSIONS ET
C                           DES PARAMETRES D'UNE FONCTION
         IF(  NETAT(LHEXPD) .EQ. 6 .OR.
     %       (NETAT(LHEXPD) .GE. 104 .AND. NETAT(LHEXPD) .LT. 200 ) .OR.
     %       (NETAT(LHEXPD) .GE. 215 ) ) THEN
C
C           ( D'UNE FONCTION
C           ----------------
            NETA0 = 7
            NETAT(LHEXPD) = NETA0
C
C           DERRIERE CETTE ( EXISTE UNE NOUVELLE <EXPRESS_D>
C           CETTE NOUVELLE EXPRESSION EST EMPILEE
            LHEXPD = LHEXPD + 1
            IF( LHEXPD .GT. LXEXPD ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
               KERR(1) =
     %        'LU: PILE DES EXPRESSIONS SATUREE.AUGMENTER LXEXPD'
            ELSE
               KERR(1) =
     %        'LU: STACK of EXPRESSIONS is SATURATED. AUGMENT LXEXPD'
            ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
C           MISE A 0 DES VALEURS POUR LE RETOUR
            NBPAOU(LHEXPD) = 0
            IF( NDERP2(LHEXPD-1) .GT. 0 ) THEN
               NDERP2(LHEXPD) = NDERP2(LHEXPD-1)
            ELSE
               NDERP2(LHEXPD) = 0
            ENDIF
            NETAT (LHEXPD) = NETA0
C
C           AJOUT D'UNE NOUVELLE EXPRESSION
            NBEXPD = NBEXPD + 1
            IF( NBEXPD .GT. MXEXPD ) GOTO 9100
C           LE DEBUT DANS NCOPER DE L'EXPRESSION
            LSEXPD(1,NBEXPD) = NOPER + 2
            LSEXPD(2,NBEXPD) = 0
C           LE NIVEAU DE L'EXPRESSION DE 1 A NBEXPD
            LNEXPD(LHEXPD) = NBEXPD
C           L'OPERATION EST AJOUTEE
            NOPER = NOPER + 1
            IF( NOPER .GE. MXNCOP ) GOTO 9200
            NCOPER(1,NOPER) = NETAT(LHEXPD)
            NCOPER(3,NOPER) = 0
C           CHAINAGE AVANT ARRIERE AVEC LA ( OU , QUI PRECEDE
            NCOPER(2,NOPER) = NDERP2(LHEXPD)
            IF( NDERP2(LHEXPD) .GT. 0 ) NCOPER(3,NDERP2(LHEXPD)) = NOPER
C           LA DERNIERE PARENTHESE VUE
            NDERP2(LHEXPD) = NOPER
C
         ELSE
C
C           ( D'UNE EXPRESSION ARITHMETIQUE
C           -------------------------------
            NETAT(  LHEXPD ) = 1
            NBPAOU( LHEXPD ) = NBPAOU( LHEXPD ) + 1
            NOPER = NOPER + 1
            IF( NOPER .GE. MXNCOP ) GOTO 9200
            NCOPER(1,NOPER) = NETAT(LHEXPD)
            NCOPER(3,NOPER) = 0
C           CHAINAGE AVANT ARRIERE AVEC LA ( OU ) QUI PRECEDE
            NCOPER(2,NOPER) = NDERP2(LHEXPD)
            IF( NDERP2(LHEXPD) .GT. 0 ) NCOPER(3,NDERP2(LHEXPD)) = NOPER
C           LA DERNIERE PARENTHESE VUE
            NDERP2(LHEXPD) = NOPER
         ENDIF
         GOTO 1000
C
      ELSE IF( CAR .EQ. ')' ) THEN
C
C        ATTENTION : 2 TYPES DE )     ( ) ET (  ,  ,  )
C                    CELLES DES EXPRESSIONS ET
C                           DES PARAMETRES D'UNE FONCTION
         IF( NBPAOU(LHEXPD) .EQ. 0 ) THEN
C
C           ) D'UNE FONCTION
C           ----------------
            NETAT(LHEXPD) = 9
C           FIN DE L'EXPRESSION DU DERNIER OPERANDE
            LSEXPD(2,LNEXPD(LHEXPD)) = NOPER
C           RECHERCHE DU NOMBRE DE PARAMETRES DE CETTE FONCTION
            NBPAR = 1
            NOPE  = NDERP2(LHEXPD)
C           LES ( ) , SONT REMONTEES JUSQU'A RETROUVER LA ( INITIALE
C           DE LA FONCTION   F( (...)  ... , ... (...)  ... )
C                             |                            |
C           ICI NOPE POINTE SUR L'OPERATION AVANT LA ) FINALE
 170        NCO = NCOPER(1,NOPE)
            IF( NCO .EQ. 8 ) THEN
C              NOPE EST UNE , DANS UNE FONCTION
C              LA ( OU , QUI PRECEDE
               NOPE  = NCOPER(2,NOPE)
C              UN PARAMETRE DE PLUS
               NBPAR = NBPAR + 1
               GOTO 170
            ELSE IF( NCO .EQ. 1 .OR. NCO .EQ. 2 ) THEN
C              NOPE EST ( OU ) D'EXPRESSION ARITHMETIQUE
               NOPE = NCOPER(2,NOPE)
               GOTO 170
            ELSE IF( NCO .EQ. 7) THEN
C              NOPE EST ( INITIALE DE LA FONCTION
               IF( NCOPER(1,NOPE-1) .EQ. 6 ) THEN
C                 OUI : STOCKAGE DU NOMBRE DE VARIABLES ET
C                       DU NUMERO DE LA FONCTION
C                       SELON LE CODAGE :  300 + NOFONC + 1000 * NBPAR
                  NCOPER(1,NOPE-1) = 300+NCOPER(3,NOPE-1) + 1000 * NBPAR
               ENDIF
            ENDIF
C
C           RETOUR A L'EXPRESSION DE DEBUT DE LA FONCTION
            NETAT(LHEXPD-1) = 9
C
         ELSE
C
C           ) D'UNE EXPRESSION ARITHMETIQUE
C           -------------------------------
            NETAT(LHEXPD) = 2
            NBPAOU( LHEXPD ) = NBPAOU( LHEXPD ) - 1
         ENDIF
C
         NOPER = NOPER + 1
         IF( NOPER .GT. MXNCOP ) GOTO 9200
         NCOPER(1,NOPER ) = NETAT(LHEXPD)
C        CHAINAGE AVANT ARRIERE AVEC LA ( OU ) OU , QUI PRECEDE
         NCOPER(2,NOPER ) = NDERP2(LHEXPD)
         NCOPER(3,NOPER ) = 0
         NCOPER(3,NDERP2(LHEXPD)) = NOPER
C        LA DERNIERE PARENTHESE VUE
         NDERP2(LHEXPD) = NOPER
C        SI FIN DE FONCTION ON DEPILE LES EXPRESSIONS
C        ARITHMETIQUES PARAMETRES DE LA FONCTION
         IF( NETAT(LHEXPD) .EQ. 9 )  LHEXPD = LHEXPD - 1
         GOTO 1000
C
      ELSE IF( CAR .EQ. ',' ) THEN
C        FIN D'UNE EXPRESSION OPERANDE ET DEBUT D'UNE AUTRE
         IF( NBPAOU(LHEXPD) .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: ( NON FERMEE OU ) EN TROP DANS'
            ELSE
               KERR(1) = 'LU: ( NOT CLOSED or TOO ) in'
            ENDIF
            KERR(2) =  KLG(NL)(1:NC)
            CALL LEREUR
            GOTO 9900
         ENDIF
C        FIN DE L'EXPRESSION EN COURS AVANT ,
         LSEXPD(2,LNEXPD(LHEXPD)) = NOPER
         NETAT (LHEXPD) = 8
C        ETAT DANS LA FONCTION
         NETAT (LHEXPD-1) = 8
C        L'OPERATION EN COURS
         NOPER = NOPER + 1
         IF( NOPER .GT. MXNCOP ) GOTO 9200
         NCOPER(1,NOPER) = NETAT(LHEXPD)
         NCOPER(2,NOPER) = NDERP2(LHEXPD)
         NCOPER(3,NOPER) = 0
         NCOPER(3,NDERP2(LHEXPD)) = NOPER
C        L'EXPRESSION SUIVANTE
C        MNOMBRE DE ( DE L'EXPRESSION
         NBPAOU(LHEXPD) = 0
C        MISE A JOUR DE LA DERNIERE PARENTHESE DE L'EXPRESSION
C        PRECEDENTE
         NDERP2(LHEXPD) = NOPER
         NBEXPD = NBEXPD + 1
         IF( NBEXPD .GT. MXEXPD ) GOTO 9100
C        LE DEBUT DANS NCOPER DE L'EXPRESSION
         LSEXPD(1,NBEXPD) = NOPER + 1
         LSEXPD(2,NBEXPD) = 0
         LNEXPD(LHEXPD) = NBEXPD
         GOTO 1000
C
      ELSE IF( CAR .EQ. '''' ) THEN
C        PARAMETRE 'CHAINE DE '' CARACTERES' D'UNE FONCTION
         IF( NETAT(LHEXPD) .NE. 7 .AND.
     %       NETAT(LHEXPD) .NE. 8 ) THEN
             NBLGRC(NRERR) = 1
             IF( LANGAG .EQ. 0 ) THEN
      KERR(1) =
     %'LU: UNE CONSTANTE CARACTERES DOIT ETRE UN PARAMETRE DE FONCTION'
             ELSE
      KERR(1) =
     %'LU: A STRING of CHARACTERS MUST BE a FUNCTION PARAMETER'
             ENDIF
             CALL LEREUR
             GOTO 9900
         ENDIF
C        RECHERCHE DE LA FIN DE CHAINE
         CALL CHPACH( NL , NC-1 , NL1 , NC1 , NLF , NCF )
         IF( NL1 .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: CHAINE NON TERMINEE PAR '' '
            ELSE
               KERR(1) = 'LU: STRING NOT FINISHED by '' '
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
         IF( NLF .GT. NL1 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =
     %     'LU: PARAMETRE CHAINE D''UNE FONCTION A CHEVAL SUR 2 LIGNES'
            ELSE
            KERR(1) =
     %     'LU: STRING PARAMETER ONTO 2 LINES is FORBIDDEN'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C
C        LA CONSTANTE CHAINE EXISTE-ELLE DANS LA TABLE DES CONSTANTES ?
C        SI ELLE N'EXISTE PAS ELLE Y EST AJOUTEE AVEC LE NUMERO NOPACH
         CALL NUPACH( KLG(NL1)(NC1+1:NCF-1) , NOPACH )
         IF( NOPACH .EQ. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =
     %     'LU: SATURATION DE LA TABLE DES CONSTANTES CHAINES'
            ELSE
            KERR(1) =
     %     'LU: ARRAY of STRINGS of CHARACTERS is SATURATED'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
         NETAT(LHEXPD) = 14
         NOPER = NOPER + 1
         IF( NOPER .GT. MXNCOP ) GOTO 9200
         NCOPER(1,NOPER) = NETAT(LHEXPD)
         NCOPER(2,NOPER) = 0
         NCOPER(3,NOPER) = NOPACH
         NCOPER(4,NOPER) = 0
C        PASSAGE EN FIN DE CONSTANTE
         NL = NLF
         NC = NCF
         GOTO 1000
C
      ELSE IF( CAR .EQ. '*' ) THEN
         IF( NC .LT. NBCALI ) THEN
            IF( KLG(NL)(NC:NC+1) .EQ. '**' ) THEN
C              OPERATEUR **
               IF( NETAT(LHEXPD) .EQ. 1     .OR.
     %             NETAT(LHEXPD) .GE. 13 ) GOTO 9000
               NETAT(LHEXPD) = 214
               NC = NC + 1
               GOTO 800
            ENDIF
         ENDIF
C        OPERATEUR *
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 212
         GOTO 800
C
      ELSE IF( CAR .EQ. '/' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 213
         GOTO 800
C
      ELSE IF( CAR .EQ. '+' ) THEN
         IF( NETAT(LHEXPD) .EQ. 102 .OR. NETAT(LHEXPD) .EQ. 103 .OR.
     %      (NETAT(LHEXPD) .GE. 210 .AND. NETAT(LHEXPD) .LE. 214 ) .OR.
     %       NETAT(LHEXPD) .EQ. 14
     %     ) GOTO 9000
         IF( NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .EQ. 101 .OR.
     %      (NETAT(LHEXPD) .GE. 104 .AND. NETAT(LHEXPD) .LE. 209 ) .OR.
     %      (NETAT(LHEXPD) .GE. 215) ) THEN
C           + UNAIRE
            NETAT(LHEXPD) = 102
         ELSE
C           + BINAIRE
            NETAT(LHEXPD) = 210
         ENDIF
         GOTO 800
C
      ELSE IF( CAR .EQ. '-' ) THEN
         IF( NETAT(LHEXPD) .EQ. 102 .OR. NETAT(LHEXPD) .EQ. 103 .OR.
     %      (NETAT(LHEXPD) .GE. 210 .AND. NETAT(LHEXPD) .LE. 214 ) .OR.
     %       NETAT(LHEXPD) .EQ. 14
     %     ) GOTO 9000
         IF( NETAT(LHEXPD) .EQ. 0   .OR.
     %       NETAT(LHEXPD) .EQ. 1   .OR.  NETAT(LHEXPD) .EQ. 7     .OR.
     %       NETAT(LHEXPD) .EQ. 8   .OR.  NETAT(LHEXPD) .EQ. 101   .OR.
     %      (NETAT(LHEXPD) .GE. 104 .AND. NETAT(LHEXPD) .LE. 209 ) .OR.
     %       NETAT(LHEXPD) .GE. 215   )   THEN
C           - UNAIRE  => AJOUT DE ZERO AVANT => 0 - BINAIRE
ccc         NETAT(LHEXPD) = 103
            IF( NC .GT. 1 ) THEN
               IF( KLG(NL)(NC-1:NC-1) .EQ. ' ' ) THEN
                  NC = NC - 1
                  KLG(NL)(NC:NC) ='0'
                  GOTO 105
               ENDIF
            ENDIF
C           DECALAGE POUR AVOIR UN '0' EN CETTE POSITION NL,NC
            NLF  = NL
            NCF  = NC
            CALL DEPKLG( NLF, NCF, NLMAX, NCMAX, 1 )
            KLG(NLF)(NCF:NCF) ='0'
            GOTO 105
         ELSE
C           - BINAIRE
            NETAT(LHEXPD) = 211
         ENDIF
         GOTO 800
C
C     LES OPERATEURS A 2 CARACTERES D'ABORD
      ELSE IF( KLG(NL)(NC:NC+1) .EQ. '<=' .OR.
     %         KLG(NL)(NC:NC+1) .EQ. '=<' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 205
         NC = NC + 1
         GOTO 800
C
      ELSE IF( KLG(NL)(NC:NC+1) .EQ. '<>' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 207
         NC = NC + 1
         GOTO 800
C
      ELSE IF( KLG(NL)(NC:NC+1) .EQ. '>=' .OR.
     %         KLG(NL)(NC:NC+1) .EQ. '=>' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 208
         NC = NC + 1
         GOTO 800
C
      ELSE IF( CAR .EQ. '<' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 204
         GOTO 800
C
      ELSE IF( CAR .EQ. '=' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 206
         GOTO 800
C
      ELSE IF( CAR .EQ. '>' ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 209
         GOTO 800
C
      ELSE IF( (NC+2 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+3).EQ.'NON '.OR. KLG(NL)(NC:NC+3).EQ.'NON(' ))
     %   .OR.  (NC+2 .EQ. NBCALI .AND.
     %   KLG(NL)(NC:NC+2) .EQ. 'NON') ) THEN
         IF( NETAT(LHEXPD) .EQ. 102 .OR.  NETAT(LHEXPD) .EQ. 103 .OR.
     %      (NETAT(LHEXPD) .GE. 210 .AND. NETAT(LHEXPD) .LE. 214 )
     %     ) GOTO 9000
         NETAT(LHEXPD) = 101
         NC = NC + 2
         GOTO 800
C
      ELSE IF( (NC+2 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+3).EQ.'NOT '.OR. KLG(NL)(NC:NC+3).EQ.'NOT(' ))
     %   .OR.  (NC+2 .EQ. NBCALI .AND.
     %   KLG(NL)(NC:NC+2) .EQ. 'NOT') ) THEN
         IF( NETAT(LHEXPD) .EQ. 102 .OR.  NETAT(LHEXPD) .EQ. 103 .OR.
     %      (NETAT(LHEXPD) .GE. 210 .AND. NETAT(LHEXPD) .LE. 214 )
     %     ) GOTO 9000
         NETAT(LHEXPD) = 101
         NC = NC + 2
         GOTO 800
C
      ELSE IF( (NC+1 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+2).EQ.'OX '.OR. KLG(NL)(NC:NC+2).EQ.'OX(' ))
     %   .OR.  (NC+1 .EQ. NBCALI .AND.
     %   KLG(NL)(NC:NC+1) .EQ. 'OX') ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 201
         NC = NC + 1
         GOTO 800
C
      ELSE IF( (NC+1 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+2).EQ.'OU '.OR. KLG(NL)(NC:NC+2).EQ.'OU(' ))
     %   .OR.  (NC+1 .EQ. NBCALI .AND.
     %   KLG(NL)(NC:NC+1) .EQ. 'OU') ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 202
         NC = NC + 1
         GOTO 800
C
      ELSE IF( (NC+1 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+2).EQ.'OR '.OR. KLG(NL)(NC:NC+2).EQ.'OR(' ))
     %   .OR.  (NC+1 .EQ. NBCALI .AND.
     %   KLG(NL)(NC:NC+1) .EQ. 'OR') ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 202
         NC = NC + 1
         GOTO 800
C
      ELSE IF( (NC+1 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+2).EQ.'ET '.OR. KLG(NL)(NC:NC+2).EQ.'ET(' ))
     %  .OR. (NC+1 .EQ. NBCALI .AND. KLG(NL)(NC:NC+1) .EQ. 'ET') ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 203
         NC = NC + 1
         GOTO 800
C
      ELSE IF( (NC+1 .LT. NBCALI .AND.
     % ( KLG(NL)(NC:NC+3).EQ.'AND '.OR. KLG(NL)(NC:NC+3).EQ.'AND(' ))
     %  .OR. (NC+1 .EQ. NBCALI .AND.KLG(NL)(NC:NC+2) .EQ. 'AND') ) THEN
         IF(NETAT(LHEXPD) .EQ. 1 .OR. NETAT(LHEXPD) .GE. 13) GOTO 9000
         NETAT(LHEXPD) = 203
         NC = NC + 2
         GOTO 800
C
      ELSE
C        RECHERCHE D'UNE FONCTION MATHEMATIQUE USUELLE
C        RECHERCHE DE LA ( SUIVANTE SUR LA MEME LIGNE
         I = INDEX( KLG(NL)(NC:NBCALI) , '(' )
         IF( I .GT. 0 ) THEN
C           ELIMINATION DES ' '
            I = NC - 2 + I
 300        IF( KLG(NL)(I:I) .EQ. ' ' ) THEN
               I = I - 1
               GOTO 300
            ENDIF
         ENDIF
C        FONCTION UNAIRE
         CALL CHOPUN( NL , NC , I , NLF , NCF , NCODE )
         IF( NCODE .GT. 0 ) THEN
            IF(NETAT(LHEXPD).GE.2 .AND. NETAT(LHEXPD).LE.6) GOTO 9000
            NETAT(LHEXPD) = NCODE
            NL    = NLF
            NC    = NCF
            GOTO 800
         ENDIF
C        FONCTION BINAIRE
         CALL CHOPBI( NL , NC , I , NLF , NCF , NCODE )
         IF( NCODE .GT. 0 ) THEN
            IF( NETAT(LHEXPD).GE.2 .AND. NETAT(LHEXPD).LE.6) GOTO 9000
            NETAT(LHEXPD) = NCODE
            NL    = NLF
            NC    = NCF
            GOTO 800
         ENDIF
C
C        LE MOT DEBUTANT EN NL,NC DOIT ETRE UNE
C        <VARIABLE>  := |<CONSTANTE_D>
C                       |<IDENT_VAR>
C                       |<VARIABLE_TMS>
C                       |<IDENT_FONC>(<EXPRESS_D> {,<EXPRESS_D> )
C
C        <CONSTANTE_D> ?
C        ---------------
         CALL CHCTE( NL , NC , NLF , NCF , I , DBLVAL )
C        SI NLF=0 ALORS <VARIABLE> NE DEBUTE PAS PAR 0...9
         IF( NLF .GT. 0 ) THEN
C           CONSTANTE RETROUVEE JUSQU'A NLF,NCF
            IF(NETAT(LHEXPD).GE.2 .AND. NETAT(LHEXPD).LE.6 ) GOTO 9000
            NLFF  = NLF
            NCFF  = NCF
            CALL CARPNB( NLFF , NCFF )
            IF( NOPER .EQ. 0 .AND.
     %          NLFF .EQ. NLMAX .AND. NCFF .GE. NCMAX ) THEN
C               IL N'EXISTE PAS D'OPERATIONS AVANT
C               L'EXPRESSION SE REDUIT A CETTE CONSTANTE
                NCODEV = 1
                RETURN
            ENDIF
            IF( NOPER .EQ. 1 .AND. NCOPER(1,1) .EQ. 103 .AND.
     %          NLFF .EQ. NLMAX .AND. NCFF .GE. NCMAX ) THEN
C               CETTE CONSTANTE NEGATIVE EST SEULE
C               L'EXPRESSION SE REDUIT A CETTE CONSTANTE
                DBLVAL = -DBLVAL
                NCODEV = 1
                RETURN
            ENDIF
C           LA CONSTANTE EXISTE-ELLE DANS LA TABLE DES CONSTANTES ?
C           SI ELLE N'EXISTE PAS ELLE Y EST AJOUTEE AVEC LE NUMERO NOCTE
            CALL NUMCTE( DBLVAL , NOCTE )
            IF( NOCTE .EQ. 0 ) THEN
C              TABLE SATUREE DES CONSTANTES
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
               KERR(1) =
     %        'LU: SATURATION DE LA TABLE DES CONSTANTES'
            ELSE
               KERR(1) =
     %        'LU: THE ARRAY of CONSTANTS is SATURATED'
            ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
            NETAT(LHEXPD) = 3
            NOPER = NOPER + 1
            IF( NOPER .GT. MXNCOP ) GOTO 9200
            NCOPER(1,NOPER) = NETAT(LHEXPD)
            NCOPER(2,NOPER) = 0
            NCOPER(3,NOPER) = NOCTE
            NCOPER(4,NOPER) = 0
C           PASSAGE EN FIN DE CONSTANTE
            NL = NLF
            NC = NCF
            GOTO 1000
         ENDIF
C
C        <IDENT_VARIABLE_UTILISATEUR OU LOCALE OU FONCTION> ?
C        ----------------------------------------------------
         CALL CHVARU( NL , NC , NLF , NCF , NOVARU )
         IF( NLF .GT. 0 ) THEN
            IF( NOVARU .GT. 0 ) THEN
C              VARIABLE RETROUVEE
               IF( NETAT(LHEXPD) .GE. 2 .AND. NETAT(LHEXPD) .LE. 6 )
     %         GOTO 9000
               IF( NOVARU .LE. NBVAR0) THEN
C                 VARIABLE UTILISATEUR
                  NETAT(LHEXPD) = 4
                  NV = NOVARU
               ELSE IF( NOVARU .EQ. NBVAR1 ) THEN
C                 NOM DE LA FONCTION
                  NETAT(LHEXPD) = 11
                  NV = NOVARU - NBVAR0
               ELSE IF( NOVARU.GT.NBVAR0 .AND. NOVARU.LT.NBVAR1 ) THEN
C                 PARAMETRE
                  NETAT(LHEXPD) = 10
                  NV = NOVARU - NBVAR0
               ELSE
C                 VARIABLE LOCALE
                  NETAT(LHEXPD) = 12
                  NV = NOVARU - NBVAR0
               ENDIF
               NOPER = NOPER + 1
               IF( NOPER .GT. MXNCOP ) GOTO 9200
               NCOPER(1,NOPER) = NETAT(LHEXPD)
               NCOPER(2,NOPER) = 0
               NCOPER(3,NOPER) = NV
               NCOPER(4,NOPER) = 0
C              PASSAGE EN FIN DE VARIABLE
               NL = NLF
               NC = NCF
               GOTO 1000
            ENDIF
         ENDIF
C
C        <VARIABLE TMS> ?
C        ----------------
         CALL CHVATS( NL , NC , NLF , NCF , NOVATS )
C        SI LA VARIABLE EST RETROUVEE, SA VALEUR EST DVATS(NOVATS)
         IF( NLF .GT. 0 ) THEN
            IF( NOVATS .GT. 0 ) THEN
C              VARIABLE RETROUVEE OU AJOUTEE
               NETAT(LHEXPD) = 5
               NOPER = NOPER + 1
               IF( NOPER .GT. MXNCOP ) GOTO 9200
               NCOPER(1,NOPER) = NETAT(LHEXPD)
               NCOPER(2,NOPER) = 0
               NCOPER(3,NOPER) = NOVATS
               NCOPER(4,NOPER) = 0
C              PASSAGE EN FIN DE VARIABLE
               NL = NLF
               NC = NCF
               GOTO 1000
            ENDIF
         ENDIF
C
C        <IDENT_FONCTION( ,  , ) ?
C        -------------------------
         CALL CHIDFO( NL , NC , NLF , NCF , NOFONC )
         IF( NLF .GT. 0 ) THEN
            IF( NOFONC .GT. 0 ) THEN
C              FONCTION RETROUVEE
               IF( NETAT(LHEXPD) .GE. 2 .AND. NETAT(LHEXPD) .LE. 6 )
     %         GOTO 9000
               NETAT(LHEXPD) = 6
               NOPER = NOPER + 1
               IF( NOPER .GT. MXNCOP ) GOTO 9200
               NCOPER(1,NOPER) = NETAT(LHEXPD)
               NCOPER(2,NOPER) = 0
               NCOPER(3,NOPER) = NOFONC
               NCOPER(4,NOPER) = 0
C              PASSAGE EN FIN D'IDENTIFICATEUR DE FONCTION
               NL = NLF
               NC = NCF
               GOTO 1000
            ENDIF
         ENDIF
C
C        ICI , CE N'EST PAS UNE VARIABLE
         IF( CAR .EQ. ';' ) GOTO 8000
      ENDIF
C
C     ERREUR
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='ALEXPD: MOT CLEF NON RECONNU DANS LA LIGNE SUIVANTE:'
      ELSE
         KERR(1) ='ALEXPD: UNKNOWN KEYWORD IN THE FOLLOWING LINE:'
      ENDIF
      IF( NL .LT. MXLGER ) THEN
         DO 799 KKK=1,NL
            KERR(1+KKK) = KLG(KKK)
 799     CONTINUE
         KKK = 1+NL
      ELSE
         KERR(2) = KLG(NL)
         KKK = 2
      ENDIF
      NBLGRC(NRERR) = KKK
      CALL LEREUR
      GOTO 9900
C
C     INITIALISATION COMMUNE
 800  NOPER = NOPER + 1
      IF( NOPER .GT. MXNCOP ) GOTO 9200
      NCOPER(1,NOPER) = NETAT(LHEXPD)
      NCOPER(2,NOPER) = 0
      NCOPER(3,NOPER) = 0
      NCOPER(4,NOPER) = 0
C
C     PASSAGE AU CARACTERE SUIVANT
1000  NC = NC + 1
      IF( NC .GT. NBCALI ) THEN
         NL = NL + 1
         IF( NL .GT. NLMAX ) THEN
C           ERREUR
            GOTO 9000
         ENDIF
      ENDIF
      IF( LHEXPD .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'ALEXPD: LHEXPD=',LHEXPD,' NC=',NC,' NL=',NL
         WRITE(IMPRIM,*) 'ALEXPD: KLG(NL)=',KLG(NL)
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: PILE DES EXPRESSIONS AU NIVEAU 0'
         ELSE
            KERR(1) = 'LU: EXPRESSION STACK AT LEVEL 0'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     TRAITEMENT DE L'ITEM SUIVANT
      GOTO 100
C
C     FIN DE L'EXPRESSION
C     ===================
 8000 IF( NBPAOU( LHEXPD ) .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: '')'' OUBLIEE'
         ELSE
            KERR(1) = 'LU: FORGOTTEN '')'' '
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF
C     LA DERNIERE OPERATION DE L'EXPRESSION
      LSEXPD(2,1) = NOPER
C
      IF( LUIMPR .GE. 10 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'LU: SORTIE DE ALEXPD avec CODES OPERATIONS'
         ELSE
            WRITE(IMPRIM,*) 'LU: EXIT of ALEXPD with OPERATION CODES'
         ENDIF
         WRITE(IMPRIM,18000) (J,(NCOPER(I,J),I=1,4),J=1,LSEXPD(2,1))
      ENDIF
18000 FORMAT(' NCOPER(',I3,')=',4I4)
      RETURN
C
C     ERREUR
C     ======
 9000 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: SYNTAXE INCORRECTE DANS'
      ELSE
         KERR(1) = 'LU: INCORRECT SYNTAX IN'
      ENDIF
      NC = NCD
      IF( NLD .EQ. NLMAX ) GOTO 9030
      KERR(2) = KLG(NLD)(NCD:NBCALI)
      DO 9010 NL=NLD+1,NLMAX-1
         NBLGRC(NRERR) = 1
         KERR(2+NL-NLD-1) = KLG(NL)
 9010 CONTINUE
      NC = 1
 9030 NBLGRC(NRERR) = NLMAX - NLD + 1
      KERR(NBLGRC(NRERR)) = KLG(NLMAX)(NC:NCMAX)
      CALL LEREUR
      GOTO 9900
C
 9100 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
      KERR(1) =
     %  'LU: TABLEAU DES EXPRESSIONS SATUREE.AUGMENTER MXEXPD'
      ELSE
      KERR(1)='LU: SATURATED ARRAY OF EXPRESSIONS. AUGMENT MXEXPD'
      ENDIF
      CALL LEREUR
      GOTO 9000
C
C     TABLE DES OPERATEURS SATUREE
 9200 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: TABLE DES OPERATIONS SATUREE'
         KERR(2) = 'LU: AUGMENTER MXNCOP DANS incl/lu.inc'
      ELSE
         KERR(1) = 'LU: THE ARRAY OF OPERATORS IS SATURATED'
         KERR(2) = 'LU: AUGMENT MXNCOP IN incl/lu.inc'
      ENDIF
      CALL LEREUR
C
9900  NRETOU = 1
      NCODEV = 0
      RETURN
      END
