      SUBROUTINE AFTSTD( KNOMTS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LES VALEURS D'UN TABLEAU TMS A PARTIR DE SON
C ----- TABLEAU DESCRIPTEUR
C
C ENTREE :
C --------
C KNOMTS : NOM DU TABLEAU MS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C23456---------------------------------------------------------------012
      PARAMETER (NF=91)
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/td.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON /UNITES/LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NOAFTS,NUNIT(26)
      CHARACTER*160     KNOM
      CHARACTER*4       KTYPE
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOMTS
C
C     INITIALISATIONS
C     LA DERNIERE LIGNE DU TABLEAU KTD OCCUPEE EST SAUVEGARDEE
      LHKTD0 = LHKTD
C     AFFICHAGE FORCE
      NOAFTS = 1
C     LA PILE DES OPERATIONS
      LHPILE = 1
C     LA PILE DES TMS
      LHTMS = 1
C
C     CODE OPERATION DE RECHERCHE D'UN TMS
      NOMTS(LHTMS)= KNOMTS
      LAPILE(0,1) = LHTMS
      LAPILE(1,1) = 1
      LAPILE(2,1) = 0
      LAPILE(3,1) = 0
      LAPILE(4,1) = 0
      LAPILE(5,1) = 0
C     IMPRESSION DEMANDEE
      LAPILE(6,1) = 1
C
C     OUVERTURE DU TABLEAU MS
      CALL TNOUVR( KNOMTS , KTYPE , NOTAMS(LHTMS) , MNTAMS(LHTMS) )
      IF( NOTAMS(LHTMS) .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NON EXISTENCE du TMS: ' // KNOMTS
         ELSE
            KERR(1) = 'UNKNOWN TMS: ' // KNOMTS
         ENDIF
         CALL LEREUR
         GOTO 9999
      ENDIF
      IF( KTYPE(1:4) .EQ. 'LEXI' ) THEN
C        LE TABLEAU EST UN LEXIQUE
         CALL LXIM( NOTAMS(LHTMS) )
         GOTO 9999
      ENDIF
C
C     ADRESSE DANS MCN DU TMS
      MCTAMS(LHTMS) = MNTAMS(LHTMS)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMTS
      ELSE
         WRITE(IMPRIM,20000) KNOMTS
      ENDIF
10000 FORMAT(/' CONTENU du TMS ',A)
20000 FORMAT(/' CONTENTS of TMS ',A)
C
C     LE NOMBRE DE VARIABLE DU TMS A RECHERCHER
      NBVATA(LHTMS) = 1
C     LE NUMERO DANS DICOTD DU NOM DU TABLEAU DESCRIPTEUR
C     C'EST LE MOT DERRIERE LA DATE DANS LE TABLEAU TMS
      NOTADS(LHTMS) = MCN( MNTAMS(LHTMS) + MOTVAR(6) )
      IF( NOTADS(LHTMS) .LE. 0 .OR.
     %    NOTADS(LHTMS) .GT. NBTD ) THEN
          NBLGRC(NRERR) = 1
          WRITE(KERR(2)(1:15),'(I15)') NOTADS(LHTMS)
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) =
     %      'AFTSTD: NUMERO INCORRECT DE TABLEAU DESCRIPTEUR'
     %       //KERR(2)(1:15)
          ELSE
             KERR(1) =
     %      'AFTSTD: INCORRECT NUMBER of DESCRIPTOR ARRAY'
     %       //KERR(2)(1:15)
          ENDIF
          CALL LEREUR
          GOTO 9999
      ENDIF
C     LE DECALAGE DANS LE TS POUR ATTEINDRE LA VALEUR NUMERIQUE TRAITEE
      LDTS(LHTMS) = 0
C
C     LE PARAMETRE D'IMPRESSION
      IMPRES = 1
      NOMUET = 1
C     LE POINTEUR SUR LE DERNIER CARACTERE STOCKE DANS LA LIGNE BUFFER
      LCLIGN = 0
C     LA LIGNE BUFFER EST BLANCHE
      KLIGNE = ' '
C     LE NOMBRE D'IDENTIFICATEURS STOCKES
      NBIDEN(LHTMS) = 0
C     LA HAUTEUR DU TAS
      LHTAS  = 0
C
C     ==================================================================
C     TRAITEMENT D'UN TMS
C     ==================================================================
C     TANT QUE LA PILE EST NON VIDE FAIRE
 100  CALL AFLIGN
      IF( LHPILE .GT. 0 ) THEN
C
C        LE NUMERO LOCAL DU TMS DANS LES TABLEAUX NBVATA , ...
         LHTMS = LAPILE(0,LHPILE)
C
         IF( LAPILE(1,LHPILE) .EQ. 2 ) THEN
C
C           TRAITEMENT D'UN CAS
            GOTO 200
C
         ELSE IF( LAPILE(1,LHPILE) .EQ. 1 ) THEN
C
C           TRAITEMENT D'UN TMS
            IF( NBVATA(LHTMS) .LE. 0 ) THEN
C              ON DEPILE LE TMS EPUISE
C              ON REND LA PLACE DANS KTD DU TMS
               LHKTD  = LAPILE(2,LHPILE) - 1
C              LE CARACTERE DE RETOUR
               NL     = LAPILE(4,LHPILE)
               NC     = LAPILE(5,LHPILE)
               LHPILE = LHPILE - 1
               GOTO 100
            ELSE
C              UN TMS DE PLUS EST TRAITE
               NBVATA(LHTMS) = NBVATA(LHTMS) - 1
            ENDIF
C
C           deftms  NOM_LEXIQUE [CHAINE] TYPE  fintms
C           =========================================
C
C           OUVERTURE DU TABLEAU MS
C           -----------------------
C           RECHERCHE DU TABLEAU DESCRIPTEUR DU TS PAR FUSION DU
C           NOM DU TD ET DU TS
            CALL FUSNOM( DICOTD(NOTADS(LHTMS)) , KNOMTS , NOMTS(LHTMS)
     %                 , NB )
            IF( NB .EQ. 0 ) THEN
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = ' AFTSTD: FUSION DES NOMS ' //
     %            DICOTD(NOTADS(LHTMS))
               ELSE
                  KERR(1) = ' AFTSTD: FUSION of NAMES ' //
     %            DICOTD(NOTADS(LHTMS))
               ENDIF
               KERR(2) = KNOMTS // ' IMPOSSIBLE'
               CALL LEREUR
               GOTO 9000
            ENDIF
C
            IF( NOTAMS(LHTMS) .GT. 0 ) THEN
C              OUVERTURE
               CALL TAMSOU( NOTAMS(LHTMS) , MNTAMS(LHTMS) )
               IF( MNTAMS(LHTMS) .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'AFTSTD: TMS NON OUVRABLE'
                  ELSE
                     KERR(1) = 'AFTSTD: NOT OPENABLE TMS'
                  ENDIF
                  KERR(2) =  NOMTS(LHTMS)
                  CALL LEREUR
                  GOTO 9999
               ENDIF
               MCTAMS(LHTMS) = MNTAMS(LHTMS)
C
C              LE TABLEAU EXISTE .  EST-IL UN NOM_LEXIQUE ?
               IF( KTYPE .EQ. 'LEXI' ) THEN
C                 AFFICHAGE DU NOM_LEXIQUE
                  CALL LXIM0( MNTAMS(LHTMS) )
                  GOTO 9000
               ENDIF
            ENDIF
C
C           NON . C'EST UN TMS
C
C           SI LE TABLEAU DESCRIPTEUR N'EST PAS DANS KTD
C           ALORS IL EST CHARGE A PARTIR DE LA LIGNE LHKTD+1
C           ------------------------------------------------
            IF( LAPILE(2,LHPILE) .LE. 0 ) THEN
C
C              LE TABLEAU TD N'EST PAS DANS KTD.IL EST CHARGE
C              FORMATION DU NOM DU FICHIER SUPPORT DU TD
               CALL NTDFIC( DICOTD(NOTADS(LHTMS)) , KNOM )
C              OUVERTURE DU FICHIER SUPPORT DU TD
               OPEN( FILE=KNOM , UNIT=NF , IOSTAT=I )
               IF( I .NE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='AFTSTD: PB pour OUVRIR le FICHIER'
                  ELSE
                     KERR(1)='AFTSTD: PROBLEM to OPEN the FILE'
                  ENDIF
                  NK = NUDCNB( KNOM )
                  KERR(2) = KNOM(1:NK)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C              LE NUMERO DE LIGNE DE DEBUT DU TD
               NL = LHKTD + 1
               NC = 0
C
 110           LHKTD = LHKTD + 1
               IF( LHKTD .GT. MXLITD ) THEN
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(2)(1:7),'(I7)') MXLITD
                  IF( LANGAG .EQ. 0 ) THEN
                  KERR(1)='AFTSTD: TABLEAU KTD SATURE' // KERR(2)(1:7)
                  ELSE
                  KERR(1)='AFTSTD: SATURATED ARRAY KTD. AUGMENT MXLITD'
                  ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              LECTURE DE LA LIGNE
               READ(NF,'(A)',END=120,IOSTAT=I) KTD(LHKTD)
               IF( I .EQ. 0 ) GOTO 110
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1)='AFTSTD: ERREUR en LECTURE du FICHIER'
               ELSE
                  KERR(1) ='AFTSTD: READING ERROR on FILE'
               ENDIF
               NK = NUDCNB( KNOM )
               KERR(2) = KNOM(1:NK)
               CALL LEREUR
               GOTO 9000
C
C              FIN DE LECTURE DU FICHIER
 120           LHKTD = LHKTD - 1
C              CARACTERE DE DEBUT DU TMS
               LAPILE(2,LHPILE) = NL
               LAPILE(3,LHPILE) = NC
C              LE CARACTERE DE RETOUR EST DEJA INITIALISE
C              IMPRESSION DEMANDEE
               LAPILE(6,LHPILE) = 1
C              FERMETURE DU FICHIER SUPPORT DU TD
               CLOSE( UNIT=NF )
C
            ENDIF
C
C           LE PREMIER CARACTERE DU TMS
            NL = LAPILE(2,LHPILE)
            NC = LAPILE(3,LHPILE)
C
C           deftms NOM_LEXIQUE [CHAINE] TYPE fintms
C           ===========================================
            CALL CHMOT( 'deftms' , NL,NC , NLD,NCD, NLF,NCF )
C           EN SORTIE NLD,NCD ET NLF,NCF DEBUT ET FIN DU MOT DANS KTD
C
C           RECHERCHE DU NOM_LEXIQUE DU TMS ET AFFICHAGE
            NL = NLF
            NC = NCF
            CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLD .LE. 0 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'AFTSTD: NOM_LEXIQUE INCORRECT '
               ELSE
                  KERR(1) = 'AFTSTD: INCORRECT NAME_LEXIQUE '
               ENDIF
               KERR(2) = KTD(NL)
               CALL LEREUR
               GOTO 9000
            ENDIF
C
C           AFFICHAGE DU TABLEAU TMS
            CALL AFCAR( 'tms :' )
            CALL AFNOM( NLD,NCD, NLF,NCF )
C
C           RECHERCHE DE [CHAINE] ET AFFICHAGE
            NL = NLF
            NC = NCF
            CALL CHCHAI( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLD .GT. 0 ) THEN
               CALL AFCHAI( NLD,NCD, NLF,NCF )
               NL = NLF
               NC = NCF
            ENDIF
            CALL AFLIGN
C
C           TRAITEMENT D'UN TYPE
c           --------------------
C           RECHERCHE DES MOTS CLES variable tableau cas
C                                   fintms fincas
 200        IF( LAPILE(1,LHPILE) .EQ. 2 ) THEN
               IMPRES = LAPILE(6,LHPILE)
            ENDIF
            CALL CHNOM( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLD .LE. 0 ) THEN
C              LE 1-ER NOM NE COMMENCE PAS PAR UNE LETTRE
C              ES CE UN <INTERVALLE_E> ?
               GOTO 500
            ENDIF
C
            IF( KTD(NLD)(NCD:NCF) .EQ. 'variable' ) THEN
C
C              variable [muet] IDENT [CHAINE] TYPE_VAR ;
C              -----------------------------------------
C              RECHERCHE DE [muet]
               NL = NLF
               NC = NCF
               CALL CHNOM( NL,NC, NLD,NCD, NLF,NCF )
               IF( KTD(NLD)(NCD:NCF) .EQ. 'muet' ) THEN
C                 LA VARIABLE EST MUETTE
c                 LA VARIABLE NOMUET A MEME VOCATION QUE IMPRES
                  NOMUET = 0
                  NL     = NLF
                  NC     = NCF
               ELSE
C                 LA VARIABLE N'EST PAS MUETTE
                  NOMUET = 1
               ENDIF
C
C              RECHERCHE DE IDENT
               CALL CHIDEN( NL,NC, NLD,NCD, NLF,NCF , NOIDEN )
C              NOIDEN NUMERO DE L'IDENTICATEUR, 0 SI NON RETROUVE
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'AFTSTD: IDENTIFICATEUR INCORRECT'
                  ELSE
                     KERR(1) = 'AFTSTD: INCORRECT IDENTIFICATOR'
                  ENDIF
                  KERR(2) =  KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOMUET .GT. 0 .AND. NOIDEN .GT. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                    KERR(1) = 'AFTSTD: '//KIDENT(NOIDEN)//' EXISTE DEJA'
                  ELSE
                   KERR(1)='AFTSTD: '//KIDENT(NOIDEN)//' ALREADY EXISTS'
                  ENDIF
                  CALL LEREUR
               ENDIF
C              AFFICHAGE DE L'IDENT
               CALL AFNOM( NLD,NCD, NLF,NCF )
C              L'IDENTIFICATEUR NON RETROUVE EST AJOUTE
               IF( NBIDEN(LHTMS) .GE. MXIDEN ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='AFTSTD: TABLE des IDENTIFICATEURS SATUREE'
                  ELSE
                     KERR(1)='AFTSTD: SATURATED TABLE of IDENTIFICATORS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               NBIDEN(LHTMS) = NBIDEN(LHTMS) + 1
               KIDENT( NBIDEN(LHTMS) ) = KTD(NLD)(NCD:NCF)
               IDENT(0,NBIDEN(LHTMS) ) = LHPILE
               IDENT(2,NBIDEN(LHTMS) ) = 0
               IDENT(3,NBIDEN(LHTMS) ) = LDTS(LHTMS)
               IDENT(4,NBIDEN(LHTMS) ) = NLD
               IDENT(5,NBIDEN(LHTMS) ) = NCD
C
C              RECHERCHE DE [CHAINE] ET AFFICHAGE
               NL = NLF
               NC = NCF
               CALL CHCHAI( NL,NC, NLD,NCD, NLF,NCF )
               IF( NLD .GT. 0 ) THEN
                  CALL AFCHAI( NLD,NCD, NLF,NCF )
                  NL = NLF
                  NC = NCF
               ENDIF
C
C              RECHERCHE DE TYPE_VAR POUR UNE VARIABLE
               NBV = 1
C              LA VARIABLE EST TRAITEE COMME UN TABLEAU D'UNE VARIABLE
               GOTO 275
C
            ELSE IF( KTD(NLD)(NCD:NCF) .EQ. 'tableau' ) THEN
C
C              tableau [muet] IDENT ( INTERVAL_ENT { , INTERVAL_ENT } )
C                      [CHAINE] TYPE_VAR ;
C              --------------------------------------------------------
C              RECHERCHE DE [muet]
               NL = NLF
               NC = NCF
               CALL CHNOM( NL,NC, NLD,NCD, NLF,NCF )
               IF( KTD(NLD)(NCD:NCF) .EQ. 'muet' ) THEN
C                 LA VARIABLE EST MUETTE
                  NOMUET = 0
                  NL     = NLF
                  NC     = NCF
               ELSE
C                 LA VARIABLE N'EST PAS MUETTE
                  NOMUET = 1
               ENDIF
C
C              RECHERCHE DE IDENT
               CALL CHIDEN( NL,NC, NLD,NCD, NLF,NCF , NOIDEN )
C              NOIDEN NUMERO DE L'IDENTICATEUR, 0 SI NON RETROUVE
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'AFTSTD: IDENTIFICATEUR INCORRECT'
                  ELSE
                     KERR(1) = 'AFTSTD: INCORRECT IDENTIFICATOR'
                  ENDIF
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOIDEN .GT. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                    KERR(1) = 'AFTSTD: '//KIDENT(NOIDEN)//' EXISTE DEJA'
                  ELSE
                   KERR(1)='AFTSTD: '//KIDENT(NOIDEN)//' ALREADY EXISTS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C              AFFICHAGE DE L'IDENT
               CALL AFNOM( NLD,NCD, NLF,NCF )
C              L'IDENTIFICATEUR NON RETROUVE EST AJOUTE
               IF( NBIDEN(LHTMS) .GE. MXIDEN ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='AFTSTD: TABLE DES IDENTIFICATEURS SATUREE'
                  ELSE
                     KERR(1)='AFTSTD: SATURATED TABLE of IDENTIFICATORS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               NBIDEN(LHTMS) = NBIDEN(LHTMS) + 1
               KIDENT( NBIDEN(LHTMS) ) = KTD(NLD)(NCD:NCF)
               IDENT(0,NBIDEN(LHTMS) ) = LHPILE
               IDENT(3,NBIDEN(LHTMS) ) = LDTS(LHTMS)
               IDENT(4,NBIDEN(LHTMS) ) = NLD
               IDENT(5,NBIDEN(LHTMS) ) = NCD
C
C              RECHERCHE DE ( INTERVAL_ENT {, INTERVAL_ENT } )
C              RECHERCHE DES ( )
               NL = NLF
               NC = NCF
               CALL CHPAOF( NL,NC , NLD,NCD , NLF,NCF )
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = KTD(NL)
                  KERR(2) = 'INCORRECT ( )'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              RECHERCHE DU NOMBRE D'INDICES DU TABLEAU
               NL = NLD
               NC = NCD
               CALL CARAVA( NL , NC )
C              NB LE NOMBRE D'INDICES DU TABLEAU
               CALL CHINTA( NL,NC , NLD,NCD , NLF,NCF , NB )
C              STOCKAGE DE NB , (VAL_MIN,VAL_MAX) DANS LE TAS
C                       DE CHAQUE INDICE
               IF( LHTAS + 1 + NB + NB .GT. MXTAS ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'AFTSTD: TAS SATURE'
                  ELSE
                     KERR(1) = 'AFTSTD: SATURATED HEAP'
                  ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              MISE A JOUR DES INDICES DU TABLEAU
               LHTAS = LHTAS + 1
               LETAS(  LHTAS  ) = NB
               IDENT(2,NBIDEN(LHTMS) ) = LHTAS
               CALL AFCAR( '(' )
C              NBV LE NOMBRE DE VARIABLES DU TABLEAU
               NBV = 1
               NL  = NLD
               NC  = NCD
               DO 250 I=1,NB
C                 RECHERCHE DE INTERVAL_ENT ET MISE A JOUR DANS LE TAS
                  CALL CHINEN(  NL,NC , NLD,NCD , NLF,NCF , N1 , N2 )
                  LHTAS = LHTAS + 2
                  LETAS(LHTAS-1 ) = N1
                  LETAS(LHTAS)    = N2
                  NBV = NBV * ( N2 - N1 + 1 )
C                 AFFICHAGE DE L'INTERVAL_ENT
                  CALL AFENTI( N1 )
                  CALL AFCAR( ' .. ' )
                  CALL AFENTI( N2 )
                  IF( I .NE. NB ) THEN
                     CALL AFCAR( ',' )
                  ELSE
                     CALL AFCAR( ')' )
                  ENDIF
C                 PASSAGE EN FIN INTERVAL_ENT
                  NL = NLF
                  NC = NCF
C                 SAUT DE , OU )
                  CALL CAR1NB( NL , NC )
 250           CONTINUE
C
C              RECHERCHE DE [CHAINE] AU DELA DE ) ET AFFICHAGE
               CALL CHCHAI( NL,NC , NLD,NCD , NLF,NCF )
               IF( NLD .GT. 0 ) THEN
                  CALL AFCHAI( NLD,NCD , NLF,NCF )
                  NL = NLF
                  NC = NCF
               ENDIF
C
C              TRAITEMENT DE TYPE_VAR :=
C             |entier {(<INTERVAL_ENT>:<CHAINE> {,<INTERVAL_ENT> : <CHAINE>})}
C             |^<NOM_LEXIQUE>    $ pointe sur un numero de nom dans ce lexique
C             |caractere
C             |reel
C             |reel2
C             |xyz
C             |tms <NOM_TMS>
C
  275          CALL CHTYPV( NL,NC , NLD,NCD , NLF,NCF , NOTYPE )
C              EN RETOUR SELON LE TYPE DE LA VARIABLE
C              LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C              REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C              COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C              TYPEOBJET=>13
C              TMS      =>21
C
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'INCORRECT TYPE_VAR'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
               IDENT(1,NBIDEN(LHTMS)) = NOTYPE
               NL = NLF
               NC = NCF
C
               IF( NOTYPE .GE. 1 .AND. NOTYPE .LE. NBTYPV ) THEN
C
C                 TYPE_VAR:= |entier {(<INTERVAL_ENT>:<CHAINE>
C                                    {,<INTERVAL_ENT> : <CHAINE>})}
C                            |^<NOM_LEXIQUE>
C                               $ pointe sur un numero de nom dans ce lexique
C                            |caractere
C                            |reel
C                            |reel2
C                            |xyz
C                 ...........................................................
C                 AFFICHAGE DU TYPE
                  CALL AFTYVA( NBIDEN(LHTMS) , NRETOU )
                  IF( NRETOU .GT. 0 ) GOTO 9000
C
               ELSE IF( NOTYPE .EQ. 21 ) THEN
C
C                 TYPE_VAR = tms NOM_TMS
C                 ......................
                  CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
                  IF( NLD .LE. 0 ) THEN
                     NBLGRC(NRERR) = 2
                     KERR(1) =  KTD(NL)
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(2) = 'NOM_LEXIQUE INCONNU'
                     ELSE
                        KERR(2) = 'UNKNOWN NAME_LEXIQUE'
                     ENDIF
                     CALL LEREUR
                     GOTO 9000
                  ENDIF
C                 AFFICHAGE
                  CALL AFCAR( ' tms' )
                  CALL AFNOM( NLD,NCD, NLF,NCF )
                  CALL AFCAR( '=' )
C
C                 NOMBRE DE VARIABLES ENTIERES (NO DU TMS )
                  CALL NVARID( NBIDEN(LHTMS) , NBV , NOTYPE )
C
C                 AFFICHAGE DE CES ENTIERS
                  DO 280 I=1,NBV
                     CALL AFENTI( MCN( MNTAMS(LHTMS) + LDTS(LHTMS) ) )
C                    DECALAGE D'UN ENTIER
                     LDTS(LHTMS) = LDTS(LHTMS) + 1
 280              CONTINUE
C
C                 PASSAGE AU ; SUIVANT DERRIERE LE NOM DU TMS
                  NL = NLF
                  NC = NCF
C
               ELSE
                  NBLGRC(NRERR) = 1
                  KERR(1) = ' AFTSTD: TYPE INCORRECT'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              FIN <TYPE_VAR> PASSAGE AU CARACTERE ;
               CALL CHCAR( ';' , NL , NC )
               CALL AFLIGN
               GOTO 200
C
            ELSE IF( KTD(NLD)(NCD:NCF) .EQ. 'fintms' ) THEN
C
               IF( LHTMS .GT. 1 ) THEN
C                 DECALAGE D'UN ENTIER DANS CE TABLEAU
                  LDTS(LHTMS-1) = LDTS(LHTMS-1) + 1
               ENDIF
               GOTO 100
C
            ELSE IF( KTD(NLD)(NCD:NCF) .EQ. 'cas' ) THEN
C              cas IDENT { INTERVALLE_E : TYPE } fincas
C              ----------------------------------------
C              LE CAS DOIT ETRE EMPILE
               IF( LHPILE .GE. MXPILE ) THEN
                  NBLGRC(NRERR) = 1
                   IF( LANGAG .EQ. 0 ) THEN
                      KERR(1) = 'AFTSTD: PILE SATUREE'
                   ELSE
                      KERR(1) = 'AFTSTD: SATURATED STACK'
                   ENDIF
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               LHPILE = LHPILE + 1
               LAPILE(0,LHPILE) = LHTMS
               LAPILE(1,LHPILE) = 2
C              CARACTERE DE IDENT DE cas
               LAPILE(2,LHPILE) = NLD
               LAPILE(3,LHPILE) = NCD
C              LE PARAMETRE D'IMPRESSION
               LAPILE(6,LHPILE) = IMPRES * NOMUET
C
C              RECHERCHE DE IDENT
               NL = NLF
               NC = NCF
               CALL CHIDEN( NL,NC, NLD,NCD, NLF,NCF , NOIDEN )
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = KTD(NL)
                  KERR(2) = 'INCORRECT  INTERVALLE_E'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOIDEN .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                   KERR(1)='AFTSTD: IDENTIFICATEUR '// KTD(NLD)(NCD:NCF)
     %                      //' INCONNU'
                  ELSE
                    KERR(1)='AFTSTD: IDENTIFICATOR '// KTD(NLD)(NCD:NCF)
     %                     //' UNKNOWN'
                 ENDIF
                 CALL LEREUR
                 GOTO 9000
               ENDIF
C
C              LE NUMERO DE L'IDENTIFICATEUR DU CAS
               LAPILE(4,LHPILE) = NOIDEN
C              LE NOMBRE DE CAS RENCONTRES
               LAPILE(5,LHPILE) = 1
C
C              RECHERCHE DE CETTE VALEUR PARMI LES INTERVALLE_E
               NL    = NLF
               NC    = NCF
C
C              TRAITEMENT D'UN CAS
 500           NOIDEN = LAPILE(4,LHPILE)
               IMPRES = LAPILE(6,LHPILE)
C              LA VALEUR ENTIERE DE L'IDENTIFICATEUR
               CALL VAIDEN( NOIDEN , NVALEN )
C              LA VALEUR DES BORNES DE L'INTERVALLE_E TRAITE
               CALL CHINTR( NL,NC, NLD,NCD, NLF,NCF , N1 , N2 )
C
C              LA VALEUR DU CAS EST-ELLE DANS L'INTERVALLE_E ?
               IF( NVALEN .LT. N1 .OR. NVALEN .GT. N2 ) THEN
C
C                 LA VALEUR N'EST PAS DANS L'INTERVALLE N1 N2
C                 IMPRESSION ALIENEE
                  IMPRES = 0
C                 RECHERCHE DE L'INTERVALLE SUIVANT EN SAUTANT
C                 LES EVENTUELS CAS INTERMEDIAIRES
C                 NLF NCF POINTE SUR LE DERNIER CARACTERE DE INTERVALLE_E
C
C                 cas <IDENT> {<INTERVALLE_E> : <TYPE>} fincas ;
C
C                 <TYPE> SE TERMINE PAR ; . RECHERCHE DU PREMIER ;
                  NLMAX = NLF
                  NCMAX = NCF
C
C                 RECHERCHE DE ;
 510              CALL CHMOT( ';' , NLMAX,NCMAX , NLD,NCD,
     %                              NLMAX,NCMAX )
C                 NLMAX NCMAX POSITION DU PROCHAIN ;
C
C                 EXISTE-T-IL UN OU PLUSIEURS cas AVANT LE ; ?
 520              CALL CHMOSC( 'cas' , NLMAX,NCMAX, NL,NC,
     %                         NLD,NCD, NLF,NCF )
C                 NLD,NCD DEBUT DE cas , NLF,NCF FIN DE cas
                  IF( NLD .GT. 0 ) THEN
C                    IL EXISTE cas . ES CE fincas ?
                     IF( NCD .LE. 3 ) GOTO 530
                     IF( KTD(NLD)(NCD-3:NCD-1) .EQ. 'fin' ) THEN
C                       UN fincas DE PLUS. SUIVI DE ;
                        LAPILE(5,LHPILE) = LAPILE(5,LHPILE) - 1
C
                        IF( LAPILE(5,LHPILE) .LE. 0 ) THEN
C                          LA FIN DU CAS INITIAL EST ATTEINTE.
                           IMPRES = LAPILE(6,LHPILE)
                           LHPILE = LHPILE - 1
C                          DEPLACEMENT AU ; QUI SUIT fincas
                           NL = NLMAX
                           NC = NCMAX
                           GOTO 200
C
                        ELSE IF( LAPILE(5,LHPILE) .EQ. 1 ) THEN
C                          RETOUR A L'INTERIEUR DU CAS INITIAL
C                          IL PEUT ETRE SUIVI DE
C                          <INTERVALLE_E> OU <TYPE> . PASSAGE AU ;
                           NL1 = NLMAX
                           NC1 = NCMAX
                           CALL CHINTR( NL1,NC1, NLD1,NCD1, NLF1,NCF1,
     %                                  N1,N2)
                           IF( NLD1 .GT. 0 ) THEN
C                             IL S'AGIT D'UN <INTERVALLE_E>
                              NL = NLMAX
                              NC = NCMAX
                              GOTO 500
                           ELSE
C                             IL S'AGIT D'UN T<YPE>
                              NL = NLMAX
                              NC = NCMAX
C                             RECHERCHE DU ; SUIVANT
                              GOTO 510
                           ENDIF
C
                        ELSE
C                          CAS DE CAS . PASSAGE AU ;
                           NL = NLMAX
                           NC = NCMAX
C                          RECHERCHE DU ; SUIVANT
                           GOTO 510
C
                        ENDIF
                     ENDIF
C
C                    cas ET NON PAS fincas
 530                 LAPILE(5,LHPILE) = LAPILE(5,LHPILE) + 1
C
C                    IL PEUT EXISTER UN AUTRE CAS AVANT LE PROCHAIN;
                     NL = NLF
                     NC = NCF
C                    NL,NC DERNIER CARACTERE DU CAS AJOUTE
                     GOTO 520
                  ENDIF
C
C                 IL N'EXISTE PLUS DE cas OU fincas ENTRE NL,NC ET
C                 NLMAX,NCMAX DU CARACTERE ;
                  NL = NLMAX
                  NC = NCMAX
C
                  IF( LAPILE(5,LHPILE) .GT. 1 ) GOTO 510
C
C                 SUIVENT <INTERVALLE_E> OU <TYPE>
                  CALL CHINTR( NL,NC, NLD1,NCD1, NLF1,NCF1, N1,N2 )
                  IF( NLD1 .GT. 0 ) THEN
C                    <INTERVALLE_E>
                     GOTO 500
                  ELSE
C                    <TYPE>
                     GOTO 510
                  ENDIF
C
               ELSE
C
C                 LA VALEUR DU CAS EST DANS LE BON INTERVALLE
C                 IMPRESSION REGENEREE
                  IMPRES = LAPILE(6,LHPILE)
C                 PASSAGE AU :
                  NL = NLF
                  NC = NCF
                  CALL CHCAR( ':' , NL , NC )
                  GOTO 200
C
               ENDIF
C
            ELSE IF( KTD(NLD)(NCD:NCF) .EQ. 'fincas' ) THEN
C
C              LE CAS EST DEPILE
               LHPILE = LHPILE - 1
               NL     = NLD
               NC     = NCF
               CALL CHCAR( ';' , NL , NC )
               IF( NL .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'AFTSTD: ; INCORRECT '
                  KERR(2) =  KTD(NLD)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               GOTO 200
C
            ELSE
C
C              ERREUR
               NBLGRC(NRERR) = 2
               KERR(1) = 'AFTSTD: TYPE INCORRECT '
               KERR(2) = KTD(NLD)
               CALL LEREUR
               GOTO 9000
            ENDIF
         ENDIF
      ENDIF
C
C     ON VIDE LE BUFFER LIGNE
      CALL AFLIGN
      NRETOU = 0
      GOTO 9999
C
C     ERREUR
C     ======
 9000 NRETOU = 1
C
C     AFFICHAGE DES TMS DESAMORCE
 9999 NOAFTS = 0
      LHKTD  = LHKTD0
      END
