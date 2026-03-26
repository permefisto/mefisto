      SUBROUTINE VATSTD( KNOMVA , NOTYPE , MNTMS , LDTMS , NCFVTS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETROUVER L'ADRESSE MCN DU TABLEAU TMS ET DU DECALAGE POUR
C ----- ATTEINDRE UNE VARIABLE D'UN TMS GRACE AU TABLEAU DESCRIPTEUR
C
C ENTREE :
C --------
C KNOMVA : NOM DE LA VARIABLE   :=  NOM TMS( VARIABLE INDICEE OU NON )
C
C          EXEMPLES : ~>POINT>P1>DEFINITION(NOTYPO)      VARIABLE
C                     ~>POINT>P1>XYZSOMMET(XYZSOM(2))    TABLEAU
C                     ~>POINT>P1>XYZSOMMET(XYZSOM(2)(3)) TABLEAU DE XYZ
C          POUR UN TABLEAU L'ABSENCE DE () RETOURNE
C          L'ADRESSE DU DEBUT DU TABLEAU
C
C SORTIES :
C ---------
C NOTYPE  : NUMERO DU TYPE DE LA VARIABLE
C           LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C           REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C           COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C           TMS      =>21
C MNTMS   : ADRESSE MCN DU TABLEAU TMS CONTENANT LA VARIABLE
C           0 SI LA VARIABLE N'EXISTE PAS
C LDTMS   : DECALAGE DANS LE TMS POUR ATTEINDRE LE PREMIER MOT DE
C           LA VARIABLE
C NCFVTS  : POSITION DU DERNIER CARACTERE DANS KNOMVA DU NOM DE LA
C           VARIABLE DU TMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE  UPMC  PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
      PARAMETER (NF=91)
C.......................................................................
      include"./incl/td.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
      include"./incl/impres.inc"
C.......................................................................
      CHARACTER*160     KNOM, KNOMTS
      CHARACTER*8       KIDENV
      CHARACTER*4       KTYPE
      CHARACTER*9       TYPNUM
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOMVA
C
C     INITIALISATIONS
      NOTYPE = 0
      MNTMS  = 0
      LDTMS  = 0
      NCFVTS = 0
C     LA DERNIERE LIGNE DU TABLEAU KTD OCCUPEE EST SAUVEGARDEE
      LHKTD0 = LHKTD
C     LA PILE DES OPERATIONS
      LHPILE = 1
C     LA PILE DES TMS
      LHTMS0 = LHTMS
      IF( LHTMS .GE. MXTMS ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VATSTD: PILE SATUREE DES TMS'
         CALL LEREUR
         GOTO 9000
      ENDIF
      LHTMS  = LHTMS + 1
C
C     LE NOMBRE DE CARACTERES DE KNOMVA
      L = LEN( KNOMVA )
C
C     RECHERCHE DU NOM DU TMS ASSOCIE A LA VARIABLE
      I = INDEX( KNOMVA , '(' )
      IF( I .GT. 0 ) THEN
C
C        LE NOM DU TMS
         KNOMTS = KNOMVA(1:I-1)
C
C        RECHERCHE VARIABLE OU TABLEAU ?
         I1 = INDEX( KNOMVA(I+1:L) , '(' )
         IF( I1 .LE. 0 ) THEN
C           IL N'EXISTE PAS DE ( (  => VARIABLE SIMPLE
            I1 = INDEX( KNOMVA(I+1:L) , ')' )
            IF( I1 .LE. 0 ) THEN
               IF( IMNMTS .NE. 0 ) THEN
                   NBLGRC(NRERR) = 1
                   KERR(1) = 'NOM INCORRECT '//KNOMVA
                   CALL LEREUR
               ENDIF
               GOTO 9000
            ENDIF
C           IL EXISTE ( )  => VARIABLE
C           ELSE
C           IL EXISTE ( (  => TABLEAU
         ENDIF
C
C        L'IDENTIFICATEUR DE LA VARIABLE OU DU TABLEAU
         NCFVTS = I + I1 - 1
         KIDENV = KNOMVA( I+1 : NCFVTS )
         NCFVTS = NCFVTS + 1
C
      ELSE
C
C        TMS SIMPLE . RETOUR AVEC MNTMS ET LDTMS = 0
         CALL TNOUVR( KNOMVA , KTYPE , NTMS , MNTMS )
         LDTMS  = 0
C        TYPE MOT
         NOTYPE = 10
         GOTO 9999
C
      ENDIF
C
C     CODE OPERATION DE RECHERCHE D'UN TMS
      M = NUDCNB(KNOMTS)
      NOMTS(LHTMS)= KNOMTS(1:M)
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
      CALL TNOUVR( KNOMTS(1:M), KTYPE, NOTAMS(LHTMS), MNTAMS(LHTMS) )
      IF( NOTAMS(LHTMS) .LE. 0 ) THEN
         IF( IMNMTS .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'NON EXISTENCE DU TMS'
            KERR(2) = KNOMTS(1:NBCAER)
            CALL LEREUR
         ENDIF
         GOTO 9999
      ENDIF
      MCTAMS(LHTMS) = MNTAMS(LHTMS)
C     LE NOMBRE DE VARIABLES DECLAREES DU TABLEAU MS
      CALL TAMSTV( NOTAMS(LHTMS) , TYPNUM , NBVARI )
      NBVATC(LHTMS) = NBVARI
C
C     LE NOMBRE DE VARIABLE DU TMS A RECHERCHER
      NBVATA(LHTMS) = 1
C     LE NUMERO DANS DICOTD DU NOM DU TABLEAU DESCRIPTEUR
C     C'EST L'ENTIER APRES LE DOUBLE PRECISION DATE
      NOTADS(LHTMS) = MCN( MNTAMS(LHTMS) + MOTVAR(6) )
      IF( NOTADS(LHTMS) .LE. 0 .OR.
     %    NOTADS(LHTMS) .GT. NBTD ) THEN
          NBLGRC(NRERR) = 1
          WRITE(KERR(MXLGER)(1:4),'(I4)') NOTADS(LHTMS)
          KERR(1) =
     %   'VATSTD: NUMERO INCORRECT DE TABLEAU DESCRIPTEUR'
     %           //KERR(MXLGER)(1:4)
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
 100  IF( LHPILE .GT. 0 ) THEN
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
            CALL FUSNOM( DICOTD(NOTADS(LHTMS)), KNOMTS, NOMTS(LHTMS),
     %                   NB )
            IF( NB .EQ. 0 ) THEN
               NBLGRC(NRERR) = 3
               KERR(1) = 'VATSTD: FUSION IMPOSSIBLE DES NOMS '
               KERR(2) = DICOTD(NOTADS(LHTMS))
               NK = NUDCNB( KNOMTS )
               KERR(3) = KNOMTS(1:NK)
               CALL LEREUR
               GOTO 9000
            ENDIF
C
            IF( NOTAMS(LHTMS) .GT. 0 ) THEN
C              OUVERTURE
               CALL TAMSOU( NOTAMS(LHTMS) , MNTAMS(LHTMS) )
               IF( MNTAMS(LHTMS) .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'VATSTD: TMS NON OUVRABLE'
                  KERR(2) = NOMTS(LHTMS)
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
C
C              LE NOMBRE DE VARIABLES DECLAREES DU TABLEAU MS
               CALL TAMSTV( NOTAMS(LHTMS) , TYPNUM , NBVARI )
               NBVATC(LHTMS) = NBVARI
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
                  KERR(1) = 'PB POUR OUVRIR LE FICHIER'
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
                  WRITE(KERR(MXLGER)(1:12),'(I12)') MXLITD
                  KERR(1) = 'TABLEAU KTD SATURE'//KERR(MXLGER)(1:12)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              LECTURE DE LA LIGNE
               READ(NF,'(A)',END=120,IOSTAT=I) KTD(LHKTD)
               IF( I .EQ. 0 ) GOTO 110
               NBLGRC(NRERR) = 2
               KERR(1) ='ERREUR dans la LECTURE DU FICHIER de'
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
               NBLGRC(NRERR) = 2
               KERR(1) ='NOM_LEXIQUE INCORRECT '
               KERR(2) = KTD(NL)
               CALL LEREUR
               GOTO 9000
            ENDIF
C
C           RECHERCHE DE [CHAINE] ET AFFICHAGE
            NL = NLF
            NC = NCF
            CALL CHCHAI( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLD .GT. 0 ) THEN
               NL = NLF
               NC = NCF
            ENDIF
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
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'IDENTIFICATEUR INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOMUET .GT. 0 .AND. NOIDEN .GT. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = KIDENT(NOIDEN)//' EXISTE DEJA'
                  CALL LEREUR
               ENDIF
C
C              L'IDENTIFICATEUR NON RETROUVE EST AJOUTE
               IF( NBIDEN(LHTMS) .GE. MXIDEN ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) =  ' TABLE DES IDENTIFICATEURS SATUREE'
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
                  NL = NLF
                  NC = NCF
               ENDIF
C
C              RECHERCHE DE TYPE_VAR POUR UNE VARIABLE
               NBV    = 1
C              LA VARIABLE EST TRAITEE COMME UN TABLEAU D'UNE VARIABLE
               INTAVA = 0
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
                  KERR(1) ='IDENTIFICATEUR INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOIDEN .GT. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) =  KIDENT(NOIDEN)//' EXISTE DEJA'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C              L'IDENTIFICATEUR NON RETROUVE EST AJOUTE
               IF( NBIDEN(LHTMS) .GE. MXIDEN ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) ='TABLE DES IDENTIFICATEURS SATUREE'
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
                  KERR(1) = ' ( ) INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              RECHERCHE DU NOMBRE D'INDICES DU TABLEAU
               INTAVA = 1
               NL = NLD
               NC = NCD
               CALL CARAVA( NL , NC )
C              NB LE NOMBRE D'INDICES DU TABLEAU
               CALL CHINTA( NL,NC , NLD,NCD , NLF,NCF , NB )
C              STOCKAGE DE NB , (VAL_MIN,VAL_MAX) DANS LE TAS
C                       DE CHAQUE INDICE
               IF( LHTAS + 1 + NB + NB .GT. MXTAS ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'TAS SATURE'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              MISE A JOUR DES INDICES DU TABLEAU
               LHTAS = LHTAS + 1
               LETAS(  LHTAS  ) = NB
               IDENT(2,NBIDEN(LHTMS) ) = LHTAS
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
                  NL = NLF
                  NC = NCF
               ENDIF
C
C              TRAITEMENT DE TYPE_VAR :=
C             |entier {(<INTERVAL_ENT>:<CHAINE> {,<INTERVAL_ENT> : <CHAINE>})}
C             |^<NOM_LEXIQUE>    $ pointe sur un numero de nom dans ce lexique
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
C              TMS      =>21
C
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) ='TYPE_VAR INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C              MISE A JOUR DU TYPE DE L'IDENTIFICATEUR
               IDENT(1,NBIDEN(LHTMS)) = NOTYPE
C
C              ES CE LE BON IDENTIFICATEUR ?
               NOIDEN = NBIDEN(LHTMS)
               IF( KIDENV .EQ. KIDENT( NOIDEN ) ) THEN
C
C                 L'IDENTIFICATEUR EST RETROUVE
C                 -----------------------------
                  MNTMS = MCTAMS( LHTMS )
C
C                 RECHERCHE DU BON DECALAGE A PARTIR DU 1-ER MOT
C                 DE CE TABLEAU
                  LDTMS = LDTS( LHTMS )
C
C                 SON TYPE
                  NOTYPE = IDENT(1,NOIDEN)
C
C                 S'IL S'AGIT D'UNE VARIABLE RETOUR
                  IF( INTAVA .EQ. 0 ) GOTO 9999
C
C                 LE NOMBRE NB D'INDICES DU TABLEAU
                  I1 = IDENT(2,NOIDEN)
                  IF( I1 .LE. 0 ) THEN
                     NBLGRC(NRERR) = 1
                     KERR(1) = 'VATSTD:VARIABLE ET NON TABLEAU '//
     %                               KIDENV
                     CALL LEREUR
                     GOTO 9999
                  ENDIF
                  NB = LETAS( I1 )
C
C                 LE NOM DE LA VARIABLE EST COPIEE DANS LES DERNIERES
C                 LIGNES DE KTD POUR UTILISER LES SP DE ~/td/F
                  NLV = MXLITD-2
                  IF( LHKTD .GE. NLV ) THEN
                     NBLGRC(NRERR) = 1
                     KERR(1) = 'VATSTD: TABLE KTD SATUREE'
                     CALL LEREUR
                     GOTO 9000
                  ENDIF
                  KTD( NLV ) = KNOMVA
C                 LA POSITION DE LA 1-ERE ( DANS KTD(NLV) )
                  NCV  = INDEX( KTD(NLV) , '(' )
C                 LA POSITION DE LA 2-EME (
                  N1 = INDEX( KTD(NLV)(NCV+1:NCKTD) , '(' )
                  IF( N1 .LE. 0 ) THEN
                     NBLGRC(NRERR) = 1
                     KERR(1) =
     %              'VATSTD:VARIABLE D''UN TABLEAU SANS ( ( :'//KNOMVA
                     CALL LEREUR
                     GOTO 9000
                  ENDIF
                  NCV = NCV + N1
C
C                 LE NOM DE LA VARIABLE EST COPIEE DANS LES DERNIERES
C                 LIGNES DE KLG POUR UTILISER LES SP DE ~/LU/F
                  NLVU = MXKLG - 1
                  KLG(NLVU) = KNOMVA
C                 LA POSITION DE LA 2-EME (  DANS KLG(NLVU) )
                  NCVU = NCV
C
C                 LE NOMBRE DE VARIABLES POUR LES INDICES QUI PRECEDENT
                  NBV0 = 1
C                 LE NOMBRE DE VARIABLES QUI PRECEDENT DANS LE TABLEAU
                  L    = 0
C
C                 LA BOUCLE SUR LES INDICES DU TABLEAU
                  DO 230 I = 1 , NB
C
C                    L'INTERVALLE N1 N2 DE L'INDICE I DU TABLEAU
                     I1 = I1 + 2
                     N1 = LETAS( I1 - 1 )
                     N2 = LETAS( I1     )
                     NOVARU = 0
C
C                    LA VALEUR DEMANDEE DE L'INDICE
                     CALL VAENTI( NLV,NCV, NLDV,NCDV, NLFV,NCFV, NVAL )
                     IF( NLDV .LE. 0 ) THEN
C                       CE N'EST PAS UN ENTIER. EST CE UNE VARIABLE UTILISATEUR?
                        IF( NCVU .GT. 0 ) THEN
                           CALL CARPNB( NLVU , NCVU )
                           CALL CHVARU( NLVU,NCVU,NLFVU,NCFVU, NOVARU )
                           IF( NOVARU .GT. 0 ) THEN
C                             VARIABLE UTILISATEUR RETROUVEE
                              NVAL = NINT( DVARU(NOVARU) )
                              GOTO 280
                           ENDIF
                        ENDIF
C                       ERREUR
                        NBLGRC(NRERR) = 2
                        WRITE(KERR(MXLGER)(1:4),'(I4)') NB
                        KERR(1) ='VATSTD: VARIABLE '//KNOMVA
     %                          //' AVEC MOINS DE '//
     %                          KERR(MXLGER)(1:12) // ' INDICES'
                        KERR(2) = 'VARIABLE:'//KNOMVA
                        CALL LEREUR
                        GOTO 9000
                     ENDIF
C
C                    L'INDICE EST IL DANS LE BON INTERVALLE
 280                 IF( NVAL .LT. N1 .OR. NVAL .GT. N2 ) THEN
C                       NON
                        NBLGRC(NRERR) =3
                        WRITE(KERR(MXLGER)(1:12),'(I12)') I
                        KERR(1) ='VATSTD: INDICE'//KERR(MXLGER)(1:12)
     %                          //' DE LA VARIABLE '//KNOMVA
                        WRITE(KERR(MXLGER)(1:12),'(I12)') NVAL
                        KERR(2) =' A POUR VALEUR'//KERR(MXLGER)(1:12)
                        WRITE(KERR(MXLGER)( 1:12),'(I12)') N1
                        WRITE(KERR(MXLGER)(21:32),'(I12)') N2
                        KERR(3) =' NON DANS L''INTERVALLE '
     %                          //KERR(MXLGER)( 1:12)
     %                          //KERR(MXLGER)(21:32)
                        CALL LEREUR
                        GOTO 9000
                     ENDIF
C
C                    LA VALEUR DE L'INDICE EST CORRECTE
C                    LE DECALAGE EST ESTIME
                     L = L + NBV0 * ( NVAL - N1 )
C                    LE NOMBRE DE VARIABLES QUI PRECEDENT
                     NBV0 = NBV0 * ( N2 - N1 + 1 )
C
C                    PASSAGE A L'INDICE SUIVANT
                     IF( NOVARU .EQ. 0 ) THEN
                        NLV = NLFV
                        NCV = NCFV
C                       PASSAGE AU CARACTERE SUIVANT QUI DOIT ETRE ','
C                       OU BIEN )
                        CALL CAR1NB( NLV , NCV )
                        NCVU = NCV
                     ELSE
                        NLVU  = NLFVU
                        NCVU  = NCFVU
C                       PASSAGE AU CARACTERE SUIVANT QUI DOIT ETRE ','
C                       OU BIEN )
                        CALL CARPNB( NLVU , NCVU )
                        NCV = NCVU
                     ENDIF
                     NCFVTS = NCV
C
 230              CONTINUE
C
C                 CONVERSION EN MOTS ET DECALAGE
                  LDTMS = LDTMS + L * MOTVAR( NOTYPE )
C
                  IF( NOTYPE .EQ. 12 ) THEN
C                    EXTRACTION DE L'UNE DES COMPOSANTES DE XYZ
C                    EXISTE-T-IL '(1)' OU '(2)' OU '(3)' APRES LE XYZ ?
                     I = NCV
                     CALL CAR1NB( NLV , NCV )
                     IF( KTD(NLV)(NCV:NCV) .EQ. '(' ) THEN
                        CALL VAENTI(NLV,NCV, NLDV,NCDV, NLFV,NCFV, NVAL)
                        IF( NVAL .GT. 0 .AND. NVAL .LE. 3 ) THEN
                           NOTYPE = 5
                           LDTMS  = LDTMS + NVAL - 1
                        ENDIF
C                       PASSAGE A LA PARENTHESE FERMANTE
                        CALL CAR1NB( NLFV , NCFV )
                        CALL CAR1NB( NLFV , NCFV )
                        IF( KTD(NLFV)(NCFV:NCFV) .EQ. ')' ) THEN
                           NLV = NLFV
                           I   = NCFV
                        ELSE
                           I = 0
                        ENDIF
                     ELSE
C                       PAS DE (no de 1 a 3) du XYZ  le 29/9/97
C                       PASSAGE A LA PARENTHESE FERMANTE
                        CALL CAR1NB( NLFV , I )
                     ENDIF
                     NCFVTS = I
                     IF( KTD(NLV)(NCFVTS:NCFVTS) .NE. ')' ) THEN
C                       XYZ NE SE TERMINE PAS AVEC UNE PARENTHESE FERMANTE
                        NBLGRC(NRERR) = 2
                        KERR(1) = 'VATSTD: VARIABLE TMS '//KNOMVA
                        KERR(2) = ' ) FINALE OUBLIEE'
                        CALL LEREUR
                        GOTO 9000
                     ENDIF
                  ENDIF
                  GOTO 9999
               ENDIF
C
C              LE CARACTERE DE FIN DE TYPE
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
C                 CALCUL DU DECALAGE SANS MODIFICATION PAR ALIENATION
C                 DU PARAMETRE NOMUET
                  NOMUE0 = NOMUET
                  NOMUET = 0
                  CALL MOTYVA( NBIDEN(LHTMS) , NRETOU )
C                 REGENERATION DE NOMUET A PARTIR DE SA SAUVEGARDE
                  NOMUET = NOMUE0
                  IF( NRETOU .GT. 0 ) GOTO 9000
C
               ELSE IF( NOTYPE .EQ. 21 ) THEN
C
C                 TYPE_VAR = tms NOM_TMS
C                 ......................
                  CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
                  IF( NLD .LE. 0 ) THEN
                     NBLGRC(NRERR) = 2
                     KERR(1) = 'NOM_LEXIQUE INCONNU'
                     KERR(2) = KTD(NL)
                     CALL LEREUR
                     GOTO 9000
                  ENDIF
C
C                 NOMBRE DE VARIABLES ENTIERES (NO DU TMS )
                  CALL NVARID( NBIDEN(LHTMS) , NBV , NOTYPE )
C
C                 DECALAGE DE NBV ENTIERS
                  LDTS(LHTMS) = LDTS(LHTMS) + NBV
C
C                 PASSAGE AU ; SUIVANT DERRIERE LE NOM DU TMS
                  NL = NLF
                  NC = NCF
C
               ELSE
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
                  KERR(1) = 'VATSTD:TYPE INCORRECT '//KERR(MXLGER)(1:4)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              FIN <TYPE_VAR> PASSAGE AU CARACTERE ;
               CALL CHCAR( ';' , NL , NC )
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
                  KERR(1) = 'PILE SATUREE'
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
                  KERR(1) = 'INTERVALLE_E INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOIDEN .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'IDENTIFICATEUR INCONNU '//KTD(NLD)(NCD:NCF)
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
                  KERR(1) = 'VATSTD: ; INCORRECT'
                  KERR(2) = KTD(NLD)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               GOTO 200
C
            ELSE
C
C              ERREUR
               NBLGRC(NRERR) = 2
               KERR(1) = 'VATSTD: TYPE INCORRECT'
               KERR(2) = KTD(NLD)
               CALL LEREUR
               GOTO 9000
            ENDIF
         ENDIF
      ENDIF
      GOTO 9999
C
C     ERREUR
C     ======
 9000 NOTYPE = 0
      MNTMS  = 0
      LDTMS  = 0
C
C     MISE A JOUR DANS L'ETAT INITIAL
 9999 LHTMS  = LHTMS0
      LHKTD  = LHKTD0
      END
