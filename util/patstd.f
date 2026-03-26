      SUBROUTINE PATSTD( KNOMTD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  GENERER LES PARAMETER D'UN TMS    ( PAS DE RECURSIVITE DE TMS )
C -----
C
C ENTREE :
C --------
C KNOMTD : NOM DU TABLEAU DESCRIPTEUR DU TMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C23456---------------------------------------------------------------012
      PARAMETER (NF=91)
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
C.......................................................................
      CHARACTER*160     KNOM
      CHARACTER*24      KNOMPA
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOMTD
      CHARACTER*96      KMOT
C
C     INITIALISATIONS
C     LA PILE DES OPERATIONS
      LHPILE = 1
C     LA PILE DES TMS
      LHTMS  = 1
C     LE NUMERO DE LIGNE DU DERNIER TD RANGE DANS KTD
      LHKTD0 = LHKTD
C     IMPOSSIBILITE DE GENERER LE PARAMETER ( 1 SINON 0)
      IMPOSS = 0
      IMPOS1 = 0
C
C     CODE OPERATION DE TRAITEMENT D'UN TMS
      NOMTS(LHTMS)= KNOMTD
      LAPILE(0,1) = LHTMS
      LAPILE(1,1) = 1
      LAPILE(2,1) = 0
      LAPILE(3,1) = 0
      LAPILE(4,1) = 0
      LAPILE(5,1) = 0
C     IMPRESSION DEMANDEE
      LAPILE(6,1) = 1
C
      NBLGRC(NRERR) = 1
      KERR(1) = 'PARAMETRES DU TABLEAU: '//KNOMTD
      CALL LERESU
C
C     LE NOMBRE DE VARIABLE DU TMS A RECHERCHER
      NBVATA(LHTMS) = 1
C     LE NUMERO DANS DICOTD DU NOM DU TABLEAU DESCRIPTEUR
      NOTADS(LHTMS) = NONMTD( KNOMTD )
      IF( NOTADS(LHTMS) .LE. 0 .OR.
     %    NOTADS(LHTMS) .GT. NBTD ) THEN
          NBLGRC(NRERR) = 1
          WRITE(KERR(MXLGER)(1:4),'(I4)') NOTADS(LHTMS)
          KERR(1) =
     %   'PATSTD: NUMERO INCORRECT DE TABLEAU DESCRIPTEUR'
     %    //KERR(MXLGER)(1:4)
          CALL LEREUR
          RETURN
      ENDIF
C
C     OUVERTURE DU FICHIER SUPPORT DU TD A LIRE DANS KTD
C     --------------------------------------------------
C     FORMATION DU NOM DU FICHIER SUPPORT DU TD
      CALL NTDFIC( DICOTD(NOTADS(LHTMS)) , KNOM )
      OPEN( FILE=KNOM , UNIT=NF , IOSTAT=I )
      IF( I .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) ='POUR OUVRIR LE FICHIER'
         KERR(2) = KNOM(1:NBCAER)
         CALL LEREUR
         GOTO 9000
      ENDIF
C     LE NUMERO DE LIGNE DE DEBUT DU TD
      NL = LHKTD + 1
      NC = 0
C
 10   LHKTD = LHKTD + 1
      IF( LHKTD .GT. MXLITD ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:12),'(I12)') MXLITD
         KERR(1) ='TABLEAU KTD SATURE '//KERR(MXLGER)(1:12)
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     LECTURE DE LA LIGNE
      READ(NF,'(A)',END=20,IOSTAT=I) KTD(LHKTD)
      IF( I .EQ. 0 ) GOTO 10
      NBLGRC(NRERR) = 2
      KERR(1) = 'ERREUR LECTURE DU FICHIER'
      KERR(2) = KNOM(1:NBCAER)
      CALL LEREUR
      GOTO 9000
C
C     FIN DE LECTURE DU FICHIER
 20   LHKTD = LHKTD - 1
C     CARACTERE DE DEBUT DU TMS
      LAPILE(2,LHPILE) = NL
      LAPILE(3,LHPILE) = NC
C     LE CARACTERE DE RETOUR EST DEJA INITIALISE
C     IMPRESSION DEMANDEE
      LAPILE(6,LHPILE) = 1
C     FERMETURE DU FICHIER SUPPORT DU TD
      CLOSE( UNIT=NF )
C
C     CREATION DU FICHIER SUPPORT DES PARAMETER DU TABLEAU DESCRIPTEUR
C     ----------------------------------------------------------------
      I    = INDEX( KNOM , ' ' )
C     ON AJOUTE .INC DERRIERE LE NOM DU FICHIER
      KNOM(I:I+3) = '.inc'
C     ON REMPLACE  /td/d/  par  /incl/  ( CF LE SP NTDFIC )
      I    = INDEX( KNOM , '/td/d/' )
      KNOM(I:I+5) = '/incl/'
C
      OPEN( FILE=KNOM , UNIT=NF , IOSTAT=I )
      IF( I .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'PROBLEME POUR CREER LE FICHIER'
         KERR(2) = KNOM(1:NBCAER)
         CALL LEREUR
         GOTO 9000
      ENDIF
CCCC
CCCC     GENERATION DE L'EN TETE DU FICHIER PARAMETER
CCC      WRITE( NF , 10020 )
CCC      WRITE( NF , 10021 ) KNOM
CCC      WRITE( NF , 10020 )
CCC10020 FORMAT('C')
CCC10021 FORMAT('C     FICHIER ',A)
C
C     GENERATION DE LA 1-ERE LIGNE DES PARAMETER
      KNOMPA = '      PARAMETER ('
      WRITE( NF  , '(A)' ) KNOMPA
C
C     LE DECALAGE DANS LE TS POUR ATTEINDRE LA VARIABLE TRAITEE
      LDTS(LHTMS) = 0
C
C     LE PARAMETRE D'IMPRESSION
      IMPRES = 1
      NOMUET = 1
C     LE POINTEUR SUR LE DERNIER CARACTERE STOCKE DANS LA LIGNE BUFFER
      LCLIGN = 0
C     LA LIGNE BUFFER EST BLANCHE
      KLIGNE = ' '
C     LE NOMBRE D'IDENTIFICATEURS STOCKES DU CAS
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

         ELSE IF( LAPILE(1,LHPILE) .EQ. 1 ) THEN

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

C           LE PREMIER CARACTERE DU TMS
            NL = LAPILE(2,LHPILE)
            NC = LAPILE(3,LHPILE)

C           deftms NOM_LEXIQUE [CHAINE] TYPE fintms
C           =======================================
            CALL CHMOT( 'deftms' , NL,NC , NLD,NCD, NLF,NCF )
C           EN SORTIE NLD,NCD ET NLF,NCF DEBUT ET FIN DU MOT DANS KTD

C           RECHERCHE DU NOM_LEXIQUE DU TMS
            NL = NLF
            NC = NCF
            CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLD .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) ='patstd: NOM_LEXIQUE INCORRECT'
               KERR(2) = KTD(NL)
               CALL LEREUR
               GOTO 9000
            ENDIF

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

            KMOT = KTD(NLD)(NCD:NCF)
            NBCMOT = NCF-NCD+1
            IF( KMOT(1:NBCMOT) .EQ. 'variable' ) THEN
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
                  KERR(1) = 'patstd:IDENTIFICATEUR INCORRECT'
                  KERR(2) =  KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IDEM = 0
               IF( NOMUET .GT. 0 .AND. NOIDEN .GT. 0 ) THEN
C                 L'IDENTIFICATEUR EXISTE DEJA.
C                 VERIFICATION DU MEME DECALAGE DANS LE TMS
                  IF( LDTS(LHTMS) .NE. IDENT(3,NOIDEN) ) THEN
C                    LES DECALAGES SONT DIFFERENTS
                     NBLGRC(NRERR) = 1
                     KERR(1) ='patstd:'//KIDENT(NOIDEN)//
     %              ' EXISTE DEJA AVEC DES DECALAGES DIFFERENTS'
                     CALL LEREUR
                  ELSE
C                    LES DECALAGES SONT IDENTIQUES
                     IDEM = 1
                  ENDIF
               ENDIF
C              L'IDENTIFICATEUR NON RETROUVE EST AJOUTE
               IF( NBIDEN(LHTMS) .GE. MXIDEN ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'TABLE DES IDENTIFICATEURS SATUREE'
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
C              LA VARIABLE EST TRAITEE COMME UN TABLEAU D'UNE VARIABLE
               GOTO 275
C
            ELSE IF( KMOT(1:NBCMOT) .EQ. 'tableau' ) THEN
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
                  KERR(1) ='patstd:IDENTIFICATEUR INCORRECT'
                  KERR(2) = KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IDEM = 0
               IF( NOIDEN .GT. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'patstd:TABLEAU '//KIDENT(NOIDEN)//
     %          ' EXISTE DEJA'
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
               CALL CHPAOF( NL,NC , NLD,NCD , NLFPA,NCFPA )
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'patstd: ( ) INCORRECT'
                  KERR(2) = KTD(NL)
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
                  KERR(1) = 'TAS SATURE'
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              MISE A JOUR DES INDICES DU TABLEAU
               LHTAS = LHTAS + 1
               LETAS( LHTAS ) = NB
               IDENT(2, NBIDEN(LHTMS) ) = LHTAS
               NL  = NLD
               NC  = NCD
               DO 250 I=1,NB
C                 RECHERCHE DE INTERVAL_ENT ET MISE A JOUR DANS LE TAS
                  CALL PAINEN(  NL,NC , NLD,NCD , NLF,NCF ,
     %                          N1,N2 , IMPOS )
C                 SI UN DES INDICES EST UN IDENTIFICATEUR  ALORS
C                 IMPOS>0 ET N1=N2=0
                  IMPOSS = MAX( IMPOSS , IMPOS )
                  IF( IMPOS .NE. 0 ) GOTO 260
C
                  LHTAS = LHTAS + 2
                  LETAS(LHTAS-1) = N1
                  LETAS(LHTAS  ) = N2
C                 PASSAGE EN FIN INTERVAL_ENT
                  NL = NLF
                  NC = NCF
C                 SAUT DE , OU )
                  CALL CAR1NB( NL , NC )
 250           CONTINUE
C
C              RECHERCHE DE [CHAINE] AU DELA DE ) ET AFFICHAGE
 260           NL = NLFPA
               NC = NCFPA
               CALL CHCHAI( NL,NC , NLD,NCD , NLF,NCF )
               IF( NLD .GT. 0 ) THEN
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
C              EN RETOUR NOTYPE LE NUMERO DU TYPE DE LA VARIABLE
C              LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C              REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C              COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C              TYPEOBJET=>13 TMS      => 21
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
                  KERR(1) = 'patstd:TYPE_VARIABLE INCORRECT '//
     %                       KERR(MXLGER)(1:4)
                  KERR(2) =  KTD(NL)
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
C                 TYPE_VAR:= |entier {(<INTERVAL_ENT> : <CHAINE>
C                                    {,<INTERVAL_ENT> : <CHAINE>})}
C                            |^<NOM_LEXIQUE>
C                               $ pointe sur un numero de nom dans ce lexique
C                            |caractere
C                            |reel
C                            |reel2
C                            |xyz
C                            |typeobjet
C                 ...........................................................
C                 SAUVEGARDE DE NOMUET
                  NN = NOMUET
                  IF( IDEM .NE. 0 ) THEN
C                    L'IDENTIFICATEUR A DEJA ETE VU .
C                    PROTECTION DE NOMUET
                     NN     = NOMUET
C                    NON GENERATION DE LA VARIABLE
                     NOMUET = 0
                  ENDIF
C                 GENERATION DU PARAMETER DE L'IDENTIFICATEUR
 280              CALL PATYVA( NF , NBIDEN(LHTMS) , IMPOSS )
C                 SAUVEGARDE DE IMPOSS EN SORTIE
                  IMPOS1 = MAX( IMPOS1 , IMPOSS )
C                 RESTAURATION DE NOMUET
                  NOMUET = NN
C
               ELSE IF( NOTYPE .EQ. 21 ) THEN
C
C                 TYPE_VAR = tms NOM_TMS
C                 ......................
                  CALL CHLEXI( NL,NC, NLD,NCD, NLF,NCF )
                  IF( NLD .LE. 0 ) THEN
                     NBLGRC(NRERR) = 2
                     KERR(1) = KTD(NL)
                     KERR(2) = 'NOM_LEXIQUE INCONNU'
                     CALL LEREUR
                     GOTO 9000
                  ENDIF
C
C                 PASSAGE AU ; SUIVANT DERRIERE LE NOM DU TMS
                  NL = NLF
                  NC = NCF
                  NN = NOMUET
C
C                 LE PARAMETER DU NUMERO DU TMS
                  GOTO 280
C
               ELSE
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
                  KERR(1) ='patstd: TYPE INCORRECT '// KERR(MXLGER)(1:4)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
C
C              FIN <TYPE_VAR> PASSAGE AU CARACTERE ;
               CALL CHCAR( ';' , NL , NC )
               GOTO 200
C
            ELSE IF( KMOT(1:NBCMOT) .EQ. 'fintms' ) THEN
C
               IF( LHTMS .GT. 1 ) THEN
C                 DECALAGE D'UN ENTIER DANS CE TABLEAU
                  LDTS(LHTMS-1) = LDTS(LHTMS-1) + 1
               ENDIF
               GOTO 100
C
            ELSE IF( KMOT(1:NBCMOT) .EQ. 'cas' ) THEN
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
C
C              POUR GENERER LES PARAMETER SAUVEGARDE
C              DU DEBUT DU DECALAGE LDTS POUR LE RESTAURER
C              AU DEBUT DE TOUT <INTERVALLE_E> RENCONTRE
               LAPILE(2,LHPILE) = LDTS(LHTMS)
C              SAUVEGARDE DE IMPOSS POUR CHAQUE DEBUT DE INTERVALLE_E
               LAPILE(3,LHPILE) = IMPOSS
C
C              LE PARAMETRE D'IMPRESSION
               LAPILE(6,LHPILE) = IMPRES * NOMUET
C
C              RECHERCHE DE IDENT
               NL = NLF
               NC = NCF
               CALL CHIDEN( NL,NC, NLD,NCD, NLF,NCF , NOIDEN )
               IF( NLD .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'patstd:INTERVALLE_E INCORRECT'
                  KERR(2) =  KTD(NL)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               IF( NOIDEN .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'IDENTIFICATEUR '//KTD(NLD)(NCD:NCF)
     %                    //' INCONNU'
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
C              LA VALEUR DES BORNES DE L'<INTERVALLE_E> TRAITE
               CALL CHINT2( NL,NC, NLD,NCD, NLF,NCF )
C
C              RESTAURATION DU DECALAGE DANS TMS DU DEBUT DE CAS
               LDTS(LHTMS) = LAPILE(2,LHPILE)
C              RESTAURATION DE IMPOSS DU DEBUT DU CAS
               IMPOSS = LAPILE(3,LHPILE)
C
C              PASSAGE AU : QUI SUIT <INTERVALLE_E>
               NL = NLF
               NC = NCF
               CALL CHCAR( ':' , NL , NC )
               GOTO 200
C
            ELSE IF( KMOT(1:NBCMOT) .EQ. 'fincas' ) THEN
C
C              LE PARAMETRE IMPOSS EST MIS A JOUR PAR SA VALEUR MAXIMALE
C              ATTEINTE AU COURS DU CAS
               IMPOSS = IMPOS1
C
C              LE CAS EST DEPILE
               LHPILE = LHPILE - 1
               NL     = NLD
               NC     = NCF
               CALL CHCAR( ';' , NL , NC )
               IF( NL .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'patstd: ; INCORRECT'
                  KERR(2) = KTD(NLD)
                  CALL LEREUR
                  GOTO 9000
               ENDIF
               GOTO 200

            ELSE

C              ERREUR
               NBLGRC(NRERR) = 3
               KERR(1) = 'patstd: TYPE INCORRECT'
               KERR(2) = KTD(NLD)
               KERR(3) = KMOT
               PRINT*,'patstd: KMOT=',KMOT(1:NBCMOT) //'??????? '
               PRINT*,'patstd: KTD(NLD=',NLD,')=',KTD(NLD)
               do kl=1,LHKTD
                  print*,'patstd: ktd(',kl,')=',ktd(kl)
               enddo
               CALL LEREUR
               GOTO 9000

            ENDIF
         ENDIF
      ENDIF

C     GENERATION DE LA DERNIERE LIGNE DES PARAMETER
      KNOMPA = '     % )'
      WRITE( NF, '(A)' ) KNOMPA

C     FERMETURE DU FICHIER SUPPORT DES PARAMETER DU TD
      CLOSE( UNIT=NF )

C     CONDENSATION DE L'INSTRUCTION PARAMETER AFIN DE REDUIRE
C     LE NOMBRE DE CARTES SUITES
      CALL COPAFI( KNOM , NF )

C     ERREUR
C     ======
 9000 LHKTD = LHKTD0
      RETURN
      END
