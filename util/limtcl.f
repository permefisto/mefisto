      SUBROUTINE LIMTCL( NMMTCL, NOMTCL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CHOISIR UN MOT CLE PARMI UNE LISTE
C -----
C
C ENTREES :
C ---------
C NMMTCL  : LE NOM DE LA LISTE DES MOTS CLES
C
C SORTIES :
C ---------
C NOMTCL  : >= 0 LE NUMERO DU MOT CLE CHOISI
C           = -2 UN PROBLEME  LORS DE LA LECTURE DU FICHIER ASSOCIE
C           = -1 UN ABANDON ( @ ) EST DEMANDE DANS L'ARBORESCENCE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/msvaau.inc"
      include"./incl/homdir.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/xvfontes.inc"
C.......................................................................
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NFDOCU, NFFRAP,
     %                  NUNITE(27)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NMMTCL
      CHARACTER*160     KNOM
      CHARACTER*64      KNM

      IF( INTERA .GT. 0 ) CALL CHOIXFONTE( NPHFCO )
C
C     INITIALISATIONS POUR FAIRE CROIRE AU TRAITEMENT D'UN TD
      NOMTCL = 0
      NBLGRC(NRMENU) = 0
      KMENU(1) = ' '
      LHTMS  = 1
      LHPILE = 1
      NOMUET = 1
      IMPRES = 1
C     LE NOMBRE DE MOTS DU RESULTAT
      NBVATC( LHTMS ) = 2
C     L'ADRESSE MCN DU TABLEAU OU STOCKER LE RESULTAT
      MCTAMS( LHTMS ) = MNMTCL
C     LE DECALAGE DU RESULTAT DANS LE TABLEAU
      LDTS  ( LHTMS ) = 0
C
C     LA LECTURE DE LA LISTE DES MOTS CLES DANS KTD
C     =============================================
      KNM = NMMTCL
      CALL MINUSC( KNM )
      KNOM   = HOMDIR // '/td/m/' // KNM
      OPEN( FILE=KNOM, UNIT=NFDOCU , STATUS='OLD' , IOSTAT=NOMTCL )
      IF( NOMTCL .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  NMMTCL//' LISTE INCONNUE de MOTS CLES'
         ELSE
            KERR(1) =  NMMTCL//' UNKNOWN LIST of KEYWORDS'
         ENDIF
         CALL LEREUR
         GOTO 9990
      ENDIF
      LHKTD = 0
C
C     LECTURE DANS KTD DE LA LISTE DE MOTS CLES
 10   READ( UNIT=NFDOCU , FMT='(A)' , IOSTAT=NOMTCL ) KTD(LHKTD+1)
      IF( NOMTCL .NE. 0 ) GOTO 100
      LHKTD = LHKTD + 1
      GOTO 10
C
C     FIN DE LECTURE DE LA LISTE DES MOTS CLES
 100  CLOSE( UNIT=NFDOCU )
C
C     variable IDENT [CHAINE] TYPE_VAR ;
C     ----------------------------------
C     RECHERCHE DE variable
      NL = 1
      NC = 0
      CALL CHMOT( 'variable' , NL,NC , NLD,NCD, NLF,NCF )
      IF( NLD .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='LU: FICHIER SANS MOT CLE variable'
         ELSE
            KERR(1) ='LU: FILE WITHOUT KEYWORD variable'
         ENDIF
         KERR(2) = KNOM(1:NUDCNB(KNOM))
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     RECHERCHE DE IDENT
      NL = NLF
      NC = NCF
      CALL CHIDEN( NL,NC, NLD,NCD, NLF,NCF , NOIDEN )
C     NOIDEN NUMERO DE L'IDENTICATEUR, 0 SI NON RETROUVE
      IF( NLD .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='LIMTCL: IDENTIFICATEUR INCORRECT'
         ELSE
            KERR(1) ='LIMTCL: INCORRECT IDENTIFICATOR'
         ENDIF
         KERR(2) =  KTD(NL)
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     AFFICHAGE DE L'IDENT
      CALL AFNOM( NLD,NCD, NLF,NCF )
      NOIDEN = 1
      NBIDEN( LHTMS  ) = NOIDEN
      KIDENT( NOIDEN ) = KTD(NLD)(NCD:NCF)
C
C     L'IDENTIFICATEUR DE CETTE VARIABLE
      IDENT( 0 , NOIDEN ) = LHPILE
      IDENT( 2 , NOIDEN ) = 0
      IDENT( 3 , NOIDEN ) = LDTS( LHTMS )
      IDENT( 4 , NOIDEN ) = NLD
      IDENT( 0 , NOIDEN ) = NCD
C
C     NUIDEN SAUVEGARDE LE NUMERO DE L'IDENTIFICATEUR DANS TD.INC
      NUIDEN = NOIDEN
C
C     RECHERCHE DE [CHAINE] ET MODIFICATION
      NL = NLF
      NC = NCF
      CALL CHCHAI( NL,NC, NLD,NCD, NLF,NCF )
      IF( NLD .GT. 0 ) THEN
         CALL AFCHAI( NLD,NCD, NLF,NCF )
         NL = NLF
         NC = NCF
      ENDIF
C
C     RECHERCHE DE TYPE_VAR POUR CETTE VARIABLE
C     TRAITEMENT DE TYPE_VAR :=  ( DOIT ETRE REDUIT ICI A entier )
C     |entier {(<INTERVAL_ENT>:<CHAINE> {,<INTERVAL_ENT> : <CHAINE>})}
C     |^<NOM_LEXIQUE>    $ pointe sur un numero de nom dans ce lexique
C     |reel
C     |reel2
C     |xyz
C     |tms <NOM_TMS>
C
      CALL CHTYPV( NL,NC , NLD,NCD , NLF,NCF , NOTYPE )
C     LE TYPE DE LA VARIABLE   NOTYPE  VAUT
C     LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C     REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C     COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C     TMS      => 21
      IF( NLD .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) ='LIMTCL: TYPE_VAR INCORRECT'
         KERR(2) = KTD(NL)
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     LE MOT QUI SUIT
      N2 = INDEX( KTD(NLD)(NCD:NCKTD) , ' ' )
      IF( N2 .LE. 0 ) THEN
         N2 = NCKTD
      ELSE
         N2 = NCD + N2 - 1
      ENDIF
C
C     S'IL EXISTE UN COUPLE (  ) AFFICHAGE DU TYPE D'ABORD
      IF( KTD(NLF)(NCF:NCF) .EQ. ')' ) THEN
CCC         IF( INTERA .LE. 1 ) THEN
CCC            CALL AFNOM( NLD,NCD, NLD,N2 )
CCC            CALL AFLIGN
CCC         ENDIF
         NCD = N2
         CALL CAR1NB( NLD , NCD )
      ENDIF
C
      IDENT(1,NBIDEN(LHTMS)) = NOTYPE
      CALL AFLIGN
C
      IF( NOTYPE .EQ. 4 ) THEN
C
C        TYPE_VAR:= |entier {(<INTERVAL_ENT>:<CHAINE>
C                           {,<INTERVAL_ENT> : <CHAINE>})}
         IF( INTERA .LE. 1 ) THEN
            CALL MOTYVA( NBIDEN(LHTMS) , NRETOU )
         ELSE
            CALL MOTYVM( NBIDEN(LHTMS) , NRETOU )
         ENDIF
         IF( NRETOU .GT. 0 ) THEN
C           ABANDON DE LA LECTURE
            NOMTCL = -1
            CALL RECTEF( NRMENU )
            CALL RECTEF( NRINVI )
            GOTO 9999
         ENDIF
C
      ELSE
C
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='LIMTCL: TYPE INCORRECT DIFFERENT DE entier'
         ELSE
            KERR(1) ='LIMTCL: INCORRECT TYPE DIFFERENT OF entier'
         ENDIF
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     SI RETOUR CORRECT LE NUMERO DU MOT CLE EST DANS
      NOMTCL = MCN( MNMTCL )
      GOTO 9999
C
C     ERREUR
 9990 NOMTCL = -2
C     REINITIALISATION COMPLETE DU TABLEAU KTD
C
 9999 LHKTD = 0
      RETURN
      END
