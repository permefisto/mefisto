      SUBROUTINE AREXPD( NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERATION DE L'ARBRE DES NBEXPD EXPRESSIONS
C -----
C
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES DANS LE COMMON / ILANUT /:
C ----------------------------------
C NBEXPD : NOMBRE D'EXPRESSIONS DEFINIES DANS LE TABLEAU NCOPER
C
C SORTIES :
C ---------
C NRETOU : 1 L'ARBRE NE PEUT ETRE FORME A CAUSE D'UNE ERREUR
C          0 L'ARBRE EST FORME ET PRET A ETRE CALCULE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     LE CHAINAGE DES OPERATIONS
      DO 10 I=LSEXPD(1,1),LSEXPD(2,1)
C        LA PRECEDENTE
         NCOPER(-1,I) = I - 1
C        LA SUIVANTE
         NCOPER( 0,I) = I + 1
 10   CONTINUE
C     LE PRECEDENT DU PREMIER OPERATEUR
      NCOPER( -1 , LSEXPD(1,1) ) = 0
C     SON PRECEDENT
      NCOPER( -1 , 0 ) = -1
C     SON SUIVANT
      NCOPER( 0 , 0 ) = LSEXPD(1,1)
C     SON CODE OPERATION
      NCOPER( 1 , 0 ) = -1
      IF( NOPER .GE. MXNCOP ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: TABLE DES OPERATEURS SATUREE'
            KERR(2) = 'LU: AUGMENTER MXNCOP'
         ELSE
            KERR(1) = 'LU: SATURATED ARRAY OF OPERATORS'
            KERR(2) = 'LU: AUGMENT MXNCOP'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C     LE SUIVANT DU DERNIER
      NCOPER( 0 , LSEXPD(2,1) ) = MXNCOP
C     SON PRECEDENT
      NCOPER( -1 , MXNCOP ) = LSEXPD(2,1)
C     SON SUIVANT
      NCOPER(  0 , MXNCOP ) = -MXNCOP
C     SON CODE OPERATION
      NCOPER(  1 , MXNCOP ) = -2
C
C     LA DERNIERE PARENTHESE FERMANTE A POUR SUIVANTE MXNCOP
      NOP1 = 0
      DO 20 I=LSEXPD(2,1),LSEXPD(1,1),-1
         IF( NCOPER(1,I) .EQ. 2 .OR. NCOPER(1,I) .EQ. 9 ) THEN
            NOP1 = I
            GOTO 25
         ENDIF
 20   CONTINUE
 25   NCOPER(3,NOP1) = MXNCOP
C
C     INITIALISATION DU CHAINAGE DES RACINES VIDES
C     POINTEUR SUR LA 1-ERE RACINE VIDE
      LSRACI( 0 ) = 1
      DO 30 I=1,MXRACI
         LSRACI( I ) = - ( I + 1 )
 30   CONTINUE
      LSRACI( MXRACI ) = 0
C
C     LA BOUCLE DE CHAINAGE DES NBEXPD EXPRESSIONS
      DO 1000 LEXPD = NBEXPD , 1 , -1
C        LE PREMIER ET DERNIER CODE OPERATION DE L'EXPRESSION
C        ( ,  OU ,  ,  OU   ,  )
         N1EXPD = LSEXPD( 1 , LEXPD )
         N2EXPD = LSEXPD( 2 , LEXPD )
C        LES EXTREMES DE L'EXPRESSION A TRAITER
C        PAS DE PASSAGE AU DELA DANS UN PREMIER TEMPS
C        NDEBUT ET NFIN LES POINTEURS SUR ( )  OU  ( ,  OU  , ,  OU  , )
         NDEBUT = NCOPER( -1 , N1EXPD )
         NFIN   = NCOPER(  0 , N2EXPD )
C        PAS DE RACINE POUR CETTE EXPRESSION
         NP1    = -1
C
C        LE TRAITEMENT DES ( ) EN COMMENCANT PAR LES PLUS INTERNES
C        ---------------------------------------------------------
C        RECHERCHE DE LA 1-ERE (
 80      NOP1 = NCOPER( 0 , NDEBUT )
         NOP2 = NCOPER( 0 , NOP1 )
         IF( NOP2 .EQ. NFIN .AND. NP1 .GT. 0 ) THEN
C           EXPRESSION TOTALEMENT TRAITEE
C           ON DOIT AVOIR NDEBUT -> RACINE -> NFIN
            NOP = LSRACI( NP1 )
            IF( NCOPER(0,NDEBUT) .NE. NOP .OR.
     %          NCOPER(-1,NFIN ) .NE. NOP ) THEN
                NBLGRC(NRERR) = 1
                WRITE(IMPRIM,*) 'LU: ',NOP,NFIN,
     %        ' DIFFERENT DE ',NCOPER(-1,NOP),
     %                         NCOPER(0,NOP)
                GOTO 9900
            ENDIF
            GOTO 800
         ENDIF
C
C        NOP1 EST LE DEBUT DE L'EXPRESSION SOIT PAR EXEMPLE (
C        NOP2 EST LA FIN   DE L'EXPRESSION SOIT PAR EXEMPLE )
 90      IF( NCOPER(1,NOP1) .NE. 1 ) THEN
C           L'OPERATION SUIVANTE
            NOP1 = NCOPER(0,NOP1)
            IF( NOP1 .EQ. NFIN ) THEN
C              IL N'EXISTE PLUS DE (
               NOP1 = NCOPER(0,NDEBUT)
               NOP2 = NCOPER(-1,NFIN)
               GOTO 100
            ELSE
C              PASSAGE A L'OPERATION SUIVANTE
               GOTO 90
            ENDIF
         ENDIF
C
C        NOP1 EST UNE (
C        RECHERCHE DE LA PREMIERE ) A PARTIR DE NOP1
         NOP2 = NOP1
C
 95      NOP2 = NCOPER(3,NOP2)
         IF( NCOPER(1,NOP2) .NE. 2 ) GOTO 95
C
C        NOP2 EST LA PREMIERE )
C        RECHERCHE DE LA PREMIERE ( QUI PRECEDE
         NOP1 = NOP2
C
 98      NOP1 = NCOPER(2,NOP1)
         IF( NCOPER(1,NOP1) .NE. 1 ) GOTO 98
C        NOP1 EST LA ( LA PLUS INTERNE
C        NOP2 EST LA ) CORRESPONDANTE
C
C        LE SUIVANT DE NOP1 (  , LE PREDECESSEUR DE NOP2 )
         NOP1 = NCOPER( 0,NOP1)
         NOP2 = NCOPER(-1,NOP2)
C        NOP1 SUIT ( , NOP2 PRECEDE )
C
C        FORMATION DE L'ARBRE ENTRE NOP1 ET NOP2 NON PARENTHESE
C        ======================================================
C
C        RECHERCHE DES OPERATEURS UNAIRES NON + ET - SANS ( )
 100     LPASS  = 0
         NBCRAC = 0
C        NOMBRE DE CREATION DE RACINES
C        CE QUI PRECEDE NOP1 ET SUIT NOP2  C-A-D LES ( )
         NPDBT = NCOPER(-1,NOP1)
         NPFIN = NCOPER( 0,NOP2)
C
C        TRAITEMENT DES OPERATEURS BINAIRES NON ( )
 105     LPASS = LPASS + 1
         GOTO ( 108, 110, 120, 130, 140, 145, 150, 160, 170, 300 ),LPASS
C
C        + - UNAIRE
 108     NCOPT1 = 102
         NCOPT2 = 103
         GOTO 200
C
C        **
 110     NCOPT1 = 214
         NCOPT2 = 214
         GOTO 200
C
C        */
 120     NCOPT1 = 212
         NCOPT2 = 213
         GOTO 200
C
C        +-  BINAIRE
 130     NCOPT1 = 210
         NCOPT2 = 211
         GOTO 200
C
C        <  <=  =  <>  >=  >
 140     NCOPT1 = 204
         NCOPT2 = 209
         GOTO 200
C
C        ET LOGIQUE
 145     NCOPT1 = 101
         NCOPT2 = 101
         GOTO 200
C
C        ET LOGIQUE
 150     NCOPT1 = 203
         NCOPT2 = 203
         GOTO 200
C
C        OU LOGIQUE
 160     NCOPT1 = 202
         NCOPT2 = 202
         GOTO 200
C
C        OX LOGIQUE
 170     NCOPT1 = 201
         NCOPT2 = 201
C
 200     NOPET = NCOPER(0,NPDBT)
         IF( NCOPER(1,NOPET) .NE. 1 ) GOTO 220
C
C        LA 1-ERE OPERATION
 210     NOPET = NCOPER(0,NOPET)
 220     IF( NOPET .GE. NPFIN ) THEN
C           TOUTES LES OPERATIONS ONT ETE VUES . PASSAGE AU SUIVANT
            GOTO 105
         ENDIF
C        LE TYPE DE L'OPERATION EST IL CELUI RECHERCHE ?
         IF( NCOPER(1,NOPET) .LT. NCOPT1 .OR.
     %       NCOPER(1,NOPET) .GT. NCOPT2 ) GOTO 210
C        OUI . IL EST AJOUTE A L'ARBRE
         IF( NCOPT1 .NE. 101 .AND. NCOPT1 .NE. 102 .AND.
     %       NCOPT1 .NE. 103 ) THEN
C           OPERATEUR BINAIRE : LES OPERANDES 1 ET 2
            NCOPD1 = NCOPER(-1,NOPET)
            NCOPD2 = NCOPER( 0,NOPET)
C           CES OPERANDES SONT ILS RACINES ?
            NP1 = NCOPER(2,NCOPD1)
            NP2 = NCOPER(2,NCOPD2)
            IF( NP1 .LT. 0 .AND. NP2 .LT. 0 ) THEN
C              LES 2 SONT DES RACINES > LA SECONDE EST RENDUE VIDE
               NP2 = - NP2
               LSRACI( NP2 ) = - LSRACI( 0 )
               LSRACI(  0  ) = NP2
C              LA RACINE NP1
               NP1 = - NP1
            ELSE IF( NP1 .GE. 0 .AND. NP2 .GE. 0 ) THEN
C              LES 2 OPERANDES NE SONT PAS DES RACINES. UNE RACINE EST CREE
               NP1 = LSRACI( 0 )
               IF( NP1 .EQ. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: TABLE DES RACINES SATUREE (AREXPD)'
                     KERR(2) = 'LU: TROP DE PARAMETRES'
                  ELSE
                     KERR(1) = 'LU: SATURATED TABLE OF ROOTS (AREXPD)'
                     KERR(2) = 'LU: TOO PARAMETERS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
               LSRACI( 0 ) = -LSRACI( NP1 )
            ELSE IF( NP2 .LT. 0 ) THEN
               NP1 = -NP2
            ELSE
               NP1 = -NP1
            ENDIF
C           L'OPERATION AVANT
            NV = NCOPER( -1 , NCOPD1 )
C           L'OPERATION APRES
            NP = NCOPER(  0 , NCOPD2 )
C           LE HAUT DE L'ARBRE LOCAL
            NBCRAC = NBCRAC + 1
            LSRACI( NP1 )    = NOPET
            NCOPER(-1,NOPET) = NV
            NCOPER(0 ,NV)    = NOPET
            NCOPER(0 ,NOPET) = NP
            NCOPER(-1,NP)    = NOPET
C           LA RACINE
            NCOPER(2,NOPET)  = -NP1
            NCOPER(3,NOPET)  = NCOPD1
            NCOPER(4,NOPET)  = NCOPD2
            NCOPER(2,NCOPD1) = NOPET
            NCOPER(2,NCOPD2) = NOPET
C
         ELSE
C
C           NON LOGIQUE OU + OU - EST UNAIRE D'OPERANDE 2
            NCOPD2 = NCOPER( 0,NOPET)
            NP2    = NCOPER(2,NCOPD2)
C           L'OPERANDE 2 EST IL UNE RACINE ?
            IF( NP2 .GE. 0 ) THEN
C              NON: UNE RACINE EST CREE
               NP1 = LSRACI( 0 )
               IF( NP1 .EQ. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: TABLE DES RACINES SATUREE (AREXPD)'
                     KERR(2) = 'LU: TROP DE PARAMETRES'
                  ELSE
                     KERR(1) = 'LU: SATURATED TABLE OF ROOTS (AREXPD)'
                     KERR(2) = 'LU: TOO PARAMETERS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
               LSRACI( 0 ) = -LSRACI( NP1 )
            ELSE
C              OUI:
               NP1 = -NP2
            ENDIF
C           LE HAUT DE L'ARBRE LOCAL
            NBCRAC = NBCRAC + 1
            LSRACI( NP1 )    = NOPET
C           L'OPERATION APRES
            NP = NCOPER(  0 , NCOPD2 )
            NCOPER(0 ,NOPET) = NP
            NCOPER(-1,NP)    = NOPET
C           LA RACINE
            NCOPER(2,NOPET)  = -NP1
            NCOPER(3,NOPET)  =  0
            NCOPER(4,NOPET)  = NCOPD2
            NCOPER(2,NCOPD2) = NOPET
         ENDIF
C        L'OPERATEUR EST RENDU NEGATIF POUR NE PAS LE RECONNAITRE
C        DANS LA RECHERCHE DES OPERATEURS
         NCOPER(1,NOPET) = -NCOPER(1,NOPET)
         GOTO 210
C
C        SI AUCUNE CREATION DE RACINE, IL N'Y A PAS D'OPERATEUR DANS
C        L'EXPRESSION QUI SE REDUIT A 1 OPERANDE
C        ===========================================================
 300     IF( NBCRAC .EQ. 0 ) THEN
C           CETTE OPERATION EST ELLE UNE RACINE?
            NOPET = NCOPER(0,NPDBT)
            IF( NCOPER(2,NOPET) .LT. 0 ) THEN
C              C'EST UNE RACINE . ELLE EST CONSERVEE
               NP1 = -NCOPER(2,NOPET)
            ELSE
C              CREATION D'UNE RACINE
               NP1 = LSRACI( 0 )
               IF( NP1 .LE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: RACINES SATUREES (AREXPD)'
                     KERR(2) = 'LU: TROP DE PARAMETRES'
                  ELSE
                     KERR(1) = 'LU: SATURATED ROOTS (AREXPD)'
                     KERR(2) = 'LU: TOO PARAMETERS'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
               LSRACI( 0 ) = - LSRACI( NP1 )
C              PAS DE PERE
               NCOPER( 2 , NOPET ) = - NP1
               LSRACI( NP1 ) = NOPET
            ENDIF
         ENDIF
C
C        L'EXPRESSION ( NOP1 ...  NOP2 ) EST DEVENUE UN ARBRE
C        DE RACINE NP1 . LES (  ) SONT SUPPRIMEES
C        ====================================================
         NOP = LSRACI( NP1 )
C        VERIFICATION QUE NPDBT EST LE PREDECESSEUR DE NOP
C                         NPFIN EST LE SUCCESSEUR   DE NOP
         IF( NCOPER(-1,NOP) .NE. NPDBT .OR.
     %       NCOPER( 0,NOP) .NE. NPFIN ) THEN
            WRITE(IMPRIM,*)
     %     'LU: AREXPD:PREDECESSEUR OU SUIVANT DE RACINE',
     %      NOP,NCOPER(-1,NOP),NCOPER(0,NOP),
     %    ' DIFFERENT DE ',NPDBT,NPFIN
            GOTO 9900
         ENDIF
C
C        EXISTE T IL DES ( AU DEBUT ET ) A LA FIN ?
         IF( NCOPER(1,NPDBT) .EQ. 1 .AND. NCOPER(1,NPFIN) .EQ. 2 ) THEN
C           OUI
C           LA ( PRECEDENTE DE ( NPDBT
            NPDBT0 = NCOPER(2,NPDBT)
C           LA ) SUIVANTE DE NPFIN )
            NPFIN1 = NCOPER(3,NPFIN)
C           ELIMINATION DE (NPDBT  NPFIN)
            NCOPER(3,NPDBT0) = NPFIN1
            NCOPER(2,NPFIN1) = NPDBT0
C           L'OPERATION PRECEDANT LE RESULTAT DE L'EXPRESSION
            NPDBT0 = NCOPER(-1,NPDBT)
            NCOPER( 0,NPDBT0) = NOP
            NCOPER(-1,NOP   ) = NPDBT0
C           L'OPERATION SUIVANT LE RESULTAT DE L'EXPRESSION
            NPFIN1 = NCOPER(0,NPFIN)
            NCOPER( 0,NOP   ) = NPFIN1
            NCOPER(-1,NPFIN1) = NOP
C           DESTRUCTION DU CODE DE ( )
            NCOPER(1,NPDBT) = 0
            NCOPER(1,NPFIN) = 0
         ENDIF
         GOTO 80
C
C        L'EXPRESSION EST CALCULEE
C        ESSAI DE LA CHAINER AVEC CE QUI PRECEDE OU SUIT
C
C        RECHERCHE DE FONCTION( PARA1, PARA2 , PARA3, ... PARA N)
 800     IF( NCOPER(1,NDEBUT) .EQ. 7 ) THEN
C           CETTE EXPRESSION EST LE 1-ER PARAMETRE D'UNE FONCTION DONT
C           TOUS LES AUTRES PARAMETRES ONT ETE CALCULES A CAUSE DU
C           SENS RETROGRADE DU CALCUL DES EXPRESSIONS
c
C           L'APPEL DE LA FONCTION EST MAINTENANT A ELIMINER
C           LES PARAMETRES SONT CHAINES ENTRE EUX COMME SUIT DANS L'ARBRE
C
C FONCTION UNAIRE  FONC1 --> PARA 1
C                        \-> 0
C FONCTION BINAIRE FONC2 --> PARA 2
C                        \>- PARA 1
C FONCTION N-AIRE  FONCN --> ENTPARA1 --> ENTPARA 2 ... ENTPARA N-2 -->PARA N
C                        \>- PARA 1   \-> PARA2                     \->PARA N-1
C
C           LA FONCTION
            NOPF  = NCOPER(-1,NDEBUT)
            NCOPF = NCOPER( 1,NOPF  )
            IF( NCOPF .NE. 6 .AND. NCOPF .LT. 100 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LU: L''OPERATEUR DEVRAIT ETRE UNE FONCTION'
               ELSE
                  KERR(1) = 'LU: THE OPERATOR SHOULD BE A FUNCTION'
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
C           LE NOMBRE DE PARAMETRES DECLARES DE LA FONCTION
            IF( NCOPF .LT. 200 ) THEN
               NBPARD = 1
            ELSE IF( NCOPF .LT. 300 ) THEN
               NBPARD = 2
            ELSE
C              FONCTION UTILISATEUR
               NBPARD = NCOPF / 1000
C              POUR MEMOIRE LE NUMERO DE LA FONCTION UTILISATEUR
C              NOFONC = MOD( NCOPF , 1000 ) - 300
            ENDIF
C           PROTECTION DE LA ( AVANT L'APPEL DE LA FONCTION
            NPAVFO = NCOPER(2,NDEBUT)
C
C           LE 1-ER PARAMETRE DE LA FONCTION EST LA DERNIERE RACINE NP1
C           INTRODUITE DANS LSRACI
            NOPAR1 = LSRACI( NP1 )
C           LE PERE DE LA FONCTION EST LA RACINE NP1
            NCOPER(2,NOPF) = -NP1
C           LA RACINE EST MISE A JOUR
            LSRACI(  NP1 ) = NOPF
C           LE PERE DU 1-ER PARAMETRE EST LA FONCTION
            NCOPER(2,NOPAR1) = NOPF
C
C           AFFECTATION DES FILS SELON LE NOMBRE DE PARAMETRES
            IF( NBPARD .EQ. 1 ) THEN
C
C              1 PARAMETRE = FILS DROIT DE LA FONCTION
C              =======================================
               NCOPER(3,NOPF) = 0
               NCOPER(4,NOPF) = NOPAR1
C              LA PARENTHESE FINALE
               NOVIR2 = NCOPER(3,NDEBUT)
C
            ELSE IF( NBPARD .EQ. 2 ) THEN
C
C              2 PARAMETRES = FILS GAUCHE ET DROIT DE LA FONCTION
C              ==================================================
C              RECHERCHE DU 2-EME PARAMETRE
C              LA VIRGULE ENTRE LES 2 PARAMETRES
               NOVIR1 = NCOPER(3,NDEBUT)
               IF( NCOPER(1,NOVIR1) .NE. 8 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: PAS DE '','' '
                  ELSE
                     KERR(1) = 'LU: FORGOTTEN '','' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LE 2-EME PARAMETRE SUIT LA ,
               NOPAR2 = NCOPER(0,NOVIR1)
               IF( NOPAR2 .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU:PAS D''OPERANDE DERRIERE '','' '
                  ELSE
                     KERR(1) = 'LU:NO OPERAND AFTER '','' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA RACINE DE CE PARAMETRE
               NP2 = NCOPER(2,NOPAR2)
               IF( NP2 .GE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU:LE PARAMETRE DEVRAIT ETRE UNE RACINE'
                     KERR(2) = 'LU:REVOIR LES PARAMETRES D''APPEL'
                  ELSE
                     KERR(1) = 'LU:THE PARAMETER SHOULD BE A ROOT'
                     KERR(2) = 'LU:SEE THE PARAMETERS OF THE CALL'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA RACINE DE NOPAR2 EST RENDUE A LSRACI
               NP2 = -NP2
               LSRACI(NP2) = -LSRACI(0)
               LSRACI( 0 ) = NP2
C              LE CHAINAGE PERE FILS GAUCHE FILS DROIT
               NCOPER(2,NOPAR2) = NOPF
               NCOPER(3,NOPF  ) = NOPAR1
               NCOPER(4,NOPF  ) = NOPAR2
C              LA PARENTHESE FINALE DE L'APPEL DE LA FONCTION
               NOVIR2 = NCOPER(3,NOVIR1)
C
            ELSE
C
C              FONCTION AVEC PLUS DE 2 PARAMETRES
C              ==================================
C              LE FILS GAUCHE DE LA FONCTION
               NCOPER(3,NOPF) = NOPAR1
C              LA FONCTION JOUE LE ROLE D'UNE VIRGULE
               NOVIR0 = NOPF
C              LA 1-ERE VIRGULE
               NOVIR1 = NCOPER(3,NDEBUT)
               IF( NCOPER(1,NOVIR1) .NE. 8 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: '','' OUBLIEE'
                  ELSE
                     KERR(1) = 'LU: FORGOTTEN '','' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA 2-EME VIRGULE
               NOVIR1 = NCOPER(3,NOVIR1)
               IF( NCOPER(1,NOVIR1) .NE. 8 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU: '','' OUBLIEE'
                  ELSE
                     KERR(1) = 'LU: FORGOTTEN '','' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LE NOMBRE DE PARAMETRES
               NBPARA = 2
C
C              LA BOUCLE SUR LES PARAMETRES 2 A NBPARD-1
C              -----------------------------------------
 900           NBPARA = NBPARA + 1
               IF( NBPARA .GT. NBPARD ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                   KERR(1)='LU:FONCTION APPELEE AVEC TROP DE PARAMETRES'
                  ELSE
                   KERR(1)='LU:TOO ARGUMENTS FOR THE CALLED FUNCTION'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LE PARAMETRE AVANT LA ,
               NOPAR2 = NCOPER(-1,NOVIR1)
               IF( NOPAR2 .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) ='LU:PAS D''OPERANDE DERRIERE '','' '
                  ELSE
                     KERR(1) ='LU:NO OPERAND AFTER '','' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA RACINE DE CE PARAMETRE
               NP2 = NCOPER(2,NOPAR2)
               IF( NP2 .GE. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'LU:LE PARAMETRE DEVRAIT ETRE UNE RACINE'
                     KERR(2) = 'LU:REVOIR LES PARAMETRES D''APPEL'
                  ELSE
                     KERR(1) ='LU:THE PARAMETER SHOULD BE A ROOT'
                     KERR(2) ='LU:SEE THE PARAMETERS OF THE CALL'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA RACINE DE NOPAR2 EST RENDUE A LSRACI
               NP2 = -NP2
               LSRACI(NP2) = -LSRACI(0)
               LSRACI( 0 ) = NP2
C
C              PROTECTION DE L'EVENTUELLE VIRGULE SUIVANTE
               NOVIR2 = NCOPER(3,NOVIR1)
C
C              LA , DEVIENT ENTRE_PARAMETRES CODE 200
               NCOPER(1,NOVIR1) = 200
C              SON PERE
               NCOPER(2,NOVIR1) = NOVIR0
C              LE FILS DROIT DE NOVIR0 EST L'ENTRE_PARAMETRES PRECEDENT
               NCOPER(4,NOVIR0) = NOVIR1
C              LE FILS GAUCHE DE L'ENTRE_PARAMETRES 1 EST LE PARAMETRE 2
               NCOPER(3,NOVIR1) = NOPAR2
C              LE PERE DU PARAMETRE 2 EST L'ENTRE_PARAMETRES 1
               NCOPER(2,NOPAR2) = NOVIR1
C              EXISTE T IL UNE AUTRE , ET PARAMETRE ?
               IF( NCOPER(1,NOVIR2) .EQ. 8 ) THEN
C                 OUI: PASSAGE AU PARAMETRE SUIVANT
                  NOVIR0 = NOVIR1
                  NOVIR1 = NOVIR2
                  GOTO 900
               ENDIF
C
C              FIN DE LA BOUCLE SUR LES PARAMETRES
C              -----------------------------------
C              NOVIR2 DOIT ETRE LA ) DE L'APPEL DE LA FONCTION
               IF( NCOPER(1,NOVIR2) .NE. 9 .AND.
     %             NCOPER(1,NOVIR2) .NE. 8 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) ='LU: PAS DE ) DANS APPEL DE FONCTION'
                  ELSE
                     KERR(1) ='LU: NO '')'' IN THE FUNCTION CALL'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C
C              LE DERNIER PARAMETRE
               NOPAR3 = NCOPER(-1,NOVIR2)
C
C              LA RACINE DE CE PARAMETRE
               NP2 = NCOPER(2,NOPAR3)
               IF( NP2 .GE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) ='LU:LE PARAMETRE DEVRAIT ETRE UNE RACINE'
                  ELSE
                     KERR(1) ='LU:THE PARAMETER SHOULD BE A ROOT'
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              LA RACINE DE NOPAR2 EST RENDUE A LSRACI
               NP2 = -NP2
               LSRACI(NP2) = -LSRACI(0)
               LSRACI( 0 ) = NP2
C
C              SON PERE EST L'ENTRE_PARAMETRES PRECEDENT
               NCOPER(2,NOPAR3) = NOVIR1
C              LE FILS DROIT DE L'ENTRE_PARAMETRES EST LE DERNIER PARAMETRE
               NCOPER(4,NOVIR1) = NOPAR3
            ENDIF
C
C           NOVIR2 DOIT ETRE LA ) DE L'APPEL DE LA FONCTION
            IF( NCOPER(1,NOVIR2) .NE. 9 .AND.
     %          NCOPER(1,NOVIR2) .NE. 8 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='LU: PAS DE ) DANS L''APPEL DE LA FONCTION'
               ELSE
                  KERR(1) ='LU: NO '')'' IN THE FUNCTION CALL'
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
C
C           LA ( OU , OU  OU ) APRES ) DE LA FONCTION
            NPAPFO = NCOPER(3,NOVIR2)
            IF( NPAPFO .GT. 0 .AND. NPAVFO .GT. 0 ) THEN
C              TOUTES LES ( ) INTERNES A L'APPEL DE LA FONCTION
C              SONT ELIMINEES
               NCOPER(3,NPAVFO) = NPAPFO
               NCOPER(2,NPAPFO) = NPAVFO
            ENDIF
C
C           LE SUIVANT DERRIERE LA ) FINALE DE LA FONCTION
            NOP2 = NCOPER(0,NOVIR2)
C           LE SUIVANT DE LA FONCTION EST LE SUIVANT DE LA ) FINALE
            NCOPER( 0,NOPF) = NOP2
            NCOPER(-1,NOP2) = NOPF
C           FIN DE L'EXPRESSION  LEXPD
         ENDIF
 1000 CONTINUE
C
C     VERIFICATION : RESTE T IL UNE SEULE RACINE  ?
C                    ET MISE DANS LSRACI(1) DE CELLE CI
      NBCRAC = 0
      DO 2000 I=1,MXRACI
         IF( LSRACI(I) .GT. 0 ) THEN
            NBCRAC = NBCRAC + 1
            NOP    = I
         ENDIF
 2000 CONTINUE
C
      IF( NBCRAC .NE. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: VOIR les PARAMETRES D''APPEL'
         ELSE
            KERR(1) = 'LU: SEE the INCORRECT PARAMETERS'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     LE CODE OPERATION REDEVIENT POSITIF
      DO 3000 I=LSEXPD(1,1),LSEXPD(2,1)
         NCOPER(1,I) = ABS( NCOPER(1,I) )
 3000 CONTINUE
C
C     PASSAGE DE LA RACINE EN LSRACI(1)
      LSRACI( 1 ) = LSRACI( NOP )
      NCOPER( 2 , LSRACI(NOP) ) = -1
      NRETOU = 0
C
      IF( LUIMPR .GE. 10 ) THEN
         WRITE(IMPRIM,*)'LU: SORTIE DE AREXPD ============='
         WRITE(IMPRIM,18000)(J,(NCOPER(I,J),I=1,4),J=1,LSEXPD(2,1))
      ENDIF
18000 FORMAT(' NCOPER(',I3,')=',4I4)
      RETURN
C
C     ERREUR
 9900 NRETOU = 1
      RETURN
      END
