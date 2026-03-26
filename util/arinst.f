      SUBROUTINE ARINST( NLMIN , NCMIN , NLMAX , NCMAX ,
     %                   MXINST, NCINST, NOINST, NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : ANALYSE LEXICALE ET GENERATION DE L'ARBRE DES INSTRUCTIONS
C ----- D'UNE FONCTION C-A-D LA TRADUCTION DES CARACTERES LUS
C       EN UNE SUITE CHAINEE D'OPERATIONS
C
C        {ETIQ :} INSTRUCTION
C       {{ETIQ :} INSTRUCTION }
C        {ETIQ :} FINFONC ;
C
C        INSTRUCTION := | NOM_FONC = EXP_ARITH ;
C                       | NOM_VAR = EXP_ARITH ;
C                       | SI EXP_LOGIQ ALORS INSTRUCTION
C                                     {SINON INSTRUCTION}
C                                      FINSI ;
C                       | POUR NOM_VAR = EXP_ARITH A EXP_ARITH ;
C                            INSTRUCTION;
C                         FINPOUR ;
C                       | REPETER ;
C                           INSTRUCTION
C                         JUSQUE  EXP_LOGIQ ;
C                       | TANTQUE EXP_LOGIQ ;
C                            INSTRUCTION
C                         FINTANT ;
C                       | RETOUR ;
C                       | ALLER ETIQ ;
C
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NLMIN,NCMIN : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C NLMAX,NCMAX : POSITION DANS KLG DU DERNIER CARACTERE A TRAITER
C MXINST      : NOMBRE MAXIMAL D'INSTRUCTIONS
C
C MODIFIES :
C ----------
C NCINST : LES INSTRUCTIONS
C
C SORTIES :
C ---------
C NOINST : NUMERO DE LA DERNIERE INSTRUCTION  DANS NCINST
C NRETOU : =0 SI PAS D'ERREUR
C          >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      PARAMETER        (MXPISI=64)
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NCINST(1:4,MXINST),LAPISI(MXPISI)
      INTEGER           NCETIQ(1:MXETIQ)
C     LA POSITION DANS NCOPER DE CHAQUE ETIQUETTE LORS DE L'ANALYSE SYNTAXIQUE
      DOUBLE PRECISION  DBLVAL
C
C     LES INITIALISATIONS
C     ===================
      NBSI   = 0
      NFINSI = 0
C     LA HAUTEUR DE LA PILE DES INSTRUCTIONS
      LHPISI = 1
      NOINST = 0
C     LA PREMIERE INSTRUCTION DE LA FONCTION EST EN_REMONTANT ALLER A
      LAPISI( LHPISI ) = 0
C     TOUS LES PERES SONT MIS A ZERO
      DO 5 I=1,MXINST
         NCINST(2,I) = 0
 5    CONTINUE
C     LE NUMERO D'INSTRUCTION DES ETIQUETTES EST MIS A ZERO
      DO 10 I=1,MXETIQ
         NCETIQ( I ) = 0
 10   CONTINUE
C     NLMIN,NCMIN CARACTERE AVANT L'INSTRUCTION A ANALYSER
C
C     LES INSTRUCTIONS 'POUR' 'REPETER' 'TANTQUE' SONT REMPLACEES
C     PAR LEUR EQUIVALENT A L'AIDE DE SI ET =
C     =================================================================
      CALL REPETE( NLMIN, NCMIN, NLMAX, NCMAX, NRETOU )
      IF( NRETOU .NE. 0 ) GOTO 9900
      CALL TANTQU( NLMIN, NCMIN, NLMAX, NCMAX, NRETOU )
      IF( NRETOU .NE. 0 ) GOTO 9900
      CALL POURVA( NLMIN, NCMIN, NLMAX, NCMAX, NRETOU )
      IF( NRETOU .NE. 0 ) GOTO 9900
C
      NL = NLMIN
      NC = NCMIN - 1
C
C     =================================
C     BOUCLE D'ANALYSE DES INSTRUCTIONS
C     =================================
C     SI CARACTERE BLANC PASSAGE AU PROCHAIN CARACTERE NON BLANC
 100  CALL CARPNB( NL, NC )
C
C     L'INSTRUCTION SUIVANTE DEBUTE T ELLE PAR ETIQ : ?
C     =================================================================
      CALL CHETIQ( NL, NC, NLF, NCF, NOETIQ )
      IF( NOETIQ .GT. 0 ) THEN
C        RECHERCHE DU : SUIVANT L'ETIQUETTE
         CALL CARPNB( NLF, NCF )
         IF(  KLG(NLF)(NCF:NCF) .EQ. ':' ) THEN
C           LE : EXISTE PASSAGE AU CARACTERE SUIVANT
            CALL CARPNB( NLF, NCF )
         ELSE
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: ETIQUETTE ' //KETIQ(NOETIQ)//
     %                   ' SANS : FINAL'
            ELSE
               KERR(1) = 'LU: LABEL ' //KETIQ(NOETIQ)//
     %                   ' WITHOUT FINAL :'
            ENDIF
            CALL LEREUR
         ENDIF
C        L'ETIQUETTE EXISTE T ELLE PAR AILLEURS ?
         IF( NCETIQ(NOETIQ) .NE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: ETIQUETTE DOUBLEE'
            ELSE
               KERR(1) = 'LU: DOUBLE LABEL'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C        LE NUMERO NCINST DE L'INSTRUCTION SUIVANT CETTE ETIQUETTE
         NCETIQ( NOETIQ ) = NOINST + 1
C        SAUT DE ETIQ:
         NL = NLF
         NC = NCF
      ENDIF
C
C     L'INSTRUCTION SANS ETIQUETTE
C     =================================================================
      IF( KLG(NL)(NC:NC+6) .EQ. 'FINFONC' .OR.
     %    KLG(NL)(NC:NC+6) .EQ. 'ENDFUNC' ) THEN
C        ==============================================================
C        FIN DU TRAITEMENT
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC+5) .EQ. 'RETOUR'  .OR.
     %         KLG(NL)(NC:NC+5) .EQ. 'RETURN') THEN
C        ==============================================================
C        RETOUR <=> ALLER  ETIQ: FINFONC  AVEC L'ETIQUETTE DE NUMERO:-1
C        RECHERCHE DE ; DERRIERE RETOUR
         NLF = NL
         NC  = NC + 5
         NCF = NC
         CALL CARPNB( NLF, NCF )
         IF( KLG(NLF)(NCF:NCF) .EQ. ';' ) THEN
            NL = NLF
            NC = NCF
         ENDIF
C        L'OPERATION EST ENREGISTREE
         IF( LAPISI(LHPISI) .EQ. 0 ) THEN
C           UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
            IF( NOINST .GE. MXINST ) GOTO 9200
            NOINST = NOINST + 1
            NCINST(1,NOINST) = 174
C           PAS DE PERE
            NCINST(2,NOINST) = 0
C           L'INSTRUCTION A EXECUTER AVANT REMONTEE
            NCINST(3,NOINST) = 0
C           L'INSTRUCTION ALLER A EN REMONTANT OU INSTRUCTION SUIVANTE
            NCINST(4,NOINST) = 0
C           CETTE INSTRUCTION EST L'INSTRUCTION SUITE
            LAPISI(LHPISI) = NOINST
         ENDIF
C
         IF( NOINST .GE. MXINST ) GOTO 9200
         NOINST = NOINST + 1
C        RETOUR <=> EN_REMONTANT_ALLER_A
         NCINST(1,NOINST) = 13
C        L'INSTRUCTION SUIVANTE EST NOINST
         NCINST(3,LAPISI(LHPISI)) = NOINST
C        LE PERE DE ALLER
         NCINST(2,NOINST) = LAPISI(LHPISI)
C        IMPASSE
         LAPISI(LHPISI) = 0
C        LE NUMERO DE L'ETIQUETTE QUI DEVIENDRA L'ADRESSE DE FINFONC
         NCINST(3,NOINST) = -1
C        PAS DE FILS DROIT
         NCINST(4,NOINST) = 0
C        PASSAGE A L'INSTRUCTION SUIVANTE
         GOTO 100
C
      ELSE IF( KLG(NL)(NC:NC+4) .EQ. 'ALLER' .OR.
     %         KLG(NL)(NC:NC+4) .EQ. 'GOTO ' ) THEN
C        ===============================================================
C        RECHERCHE DE L'ETIQUETTE QUI SUIT
         NLF = NL
         NCF = NC + 4
         CALL CARPNB( NLF, NCF )
         NL = NLF
         NC = NCF
         CALL CHETIQ( NL, NC, NLF, NCF, NOETIQ )
         IF( NOETIQ .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: ALLER SANS ETIQUETTE DECLAREE'
            ELSE
               KERR(1) = 'LU: GOTO WITHOUT DECLARED LABEL'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C        RECHERCHE DE ; DERRIERE   ALLER ETIQUETTE
         NL = NLF
         NC = NCF
         CALL CARPNB( NL, NC )
         IF( KLG(NL)(NC:NC) .EQ. ';' ) THEN
            NLF = NL
            NCF = NC
         ENDIF
C
C        L'OPERATION EST ENREGISTREE
         IF( LAPISI(LHPISI) .EQ. 0 ) THEN
C           UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
            IF( NOINST .GE. MXINST ) GOTO 9200
            NOINST = NOINST + 1
            NCINST(1,NOINST) = 174
C           PAS DE PERE
            NCINST(2,NOINST) = 0
C           L'INSTRUCTION A EXECUTER AVANT REMONTEE
            NCINST(3,NOINST) = 0
C           L'INSTRUCTION ALLER A EN REMONTANT OU INSTRUCTION SUIVANTE
            NCINST(4,NOINST) = 0
C           CETTE INSTRUCTION EST L'INSTRUCTION SUITE
            LAPISI(LHPISI) = NOINST
         ENDIF
         IF( NOINST .GE. MXINST ) GOTO 9200
         NOINST = NOINST + 1
C        EN_REMONTANT_ALLER_A
         NCINST(1,NOINST) = 13
C        L'INSTRUCTION SUIVANTE EST NOINST
         NCINST(4,LAPISI(LHPISI)) = NOINST
C        LE PERE DE ALLER
         NCINST(2,NOINST) = LAPISI(LHPISI)
C        IMPASSE
         LAPISI(LHPISI) = 0
C        LE NUMERO DE L'ETIQUETTE
         NCINST(3,NOINST) = NOETIQ
C        PAS DE FILS DROIT
         NCINST(4,NOINST) = 0
C        PASSAGE A L'INSTRUCTION SUIVANTE
         NL = NLF
         NC = NCF
         GOTO 100
C
      ELSE IF( KLG(NL)(NC:NC+7) .EQ. 'AFFICHER' .OR.
     %         KLG(NL)(NC:NC+7) .EQ. 'DISPLAY ' ) THEN
C        ===============================================================
C        RECHERCHE DU ; FINAL DERRIERE LES PARAMETRES A AFFICHER
         CALL CHMTLU( ';', NL,NC+7, NLMAX,NCMAX, NLPTVI,NCPTVI )
         IF( NLPTVI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: AFFICHER SANS '';'' FINAL'
            ELSE
               KERR(1) = 'LU: DISPLAY WITHOUT FINAL '';'' '
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C        AFFICHER PARAM {,PARAM} ; EST TRAITE COMME UN APPEL DE FONCTION
         CALL ALEXPD( NL,NC, NLPTVI,NCPTVI, NCODEV, DBLVAL, NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9900
C        GENERATION DE L'ARBORESCENCE
         CALL AREXPD( NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9900
C
C        ICI AFFICHER ... ; EST CORRECT : L'OPERATION EST ENREGISTREE
         IF( LAPISI(LHPISI) .EQ. 0 ) THEN
C           UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
            IF( NOINST .GE. MXINST ) GOTO 9200
            NOINST = NOINST + 1
            NCINST(1,NOINST) = 174
C           PAS DE PERE
            NCINST(2,NOINST) = 0
C           L'INSTRUCTION A EXECUTER AVANT REMONTEE
            NCINST(3,NOINST) = 0
C           L'INSTRUCTION ALLER A EN REMONTANT OU INSTRUCTION SUIVANTE
            NCINST(4,NOINST) = 0
C           CETTE INSTRUCTION EST L'INSTRUCTION SUITE
            LAPISI(LHPISI) = NOINST
         ENDIF
C
C        LE PERE DE AFFICHER
         NI0 = LAPISI(LHPISI)
C        LA DERNIERE INSTRUCTION DANS NCOPER DE L'EXPRESSION
         I   = LSEXPD(2,1)
         IF( NOINST+I+1 .GT. MXINST ) GOTO 9200
C        UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
         IF( NOINST .GE. MXINST ) GOTO 9200
         NOINST = NOINST + 1
         NCINST(1,NOINST) = 174
C        LE PERE ET LE FILS
         NCINST(2,NOINST) = NI0
         NCINST(3,NI0)    = NOINST
C        L'INSTRUCTION A EXECUTER AVANT REMONTEE
         NCINST(3,NOINST) = 0
C        L'INSTRUCTION ALLER A EN REMONTANT OU INSTRUCTION SUIVANTE
         NCINST(4,NOINST) = 0
C        L'INSTRUCTION SUIVANTE APRES AFFICHER
         LAPISI(LHPISI)   = NOINST
C
C        COPIE DE L'EXPRESSION ARITHMETIQUE DANS NCINST
         NI1 = NOINST + LSRACI(1)
         CALL COEXPD( I, NOINST, NCINST )
C        LE PERE DE L'EXPRESSION ARITHMETIQUE NI1 EST NI0
         NCINST(2,NI1) = NI0
C        AFFICHER A EXECUTER
         NCINST(4,NI0) = NI1
C
C        PASSAGE A L'INSTRUCTION SUIVANTE
         NL = NLPTVI
         NC = NCPTVI
         GOTO 100
C
      ELSE IF( KLG(NL)(NC:NC+2) .EQ. 'SI ' .OR.
     %         KLG(NL)(NC:NC+2) .EQ. 'IF ' ) THEN
C        ===============================================================
         NBSI = NBSI + 1
C        RECHERCHE DE ALORS ( OU FIN EXP_LOGIQ )
         CALL CHMTL2( 'ALORS', 'THEN', NL, NC+2, NLMAX, NCMAX,
     %                NLALOR, NCALOR )
         IF( NLALOR .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: SI SANS ALORS'
            ELSE
               KERR(1) = 'LU: IF WITHOUT THEN'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C
C        L'EXPRESSION LOGIQUE EST ANALYSEE ET SON ARBRE GENERE
C        -----------------------------------------------------
         NCALOR = NCALOR - 1
         CALL ALEXPD( NL    , NC+2  , NLALOR, NCALOR ,
     %                NCODEV, DBLVAL, NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9900
         CALL CARPNB( NLALOR, NCALOR )
C        A SUPPRIMER APRES MISE AU POINT
         IF( KLG(NLALOR)(NCALOR:NCALOR+4) .NE. 'ALORS' .AND.
     %       KLG(NLALOR)(NCALOR:NCALOR+3) .NE. 'THEN' ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LU: PB AVEC ALORS ET SON DECALAGE'
            ELSE
               KERR(1) = 'LU: PB WITH ''THEN'' '
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
         IF( NCODEV .GT. 0 ) THEN
C           L'EXPRESSION SE REDUIT A UNE CONSTANTE
C           GENERATION D'UN ARBRE COMPATIBLE POUR L'EXECUTION
            LSRACI(1) = 1
C           LA CONSTANTE EXISTE-ELLE DANS LA TABLE DES CONSTANTES?
C           SI ELLE N'EXISTE PAS ELLE Y EST AJOUTEE AVEC
C           LE NUMERO NOCTE
            CALL NUMCTE( DBLVAL, NOCTE )
            IF( NOCTE .EQ. 0 ) GOTO 9300
            NCOPER(1,1) = 3
            NCOPER(2,1) = 0
            NCOPER(3,1) = NOCTE
            NCOPER(4,1) = 0
            LSEXPD(2,1) = 1
         ENDIF
C        GENERATION DE L'ARBORESCENCE DE L'EXPRESSION LOGIQUE
         CALL AREXPD( NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9900
C
C        GENERATION DU SI ...
C        ---------------------
         IF( LAPISI(LHPISI) .EQ. 0 ) THEN
C           UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
            IF( NOINST .GE. MXINST ) GOTO 9200
            NOINST = NOINST + 1
            NCINST(1,NOINST) = 174
C           PAS DE PERE
            NCINST(2,NOINST) = 0
C           L'INSTRUCTION A EXECUTER AVANT REMONTEE
            NCINST(3,NOINST) = 0
C           L'INSTRUCTION ALLER A EN REMONTANT OU INSTRUCTION SUIVANTE
            NCINST(4,NOINST) = 0
C           CETTE INSTRUCTION EST L'INSTRUCTION SUITE
            LAPISI(LHPISI) = NOINST
         ENDIF
C
         NOINST = NOINST + 4
         IF( NOINST .GT. MXINST ) GOTO 9200
C        L'INSTRUCTION GRAND PERE
         NI0 = LAPISI( LHPISI )
C        L'INSTRUCTION PERE FINSI
         NI1 = NOINST - 3
C        L'INSTRUCTION SI
         NI2 = NOINST - 2
C        L'INSTRUCTION EXP_LOGIQ
         NI3 = LSRACI(1) + NOINST
C        L'INSTRUCTION ALORS_SINON
         NI4 = NOINST - 1
C        L'INSTRUCTION DEBUT DU ALORS
         NI5 = NOINST
C        COPIE DE L'EXPRESSION LOGIQUE DANS NCINST
         I = LSEXPD(2,1)
         IF( NOINST+I .GT. MXINST ) GOTO 9200
         CALL COEXPD( I, NOINST, NCINST )
C
C        LE CHAINAGE AVEC LE GRAND PERE
         NCINST(4,NI0) = NI1
C
C        FINSI OU EN _REMONTANT_ALLER
         NCINST(1,NI1) = 174
         NCINST(2,NI1) = NI0
         NCINST(3,NI1) = 0
         NCINST(4,NI1) = NI2
C
C        SI
         NCINST(1,NI2) = 175
         NCINST(2,NI2) = NI1
         NCINST(3,NI2) = NI4
         NCINST(4,NI2) = NI3
C
C        EXPRESSION_LOGIQUE
         NCINST(2,NI3) = NI2
C
C        ALORS_SINON
         NCINST(1,NI4) = 176
         NCINST(2,NI4) = NI2
C        L'INSTRUCTION SINON A PRIORI INEXISTANTE
         NCINST(3,NI4) = 0
C        L'INSTRUCTION ALORS
         NCINST(4,NI4) = NI5
C
C        L'INSTRUCTION SUIVANTE DU ALORS
         NCINST(1,NI5) = 174
         NCINST(2,NI5) = NI4
         NCINST(3,NI5) = 0
         NCINST(4,NI5) = 0
C
C        LAPISI(LHPISI) POINTE SUR FINSI ( CODE 174 )
         LAPISI(LHPISI) = NI1
C
C        TRAITEMENT DU ALORS
         LHPISI = LHPISI + 2
         IF( LHPISI .GT. MXPISI ) GOTO 9400
C
C        SONT EMPILES : L'INSTRUCTION ALORS_SINON
C                       LE NUMERO DE LA PREMIERE INSTRUCTION DU ALORS
C
C        L'INSTRUCTION SUIVANTE DU ALORS TERMINE EST CELLE DU ALORS_SINON
         LAPISI( LHPISI - 1 ) = NI4
C
C        LA PREMIERE INSTRUCTION APRES ALORS
         LAPISI( LHPISI ) = NI5
C
C        PASSAGE A L'INSTRUCTION SUIVANTE
         NL = NLALOR
         NC = NCALOR + 4
         GOTO 100
C
      ELSE IF( KLG(NL)(NC:NC+4) .EQ. 'SINON' .OR.
     %         KLG(NL)(NC:NC+4) .EQ. 'ELSE '  ) THEN
C        ===============================================================
C        FIN DES INSTRUCTIONS DU ALORS  . PAS D'INSTRUCTION SUIVANTE
C        LE SINON DU ALORS_SINON
         NOINST = NOINST + 1
         IF( NOINST .GT. MXINST ) GOTO 9200
C        L'INSTRUCTION GRAND PERE
         NI0 = LAPISI( LHPISI - 1 )
C        L'INSTRUCTION SUITE OU PREMIERE INSTRUCTION DU SINON
         NI1 = NOINST
C
C        ALORS_SINON
         NCINST(3,NI0) = NI1
C
C        L'INSTRUCTION SUIVANTE DU SINON
         NCINST(1,NI1) = 174
         NCINST(2,NI1) = NI0
         NCINST(3,NI1) = 0
         NCINST(4,NI1) = 0
C
C        LA PREMIERE INSTRUCTION DU SINON
         LAPISI( LHPISI ) = NI1
C
C        PASSAGE A L'INSTRUCTION SUIVANTE
         NC = NC + 4
         GOTO 100
C
C
      ELSE IF( KLG(NL)(NC:NC+4) .EQ. 'FINSI' .OR.
     %         KLG(NL)(NC:NC+4) .EQ. 'ENDIF' ) THEN
C        ===============================================================
         NFINSI = NFINSI + 1
         NOINST = NOINST + 1
         IF( NOINST .GT. MXINST ) GOTO 9200
C        L'INSTRUCTION GRAND PERE AVANT SI
         LHPISI = LHPISI - 2
         NI0 = LAPISI( LHPISI )
C        L'INSTRUCTION SUITE DU FINSI
         NI1 = NOINST
C
C        FINSI
         NCINST(3,NI0) = NI1
C
C        L'INSTRUCTION SUIVANTE DU SINON
         NCINST(1,NI1) = 174
         NCINST(2,NI1) = NI0
         NCINST(3,NI1) = 0
         NCINST(4,NI1) = 0
C
C        L'INSTRUCTION SUIVANTE DERRIERE FINSI
         LAPISI( LHPISI ) = NI1
C
C        PASSAGE A L'INSTRUCTION SUIVANTE
         NC = NC + 4
C        RECHERCHE DU ; SUIVANT
         NLF = NL
         NCF = NC
         CALL CARPNB( NLF, NCF )
         IF( KLG(NLF)(NCF:NCF) .EQ. ';' ) THEN
            NL = NLF
            NC = NCF
         ENDIF
         GOTO 100
C
      ELSE
C
C        NOM_VAR LOCALE OU GLOBALE OU NOM_FONC SUIVI DE = EXP_ARITH ?
C        ============================================================
         CALL CHVARU( NL, NC, NLF, NCF, NOVARU )
         IF( NOVARU .GT. 0 ) THEN
C           LE CARACTERE SUIVANT EST IL = ?
            CALL CARPNB( NLF, NCF )
            IF( KLG(NLF)(NCF:NCF) .EQ. '=' ) THEN
C
C              ICI VARIABLE = .  RECHERCHE D'UNE EXPRESSION ARITHMETIQUE
               CALL CARPNB( NLF, NCF )
C              RECHERCHE DU PROCHAIN ;
               CALL CHMTLU( ';', NLF, NCF, NLMAX, NCMAX ,
     %                      NLPTVI, NCPTVI )
               IF( NLPTVI .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) =
     %           'LU: EXPRESSION ARITHMETIQUE NON TERMINEE PAR ;'
                  ELSE
                     KERR(1)='LU: ARITHMETICAL EXPRESSION NOT FINISHED b
     %y '' '
                  ENDIF
                  CALL LEREUR
                  GOTO 9900
               ENDIF
               CALL ALEXPD( NLF   , NCF   , NLPTVI , NCPTVI ,
     %                      NCODEV, DBLVAL, NRETOU )
               IF( NRETOU .NE. 0 ) GOTO 9900
C
               IF( NCODEV .GT. 0 ) THEN
C                 L'EXPRESSION SE REDUIT A UNE CONSTANTE
C                 GENERATION D'UN ARBRE COMPATIBLE POUR L'EXECUTION
                  LSRACI(1) = 1
C                 LA CONSTANTE EXISTE-ELLE DANS LA TABLE DES CONSTANTES?
C                 SI ELLE N'EXISTE PAS ELLE Y EST AJOUTEE AVEC
C                 LE NUMERO NOCTE
                  CALL NUMCTE( DBLVAL, NOCTE )
                  IF( NOCTE .EQ. 0 ) GOTO 9300
                  NCOPER(1,1) = 3
                  NCOPER(2,1) = 0
                  NCOPER(3,1) = NOCTE
                  NCOPER(4,1) = 0
                  LSEXPD(2,1) = 1
               ENDIF
C
C              GENERATION DE L'ARBORESCENCE
               CALL AREXPD( NRETOU )
               IF( NRETOU .NE. 0 ) GOTO 9900
C
C              ICI NOM_VAR OU NOM_FONC = EXP_ARITH EST CORRECT
C              -----------------------------------------------
C              GENERATION DE L'OPERATION
               IF( LAPISI(LHPISI) .EQ. 0 ) THEN
C                 UNE INSTRUCTION 174 EN_REMONTANT_ALLER_A EST CREEE
                  IF( NOINST .GE. MXINST ) GOTO 9200
                  NOINST = NOINST + 1
                  NCINST(1,NOINST) = 174
C                 PAS DE PERE
                  NCINST(2,NOINST) = 0
C                 L'INSTRUCTION A EXECUTER AVANT REMONTEE
                  NCINST(3,NOINST) = 0
C                 L'INSTRUCTION ALLER A EN REMONTANT
C                 OU INSTRUCTION SUIVANTE
                  NCINST(4,NOINST) = 0
C                 CETTE INSTRUCTION EST L'INSTRUCTION SUITE
                  LAPISI(LHPISI) = NOINST
               ENDIF
C
               NOINST = NOINST + 3
               IF( NOINST .GT. MXINST ) GOTO 9200
C
C              L'INSTRUCTION GRAND PERE
               NI0 = LAPISI( LHPISI )
C              L'INSTRUCTION PERE =
               NI1 = NOINST - 2
C              L'INSTRUCTION FILS GAUCHE  VARIABLE DU =
               NI2 = NOINST - 1
C              L'INSTRUCTION SUIVANTE
               NI3 = NOINST
C
C              COPIE DE L'EXPRESSION ARITHMETIQUE DANS NCINST
               I = LSEXPD(2,1)
               IF( NOINST+I .GT. MXINST ) GOTO 9200
               CALL COEXPD( I, NOINST, NCINST )
C              LE PERE DE L'EXPRESSION ARITHMETIQUE EST =
               NI4 = NI3 + LSRACI(1)
               NCINST(2,NI4) = NI1
C
C              LE GRAND PERE EST COMPLETE
               NCINST(3,NI0) = NI3
               NCINST(4,NI0) = NI1
C
C              L'AFFECTATION  =  ( LE PERE )
               NCINST(1,NI1) = 177
C              LE PERE DE PERE
               NCINST(2,NI1) = NI0
C              LA VARIABLE A CHARGER
               NCINST(3,NI1) = NI2
C              LE CHAINAGE SUR L'EXPRESSION ARITHMETIQUE
               NCINST(4,NI1) = NI4
C
C              LA VARIABLE A CHARGER DE LA VALEUR DE L'EXP_ARITH
               IF( NOVARU .GT. NBVAR0 ) THEN
                  IF( NOVARU .GE. NBVAR0+1 .AND.
     %                NOVARU .LT. NBVAR1 ) THEN
C                    PARAMETRE DE LA FONCTION
                     NCINST(1,NI2) = 10
                  ELSE IF( NOVARU .EQ. NBVAR1 ) THEN
C                    NOM DE LA FONCTION
                     NCINST(1,NI2) = 11
                  ELSE
C                    VARIABLE LOCALE
                     NCINST(1,NI2) = 12
                  ENDIF
C                 LE DECALAGE PAR RAPPORT A NOM_FONC
                  NCINST(3,NI2) = NOVARU - NBVAR0
               ELSE
C                 VARIABLE GLOBALE
                  NCINST(1,NI2) = 4
C                 LE NUMERO DE LA VARIABLE GLOBALE
                  NCINST(3,NI2) = NOVARU
               ENDIF
               NCINST(2,NI2) = NI1
               NCINST(4,NI2) = 0
C
C              EN_REMONTANT_ALLER A L'INSTRUCTION SUIVANTE
C              L'INSTRUCTION SUIVANTE A POUR PERE CETTE INSTRUCTION
               NCINST(1,NI3) = 174
               NCINST(2,NI3) = NI0
               NCINST(3,NI3) = 0
               NCINST(4,NI3) = 0
C              L'INSTRUCTION SUIVANTE EST MISE DANS LA PILE
               LAPISI( LHPISI ) = NI3
C
C              PASSAGE A L'INSTRUCTION SUIVANTE
               NL = NLPTVI
               NC = NCPTVI
               GOTO 100
            ELSE
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LU: PAS DE = APRES UNE VARIABLE'
               ELSE
                  KERR(1) = 'LU: ABSENCE of = AFTER a VARIABLE'
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
         ENDIF
C
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  'LU: INSTRUCTION INCONNUE'
         ELSE
            KERR(1) = 'LU: UNKNOWN STATEMENT'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     FIN DE LA FONCTION   FINFONC ET RETOUR
C     ==================
 8000 NOINST = NOINST + 1
      IF( NOINST .GT. MXINST ) GOTO 9200
      NI0 = LAPISI(LHPISI)
      NCINST(4,NI0) = NOINST
      NCINST(1,NOINST) =  99
      NCINST(2,NOINST) =  NI0
      NCINST(3,NOINST) =  0
      NCINST(4,NOINST) =  0
C
C     LES NUMEROS DES INSTRUCTIONS ALLER A DES ETIQUETTES DOIVENT ETRE
C     MISES A JOUR
      DO 8020 I=1,NOINST
         IF( NCINST(1,I) .EQ. 13 ) THEN
C           LE NUMERO DE L'ETIQUETTE
            NOETIQ = NCINST(3,I)
            IF( NOETIQ .EQ. -1 ) THEN
C              RETOUR  SUR LA DERNIERE INSTRUCTION
               NCINST(3,I) = NOINST
               GOTO 8020
            ELSE IF( NOETIQ .LE. 0 .OR. NOETIQ .GT. NBETIQ ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LU: ETIQUETTE INCORRECTE'
               ELSE
                  KERR(1) = 'LU: WRONG LABEL'
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
            IF( NCETIQ(NOETIQ) .LE. 0 .OR.
     %          NCETIQ(NOETIQ) .GT. NOINST ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LU: ETIQUETTE OUBLIEE:'//KETIQ(NOETIQ)
               ELSE
                  KERR(1) = 'LU: FORGOTTEN LABEL:'//KETIQ(NOETIQ)
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
C           LE NUMERO DE L'ETIQUETTE EST REMPLACEE PAR LE
C           LE NUMERO DE L'INSTRUCTION
            NI1 = NCETIQ(NOETIQ)
C           LE CODE DE L'OPERATION
            NI2 = NCINST(1,NI1)
            IF( NI2 .EQ. 177 .OR. NI2 .EQ. 13 ) THEN
C              L'ETIQUETTE PORTE SUR LE PERE DE L'INSTRUCTION
               NI1 = NCINST(2,NI1)
            ENDIF
            NCINST(3,I) = NI1
         ENDIF
 8020 CONTINUE
C
C     AFFICHAGE DU CHAINAGE EN CAS DE DEBOGUAGE
      IF( LUIMPR .GE. 10 ) THEN
         WRITE(IMPRIM,*)'LU: SORTIE DE ARINST ============='
         DO 8100 I=1,NOINST
            N = NCINST(1,I)
            IF( N.EQ.1 .OR. N.EQ.2 .OR. N.EQ.7 .OR. N.EQ.8 .OR.
     %          N.EQ.9 ) GOTO 8100
            WRITE(IMPRIM,18100) I,(NCINST(J,I),J=1,4)
 8100    CONTINUE
18100 FORMAT( 'LU: NCINST(',I3,')=',4I5)
         DO 8200 I=1,LHKLG
            WRITE(IMPRIM,*) KLG(I)
 8200    CONTINUE
      ENDIF
C
      IF( NBSI .NE. NFINSI ) THEN
C        AVERTISSEMENT : NOMBRE INCORRECT DE FINSI PAR RAPPORT AUX SI
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='LU: NOMBRE DE "SI" NON EGAL AU NOMBRE DE "FINSI"'
         ELSE
         KERR(1)='LU:NUMBER OF ''IF'' NON EQUAL TO NUMBER OF ''ENDIF'' '
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
      NRETOU = 0
      RETURN
C
 9200 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: TABLE DES INSTRUCTIONS SATUREE'
         KERR(2) = 'LU: DESOLE IL FAUT AUGMENTER MXINST DANS LIDFFO'
      ELSE
         KERR(1) = 'LU: SATURATED TABLE OF INSTRUCTIONS'
         KERR(2) = 'LU: SORRY MXINST MUST BE AUGMENTED IN PROC LIDFFO'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
 9300 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: TABLE DES CONSTANTES SATUREE'
      ELSE
         KERR(1) = 'LU: SATURED TABLE OF CONSTANTES'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
 9400 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: PILE DES INSTRUCTIONS SATUREE'
      ELSE
         KERR(1) = 'LU: SATURED STACK OF INSTRUCTIONS'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
 9900 NRETOU = 1
      RETURN
      END
