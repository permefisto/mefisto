      SUBROUTINE LIDFFO( NLD , NCD , NLFIN , NCFIN , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: LIRE ET DECLARER UNE FONCTION DE SYNTAXE
C ----
C      DEFFONC <NOM_FONC>(<PARAM>{,<PARAM>});
C        {DEFVAR  NOM_VAR{,NOM_VAR} ; }
C        {DEFETIQ ETIQ{,ETIQ} ; }
C        {DEFLAB  ETIQ{,ETIQ} ; }
C        {ETIQ :} INSTRUCTION
C       {{ETIQ :} INSTRUCTION }
C        {ETIQ :} FINFONC ;
C
C        INSTRUCTION := | NOM_FONC = EXP_ARITH ;
C                       | NOM_VAR = EXP_ARITH ;
C                       | AFFICHER EXP_ARITH{,EXP_ARITH} ;
C                       | SI EXP_LOGIQ ALORS INSTRUCTION
C                                     {SINON INSTRUCTION}
C                                      FINSI ;
C                       | POUR NOM_VAR = EXP_ARITH A EXP_ARITH ;
C                            INSTRUCTION
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
C       CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NLD NCD : NUMERO DE LIGNE ET COLONNE DANS KLG DU DEBUT DE DEFFONC
C
C SORTIES :
C ---------
C NLFIN NCFIN : NUMERO DE LIGNE ET COLONNE DANS KLG DU ; FINAL
C NRETOU  : =0 SI PAS D'ERREUR
C           >0 SI UNE ERREUR A ETE RENCONTREE
C           -1 SI LE CARACTERE D'ABANDON @ A ETE FRAPPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C23456---------------------------------------------------------------012
      PARAMETER        (MXINST=3072)
      include"./incl/nbcamo.inc"
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___texte.inc"
      include"./incl/a_fonction__arbre.inc"
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      INTEGER            NCINST(1:4,1:MXINST)
      CHARACTER*1        CAR
      CHARACTER*4        NOM
      CHARACTER*24       NMFONC
C
C     SAUVEGARDE DU NOMBRE DE VARIABLES ACTUELLES
      NBVAR0 = NBVARU
      NTLXFO = 0
C
C     LE DERNIER CARACTERE DE DEFFONC
      NL = NLD
      NC = NCD + 6
C
C     LECTURE DE L'ITEM SUIVANT
      CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 1
      ENDIF
C
C     RECHERCHE DU NOM DE LA FONCTION
C     ===============================
      NCID = NC
      CAR  = KLG(NL)(NC:NC)
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE 1-ER CARACTERE N'EST PAS  UNE LETTRE LICITE
         NBLGRC(NRERR) = 1
         KERR(1) = 'LU: NOM INCORRECT DE FONCTION' //
     %              KLG(NL)(NCID:NBCALI)
         GOTO 9000
      ENDIF
C
C     LES CARACTERES SUIVANTS PEUVENT ETRE DES LETTRES OU CHIFFRES
      IF( NC .EQ. NBCALI ) GOTO 300
C     NC N'EST PAS LE DERNIER CARACTERE DE LA LIGNE
 110  NC  = NC + 1
      CAR = KLG(NL)(NC:NC)
C
      IF( CAR .EQ. ' ' .OR. CAR .EQ. '(' ) THEN
C        FIN D'IDENTIFICATEUR DE LA FONCTION
         NC = NC - 1
         GOTO 300
      ENDIF
C
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %           (CAR .GE. '0' .AND. CAR .LE. '9') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE CARACTERE EST ILLICITE
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: NOM INCORRECT DE FONCTION' //
     %                    KLG(NL)(NCID: NBCALI)
         GOTO 9000
      ENDIF
C     LE CARACTERE EST LICITE : PASSAGE AU SUIVANT
      GOTO 110
C
C     TRAITEMENT DE FIN DE NOM DE FONCTION KLG(NL)(NCID:NCF)
 300  NMFONC = KLG(NL)(NCID:NC)
C
C     CETTE FONCTION NE DOIT PAS AVOIR UN NOM DE FONCTION USUELLE
      CALL CHNMO1( NMFONC(1:5), I )
      IF( I .NE. 0 ) GOTO 305
      CALL CHNMO2( NMFONC(1:6), I )
      IF( I .NE. 0 ) GOTO 305
      GOTO 310
C
C     ERREUR
 305  NBLGRC(NRERR) = 3
      KERR(1) = 'LU: LA FONCTION ' // NMFONC
      KERR(2) = 'LU: A MEME NOM QU''UNE FONCTION USUELLE'
      KERR(3) = 'LU: CE QUI EST INTERDIT'
      GOTO 9000
C
C     LA FONCTION NE DOIT PAS AVOIR LE NOM D'UNE VARIABLE
 310  DO 320 I=1,NBVARU
         IF( NMFONC .EQ. KVARU(I) ) THEN
            NBLGRC(NRERR) = 3
            KERR(1) = 'LU: LA FONCTION ' // NMFONC
            KERR(2) = 'LU: A MEME NOM QU''UNE VARIABLE'
            KERR(3) = 'LU: CE QUI EST INTERDIT'
            GOTO 9000
         ENDIF
 320  CONTINUE
C
C     RECHERCHE DE LA (
      CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 1
      ENDIF
      IF( KLG(NL)(NC:NC) .NE. '(' ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: PAS DE ( APRES LE NOM DE LA FONCTION'
         GOTO 9000
      ENDIF
C
C     RECHERCHE DES PARAMETRES ET DE LEUR NOMBRE
C     ==========================================
      NBPARA = 0
C
 400  CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 1
      ENDIF
      NCP = NC
C
C     S'AGIT IL D'UN IDENTIFICATEUR CORRECT ?
      IF( NL .EQ. 0 ) GOTO 9000
C
C     RECHERCHE DU NOM DE LA VARIABLE
      CAR  = KLG(NL)(NC:NC)
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE 1-ER CARACTERE N'EST PAS  UNE LETTRE LICITE
         NBLGRC(NRERR) = 1
         KERR(1) ='LU: NOM INCORRECT DE PARAMETRE' //
     %             KLG(NL)(NCP:NBCALI)
         GOTO 9000
      ENDIF
C
C     LES CARACTERES SUIVANTS PEUVENT ETRE DES LETTRES OU CHIFFRES
      IF( NC .EQ. NBCALI ) GOTO 600
C     NC N'EST PAS LE DERNIER CARACTERE DE LA LIGNE
 500  NC  = NC + 1
      CAR = KLG(NL)(NC:NC)
C
      IF( CAR .EQ. ' ' .OR. CAR .EQ. ',' .OR. CAR .EQ. ')' ) THEN
C        FIN DE PARAMETRE
         NC = NC - 1
         GOTO 600
      ENDIF
C
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %           (CAR .GE. '0' .AND. CAR .LE. '9') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE CARACTERE EST ILLICITE
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: NOM INCORRECT DE VARIABLE' //
     %                    KLG(NL)(NCP:NBCALI)
         GOTO 9000
      ENDIF
C     LE CARACTERE EST LICITE : PASSAGE AU SUIVANT
      GOTO 500
C
C     TRAITEMENT DE FIN DE NOM DE PARAMETRE DE LA FONCTION
C     LE PARAMETRE EST TRAITE COMME UNE VARIABLE LOCALE
C     PLACEE DANS DVARU AU DELA DU NOM DE LA FONCTION
 600  IF( NBVARU .GE. MXVARU ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: PILE SATUREE DES VARIABLES'
         GOTO 9000
      ENDIF
C     AJOUT DE LA VARIABLE
      NBVARU = NBVARU + 1
      KVARU(NBVARU) = KLG(NL)(NCP:NC)
      NBPARA = NBPARA + 1
C
C     DOUBLE PARAMETRE DE MEME IDENTIFICATEUR ?
      DO 610 I=NBVARU-1,NBVAR0+1,-1
         IF( KVARU(I) .EQ. KVARU(NBVARU) ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'LU: 2 NOMS IDENTIQUES' // KVARU(I)
            GOTO 9000
         ENDIF
 610  CONTINUE
C
C     PASSAGE AU CARACTERE NON BLANC SUIVANT
      CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 1
      ENDIF
      CAR = KLG(NL)(NC:NC)
      IF( CAR .EQ. ',' ) THEN
         GOTO 400
      ELSE IF( CAR .NE. ')' ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LU: PAS DE ) DANS LA DEFINITION DE LA FONCTION'
         GOTO 9000
      ENDIF
C
C     LE NOM DE LA FONCTION DEVIENT LA PREMIERE VARIABLE LOCALE
C     DANS DVARU DERRIERE LES PARAMETRES DE LA FONCTION
      IF( NBVARU .GE. MXVARU ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: PILE SATUREE DES VARIABLES'
         GOTO 9000
      ENDIF
C     AJOUT DE LA VARIABLE
      NBVARU = NBVARU + 1
      KVARU(NBVARU) = NMFONC
C
C     LES NBPARA PARAMETRES DE LA FONCTION SONT
C     DANS DVARU(NBVAR0+1:NBVAR0+NBPARA)
C     LA FONCTION EST LA VARIABLE DVARU(NBVAR1)
      NBVAR1 = NBVARU
C
C     RECHERCHE DU ; DERRIERE )
C     =========================
      CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 1
      ENDIF
      IF( KLG(NL)(NC:NC) .NE. ';' ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: PAS DE ; DERRIERE NOMFONC(PARAMETRES) '
         CALL LERESU
         GOTO 630
      ENDIF
C
C     RECHERCHE DU CARACTERE SUIVANT LE ; DE LA DEFINITION DE FONCTION
 620  CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) THEN
         CALL LIRLIG( I )
         IF( I .EQ. -1 ) THEN
C           LE CARACTERE D'ABANDON @ A ETE FRAPPE
            NRETOU = -1
            GOTO 9100
         ENDIF
         IF( I .NE. 0 ) GOTO 9000
         NL = LHKLG
         NC = 0
         GOTO 620
      ENDIF
C
C     LECTURE JUSQU'A TROUVER FINFONC ou ENDFUNC ;
C     ============================================
 630  NLDEBU = NL
      NCDEBU = NC
      CALL LIJUMO( 2, 'FINFONC', 'ENDFUNC', NL, NC, NLFIN, NCFIN )
      IF( NLFIN .EQ. -1 ) THEN
C        FRAPPE DE @ => ABANDON
         NRETOU = -1
         GOTO 9100
      ELSE IF( NLFIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='LU: PAS DE ''FINFONC'' POUR TERMINER LA FONCTION'
         ELSE
            KERR(1)='LU: NO ''ENDFUNC'' TO FINISH THE FUNCTION'
         ENDIF
         GOTO 9000
      ENDIF
      IF( INTERA .GE. 1 ) THEN
         CALL RECTEF( NRERR  )
      ENDIF
C
C     LE LEXIQUE DE LA FONCTION EST CREE
C     ==================================
      CALL LXLXOU( NTFONC , NMFONC , NTLXFO , MNLXFO )
      IF( NTLXFO .GT. 0 ) THEN
         CALL LXLXDS( NTFONC , NMFONC )
      ENDIF
      CALL LXLXDC( NTFONC , NMFONC , 24 , 3 )
      CALL LXLXOU( NTFONC , NMFONC , NTLXFO , MNLXFO )
C
C     LE TABLEAU >FONCTION>>TEXTE
C     ----------------------------
      NBLITX = NLFIN - NLD + 1
      NBCLTX = NBCALI / NBCAMO
      IF( NBCALI - NBCLTX*NBCAMO .NE. 0 ) NBCLTX = NBCLTX + 1
      NBMOTS = WECATX + NBLITX * NBCLTX
      CALL LXTNDC( NTLXFO , 'TEXTE' , 'MOTS' , NBMOTS )
      CALL LXTSOU( NTLXFO , 'TEXTE' , NTTXFO , MNTXFO )
C     INITIALISATION DU TABLEAU TEXTE
      MCN(MNTXFO+WBCLTX) = NBCLTX
      MCN(MNTXFO+WBLITX) = NBLITX
      MN = MNTXFO + WECATX
      DO 710 J=NLD,NLFIN
         DO 700 I=1,NBCALI,NBCAMO
            MCN(MN) = ICHARX( KLG(J)(I:I+NBCAMO-1) )
            MN = MN + 1
 700     CONTINUE
 710  CONTINUE
C     MISE A BLANC DES CARACTERES AVANT DEFFONC
      MN = MNTXFO + WECATX
      DO 720 I=1,NCD-NBCAMO,NBCAMO
         MCN(MN) = ICHARX('    ')
         MN = MN + 1
 720  CONTINUE
      J = MOD( NCD-1 , NBCAMO )
      IF( J .GT. 0 ) THEN
         NOM(1:J) = ' '
         NOM(J+1:NBCAMO) = 'DEFFONC'
         MCN(MN) = ICHARX( NOM )
      ENDIF
C     MISE A BLANC DES CARACTERES APRES FINFONC
      MN = MNTXFO + WECATX + (NLFIN-NLD)*NBCLTX + (NCFIN+6)/NBCAMO
      J  = MOD( NCFIN+6 , NBCAMO )
      IF( J .GT. 0 ) THEN
         IF( J .EQ. 1 ) THEN
            NOM = 'C   '
         ELSE IF( J .EQ. 2 ) THEN
            NOM = 'NC  '
         ELSE IF( J .EQ. 3 ) THEN
            NOM = 'ONC '
         ENDIF
         MCN(MN) = ICHARX( NOM )
      ENDIF
      DO 730 I=(NCFIN+6)/NBCAMO+1,NBCLTX
         MN = MN + 1
         MCN(MN) = ICHARX( '    ' )
 730  CONTINUE
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNTXFO) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNTXFO + MOTVAR(6) ) = NONMTD( '~>>>TEXTE' )
C
C     --------------------------------------------
C     EFFACEMENT DES COMMENTAIRES DANS LA FONCTION
C     --------------------------------------------
      CALL EFACOM( NLDEBU , NCDEBU , NLFIN , NCFIN , NRETOU )
      IF( NRETOU .NE. 0 ) GOTO 9000
C
      NL     = NLDEBU
      NC     = NCDEBU - 1
C
C     RECHERCHE DEFVAR ET/OU DEFETIQ ET/OU DEFLAB
C     ===========================================
 800  CALL CARPNB( NL , NC )
      IF( KLG(NL)(NC:NC+5) .EQ. 'DEFVAR'  ) THEN
         CALL LIDFVA( NL , NC , NLF , NCF , NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9000
         NL = NLF
         NC = NCF
         GOTO 800
C
      ELSE IF( KLG(NL)(NC:NC+6) .EQ. 'DEFETIQ'  ) THEN
         CALL LIDFET( NL , NC , NLF , NCF , NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9000
         NL = NLF
         NC = NCF
         GOTO 800
C
      ELSE IF( KLG(NL)(NC:NC+5) .EQ. 'DEFLAB'  ) THEN
         CALL LIDFLA( NL , NC , NLF , NCF , NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9000
         NL = NLF
         NC = NCF
         GOTO 800
      ENDIF
C
C     LES INSTRUCTIONS SONT ENTRE (NLDEBU,NCDEBU) ET (NLFIN,NCFIN)
C     ------------------------------------------------------------
      NLDEBU = NL
      NCDEBU = NC
C
C     -------------------------------------------------------------
C     ANALYSE ET GENERATION ( RETROUVER ET CHAINER LES OPERATIONS )
C     -------------------------------------------------------------
      NCFIN = NCFIN - 1
      CALL ARINST( NLDEBU , NCDEBU , NLFIN , NCFIN , MXINST , NCINST ,
     %             NOINST , NRETOU )
      IF( NRETOU .NE. 0 ) GOTO 9000
      NCFIN = NCFIN + 1
C
C
C     FONCTION CORRECTEMENT COMPILEE : LE TABLEAU >FONCTION>>ARBRE
C     -------------------------------------------------------------
C     NBVALO  LE NOMBRE DE VARIABLES LOCALES DE LA FONCTION
C             >0 CAR LA VARIABLE LOCALE 1 EST LA VARIABLE FONCTION
      NBVALO = NBVARU - NBVAR1 + 1
      L1ARBR = 4
      L2ARBR = NOINST
      CALL LXTNDC( NTLXFO , 'ARBRE' , 'MOTS' , WARBRE+L1ARBR*L2ARBR )
      CALL LXTSOU( NTLXFO , 'ARBRE' , NTARFO , MNARFO )
C     INITIALISATION DU TABLEAU TEXTE
      MCN(MNARFO+W1ARBR) = L1ARBR
      MCN(MNARFO+W2ARBR) = L2ARBR
      MCN(MNARFO+WBPAFO) = NBPARA
      MCN(MNARFO+WBVALO) = NBVALO
C     LA POSITION DE LA RACINE DE L'ARBRE
      MCN(MNARFO+WRARBR) = 1
      MN = MNARFO + WARBRE
      DO 8100 I=1,L2ARBR
         DO 8000 J=1,L1ARBR
            MCN(MN) = NCINST(J,I)
            MN = MN + 1
 8000    CONTINUE
 8100 CONTINUE
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARFO) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARFO + MOTVAR(6) ) = NONMTD( '~>FONCTION>>ARBRE' )
C
C     FIN DU CALCUL
C     =============
      NBVARU = NBVAR0
      NBETIQ = 0
      NRETOU = 0
C     RECHERCHE DU ; DERRIERE FINFONC
      NCFIN  = NCFIN + 6
      NLPTV2 = NLFIN
      NCPTV2 = NCFIN
      CALL CARPNB( NLFIN , NCFIN )
      IF( NLFIN .LE. 0 ) THEN
C        PROTECTION D'UN OUBLI DU ;
         KLG(NLPTV2)(NCPTV2:NCPTV2) = ';'
         NLFIN = NLPTV2
         NCFIN = NCPTV2
      ENDIF
 8200 IF( KLG(NLFIN)(NCFIN:NCFIN) .NE. ';' ) THEN
         NCFIN = NCFIN + 1
         IF( NCFIN .GT. NBCALI ) THEN
            NLFIN = NLFIN + 1
            NCFIN = 1
         ENDIF
         GOTO 8200
      ENDIF
      IF( KLG(NLFIN)(NCFIN:NCFIN) .EQ. ';' ) THEN
         NLPTV2 = NLFIN
         NCPTV2 = NCFIN
      ENDIF
C
C     AFFICHAGE DE L'AJOUT DE LA FONCTION
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: AJOUT DE LA FONCTION ' // NMFONC
      ELSE
         KERR(1) = 'LU: ADDED FUNCTION ' // NMFONC
      ENDIF
      CALL LERESU
      RETURN
C
C     ERREUR DETECTEE
C     ===============
 9000 NRETOU = NBLGRC(NRERR) + 1
      IF( LANGAG .EQ. 0 ) THEN
      KERR(NRETOU) = ' '
      KERR(NRETOU+1)=' GRAMMAIRE DE DEFINITION D''UNE FONCTION'
      KERR(NRETOU+2)=' --------------------------------------'
      KERR(NRETOU+3)=' DEFFONC <NOM_FONC>(<PARAM>{,<PARAM>}); '
      KERR(NRETOU+4)='  {DEFVAR  NOM_VAR {,NOM_VAR} ; } '
      KERR(NRETOU+5)='  {DEFETIQ ETIQ {,ETIQ} ; }'
      KERR(NRETOU+6)='  {ETIQ :} INSTRUCTION'
      KERR(NRETOU+7)=' {{ETIQ :} INSTRUCTION }'
      KERR(NRETOU+8)='  {ETIQ :} FINFONC ;'
      KERR(NRETOU+9)='  INSTRUCTION ::='
      KERR(NRETOU+10)=' | NOM_FONC = EXP_ARITH ;'
      KERR(NRETOU+11)=' | NOM_VAR = EXP_ARITH ;'
      KERR(NRETOU+12)=' | AFFICHER EXP_ARITH{,EXP_ARITH} ;'
      KERR(NRETOU+13)=' | SI EXP_LOGIQ ALORS INSTRUCTION'
      KERR(NRETOU+14)='               {SINON INSTRUCTION}'
      KERR(NRETOU+15)='   FINSI ;'
      KERR(NRETOU+16)=' | POUR NOM_VAR = EXP_ARITH A EXP_ARITH ;'
      KERR(NRETOU+17)='      INSTRUCTION'
      KERR(NRETOU+18)='   FINPOUR ;'
      KERR(NRETOU+19)=' | REPETER ;'
      KERR(NRETOU+20)='      INSTRUCTION'
      KERR(NRETOU+21)='   JUSQUE  EXP_LOGIQ ;'
      KERR(NRETOU+22)=' | TANTQUE EXP_LOGIQ ;'
      KERR(NRETOU+23)='      INSTRUCTION'
      KERR(NRETOU+24)='   FINTANT ;'
      KERR(NRETOU+25)=' | RETOUR ;'
      KERR(NRETOU+26)=' | ALLER ETIQ ;'
      KERR(NRETOU+27)='UN NOM DE FONCTION DOIT ETRE DIFFERENT d''UN NOM'
      KERR(NRETOU+28)='de VARIABLE OU de FONCTION USUELLE (COS,EXP,..)'
      ELSE
      KERR(NRETOU) = ' '
      KERR(NRETOU+1)=' THE DEFINITION OF THE GRAMMAR OF A FUNCTION'
      KERR(NRETOU+2)=' -------------------------------------------'
      KERR(NRETOU+3)=' DEFFUNC <NAME_FONC>(<PARAM>{,<PARAM>}); '
      KERR(NRETOU+4)='  {DEFVAR  NAME_VAR {,NAME_VAR} ; } '
      KERR(NRETOU+5)='  {DEFLAB LABEL {,LABEL} ; }'
      KERR(NRETOU+6)='  {LABEL :} INSTRUCTION'
      KERR(NRETOU+7)=' {{LABEL :} INSTRUCTION }'
      KERR(NRETOU+8)='  {LABEL :} ENDFUNC ;'
      KERR(NRETOU+9)='  INSTRUCTION ::='
      KERR(NRETOU+10)=' | NAME_FONC = EXP_ARITH ;'
      KERR(NRETOU+11)=' | NAME_VAR = EXP_ARITH ;'
      KERR(NRETOU+12)=' | DISPLAY EXP_ARITH{,EXP_ARITH} ;'
      KERR(NRETOU+13)=' | IF EXP_LOGIQ THEN INSTRUCTION'
      KERR(NRETOU+14)='               {ELSE INSTRUCTION}'
      KERR(NRETOU+15)='   ENDIF ;'
      KERR(NRETOU+16)=' | FOR NAME_VAR = EXP_ARITH TO EXP_ARITH ;'
      KERR(NRETOU+17)='      INSTRUCTION'
      KERR(NRETOU+18)='   ENDFOR ;'
      KERR(NRETOU+19)=' | REPEAT ;'
      KERR(NRETOU+20)='      INSTRUCTION'
      KERR(NRETOU+21)='   UNTIL EXP_LOGIQ ;'
      KERR(NRETOU+22)=' | WHILE EXP_LOGIQ ;'
      KERR(NRETOU+23)='      INSTRUCTION'
      KERR(NRETOU+24)='   ENDWHILE ;'
      KERR(NRETOU+25)=' | RETURN ;'
      KERR(NRETOU+26)=' | GOTO LABEL ;'
      KERR(NRETOU+27)='The FUNCTION NAME MUST BE DIFFERENT of a NAME'
      KERR(NRETOU+28)='of VARIABLES OR USUAL FUNCTIONS (COS,EXP,..)'
      ENDIF
      NBLGRC(NRERR) = NRETOU + 28
      CALL LEREUR
C
 9100 NBVARU = NBVAR0
      NBETIQ = 0
      NRETOU = 1
C     DESTRUCTION DU LEXIQUE DE LA FONCTION
      IF( NTLXFO .GT. 0 ) THEN
         CALL LXLXDS( NTFONC , NMFONC )
      ENDIF
C
      RETURN
      END
