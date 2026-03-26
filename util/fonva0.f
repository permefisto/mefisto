      SUBROUTINE FONVA0( NOFONC, NBPAFO, ICPAFO,
     %                   NCODEV, DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      CALCULER LA VALEUR DE LA FONCTION DE NUMERO NOFONC DANS LE
C -----      LEXIQUE DES FONCTIONS ET DE PARAMETRES DPARA
C            CF ~/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR

C ENTREES :
C ---------
C NOFONC  : NUMERO DE LA FONCTION DANS LE LEXIQUE DES FONCTIONS
C NBPAFO  : NOMBRE DE PARAMETRES DE LA FONCTION
C ICPAFO  : INDICE DANS DCTE DE LA VALEUR DU 1-ER PARAMETRE DE
C           L'APPEL DE LA FONCTION NOFONC

C SORTIES :
C ---------
C NCODEV  : 0 DBLVAL N'EST PAS INITIALISEE EN SORTIE et DBLVAL=0D0
C           1 DBLVAL   EST     INITIALISEE EN SORTIE
C DBLVAL  : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER(W)
      PARAMETER        (MXPIOP=2048,MXPIFV=48)
C     MXPIOP  NOMBRE MAXIMAL D'OPERATIONS DANS LA PILE DES OPERATIONS
C                            DE CHAQUE FONCTION DE VALEUR A CALCULER
C     MXPIFV  NOMBRE MAXIMAL DE FONCTIONS DE FONCTIONS CALCULABLES
C                            F( G ( H ( I ( ... )))
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_fonction__arbre.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  DBLVAL
      CHARACTER*24      NMFONC
      INTEGER           NOFOLX(1:MXPIFV),
     %                  MNARBR(1:MXPIFV),MNPIL1(1:MXPIFV),
     %                  MNPIL2(1:MXPIFV),LHPIOP(1:MXPIFV),
     %                  LHCTFV(1:MXPIFV),NBPAFV(1:MXPIFV),
     %                  NBVLFV(1:MXPIFV)

C     LA PREMIERE FONCTION A TRAITER : NOFONC EST EMPILEE DANS PIFV
C     =============================================================
      NBPAF1 = NBPAFO
      LHPIFV = 1
C     LE NUMERO DE LA FONCTION DANS LE LEXIQUE DES FONCTIONS
      NOFOLX( LHPIFV ) = NOFONC
C     LE NOMBRE DE PARAMETRES DE LA FONCTION A CALCULER
      NBPAFV( LHPIFV ) = NBPAF1
C     L'INDICE DANS DCTE DU DEBUT DES PARAMETRES
      LHCTFV( LHPIFV ) = ICPAFO
C     SAUVEGARDE DE NBCTE EN CAS DE PROBLEME
      NBCT0  = NBCTE
      NBTYVA = 0

C     ******************************************************************
C     TANT QUE LA PILE DES FONCTIONS N'EST PAS VIDE :
C     TRAITER LA FONCTION EN HAUT DE LA PILE
C     ******************************************************************

  1   IF( LHPIFV .GT. 0 ) THEN

C     RESTAURATION DE L'ARBRE DES OPERATIONS DE LA FONCTION
      IF( NOFOLX(LHPIFV) .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: NOM INCORRECT DE FONCTION'
         ELSE
            KERR(1) = 'LU: INCORRECT NAME of FUNCTION'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      CALL LXNLOU( NTFONC, NOFOLX(LHPIFV), NTLXFO, MNLXFO )
      IF( NTLXFO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(3)(1:10),'(I10)') NOFOLX(LHPIFV)
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: FONCTION UTILISATEUR DE NO INCORRECT '
     %                // KERR(3)(1:10)
         ELSE
            KERR(1) = 'LU: USER FUNCTION with an INCORRECT NAME '
     %                // KERR(3)(1:10)
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF

C     LE NOM DE LA FONCTION
      CALL NMOBNU( 'FONCTION', NOFOLX(LHPIFV), NMFONC )
      CALL LXTSOU( NTLXFO, 'ARBRE', NTARBR, MNARBR(LHPIFV) )
      IF( NTARBR .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: FONCTION '//NMFONC//' NON COMPILEE'
         ELSE
            KERR(1) = 'LU: FUNCTION '//NMFONC//' NOT COMPILED'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF

      IF( NBPAF1 .NE. MCN(MNARBR(LHPIFV)+WBPAFO) ) THEN
            WRITE(KERR(4)(1:3),'(I3)') MCN(MNARBR(LHPIFV)+WBPAFO)
            WRITE(KERR(4)(4:6),'(I3)') NBPAF1
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='LU: FONCTION ' // NMFONC
               KERR(2)='DEFINIE AVEC' // KERR(4)(1:3) // ' PARAMETRES'
               KERR(3)='APPELEE AVEC' // KERR(4)(4:6) // ' PARAMETRES'
            ELSE
               KERR(1)='LU: FUNCTION ' // NMFONC
               KERR(2)='DEFINED with ' // KERR(4)(1:3) // ' PARAMETERS'
               KERR(3)='CALLED  with ' // KERR(4)(4:6) // ' PARAMETERS'
            ENDIF
            CALL LEREUR
         GOTO 9900
      ENDIF

C     MNPIL1(LHPIFV) ADRESSE MCN DE LA PILE DES OPERATEURS
C     MNPIL2(LHPIFV) ADRESSE MCN DE LA PILE DU NOMBRE DE FILS DE L'OPERATEUR
      MNPIL1(LHPIFV) = 0
      CALL TNMCDC( 'ENTIER', 2*MXPIOP, MNPIL1(LHPIFV) )
      MNPIL2(LHPIFV) = MNPIL1(LHPIFV) + MXPIOP

C     LA RACINE DE L'ARBRE
      NOPERP = MCN( MNARBR(LHPIFV) + WRARBR )

C     LE NOMBRE DE VARIABLES LOCALES DE LA FONCTION
      NBVLFV(LHPIFV) = MCN( MNARBR(LHPIFV) + WBVALO )

C     LES VARIABLES LOCALES SONT EMPILEES DANS DCTE APRES LES PARAMETRES
      IF( NBCTE+NBVLFV(LHPIFV) .GE. MXDCTE ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  'LU: dans FONVA0 la PILE DCTE est SATUREE'
         ELSE
            KERR(1) =  'LU: in FONVA0 the STACK DCTE is SATURATED'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      DCTE(NBCTE+1) = 77777777777.D0
      NBCTE = NBCTE + NBVLFV(LHPIFV)

C     L'OPERATION RACINE
      MNAP   = MNARBR(LHPIFV) + WARBRE +
     %        (NOPERP-1) * MCN(MNARBR(LHPIFV)+W1ARBR) - 1
      NCOPEP = MCN(MNAP+1)

C     ANALYSE DE L'OPERATION
      IF( NCOPEP .EQ. 3 ) THEN
         DBLVAL = DCTE( MCN(MNAP+3) )
         GOTO 2
      ELSE IF( NCOPEP .EQ. 4 ) THEN
         DBLVAL = DVARU( MCN(MNAP+3) )
         GOTO 2
      ELSE IF( NCOPEP .EQ. 5 ) THEN
       CALL VAVATS( MCN(MNAP+3), NCODEV, DBLVAL )
         GOTO 2
      ELSE IF( NCOPEP .EQ. 10 ) THEN
         DBLVAL = DCTE( LHCTFV(LHPIFV) - 1 + MCN(MNAP+3) )
         GOTO 2
      ENDIF
      GOTO 3

C     LA VALEUR RETROUVEE EST EMPILEE
  2   IF( NBCTE .GE. MXDCTE ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: dans FONVA0 la PILE DCTE est SATUREE'
         ELSE
            KERR(1) = 'LU: in FONVA0 the STACK DCTE is SATURATED'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      NBCTE = NBCTE + 1
      DCTE(NBCTE) = DBLVAL
      GOTO 8900

C     LAPIOP CONTIENT LE POINTEUR SUR L'OPERATION EMPILEE
C     LAPINF CONTIENT LE NOMBRE DE FILS A TRAITER POUR CETTE OPERATION
C     LA PREMIERE OPERATION EST EMPILEE
C     ----------------------------------------------------------------
  3   LHPIOP(LHPIFV) = 1
      MCN(MNPIL1(LHPIFV)+1) = NOPERP
C     LE NOMBRE DE PARAMETRES DE L'OPERATION NOP
      MCN(MNPIL2(LHPIFV)+1) = MIN( 2, NBPRMT( NCOPEP ) )

C     *******************************
C     LE PARCOURS PRE-FIXE DE L'ARBRE
C     *******************************
C     TANT QU'IL EXISTE UNE OPERATION A FAIRE ALORS TRAITER CETTE OPERATION
10    LHPILA = LHPIOP(LHPIFV)
      DO 100 IP = LHPILA, 1, -1
         IF( MCN(MNPIL1(LHPIFV)+IP) .GT. 0 ) THEN
C           L'OPERATION DE HAUT DE PILE
C           SA POSITION DANS LARBRE
            NOPERP = MCN(MNPIL1(LHPIFV)+IP)
            MNAP   = MNARBR(LHPIFV) + WARBRE +
     %              (NOPERP-1)*MCN(MNARBR(LHPIFV)+W1ARBR) - 1
C           SON CODE OPERATION
            NCOPEP = MCN(MNAP+1)
C           LE NOMBRE DE FILS (0 1 OU 2) A TRAITER
            NBFILS = MCN(MNPIL2(LHPIFV)+IP)
            IF( NBFILS .GT. 0 ) THEN

C              LE FILS NBFILS EST EMPILE
C              =========================
               NBFILA = MIN( 2, NBFILS )
               NOFILS = MCN(MNAP+5-NBFILA)
               IF( NOFILS .EQ. 0 ) THEN
C                 LE NOMBRE DE FILS A TRAITER DE L'OPERATION EST
C                 DECREMENTE
                  MCN(MNPIL2(LHPIFV)+IP) = MCN(MNPIL2(LHPIFV)+IP) - 1
                  GOTO 10
               ENDIF

               IF( LHPIOP(LHPIFV) .GE. MXPIOP ) GOTO 9000
               LHPIOP(LHPIFV) = LHPIOP(LHPIFV) + 1
               MCN( MNPIL1(LHPIFV) + LHPIOP(LHPIFV) ) = NOFILS
C              LE NOMBRE NBPARF DE PARAMETRES DE L'OPERATION NOFILS
               MNAF   = MNARBR(LHPIFV)+WARBRE+
     %                 (NOFILS-1)*MCN(MNARBR(LHPIFV)+W1ARBR) - 1
C              EST EMPILE DANS LAPINF
               NCOPEF = MCN(MNAF+1)
               NBPARF = NBPRMT( NCOPEF )
               IF( (NCOPEF .GE. 3  .AND. NCOPEF .LE.  5) .OR.
     %             (NCOPEF .GE. 10 .AND. NCOPEF .LE. 12) .OR.
     %              NCOPEF .EQ. 14  ) THEN

C                 CONSTANTE OU VARIABLE OU VALEUR TMS OU PARAMETRE
C                 ------------------------------------------------
C                 MISE DANS LA PILE DES CONSTANTES = DCTE
                  IF( NBCTE .GE. MXDCTE ) THEN
                     NBLGRC(NRERR) = 1
                     IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='LU: dans FONVA0 la PILE DCTE est SATUREE'
                     ELSE
                     KERR(1)='LU: in FONVA0 the STACK DCTE is SATURATED'
                     ENDIF
                     CALL LEREUR
                     GOTO 9900
                  ENDIF
                  NBCTE = NBCTE + 1
                  IF( NCOPEF .EQ. 3 ) THEN
                     DCTE(NBCTE) = DCTE ( MCN(MNAF+3) )
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 4 ) THEN
                     DCTE(NBCTE) = DVARU( MCN(MNAF+3) )
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 5 ) THEN
                     CALL VAVATS( MCN(MNAF+3), NCODEV, DBLVAL )
                     DCTE(NBCTE) = DBLVAL
                     NBTYVA = 3
                  ELSE IF( NCOPEF .GE. 10 .AND.
     %                     NCOPEF .LE. 12 ) THEN
C                    PARAMETRE OU NOM FONCTION OU VARIABLE LOCALE
                     DCTE(NBCTE) = DCTE( LHCTFV(LHPIFV)-1+MCN(MNAF+3) )
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 14 ) THEN
C                    LE NUMERO DE LA CONSTANTE CHAINE EST EMPILEE
                     DCTE(NBCTE) = MCN(MNAF+3)
                     NBTYVA = 14
                  ENDIF
C                 LA VALEUR CONSTANTE EST EMPILEE
                  MCN( MNPIL1(LHPIFV) + LHPIOP(LHPIFV) ) = - NBTYVA
               ENDIF
               MCN(MNPIL2(LHPIFV)+LHPIOP(LHPIFV)) = MIN( 2, NBPARF )
C              LE NOMBRE DE FILS A TRAITER DE L'OPERATION EST DECREMENTE
               MCN(MNPIL2(LHPIFV)+IP) = MCN(MNPIL2(LHPIFV)+IP) - 1
               GOTO 10
C
            ELSE
C
C              TOUS LES FILS DE NOPERP ONT ETE TRAITES DONC
C              LES PARAMETRES SONT TOUS CALCULES ET AU DESSUS DE LUI
C              DANS LA PILE => CALCUL DE L'OPERATION
C              ======================================================
C              LE NOMBRE NBPARF DE PARAMETRES DE L'OPERATION NOPERP
               NBPARF = NBPRMT( NCOPEP )
C
C              CALCUL DE L'OPERATION SELON LE NOMBRE DE PARAMETRES
               IF( NCOPEP .LT. 300 ) THEN
C
                  IF( NBPARF .EQ. 0 ) THEN
C
C                    0 PARAMETRE
C                    ===========
                     IF( NCOPEP .EQ. 13 ) THEN
C
C                       ALLER ETIQUETTE
C                       ---------------
C                       L'INSTRUCTION ALLER EST DEPILEE
                        LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - 1
C                       L'INSTRUCTION AVEC CETTE ETIQUETTE EST EMPILEE
                        NOINST = MCN( MNAP + 3 )
                        MCN(MNPIL1(LHPIFV)+IP-1) = NOINST
C                       LE NOMBRE DE PARAMETRES DE L'OPERATION NOINST
                        MNAF   = MNARBR(LHPIFV)+WARBRE+
     %                 (NOINST-1)*MCN(MNARBR(LHPIFV)+W1ARBR) - 1
C                       EST EMPILE DANS LAPINF
                        NCOPEF = MCN(MNAF+1)
                        NBPARF = NBPRMT( NCOPEF )
                        MCN(MNPIL2(LHPIFV)+IP-1) = MIN(2,NBPARF)
                        GOTO 10
C
                     ELSE IF( NCOPEP .EQ. 99 ) THEN
C
C                       FINFONC
C                       -------
                        GOTO 8900
                     ENDIF
C
                  ELSE IF( NBPARF .EQ. 1 ) THEN
C
C                    1 PARAMETRE
C                    ===========
                     IF( NCOPEP .LT. 174 ) THEN
C
C                       OPERATEUR UNAIRE CALCULABLE AVEC DCTE(NBCTE)
C                       --------------------------------------------
                        CALL VAOPUN( NCOPEP,DCTE(NBCTE), NCODEV,DBLVAL )
                        IF( NCODEV .EQ. 0 ) GOTO 9900
C
                     ELSE
C
                        IF( NCOPEP .EQ. 174 ) THEN
C
C                          SUITE D'UNE INSTRUCTION
C                          -----------------------
C                          L'INSTRUCTION AVEC CETTE ETIQUETTE
C                          REMPLACE SUITE
                           NOINST = MCN( MNAP + 3 )
                           IF( NOINST .EQ. 0 ) THEN
C                             PAS DE SUITE: L'INSTRUCTION EST DEPILEE
                              LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - 1
                              GOTO 10
                           ENDIF
C                          NOINST EST L'INSTRUCTION SUIVANTE A EMPILER
                           MCN(MNPIL1(LHPIFV)+IP) = NOINST
C                          LE NOMBRE DE PARAMETRES DE L'OPERATION NOINST
                           MNAF   = MNARBR(LHPIFV)+WARBRE+
     %                    (NOINST-1)*MCN(MNARBR(LHPIFV)+W1ARBR) - 1
C                          EST EMPILE DANS LAPINF
                           NCOPEF = MCN(MNAF+1)
                           NBPARF = NBPRMT( NCOPEF )
                           MCN(MNPIL2(LHPIFV)+IP) = MIN(2,NBPARF)
C
                        ELSE IF( NCOPEP .EQ. 175 ) THEN
C
C                          SI  L'EXPRESSION LOGIQUE A ETE CALCULEE
C                          --  LE RESULTAT EST DANS DCTE(NBCTE)
C                              RESULTAT FAUX SI DCTE(NBCTE)=0, VRAI SINON
C                          L'INSTRUCTION A EXECUTER ENSUITE ECRASE LE SI
                           IF( DCTE(NBCTE) .NE. 0 ) THEN
                              NOINST = 1
                           ELSE
                              NOINST = 0
                           ENDIF
C                          L'INSTRUCTION DE CETTE ETIQUETTE REMPLACE SI
C                          L'INSTRUCTION ALORS_SINON
                           MNAF = MCN( MNAP + 3 )
                           MNAF = MNARBR(LHPIFV)+WARBRE
     %                          + (MNAF-1)*MCN(MNARBR(LHPIFV)+W1ARBR)-1
C                          ALORS OU SINON SELON EXP_LOGIQ
                           NOINST = MCN( MNAF + 3 + NOINST )
                           MCN(MNPIL1(LHPIFV)+IP)= NOINST
C                          LE NOMBRE DE PARAMETRES DE L'OPERATION
C                          DE POSITION NOINST
                           MNAF = MNARBR(LHPIFV)+WARBRE
     %                          + (NOINST-1)*MCN(MNARBR(LHPIFV)+W1ARBR)
                           NOINST = MCN( MNAF )
                           MCN(MNPIL2(LHPIFV)+IP)= MIN(2,NBPRMT(NOINST))
C                          LA CONSTANTE LOGIQUE EST DEPILEE
                           NBCTE = NBCTE - 1
C                          L'EXPRESSION LOGIQUE EST DEPILEE
                           LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - 1
C
                        ELSE IF( NCOPEP .EQ. 177 ) THEN
C
C                           =   AFFECTATION  VARIABLE = VALEUR
C                          ---               -----------------
C                          L'INSTRUCTION DEFINISSANT LA VARIABLE
                           NOINST = MCN( MNAP + 3 )
C                          L'ADRESSE DE CETTE VARIABLE
                           MNAF = MNARBR(LHPIFV) + WARBRE +
     %                    (NOINST-1)*MCN(MNARBR(LHPIFV)+W1ARBR) - 1
C                          SON CODE DE VARIABLE
                           NOINST = MCN(MNAF+1)
                           IF( NOINST .EQ. 4 ) THEN
C                             VARIABLE GLOBALE
                              DVARU( MCN(MNAF+3) ) = DBLVAL
                           ELSE IF( NOINST .GE. 10 .AND.
     %                              NOINST .LE. 12 )  THEN
C                             VARIABLE LOCALE STOCKEE DERRIERE LES PARAMETRES
                              DCTE(LHCTFV(LHPIFV)-1+MCN(MNAF+3))
     %                              = DCTE(NBCTE)
                           ELSE
                              WRITE(KERR(2)(1:4),'(I4)') NOINST
                              NBLGRC(NRERR) = 1
                              IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='LU: ERREUR un CODE VARIABLE est INCORRECT'
     %                        // KERR(2)(1:4)
                              ELSE
                     KERR(1)='LU: ERROR a VARIABLE CODE is INCORRECT'
     %                        // KERR(2)(1:4)
                              ENDIF
                              CALL LEREUR
                              GOTO 9900
                           ENDIF
C                          EXP_ARITH ET = SONT DEPILES
                           LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - 2
C                          LA VALEUR DE L'EXPRESSION ARITHMETIQUE
C                          EST DEPILEE DE DCTE
                           NBCTE = NBCTE - 1
                        ENDIF
                        GOTO 10
                     ENDIF
C
                  ELSE IF( NBPARF .EQ. 2 ) THEN
C
C                    2 PARAMETRES
C                    ============
                     IF( NCOPEP .NE. 200 ) THEN
C
C                       OPERATEUR BINAIRE CALCULABLE AVEC
C                       DCTE(NBCTE-1:NBCTE)
C                       ---------------------------------
                        CALL VAOPBI( NCOPEP ,
     %                               DCTE(NBCTE-1), DCTE(NBCTE) ,
     %                               NCODEV, DBLVAL )
                        IF( NCODEV .EQ. 0 ) GOTO 9900
C
                     ELSE
C
C                       NCOPEP = 200    ENTRE_PARAMETRES
C                       --------------------------------
C                       LES PARAMETRES DU DESSUS DE LA PILE SONT ABAISSES
C                       D'UN CRAN
                        DO 5 I=IP+1,LHPIOP(LHPIFV),1
                           MCN(MNPIL1(LHPIFV)+I-1)=MCN(MNPIL1(LHPIFV)+I)
                           MCN(MNPIL2(LHPIFV)+I-1)=MCN(MNPIL2(LHPIFV)+I)
 5                      CONTINUE
                        LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - 1
                        GOTO 10
                     ENDIF
                  ENDIF
C
               ELSE

C                 FONCTION UTILISATEUR CALCULABLE AVEC
C                 DCTE(NBCTE-NBPARF+1:NBCTE)
C                 LA FONCTION EST EMPILEE
C                 ------------------------------------
                  NUFONC =  MOD(NCOPEP,1000) - 300
                  IF( NUFONC .EQ. 0 ) THEN

C                    PSEUDO FONCTION  AFFICHER
C                    AFFICHER PARAM {, PARAM }
                     IPP = NBCTE - NBPARF
                     DO 90 I=IPP+1,NBCTE
                        WRITE(KERR(MXLGER)(31:34),'(I4)') I-IPP
                        WRITE(KERR(MXLGER)(1:25),'(G25.17)') DCTE(I)
                        IF( LANGAG .EQ. 0 ) THEN
                           KERR(I-IPP) = 'LU: AFFICHAGE VALEUR'//
     %                                    KERR(MXLGER)(31:34)  //
     %                                   ' =' // KERR(MXLGER)(1:25)
                        ELSE
                           KERR(I-IPP) = 'LU: DISPLAYED VALUE'//
     %                                    KERR(MXLGER)(31:34)  //
     %                                   ' =' // KERR(MXLGER)(1:25)
                        ENDIF
 90                  CONTINUE
                     NBLGRC(NRERR) = NBCTE-IPP
                     CALL LERESU
C                    LES PARAMETRES SONT ECRASES DANS DCTE
                     NBCTE  = NBCTE - NBPARF + 1
C                    MISE A ZERO EN CAS D'IMPRESSION DU RESULTAT
                     DCTE( NBCTE ) = 0D0
C                    L'OPERATION FONCTION AFFICHER ET LES PARAMETRES
C                    SONT DEPILES DE LAPIOP
                     LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - NBPARF - 1
                     GOTO 10
                  ENDIF
C
C                 VRAIE FONCTION UTILISATEUR
C                 LA FONCTION EST DEPILEE DE LA PILE DES OPERATIONS
C                 POUR NE PAS LA REEXECUTER AU RETOUR
                  LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - NBPARF
C                 LA FONCTION EST EMPILEE DANS LA PILE DES FONCTIONS
                  LHPIFV = LHPIFV + 1
                  IF( LHPIFV .GT. MXPIFV ) THEN
                     NBLGRC(NRERR) = 2
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(1)='LU: PILE DES FONCTIONS SATUREE'
                        KERR(2)='LU: DANS FONVA0 AUGMENTER MXPIFV'
                     ELSE
                        KERR(1)='LU: SATURATED STACK of FUNCTIONS'
                        KERR(2)='LU: in FONVA0 AUGMENT MXPIFV'
                     ENDIF
                     CALL LEREUR
                     GOTO 9900
                  ENDIF
C                 LE NOMBRE DE PARAMETRES DE LA FONCTION A EMPILER
                  NBPAF1 = NBPARF
C                 LE NUMERO DE LA FONCTION DANS LE LEXIQUE DES FONCTIONS
                  NOFOLX( LHPIFV ) = NUFONC
C                 LE NOMBRE DE PARAMETRES DE LA FONCTION A CALCULER
                  NBPAFV( LHPIFV ) = NBPARF
C                 L'INDICE DANS DCTE DU DEBUT DES PARAMETRES
                  LHCTFV( LHPIFV ) = NBCTE - NBPARF + 1
                  GOTO 1
C
               ENDIF
C
C              LE RESULTAT DE L'OPERATION ECRASE L'OPERATION
               NBCTE = NBCTE - NBPARF + 1
               DCTE( NBCTE ) = DBLVAL
C              LES OPERANDES SONT DEPILES
               LHPIOP(LHPIFV) = LHPIOP(LHPIFV) - NBPARF
C              LE CODE OPERATION DANS LA PILE EST ECRASE
C              PAR LE RESULTAT
               MCN(MNPIL1(LHPIFV)+ LHPIOP(LHPIFV) ) = -3
               MCN(MNPIL2(LHPIFV)+ LHPIOP(LHPIFV) ) =  0
               GOTO 10
            ENDIF
         ENDIF
 100  CONTINUE

C     C'EST LA RACINE : FIN DU CALCUL DE CETTE FONCTION
C     =================================================
C     DESTRUCTION DU TABLEAU DES 2 PILES
 8900 CALL TNMCDS( 'ENTIER', 2*MXPIOP, MNPIL1(LHPIFV) )
C     LA VALEUR DE CETTE FONCTION EN POSITION DE 1-ERE VARIABLE LOCALE
      DBLVAL = DCTE( LHCTFV(LHPIFV) + NBPAFV(LHPIFV) )
C     LES PARAMETRES DE CETTE FONCTION SONT DEPILES
      NBCTE  = LHCTFV(LHPIFV)
C     LA FONCTION EXECUTEE EST DEPILEE
      LHPIFV = LHPIFV - 1
      IF( LHPIFV .GT. 0 ) THEN

C        LE RESULTAT DE LA FONCTION EXECUTEE EST EMPILE DANS DCTE
C        A LA PLACE DU PREMIER PARAMETRE
         DCTE( NBCTE ) = DBLVAL

C        LE RESULTAT DE LA FONCTION EXECUTEE EST EMPILE
C        DANS LA PILE DES OPERATIONS COMME ETANT UNE CONSTANTE
         MCN(MNPIL1(LHPIFV)+ LHPIOP(LHPIFV) ) = -3
         MCN(MNPIL2(LHPIFV)+ LHPIOP(LHPIFV) ) =  0
         GOTO 10
      ENDIF
C     LA PILE DES FONCTIONS EST VIDE
      ENDIF

C     ******************************
C     LA PILE DES FONCTIONS EST VIDE
C     ******************************
C     LE RESULTAT DE LA FONCTION EST DBLVAL : SORTIE AVEC VALEUR CALCULEE
      NCODEV = 1
      GOTO 9999

C     ERREUR
 9000 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='LU:(SP FONVA0) PILE SATUREE DES OPERATIONS'
         KERR(2)='AUGMENTER MXPIOP'
      ELSE
         KERR(1)='LU:(SP FONVA0) SATURATED STACK of OPERATIONS'
         KERR(2)='AUGMENT MXPIOP'
      ENDIF
      CALL LEREUR

 9900 NCODEV = 0
      DBLVAL = 0D0

C     REMISE EN L'ETAT INITIAL DE LA PILE DES CONSTANTES
 9999 NBCTE = NBCT0
      RETURN
      END
