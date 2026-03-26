      SUBROUTINE VAEXPD( NCODEV , DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXECUTION DE L'ARBRE DE L' EXPRESSION
C -----
C          CF ~/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C SORTIES :
C ---------
C NCODEV : 0 DBLVAL N'EST PAS INITIALISEE EN SORTIE
C          1 DBLVAL   EST     INITIALISEE EN SORTIE
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  DBLVAL

cccC     SAUVEGARDE DE LA VALEUR D'ENTREE DU NOMBRE DE CONSTANTES STOCKEES
cccC     2/8/2014       26/11/2015
ccc      NBCTE0 = NBCTE
C
C     LAPIOP EN <=> AVEC NDERP2
C     LAPINF EN <=> AVEC NBPAOU
      MXPILE = LXEXPD
C
C     LA RACINE EST LSRACI(1)
      NOPERP = LSRACI( 1 )
C     L'OPERATION RACINE
      NCOPEP = NCOPER(1,NOPERP)
C
C     ANALYSE DE L'OPERATION
      IF( NCOPEP .EQ. 3 ) THEN
         DBLVAL = DCTE( NCOPER(3,NOPERP) )
         NCODEV = 1
         GOTO 9000
      ELSE IF( NCOPEP .EQ. 4 ) THEN
         DBLVAL = DVARU( NCOPER(3,NOPERP) )
         NCODEV = 1
         GOTO 9000
      ELSE IF( NCOPEP .EQ. 5 ) THEN
         CALL VAVATS( NCOPER(3,NOPERP) , NCODEV , DBLVAL )
         GOTO 9000
      ENDIF

C     LAPIOP CONTIENT LE POINTEUR SUR L'OPERATION EMPILEE
C     LAPINF CONTIENT LE NOMBRE DE FILS A TRAITER POUR CETTE OPERATION
C     LA PREMIERE OPERATION EST EMPILEE
C     ----------------------------------------------------------------
      LHPILE = 1
      LAPIOP(1) = NOPERP
C     LE NOMBRE DE PARAMETRES DE L'OPERATION NOP
      LAPINF(1) = MIN( 2 , NBPRMT( NCOPEP ) )
C     CODE DE NON AFFICHAGE DE VALEURS
      NCODNA = 1
C
C     *******************************
C     LE PARCOURS PRE-FIXE DE L'ARBRE
C     *******************************
C     TANT QU'IL EXISTE UNE OPERATION A FAIRE ALORS TRAITER CETTE OPERATION
10    LHPILA = LHPILE
      DO 100 IP = LHPILA , 1 , -1
         IF( LAPIOP(IP) .GT. 0 ) THEN
C           L'OPERATION DE HAUT DE PILE
C           SA POSITION DANS NCOPER
            NOPERP = LAPIOP(IP)
C           SON CODE OPERATION
            NCOPEP = NCOPER(1,NOPERP)
C           LE NOMBRE DE FILS (0 1 OU 2) A TRAITER
            NBFILS = LAPINF(IP)
            IF( NBFILS .GT. 0 ) THEN
C
C              LE FILS NBFILS EST EMPILE
C              =========================
               NBFILA = MIN( 2 , NBFILS )
               NOFILS = NCOPER(5-NBFILA,NOPERP)
C
               IF( LHPILE .GE. MXPILE ) GOTO 9800
               LHPILE = LHPILE + 1
               LAPIOP( LHPILE ) = NOFILS
C              LE NOMBRE NBPARA DE PARAMETRES DE L'OPERATION NOFILS
C              EST EMPILE DANS LAPINF
               NCOPEF = NCOPER(1,NOFILS)
               NBPARA = NBPRMT( NCOPEF )
               IF( (NCOPEF .GE. 3 .AND. NCOPEF .LE. 5) .OR.
     %              NCOPEF .EQ. 14 ) THEN
C
C                 CONSTANTE OU VARIABLE OU VALEUR TMS
C                 -----------------------------------
C                 MISE DANS LA PILE DES CONSTANTES = DCTE
                  IF( NBCTE .GE. MXDCTE ) THEN
                     NBLGRC(NRERR) = 1
                     KERR(1) = 'LU: VAEXPD: PILE DCTE SATUREE'
                     CALL LEREUR
                     NCODEV = 0
                     GOTO 9000
                  ENDIF
                  NBCTE = NBCTE + 1
                  IF( NCOPEF .EQ. 3 ) THEN
                     DCTE(NBCTE) = DCTE ( NCOPER(3,NOFILS) )
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 4 ) THEN
                     DCTE(NBCTE) = DVARU( NCOPER(3,NOFILS) )
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 5 ) THEN
                     CALL VAVATS(  NCOPER(3,NOFILS) , NCODEV , DBLVAL )
                     IF( NCODEV .EQ. 0 ) GOTO 9000
                     DCTE(NBCTE) = DBLVAL
                     NBTYVA = 3
                  ELSE IF( NCOPEF .EQ. 14 ) THEN
C                    LE NUMERO DE LA CONSTANTE CHAINE EST EMPILEE
                     DCTE(NBCTE) = NCOPER(3,NOFILS)
                     NBTYVA = 14
                  ELSE
                     NBTYVA = 0
                  ENDIF
C                 LE TYPE DE LA CONSTANTE EST EMPILE
                  LAPIOP( LHPILE ) = - NBTYVA
               ENDIF
               LAPINF(LHPILE) = MIN( 2 , NBPARA )
C              LE NOMBRE DE FILS A TRAITER DE L'OPERATION EST DECREMENTE
               LAPINF(IP) = LAPINF(IP) - 1
               GOTO 10
C
            ELSE
C
C              TOUS LES FILS DE L'OPERATEUR NOPERP ONT ETE TRAITES DONC
C              LES PARAMETRES DE L'OPERATEUR SONT TOUS CALCULES ET
C              AU DESSUS DE LUI DANS LA PILE => CALCUL DE L'OPERATION
C              ========================================================
C              LE NOMBRE NBPARA DE PARAMETRES DE L'OPERATION NOPERP
               NBPARA = NBPRMT( NCOPEP )
               IF( NBPARA .GT. LHPILA-IP ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'VAEXPD: NOMBRE INCORRECT DE PARAMETRES'
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C
C              CALCUL DE L'OPERATION SELON LE NOMBRE DE PARAMETRES
               IF( NCOPEP .LT. 300 ) THEN
C
C                 FONCTION UNAIRE OU BINAIRE STANDARD
                  IF( NBPARA .EQ. 1 ) THEN
C
C                    OPERATEUR UNAIRE CALCULABLE AVEC DCTE(NBCTE)
C                    --------------------------------------------
                     CALL VAOPUN( NCOPEP,DCTE(NBCTE),NCODEV,DBLVAL )
                     IF( NCODEV .EQ. 0 ) GOTO 9000
C
                  ELSE IF( NBPARA .EQ. 2 ) THEN
C
                     IF( NCOPEP .NE. 200 ) THEN
C
C                       OPERATEUR BINAIRE CALCULABLE AVEC DCTE(NBCTE-1:NBCTE)
C                       -----------------------------------------------------
                        CALL VAOPBI( NCOPEP ,
     %                               DCTE(NBCTE-1) , DCTE(NBCTE) ,
     %                               NCODEV , DBLVAL )
                        IF( NCODEV .EQ. 0 ) GOTO 9000
C
                     ELSE
C
C                       NCOPEP = 200    ENTRE_PARAMETRES
C                       --------------------------------
C                       LES PARAMETRES DU DESSUS DE LA PILE SONT ABAISSES
C                       D'UN CRAN
                        DO 5 I=IP+1,LHPILE,1
                           LAPIOP(I-1) = LAPIOP(I)
                           LAPINF(I-1) = LAPINF(I)
 5                      CONTINUE
                        LHPILE = LHPILE - 1
                        GOTO 10
                     ENDIF
                  ENDIF
C
               ELSE
C
C                 FONCTION UTILISATEUR
C                 --------------------
                  NOFONC = MOD(NCOPEP,1000) - 300
C                 LES PARAMETRES SONT DCTE(NBCTE-NBPARA+1:NBCTE)
                  IF( NOFONC .GT. 0 ) THEN
C                    VERITABLE FONCTION UTILISATEUR
                     CALL FONVA0( NOFONC , NBPARA , NBCTE-NBPARA+1 ,
     %                            NCODEV , DBLVAL )
                     IF( NCODEV .EQ. 0 ) GOTO 9000
                  ELSE
C                    AFFICHER PARAM {, PARAM }
                     IPP = NBCTE - NBPARA
                     DO 90 I=IPP+1,MIN(NBCTE,IPP+MXLGER)
                        WRITE(KERR(MXLGER)(1:4),'(I4)') I-IPP
                        WRITE(KERR(MXLGER)(11:35),'(G25.17)') DCTE(I)
                        KERR(I-IPP) = 'LU: AFFICHAGE VALEUR' //
     %                               KERR(MXLGER)(1:4)// ' = ' //
     %                               KERR(MXLGER)(11:35)
CCCC                       AFFICHAGE SUR IMPRIMANTE
CCC                        IF(INTERA.GE.3) WRITE(IMPRIM,'(A)') KERR(I-IPP)
 90                  CONTINUE
                     NBLGRC(NRERR) = NBCTE-IPP
C                    AFFICHAGE SUR L'ECRAN
                     CALL LERESU
C                    LES PARAMETRES SONT ECRASES DANS DCTE
                     NBCTE  = NBCTE  - NBPARA + 1
C                    MISE A ZERO EN CAS D'IMPRESSION DU RESULTAT
                     DCTE( NBCTE ) = 0D0
C                    L'OPERATION FONCTION AFFICHER ET LES PARAMETRES
C                    SONT DEPILES DE LAPIOP
                     LHPILE = LHPILE - NBPARA - 1
C                    AFFICHAGE
                     NCODNA = 0
                     GOTO 10
                  ENDIF
C
               ENDIF
C
C              LE RESULTAT DE L'OPERATION ECRASE L'OPERATION
               NBCTE = NBCTE - NBPARA + 1
               DCTE( NBCTE ) = DBLVAL
C              LES OPERANDES SONT DEPILES
               LHPILE = LHPILE - NBPARA
C              LE CODE OPERATION DANS LA PILE EST ECRASE
C              PAR LE RESULTAT
               LAPIOP( LHPILE ) = -3
               LAPINF( LHPILE ) =  0
               GOTO 10
            ENDIF
         ENDIF
 100  CONTINUE
C
C     C'EST LA RACINE : FIN DU CALCUL
C     ===============================
      DBLVAL = DCTE( NBCTE )
      NCODEV = NCODNA

cccC     LE RESULTAT EST DEPILE DANS DCTE
ccc      NBCTE  = NBCTE - 1

C     RESTAURATION DE LA VALEUR D'ENTREE DU NOMBRE DE CONSTANTES STOCKEES
C     2/8/2014
 9000 NBCTE = NBCTE0
      GOTO 9999
C
C     ERREUR
 9800 NBLGRC(NRERR) = 1
      KERR(1) = 'LU: SP VAEXPD: LAPILE EST SATUREE. AUGMENTER MXPILE'
      CALL LEREUR
 9900 NCODEV = 0

 9999 RETURN
      END
