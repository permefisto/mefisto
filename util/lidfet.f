      SUBROUTINE LIDFET( NLD , NCD , NLF , NCF , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE ET DECLARER DES NOMS D'ETIQUETTES SEPARES PAR
C ----- DES VIRGULES ET TERMINES PAR ;
C
C       CF ~LU/GRAMMAIRE DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NLD NCD : NUMERO DE LIGNE ET COLONNE DANS KLG DE LA CHAINE DEFETIQ
C
C SORTIES :
C ---------
C NLF NCF : NUMERO DE LIGNE ET COLONNE DANS KLG DU ; FINAL
C NRETOU  : =0 SI PAS D'ERREUR
C           >0 SI UNE ERREUR A ETE RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*1        CAR
C
C     LE DERNIER CARACTERE DE DEFETIQ
      NL = NLD
      NC = NCD + 6
C
C     LECTURE DE L'ITEM SUIVANT
 100  CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) GOTO 9000
C
C     RECHERCHE DU NOM DE L' ETIQUETTE
      NCID = NC
      CAR  = KLG(NL)(NC:NC)
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %           (CAR .GE. 'a' .AND. CAR .LE. 'z') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE 1-ER CARACTERE N'EST PAS  UNE LETTRE LICITE
         NBLGRC(NRERR) = 1
         KERR(1) = 'LU: NOM INCORRECT D''ETIQUETTE'//KLG(NL)(NC:NBCALI)
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     LES CARACTERES SUIVANTS PEUVENT ETRE DES LETTRES OU CHIFFRES
      IF( NC .EQ. NBCALI ) GOTO 300
C     NC N'EST PAS LE DERNIER CARACTERE DE LA LIGNE
 110  NC  = NC + 1
      CAR = KLG(NL)(NC:NC)
C
      IF( CAR .EQ. ' ' .OR. CAR .EQ. ',' .OR. CAR .EQ. ';' ) THEN
C        FIN D'IDENTIFICATEUR
         NC = NC - 1
         GOTO 300
      ENDIF
C
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %           (CAR .GE. 'a' .AND. CAR .LE. 'z') .OR.
     %           (CAR .GE. '0' .AND. CAR .LE. '9') .OR.
     %            CAR .EQ. '_' ) ) THEN
C        LE CARACTERE EST ILLICITE
         NBLGRC(NRERR) = 1
         KERR(1) = 'LU: NOM INCORRECT D''ETIQUETTE' //
     %                  KLG(NL)(NC:NBCALI)
         CALL LEREUR
         GOTO 9000
      ENDIF
C     LE CARACTERE EST LICITE : PASSAGE AU SUIVANT
      GOTO 110
C
C     TRAITEMENT DE FIN DE NOM DE ETIQUETTE KLG(NL)(NCID:NCF)
 300  IF( NBETIQ .GE. MXETIQ ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU:PILE SATUREE DES ETIQUETTES UTILISATEUR'
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     AJOUT DE L'ETIQUETTE
      NBETIQ = NBETIQ + 1
      KETIQ(NBETIQ) = KLG(NL)(NCID:NC)
      WRITE(IMPRIM,10900) KETIQ(NBETIQ)
10900 FORMAT('LU: AJOUT DE L''ETIQUETTE ',A)
C
C     PASSAGE AU CARACTERE NON BLANC SUIVANT
      CALL CARPNB( NL , NC )
      CAR = KLG(NL)(NC:NC)
      IF( CAR .EQ. ',' ) THEN
         GOTO 100
      ELSE IF( CAR .EQ. ';' ) THEN
C        LE ; FINAL
         NLF = NL
         NCF = NC
         NRETOU = 0
         RETURN
      ENDIF
C
C     ERREUR DETECTEE
 9000 NBLGRC(NRERR) = 2
      KERR(1) =
     %'LU: DEFETIQ NOM_ETIQUETTE{,NOM_ETIQUETTE} ; INCORRECT'
      KERR(2) =  KLG(NLD)
      CALL LEREUR
      NRETOU = 1
      END
