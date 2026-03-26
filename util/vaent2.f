      SUBROUTINE VAENT2( NL , NC , NLD , NCD , NLF , NCF ,
     %                   NVAL , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR DE L'ENTIER DEFINI PAR LES CARACTERES
C ----- A PARTIR DE NL,NC DANS KTD
C       ID_ENTIER := |ENTIER  |IDENT
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DE L'ID_ENTIER
C           NLD = 0 SI ENTIER INCORRECT
C NLF,NCF : POSITION DANS KTD DU DERNIER CARACTERE DE ID_ENTIER
C NVAL    : VALEUR ENTIERE DU RESULTAT
C NRETOU  : 0 SI NVAL EST CONNU   C-A-D  ID_ENTIER = ENTIER
C           1 SI NVAL EST INCONNU C-A-D  ID_ENTIER = IDENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      CHARACTER*12  KENTIE
      CHARACTER*1   CAR1
C
C     RECHERCHE DU 1-ER CARACTERE NON BLANC A PARTIR DE (NL,NC)
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
C     LA SUITE EST ELLE LE SIGNE -?
      NBK = 0
      NLF = NLD
      NCF = NCD
      CAR1 = KTD(NLF)(NCF:NCF)
      IF( CAR1 .EQ. '-' ) THEN
         NSIGNE = -1
         NBK    = 1
C        PASSAGE AU CARACTERE SUIVANT
         CALL CARAPR( NLF , NCF )
      ELSE
         NSIGNE = 1
      ENDIF
C
C     LA SUITE EST ELLE UN ENTIER OU UN IDENTIFICATEUR ENTIER ?
      KENTIE = ' '
      CAR1   = KTD(NLF)(NCF:NCF)
      IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
C
C        ENTIER
C        ======
         NBK = NBK + 1
         KENTIE(NBK:NBK) = CAR1
C
C        LE CARACTERE SUIVANT
 10      CALL CARAPR( NLF , NCF )
         CAR1 = KTD(NLF)(NCF:NCF)
         IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
            NBK = NBK + 1
            KENTIE(NBK:NBK) = CAR1
            GOTO 10
         ENDIF
C
C        CE CARACTERE N'EST PAS UN CHIFFRE
C        CADRAGE A GAUCHE DE L'ENTIER
         J = 12 - NBK
         DO 20 I=1,NBK
            KENTIE(J+I:J+I) = KENTIE(I:I)
            KENTIE(I:I)     = ' '
 20      CONTINUE
         READ( KENTIE , '(I12)' ) NVAL
C
C        MISE A JOUR DU SIGNE
         NVAL = NVAL * NSIGNE
C
         NRETOU = 0
C        LE DERNIER CARACTERE DE L'ENTIER
         CALL CARAVA( NLF , NCF )
C
      ELSE
C
C        ID_ENTIER
C        =========
         CALL CARAVA( NLF,NCF )
         CALL CHIDEN( NLF,NCF , NLD,NCD , NLF,NCF , NOIDEN )
         IF( NLD .LE. 0 ) RETURN
         IF( NOIDEN .LE. 0 ) THEN
            NLD = 0
            RETURN
         ELSE
C           ID_ENT = IDENT . IMPOSSIBLE CALCULER NVAL
            NRETOU = 1
            NVAL   = 0
         ENDIF
      ENDIF
      END
