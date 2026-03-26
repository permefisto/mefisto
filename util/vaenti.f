      SUBROUTINE VAENTI( NL , NC , NLD , NCD , NLF , NCF , NVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR DE L'ENTIER DEFINI PAR LES CARACTERES
C ----- A PARTIR DE NL,NC DANS KTD
C       ID_ENTIER := |ENTIER  |IDENT |IDENT(ENTIER) |IDENT(IDENT)
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
C           IINFO( 'GRAND' ) SI UN DES ID_ENTIER EST ERRONE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*12      KENTIE
      CHARACTER*1       CAR1
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
C        LE DERNIER CARACTERE DE L'ENTIER
         CALL CARAVA( NLF , NCF )
C
      ELSE
C
C        ID_ENTIER
C        =========
         CALL CARAVA( NLF,NCF )
         CALL CHIDEN( NLF,NCF , NLD,NCD , NLF,NCF , NOIDET )
C        NOIDET NUMERO D'IDENTIFICATEUR DU TABLEAU
         IF( NLD .LE. 0 ) RETURN
         IF( NOIDET .LE. 0 ) THEN
            NLD = 0
            RETURN
         ENDIF
C
C        SI L'IDENTIFICATEUR EST CELUI D'UN TABLEAU
C        ALORS RECHERCHE DE SES INDICES
         I1 = IDENT(2,NOIDET)
         IF( I1 .GT. 0 ) THEN
C
C           IDENT_DE_TABLEAU( ... )
C           =======================
C           LE NOMBRE NB D'INDICES DU TABLEAU
            NBIND = LETAS( I1 )
C
C           LA POSITION DE LA 1-ERE (
            CALL CAR1NB( NLF , NCF )
            IF( KTD(NLF)(NCF:NCF) .NE. '(' ) THEN
               NBLGRC(NRERR) =2
               KERR(1) ='VAENTI:TABLEAU SANS INDICE'
               KERR(2) = KTD(NLF)
               CALL LEREUR
               NLD = 0
               RETURN
            ENDIF
C
C           LECTURE DES INDICES DE CE TABLEAU
            LD     = 0
            NDEIND = 1
            DO 90 II=1,NBIND
C              LES BORNES DE L'INTERVALLE II DU TABLEAU
               I1 = I1 + 2
               N1 = LETAS( I1 - 1 )
               N2 = LETAS( I1     )
C              LECTURE DE L'INDICE A PARTIR DE NLF NCF
               CALL CAR1NB( NLF , NCF )
               CAR1   = KTD(NLF)(NCF:NCF)
               NBK    = 0
               IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
C
C                 INDICE ENTIER
C                 -------------
                  NBK = NBK + 1
                  KENTIE(NBK:NBK) = CAR1
C
C                 LE CARACTERE SUIVANT
 30               CALL CARAPR( NLF , NCF )
                  CAR1 = KTD(NLF)(NCF:NCF)
                  IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
                     NBK = NBK + 1
                     KENTIE(NBK:NBK) = CAR1
                     GOTO 30
                  ENDIF
C
C                 CE CARACTERE N'EST PAS UN CHIFFRE
C                 CADRAGE A GAUCHE DE L'ENTIER
                  J = 12 - NBK
                  DO 40 I=1,NBK
                     KENTIE(J+I:J+I) = KENTIE(I:I)
                     KENTIE(I:I)     = ' '
 40               CONTINUE
C                 LA VALEUR DE L'INDICE
                  READ( KENTIE , '(I12)' ) NVALIN
C                 LE DERNIER CARACTERE DE L'ENTIER
                  CALL CARAVA( NLF , NCF )
C
               ELSE
C
C                 INDICE = IDENT
C                 --------------
                  CALL CARAVA( NLF,NCF )
                  CALL CHIDEN( NLF,NCF , NLD,NCD , NLF,NCF , NOIDEN )
                  IF( NLD .LE. 0 ) RETURN
                  IF( NOIDEN .LE. 0 ) THEN
                     NLD = 0
                     RETURN
                  ENDIF
                  CALL VAIDEN( NOIDEN , NVALIN )
                  IF( NVALIN .EQ. IINFO('GRAND') ) THEN
C                    VALEUR INCORRECTE
                     NLD = 0
                     RETURN
                  ENDIF
               ENDIF
C
C              NVALIN EST IL DANS L'INTERVALLE N1 N2
               IF( NVALIN .LT. N1 .OR. NVALIN .GT. N2 ) THEN
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:4),'(I4)') II
                  WRITE(KERR(MXLGER)(11:22),'(I12)') NVAL
                  KERR(1) = 'VAENTI:INDICE '//KERR(MXLGER)(1:4)//'='
     %                    // KERR(MXLGER)(11:22)
                  WRITE(KERR(MXLGER)( 1: 4),'(I4)') N1
                  WRITE(KERR(MXLGER)(11:14),'(I4)') N2
                  KERR(2) ='NON COMPRIS ENTRE '//KERR(MXLGER)(1:4)//
     %                      KERR(MXLGER)(11:14)
                  CALL LEREUR
                  NLD = 0
                  RETURN
               ENDIF
C
C              LE DECALAGE DANS LE TMS
               LD = LD + ( NVAL - N1 ) * NDEIND
               NDEIND = NDEIND * ( N2 - N1 + 1 )
C
C              ICI PRESENCE OBLIGATOIRE DE , OU )
               CALL CAR1NB( NLF , NCF )
 90         CONTINUE
C
C              IDENT_DE_TABLEAU(NVAL) . RECHERCHE DE SA VALEUR
               CALL VAIDEN( NOIDET , NVAL )
C
C           LE DECALAGE FINAL DANS LE TABLEAU
            LD = LD * MOTVAR( IDENT(1,NOIDET) )
C           LHPILE DE CET IDENTIFICATEUR
            N1 = IDENT(0,NOIDET)
C           LHTMS DE CE TABLEAU
            N1 = LAPILE( 0 , N1 )
C           ADRESSE MCN DU TABLEAU DU TMS CONTENANT CET IDENTIFICATEUR
            N2 = MCTAMS(N1)
C           LA VALEUR DE LA VARIABLE
            NVAL = MCN( N2 + IDENT(3,NOIDET) + LD )
C
         ELSE
C
C           IDENT_VARIABLE
C           ==============
            CALL VAIDEN( NOIDET , NVAL )
            IF( NVAL .EQ. IINFO('GRAND') ) THEN
C              VALEUR INCORRECTE
               NLD = 0
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
C     MISE A JOUR DU SIGNE
      NVAL = NVAL * NSIGNE
      END
