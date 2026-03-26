      SUBROUTINE LIUTEN( NL , NC , NLF , NCF , NVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR DE L'ENTIER DEFINI PAR LES CARACTERES
C ----- A PARTIR DE NL,NC DANS KLG
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KLG DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           ILS RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLF,NCF : POSITION DANS KLG DU DERNIER CARACTERE DE L'ENTIER
C           NLF =  0 SI ENTIER INCORRECT OU ERREUR
C                 -1 SI @ POUR SORTIR
C NVAL    : VALEUR ENTIERE DU RESULTAT SI NLF>0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*12      KENTIE
      CHARACTER*1       CAR1
C
C     SAUVEGARDE DE NL NC
      NLF = NL
      NCF = NC
C
C     RECHERCHE DU 1-ER CARACTERE NON BLANC A PARTIR DE (NL,NC)
 1    DO 10 I=NCF,NBCALI
         IF( KLG(NLF)(I:I) .NE. ' ' ) GOTO 20
 10   CONTINUE
C
C     LA LIGNE NE CONTIENT PAS L'ENTIER . PASSAGE A LA LIGNE SUIVANTE
      NBLGRC(NRINVI) = NBLGRC(NRINVI) + 1
      KINVI( NBLGRC(NRINVI) ) ='UN ENTIER'
      CALL LIRLIG( I )
      IF( I .EQ. -1 ) THEN
C        ECHAPPATOIRE AVEC @
         NLF = -1
         RETURN
      ENDIF
      NBLGRC(NRINVI) = NBLGRC(NRINVI) - 1
      IF( I .NE. 0 ) GOTO 9000
      NLF = LHKLG
      NCF = 1
      GOTO 1
C
C     ECHAPPATOIRE AVEC @
 20   IF( KLG(NLF)(I:I) .EQ. '@' ) THEN
         NLF = -1
         RETURN
      ENDIF
C
C     L'ENTIER DEBUTE AVEC LE CARACTERE NON BLANC  NLF,NCF
C     LA SUITE EST ELLE LE SIGNE -?
      NCF = I
      NBK = 0
      CAR1 = KLG(NLF)(NCF:NCF)
      IF( CAR1 .EQ. '-' .OR. CAR1 .EQ. '+' ) THEN
         IF ( CAR1 .EQ. '-' ) THEN
            NSIGNE = -1
         ELSE
            NSIGNE = 1
         ENDIF
         NBK    = 1
C        PASSAGE AU CARACTERE SUIVANT
         IF( NCF .GE. NBCALI ) THEN
            NCF = 1
            NLF = NLF + 1
         ELSE
            NCF = NCF + 1
         ENDIF
      ELSE
         NSIGNE = 1
      ENDIF
C
C     LA SUITE EST ELLE UN ENTIER NON SIGNE?
      KENTIE = ' '
      CAR1   = KLG(NLF)(NCF:NCF)
      IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
C
C        ENTIER
C        ======
         NBK = NBK + 1
         KENTIE(NBK:NBK) = CAR1
C
C        LE CARACTERE SUIVANT
 30      IF( NCF .GE. NBCALI ) THEN
            NCF = 1
            NLF = NLF + 1
         ELSE
            NCF = NCF + 1
         ENDIF
         CAR1 = KLG(NLF)(NCF:NCF)
         IF( CAR1 .GE. '0' .AND. CAR1 .LE. '9' ) THEN
            NBK = NBK + 1
            KENTIE(NBK:NBK) = CAR1
            GOTO 30
         ENDIF
C
C        CE CARACTERE N'EST PAS UN CHIFFRE
C        CADRAGE A GAUCHE DE L'ENTIER
         J = 12 - NBK
         DO 40 I=1,NBK
            KENTIE(J+I:J+I) = KENTIE(I:I)
            KENTIE(I:I)     = ' '
 40      CONTINUE
         READ( KENTIE , '(I12)' , IOSTAT=IERR ) NVAL
         IF( IERR .NE. 0 ) GOTO 9000
C        LE DERNIER CARACTERE DE L'ENTIER
         IF( NCF .LE. 1 ) THEN
            NCF = NBCALI
            NLF = NLF - 1
         ELSE
            NCF = NCF - 1
         ENDIF
      ELSE
         GOTO 9000
      ENDIF
C
C     MISE A JOUR DU SIGNE
      NVAL = NVAL * NSIGNE
      RETURN
C
C     VALEUR INCORRECTE
 9000 NBLGRC(NRERR) = 1
      KERR(1) = 'ENTIER INCORRECT. A REDONNER'
      CALL LEREUR
      GOTO 1
      END
