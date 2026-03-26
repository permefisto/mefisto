      SUBROUTINE LIUTRE( NL , NC , NLF , NCF , DBREEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER LA VALEUR DU REEL DEFINI PAR LES CARACTERES
C -----     A PARTIR DE NL,NC DANS KLG
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KLG DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           ILS RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLF,NCF : POSITION DANS KLG DU DERNIER CARACTERE DU REEL
C           NLF =  0 SI REEL INCORRECT OU ERREUR
C                 -1 SI @ POUR SORTIR
C DBREEL  : VALEUR REELLE DU RESULTAT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION  DBREEL
      CHARACTER*25      KDREEL
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
C     LA LIGNE NE CONTIENT PAS LE REEL . PASSAGE A LA LIGNE SUIVANTE
      NBLGRC(NRINVI) = NBLGRC(NRINVI) + 1
      KINVI( NBLGRC(NRINVI) ) ='UN REEL'
      CALL LIRLIG( I )
      IF( I .EQ. -1 ) THEN
C        ECHAPPATOIRE AVEC @
         NLF = -1
         RETURN
      ENDIF
      NBLGRC(NRINVI) = NBLGRC(NRINVI) - 1
      IF( I .NE. 0 ) THEN
         NLF = 0
         RETURN
      ENDIF
      NLF = LHKLG
      NCF = 1
      GOTO 1
C
C     ECHAPPATOIRE AVEC @
 20   IF( KLG(NLF)(NCF:NCF) .EQ. '@' ) THEN
         NLF = -1
         RETURN
      ENDIF
C
C     LE REEL DEBUTE AVEC LE CARACTERE NON BLANC  NLF,I
C     RECHERCHE DU PROCHAIN BLANC
      NCF  = INDEX( KLG(NLF)(I:NBCALI) , ' ' )
      IF( NCF .NE. 0 ) THEN
         NCF = I - 2 + NCF
      ELSE
         NCF = NBCALI
      ENDIF
C     RECHERCHE DE ;
      NCFF = INDEX( KLG(NLF)(I:NBCALI) , ';' )
      IF( NCFF .NE. 0 ) THEN
         NCF = I - 2 + NCFF
      ENDIF
C
C     LES CARACTERES DE CE FUTUR REEL
      KDREEL = KLG(NLF)(I:NCF)
C
C     TRADUCTION DE LA CHAINE DE CARACTERES EN LA VALEUR REELLE
      READ( KDREEL , '(G25.17)' , IOSTAT=IERR ) DBREEL
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'REEL INCORRECT. A REDONNER'
         CALL LEREUR
         GOTO 1
      ENDIF
      END
