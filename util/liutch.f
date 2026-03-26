      SUBROUTINE LIUTCH( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER LE NOM DEFINI PAR LES CARACTERES
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
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE DU NOM
C NLF,NCF : POSITION DANS KLG DU DERNIER CARACTERE DU NOM
C           NLF =  0  SI NOM INCORRECT
C                 -1  SI @ POUR SORTIR
C           NLD = NLF SI NOM CORRECT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     SAUVEGARDE DE NL NC
      NLD = NL
      NCD = NC
C
C     RECHERCHE DU 1-ER CARACTERE NON BLANC A PARTIR DE (NL,NC)
 1    DO 10 I=NCD,NBCALI
         IF( KLG(NLD)(I:I) .NE. ' ' ) GOTO 20
 10   CONTINUE
C
C     LA LIGNE NE CONTIENT PAS DE CARACTERE . PASSAGE A LA LIGNE SUIVANTE
      NBLGRC(NRINVI) = NBLGRC(NRINVI) + 1
      KINVI( NBLGRC(NRINVI) ) ='DES CARACTERES'
      CALL LIRLIG( I )
      IF( I .EQ. -1 ) THEN
C        ECHAPPATOIRE
         NLF = -1
         RETURN
      ENDIF
C
      NBLGRC(NRINVI) = NBLGRC(NRINVI) - 1
      IF( I .NE. 0 ) GOTO 9000
      NLD = LHKLG
      NCD = 1
      GOTO 1
C
C     LE NOM DEBUTE AVEC LE CARACTERE NON BLANC NLD,NCD
C     RECHERCHE DU PROCHAIN BLANC
 20   NCD = I
      IF( KLG(NLD)(NCD:NCD) .EQ. '@' ) THEN
C        ECHAPPATOIRE
         NLF = -1
         RETURN
      ENDIF
C
C     LE DERNIER CARACTERE EST LE CARACTERE AVANT LE PREMIER BLANC
      NLF = NLD
      I   = INDEX( KLG(NLD)(NCD:NBCALI) , ' ' )
      NCF = INDEX( KLG(NLD)(NCD:NBCALI) , ';' )
      IF( NCF .LT. I ) I = NCF
      IF( I .GT. 0 ) THEN
         NCF = NCD + I - 2
      ELSE
         NCF = NBCALI
      ENDIF
      RETURN
C
C     VALEUR INCORRECTE
 9000 NBLGRC(NRERR) = 1
      KERR(1) = 'CARACTERES INCORRECTS. A REDONNER'
      CALL LEREUR
      GOTO 1
      END
