      SUBROUTINE CHCHSC( NL  , NC  , NLMAX , NCMAX ,
     %                   NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UNE CHAINE
C       CHAINE := ' ... '{' ...  '} RECHERCHEE A PARTIR DE NL NC
C       ET AVANT NLMAX NCMAX
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           INCHANGES EN SORTIE
C NLMAX   : NUMERO DE LA DERNIERE LIGNE DE RECHERCHE
C NCMAX   : NUMERO DE LA DERNIERE COLONNE DE NLMAX DE RECHERCHE
C
C SORTIES :
C ---------
C NLD, NCD: POSITION DANS KTD DE LA PREMIERE ' DE LA CHAINE
C           NLD=0 SI PAS DE CHAINE DEBUTANT PAR '
C                    OU DEBORDEMENT DANS LA RECHERCHE
C NLF, NCF: POSITION DANS KTD DE LA DERNIERE ' DE LA CHAINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  AVRIL 1989
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      NC1 = NCKTD
C     LE 1ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
 10   NNN = INDEX( KTD(NLD)(NCD:NCKTD) , '''' )
      IF( NNN .LE. 0 ) THEN
C        CARACTERE ' NON RETROUVE PASSAGE A LA LIGNE SUIVANTE
         NLD = NLD + 1
         IF( NLD .GT. NLMAX .OR. NLD .GT. LHKTD ) THEN
C           PAS DE CHAINE
            NLD = 0
         ELSE IF( NLD .EQ. NLMAX ) THEN
            NC1 = NCMAX
         ENDIF
         NCD = 1
         GOTO 10
      ENDIF
      NCD = NCD + NNN - 1
C
C     LE CARACTERE ' SUIVANT
      NLF = NLD
      NCF = NCD
C
 20   CALL CARAPR( NLF , NCF )
C
 30   NNN = INDEX( KTD(NLF)(NCF:NC1) , '''' )
      IF( NNN .LE. 0 ) THEN
         NLF = NLF + 1
         IF( NLF .GT. NLMAX .OR. NLF .GT. LHKTD ) THEN
C           PASSAGE AU DELA DE LA DERNIERE LIGNE INITIALISEE DE KTD
            NLD = 0
            RETURN
         ELSE IF( NLF .EQ. NLMAX ) THEN
            NC1 = NCMAX
         ENDIF
         NCF = 1
         GOTO 30
      ELSE
         NCF = NCF + NNN - 1
      ENDIF
C
C     CE CARACTERE EST-IL IMMEDIATEMENT SUIVI DE ' ?
      NL0 = NLF
      NC0 = NCF
      CALL CARAPR( NL0 , NC0 )
      IF( KTD(NL0)(NC0:NC0) .EQ. '''' ) THEN
C        ''  ON RECOMMENCE
         NLF = NL0
         NCF = NC0
         GOTO 20
      ENDIF
      END
