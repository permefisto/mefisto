      SUBROUTINE CHCHAI( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UNE CHAINE
C       CHAINE := ' ... '{' ...  '}
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           INCHANGES EN SORTIE
C SORTIES :
C ---------
C NLD, NCD: POSITION DANS KTD DE LA PREMIERE ' DE LA CHAINE
C           NLD=0 SI PAS DE CHAINE DEBUTANT PAR '
C                    OU DEBORDEMENT DANS LA RECHERCHE
C NLF, NCF: POSITION DANS KTD DE LA DERNIERE ' DE LA CHAINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     LE 1ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
      NNN = INDEX( KTD(NLD)(NCD:NCKTD) , '''' )
      IF( NNN .LE. 0 ) THEN
C        LE 1-ER CARACTERE DOIT ETRE '
C        CE N'EST PAS UNE CHAINE
         NLD = 0
         RETURN
      ENDIF
      NCD = NCD + NNN - 1
C
C     LE CARACTERE ' SUIVANT
      NLF = NLD
      NCF = NCD
C
 20   CALL CARAPR( NLF , NCF )
C
 30   NNN = INDEX( KTD(NLF)(NCF:NCKTD) , '''' )
      IF( NNN .LE. 0 ) THEN
         NLF = NLF + 1
         IF( NLF .GT. LHKTD ) THEN
C           PASSAGE AU DELA DE LA DERNIERE LIGNE INITIALISEE DE KTD
            NLD = 0
            RETURN
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
