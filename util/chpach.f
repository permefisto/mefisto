      SUBROUTINE CHPACH( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KLG
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UNE CHAINE
C       CHAINE := ' ... '{' ...  '}
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KLG DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           INCHANGES EN SORTIE
C SORTIES :
C ---------
C NLD, NCD: POSITION DANS KLG DE LA PREMIERE ' DE LA CHAINE
C           NLD=0 SI PAS DE CHAINE DEBUTANT PAR '
C                    OU DEBORDEMENT DANS LA RECHERCHE
C NLF, NCF: POSITION DANS KLG DE LA DERNIERE ' DE LA CHAINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
C     LE 1ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CARPNB( NLD , NCD )
C
      NNN = INDEX( KLG(NLD)(NCD:NBCALI) , '''' )
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
 20   IF( NCF .GE. NBCALI ) THEN
         NCF = 1
         NLF = NLF + 1
      ELSE
         NCF = NCF + 1
      ENDIF
C
 30   NNN = INDEX( KLG(NLF)(NCF:NBCALI) , '''' )
      IF( NNN .LE. 0 ) THEN
         NLF = NLF + 1
         IF( NLF .GT. LHKLG ) THEN
C           PASSAGE AU DELA DE LA DERNIERE LIGNE INITIALISEE DE KLG
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
      IF( NC0 .GE. NBCALI ) THEN
         NC0 = 1
         NL0 = NL0 + 1
      ELSE
         NC0 = NC0 + 1
      ENDIF
      IF( KLG(NL0)(NC0:NC0) .EQ. '''' ) THEN
C        ''  ON RECOMMENCE
         NLF = NL0
         NCF = NC0
         GOTO 20
      ENDIF
      END
