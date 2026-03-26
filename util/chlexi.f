      SUBROUTINE CHLEXI( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UN LEXIQUE
C       DE NOM GENERE PAR   ~{>NOM}
C       EXEMPLE  ~>LIGNE
C                ~ REMPLACE LE CARACTERE a POUR LEXIQUE ADAM
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU LEXIQUE
C           NLD = 0 SI CE N'EST PAS UN LEXIQUE
C NLF,NCF : POSITION DU DERNIER CARACTERE DU LEXIQUE
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
      N = INDEX( KTD(NLD)(NCD:NCKTD) , '~' )
      IF( N .LE. 0 ) THEN
C        LE 1-ER CARACTERE DOIT ETRE '~'
C        CE N'EST PAS UN LEXIQUE
         NLD = 0
         NCD = 0
         RETURN
      ENDIF
C
C     LE CARACTERE SUIVANT DOIT ETRE ' ' OU FIN DE LIGNE OU '>'
      NCD = NCD + N - 1
      IF( NCD .LT. NCKTD ) THEN
         IF( KTD(NLD)(NCD+1:NCD+1) .NE. ' ' .AND.
     %       KTD(NLD)(NCD+1:NCD+1) .NE. '>' ) THEN
C           LE CARACTERE SUIVANT N'EST PAS >  => ERREUR
            NLD = 0
            NCD = 0
            RETURN
         ENDIF
      ENDIF
C
C     LE DERNIER CARACTERE DANS KTD
      NLF = NLD
      N   = INDEX( KTD(NLD)(NCD:NCKTD) , ' ' )
      IF( N .LE. 0 ) THEN
         NCF = NCKTD
      ELSE
         NCF = NCD + N - 2
      ENDIF
      RETURN
      END
