      SUBROUTINE CHPAOF( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ( ET DERNIER CARACTERE )
C       LES ( ) INTERMEDIAIRES SONT IGNOREES A CONDITION
C       QUE CE SOIENT DES COUPLES
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DE (
C           NLD=0 SI LE 1-ER CARACTERE NON BLANC N'EST PAS (
C NLF,NCF : POSITION DE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C
C     RECHERCHE DU 1-ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
C     LE 1ER CARACTERE EST IL ( ?
      IF( KTD(NLD)(NCD:NCD) .NE. '(' ) THEN
         NLD = 0
         RETURN
      ENDIF
C
C     RECHERCHE DES ( ) .... )
      NLF  = NLD
      NCF  = NCD
      NBPA = 1
C
 20   CALL CARAPR( NLF , NCF )
C
C     SAUT DES CHAINES ' ... '
      IF( KTD(NLF)(NCF:NCF) .EQ. '''' ) THEN
C        RECHERCHE DE ' FINAL
         CALL CARAPR( NLF , NCF )
 30      NNN = INDEX(KTD(NLF)(NCF:NCKTD) , '''' )
         IF( NCF .LE. 0 ) THEN
            NCF = 1
            NLF = NLF + 1
            GOTO 30
         ENDIF
         NCF = NCF + NNN - 1
         GOTO 20
C
      ELSE IF( KTD(NLF)(NCF:NCF) .EQ. '(' ) THEN
C        UNE PARENTHESE OUVRANTE DE PLUS A REFERMER
         NBPA = NBPA + 1
         GOTO 20
C
      ELSE IF( KTD(NLF)(NCF:NCF) .EQ. ')' ) THEN
C        UNE PARENTHESE FERMANTE
         NBPA = NBPA - 1
         IF( NBPA .GT. 0 ) GOTO 20
C
      ELSE
C        CARACTERE DIFFERENT DE ' ( ) . PASSAGE AU SUIVANT
         GOTO 20
      ENDIF
      END
