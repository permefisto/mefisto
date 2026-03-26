      SUBROUTINE CHNOM( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UN NOM DEFINI PAR
C       NOM := LETTRE{ALFANUM}  SUR UNE SEULE LIGNE
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU NOM
C           NLD = 0 SI LE 1-ER CARACTERE N'EST PAS UNE LETTRE
C NLF,NCF : POSITION DU DERNIER CARACTERE DU NOM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      CHARACTER*1 CAR
      LOGICAL     LETTRE
C
C     LE 1ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
      CAR = KTD(NLD)(NCD:NCD)
      IF( .NOT. LETTRE( CAR ) ) THEN
C        CE N'EST PAS UNE LETTRE
         NLD = 0
         RETURN
      ENDIF
C
C     RECHERCHE DU PREMIER BLANC
      NCF = INDEX( KTD(NLD)(NCD:NCKTD) , ' ' )
      IF( NCF .LE. 0 ) THEN
         NCF = NCKTD
      ELSE
         NCF = NCD + NCF - 2
      ENDIF
      NLF = NLD
      END
