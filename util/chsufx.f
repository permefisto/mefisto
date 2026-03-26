      SUBROUTINE CHSUFX( NL , NCD , NCF , NCDSUF , NCFSUF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UN SUFFIXE
C       RECHERCHE A PARTIR DE LA FIN ET
C       DONT LE SEPARATEUR EST KSUFIX ET EST SUIVI DE '>' OU ' '
C
C ENTREES :
C ---------
C NL , NCD : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C            LE PREMIER CARACTERE A ANALYSER
C NL , NCF : POSITION DANS KTD DU DERNIER CARACTERE A ANALYSER
C            NL , NCD , NCF RESTENT INCHANGES APRES EXECUTION
C
C SORTIES :
C ---------
C NCDSUF : POSITION DU PREMIER CARACTERE DU SUFFIXE APRES KSUFIX
C          NCDSUF = 0  S'IL N'EXISTE PAS DE SUFFIXE
C NCFSUF : POSITION DU DERNIER CARACTERE DU SUFFIXE >= NCDSUF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    MARS 1989
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     RECHERCHE DU SEPARATEUR ENTRE NOM ET SUFFIXE
      DO 10 NCDSUF = NCF , NCD , -1
         IF( KTD(NL)(NCDSUF:NCDSUF) .EQ. KSUFIX ) GOTO 20
 10   CONTINUE
C     PAS DE SUFFIXE
      NCDSUF = 0
      RETURN
C
C     IL EXISTE UN SUFFIXE
C     RECHERCHE DE '>' OU ' '
 20   NCFSUF = INDEX( KTD(NL)(NCDSUF:NCF) , '>' )
      IF( NCFSUF .GT. 0 ) THEN
          NCFSUF = NCDSUF - 2 + NCFSUF
      ELSE
          NCFSUF = INDEX( KTD(NL)(NCDSUF:NCF) , ' ' )
          IF( NCFSUF .GT. 0 ) THEN
             NCFSUF = NCDSUF - 2 + NCFSUF
          ELSE
             NCFSUF = NCF
          ENDIF
      ENDIF
C
C     VERIFICATION SI NCDSUF  <= NCFSUF
C     LE CARACTERE KSUFIX EST OTE
      NCDSUF = NCDSUF + 1
      IF( NCDSUF .GT. NCFSUF ) THEN
C         PAS DE SUFFIXE
          NCDSUF = 0
      ENDIF
      END
