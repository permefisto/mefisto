      SUBROUTINE CHINTA( NL  , NC  ,
     %                   NLD , NCD , NLF , NCF , NBINDI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE D'INDICES D'UN TABLEAU
C ----- IDENT ( INTERVAL_ENT { , INTERVAL_ENT } )
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DE (     0 SI 1-ER CARACTERE DIFFERENT DE (
C NLF,NCF : POSITION DE )
C NBINDI  : LE NOMBRE D'INDICES DU TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,UNITE(30)
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
      IF( KTD(NLD)(NCD:NCD) .NE. '(' ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR  PAS DE ( DANS'
         ELSE
            KERR(1) ='ERROR  NO ( in'
         ENDIF
         KERR(2) = KTD(NLD)
         CALL LEREUR
         NLD = 0
         RETURN
      ENDIF
C
C     ICI (NLD,NCD) POSITION DANS KTD DE (
C     NOMBRE DE (
      NBP    = 1
      NBINDI = 1
      NLF    = NLD
      NCF    = NCD
C
C     COMPTAGE DES , EN ELIMINANT LES ( ) INTERNES
 20   CALL CARAPR( NLF , NCF )
      IF( KTD(NLF)(NCF:NCF) .EQ. '(' ) THEN
         NBP = NBP + 1
      ELSE IF( KTD(NLF)(NCF:NCF) .EQ. ')' ) THEN
         NBP = NBP - 1
         IF( NBP .EQ. 0 ) RETURN
      ELSE IF( KTD(NLF)(NCF:NCF) .EQ. ',' ) THEN
         IF( NBP .EQ. 1 ) THEN
C           1-ER NIVEAU DE PARENTHESES
            NBINDI = NBINDI + 1
         ENDIF
      ENDIF
      GOTO 20
C
      END
