      SUBROUTINE CHTYPV( NL , NC , NLD , NCD , NLF , NCF , NOTYPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER NOTYPE LE NUMERO DU TYPE DE LA VARIABLE
C ----- A PARTIR DE (NL,NC) INITIAL DANS KTD
C
C       TYPE_VAR := |entier { ( INTERVAL_ENT : CHAINE
C                            {, INTERVAL_ENT : CHAINE } ) }
C                   |^LEXIQUE   $ ASSOCIE AU NO DANS CE LEXIQUE
C                   |caractere
C                   |reel
C                   |reel2
C                   |xyz
C                   |tms NOM_TMS
C
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU TYPE_TAB
C           NLD=0 SI CE N'EST PAS UN TYPE_VAR
C NLF,NCF : POSITION DU DERNIER CARACTERE DU MOT CLE DE TYPE_TAB
C           C-A-D ) OU r ou l ou 2 ou z ou e
C SORTIE :
C --------
C NOTYPE : LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TMS      =>21
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     RECHERCHE DU MOT CLE INITIAL
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
      IF( KTD(NLD)(NCD:NCD) .EQ. '^' ) THEN
C
C        ^LEXIQUE  . PROTECTION
         CALL CHLEXI( NLD , NCD , N1 , N2 , NLF , NCF )
         IF( N1 .LE. 0 ) THEN
            NLD = 0
            RETURN
         ENDIF
         NOTYPE = 11
         RETURN
C
      ELSE
C
C        LE MOT CLE
         CALL CARAVA( NLD , NCD )
         CALL CHNOM ( NLD , NCD , NLD , NCD , NLF , NCF )
         IF( KTD(NLF)(NCF:NCF) .EQ. ';' ) THEN
C           SUPPRESSION DU ;
            NCF = NCF - 1
         ENDIF
C
         IF( KTD(NLD)(NCD:NCF) .EQ. 'entier' ) THEN
C
C           entier { ( INTERVAL_ENT:CHAINE {, INTERVAL_ENT:CHAINE } ) }
c           RECHERCHE DES ( )
            NL1 = NLF
            NC1 = NCF
            CALL CHPAOF( NL1 , NC1 , N1 , N2 , NLF , NCF )
C           SI N1=0 PAS DE ( EN 1-ER CARACTERE
C           ET NLF NCF EST INCHANGE EN SORTIE
C           SINON NLF NCF POINTE SUR )
            NOTYPE = 4
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'reel2' ) THEN
            NOTYPE = 6
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'reel' ) THEN
            NOTYPE = 5
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'caractere' ) THEN
            NOTYPE = 2
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'xyz' ) THEN
            NOTYPE = 12
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'typeobjet' ) THEN
            NOTYPE = 13
C
         ELSE  IF( KTD(NLD)(NCD:NCF) .EQ. 'tms' ) THEN
            NOTYPE = 21
C
         ELSE
            NBLGRC(NRERR) = 1
            KERR(1) = 'CHTYPV:TYPE VARIABLE INCORRECT'
            CALL LEREUR
            IF( NLD .GT. 0 .AND. NCD.LE.NCKTD ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = KTD(NLD)(NCD:NCKTD)
               CALL LEREUR
            ENDIF
            NOTYPE = 0
            NLD    = 0
         ENDIF
      ENDIF
      END
