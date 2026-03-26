      INTEGER FUNCTION MOTTAB( NUMTYP , NBVARI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE DE MOTS OCCUPES PAR NBVARI VARIABLES DE
C ----- NUMERO DE TYPE NUMTYP
C **********************************************************************
C ATTENTION : CETTE FONCTION EST DEPENDANTE MACHINE ICI VERSION IBM
C **********************************************************************
C ENTREE :
C --------
C NUMTYP : NUMERO COMPRIS ENTRE 1 ET 9 AVEC LA CORRESPONDANCE
C          LOGIQUE  <= 1 CARACTERE<= 2 ENTIER/2 <= 3 ENTIER   <= 4
C          REEL     <= 5 REEL2    <= 6 REEL4    <= 7 COMPLEXE <= 8
C          COMPLEXE2<= 9 MOTS     <=10 ^LEXIQUE <=11 XYZ      <=12
C NBVARI : NOMBRE DE VARIABLES DE CE TYPE
C
C SORTIE :
C --------
C MOTTAB : NOMBRE DE MOTS OCCUPES PAR NBVARI VARIABLES DE TYPE NUMTYP
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     POUR LES TYPES : LOGIQUE ,CARACTER,ENTIER/2=> NOMBRE DE VARIABLES
C                                                   DANS UN MOT
C     POUR LES TYPES : ENTIER ,REEL ,REEL2,REEL4,COMPLEXE,MOTS
C                                                => NOMBRE DE MOTS
C                                                   PAR VARIABLES
C
      IF( NUMTYP .GT. 0  .AND.  NUMTYP .LE. 3 ) THEN
         MOTTAB = ( NBVARI - 1 ) / MOTYPV( NUMTYP ) + 1
      ELSEIF( NUMTYP .GE. 4  .AND.  NUMTYP .LE. NBTYPV ) THEN
         MOTTAB = NBVARI * MOTYPV( NUMTYP )
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUMTYP
         KERR(1) = 'MOTTAB: TYPE INCONNU='//KERR(MXLGER)(1:4)
         CALL LEREUR
         MOTTAB = 0
      ENDIF
      END
