      INTEGER FUNCTION NTYOBJ( KNOMOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RETOURNER LE NUMERO (OU DIMENSION) DE L'OBJET DE NOM KNOMOB
C -----
C
C ENTREES :
C ---------
C KNOMOB : NOM DE L'OBJET 'POINT' 'LIGNE' 'SURFACE' 'VOLUME' 'OBJET'
C                         'FONCTION' 'TRANSFO'
C SORTIE :
C --------
C NTYOBJ : NUMERO COORRESPONDANT  1 OU 2 OU 3 OU 4 OU 5 OU 6 OU 7
C          0 SI KNOMOB N'EST PAS UN OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC SEPTEMBRE 1988
C.....................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     KNOMOB
      CHARACTER*1       NOMOBJ(1:7)
      DATA              NOMOBJ/ 'P','L','S','V','O','F','T'  /
C
C     LE NUMERO DANS LA LISTE NOMOBJ
      DO 10 NTYOBJ = 1 , 7
         IF( NOMOBJ(NTYOBJ) .EQ. KNOMOB(1:1) ) RETURN
 10   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'NTYOBJ: TYPE D''OBJET INCONNU ' // KNOMOB
      ELSE
         KERR(1) = 'NTYOBJ: UNKNOWN TYPE of OBJECT ' // KNOMOB
      ENDIF
      CALL LEREUR
      NTYOBJ = 0
      END
