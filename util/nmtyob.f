      CHARACTER*10 FUNCTION NMTYOB( NTYOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NOM DE L'OBJET DE NUMERO TYPE NTYOB
C -----
C
C ENTREE :
C --------
C NTYOB  : NUMERO CORRESPONDANT de 1 A 7
C          SI NUMERO 15 ALORS SA VALEUR DEVIENT 7
C
C SORTIE :
C --------
C NMTYOB : NOM DE L'OBJET 'POINT'  'LIGNE'  'SURFACE'
C                         'VOLUME' 'OBJET'  'FONCTION' 'TRANSFO'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  OCTOBRE 1988
C.....................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*9       NOMOBJ(1:7)
      DATA              NOMOBJ/'POINT    ' , 'LIGNE    ' ,
     %                         'SURFACE  ' , 'VOLUME   ' , 'OBJET  ',
     %                         'FONCTION ' , 'TRANSFO  ' /
C
C     NTYOB EST IL CORRECT ?
      NTY = NTYOB
      IF( NTYOB .EQ. 15 ) NTY = 7
      IF( NTY .GT. 0 .AND. NTY .LE. 7 ) THEN
C
C        LE NOM DE L'OBJET
         NMTYOB =  NOMOBJ( NTY )
      ELSE
C
C        ERREUR
         NMTYOB = 'ERREUR '
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NTYOB
         KERR(1) = 'TYPE OBJET:'//KERR(MXLGER)(1:4)//' <0 OU >7'
         CALL LEREUR
      ENDIF

      RETURN
      END
