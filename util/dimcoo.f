      SUBROUTINE DIMCOO( MXCOOR, COOR, NDIM )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LA DIMENSION DE L'ESPACE DES POINTS DEFINIS PAR LEURS
C -----    COORDONNEES COOR(3,MXCOOR)
C
C ENTREES:
C --------
C MXCOOR : NOMBRE MAXIMAL DE POINTS
C COOR   : COORDONNEES DES MXCOOR POINTS
C
C SORTIE :
C --------
C NDIM   : DIMENSION DE L'ESPACE DES POINTS ( 1 ou 2 ou 3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS         DECEMBRE 1984
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY     JUIN 2009
C.......................................................................
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL COOR(3,MXCOOR)
C
      DO 20 NDIM=3,1,-1
         DO 10 I=1,MXCOOR
            IF( COOR(NDIM,I) .NE. 0. ) RETURN
 10      CONTINUE
 20   CONTINUE
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) NDIM,' est la DIMENSION des ',MXCOOR,' POINTS'
      ELSE
         WRITE(IMPRIM,*) NDIM,' is the DIMENSION of ',MXCOOR,' POINTS'
      ENDIF
      RETURN
      END
