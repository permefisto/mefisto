      INTEGER FUNCTION NBFACE ( NCOGEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE DE FACES DE L ELEMENT DE CODE GEOMETRIQUE
C ----- NCOGEL
C
C ENTREE :
C --------
C NCOGEL : 1 NOEUD
C          2 SEGMENT
C          3 TRIANGLE
C          4 QUADRANGLE
C          5 TETRAEDRE
C          6 PENTAEDRE
C          7 HEXAEDRE
C          8 6-CUBE
C          9 PYRAMIDE
C         >9 => ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
C
      INTEGER NFACE(9)
      DATA    NFACE   / 0 , 0 , 1 , 1 , 4 , 5 , 6 , 6 , 5 /
C                 EN REALITE LE 6-CUBE A 12 FACES MAIS TRACE D'UN HEXA
C
      IF( NCOGEL .LE. 0  .OR. NCOGEL .GT. 9 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CODE GEOMETRIQUE EF NCOGEL=' // KERR(MXLGER)(1:4)
            KERR(2) = '<=0 ou >9 INTERDIT'
         ELSE
            KERR(1) = 'GEOMETRIC FE CODE NCOGEL=' // KERR(MXLGER)(1:4)
            KERR(2) = '<=0 or >9 FORBIDDEN'
         ENDIF
         CALL LEREUR
         NBFACE = 0
      ELSE
         NBFACE = NFACE( NCOGEL )
      ENDIF
C
      RETURN
      END
