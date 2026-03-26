      FUNCTION NBSOME ( NCOGEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER LE NOMBRE DE SOMMETS DE L EF DE CODE GEOMETRIQUE
C -----     NCOGEL
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
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
C
      INTEGER NSOME(9)
      DATA    NSOME / 1, 2, 3, 4, 4, 6, 8, 64, 5 /
C
      IF( NCOGEL .LE. 0  .OR. NCOGEL .GT. 9 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='CODE GEOMETRIQUE EF  NCOGEL=' // KERR(MXLGER)(1:4)
            KERR(2) ='CODE GEOMETRIQUE <=0 ou >9 INTERDIT'
         ELSE
            KERR(1) ='FE GEOMETRIC NUMBER  NCOGEL=' // KERR(MXLGER)(1:4)
            KERR(2) ='GEOMETRIC NUMBER <=0 ou >9 FORBIDDEN'
         ENDIF
         CALL LEREUR
         NBSOME = 0
      ELSE
         NBSOME = NSOME( NCOGEL )
      ENDIF
C
      RETURN
      END
