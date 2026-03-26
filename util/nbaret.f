      INTEGER FUNCTION NBARET( NCOGEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NOMBRE D ARETES DE L EF DE CODE GEOMETRIQUE
C -----    NCOGEL
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
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
C
      INTEGER  NARET(9)
      DATA     NARET   / 0 , 1 , 3 , 4 , 6 , 9 , 12 , 12 , 8 /
C                                         FAUX pour le 6-CUBE=HEXA
C
      IF( NCOGEL .LE. 0 .OR. NCOGEL .GT. 9 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NBARET: CODE EF NCOGEL=' //
     %                 KERR(MXLGER)(1:4) // ' <=0 OU >9 INTERDIT'
         ELSE
            KERR(1) = 'NBARET: FINITE ELEMENT NCOGEL=' //
     %                 KERR(MXLGER)(1:4) // ' <=0 or >9 FORBIDDEN'
         ENDIF
         CALL LEREUR
         NBARET = 0
         GOTO 9999
      ENDIF
C
C     LE NOMBRE D'ARETES DE L'EF
      NBARET = NARET( NCOGEL )
C
 9999 RETURN
      END
