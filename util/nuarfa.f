      SUBROUTINE NUARFA( NCOGEL, NBARFA, NOARFA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE CODE GEOMETRIQUE NCOGEL (1:POINT,2:SEGMENT,..)
C -----    LE NO ELEMENTAIRE DES ARETES DES FACES DE L'EF
C
C ENTREES:
C --------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:POINT 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE
C          5:TETRAEDRE 6:PENTAEDRE 7:HEXAEDRE 8:6-CUBE 9:PYRAMIDE
C         >9:ERREUR
C
C SORTIES:
C --------
C NBARFA : NBARFA(I) NOMBRE D'ARETES DE LA FACE I DE L'EF
C          I VARIANT DE 1 AU NOMBRE DE FACES DE L'EF
C NOARFA : NOARFA(I,J) NO DE LA I-EME ARETE DE LA FACE J
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C ......................................................................
      INTEGER  NBARFA(6), NOARFA(4,6)
C
C     ******************************************************************
      GOTO ( 1 , 1 , 300 , 400 , 500 , 600 , 700 , 800 , 900, 1 ),NCOGEL
C     ******************************************************************
C
C     ==================================================================
C     POINT OU SEGMENT
C     ==================================================================
 1    NBARFA(1) = 0
      RETURN
C
C     ==================================================================
C     TRIANGLE
C     ==================================================================
  300 NBSOM = 3
C
  310 NBARFA(1) = NBSOM
      DO 320 I=1,NBSOM
         NOARFA(I,1) = I
  320 CONTINUE
      RETURN
C
C     ==================================================================
C     QUADRANGLE
C     ==================================================================
  400 NBSOM = 4
      GOTO 310
C
C     ==================================================================
C     TETRAEDRE
C     ==================================================================
C     FACE 1 S132
  500 NBARFA(1) = 3
      NOARFA(1,1) = 3
      NOARFA(2,1) = 2
      NOARFA(3,1) = 1
C
C     FACE 2 S143
      NBARFA(2) = 3
      NOARFA(1,2) = 4
      NOARFA(2,2) = 6
      NOARFA(3,2) = 3
C
C     FACE 3 S124
      NBARFA(3) = 3
      NOARFA(1,3) = 1
      NOARFA(2,3) = 5
      NOARFA(3,3) = 4
C
C     FACE 4 S234
      NBARFA(4) = 3
      NOARFA(1,4) = 2
      NOARFA(2,4) = 6
      NOARFA(3,4) = 5
      RETURN
C
C     ==================================================================
C     PENTAEDRE
C     ==================================================================
C     FACE 1 S132
  600 NBARFA(1) = 3
      NOARFA(1,1) = 3
      NOARFA(2,1) = 2
      NOARFA(3,1) = 1
C
C     FACE 2 S1463
      NBARFA(2) = 4
      NOARFA(1,2) = 4
      NOARFA(2,2) = 9
      NOARFA(3,2) = 6
      NOARFA(4,2) = 3
C
C     FACE 3 S1254
      NBARFA(3) = 4
      NOARFA(1,3) = 1
      NOARFA(2,3) = 5
      NOARFA(3,3) = 7
      NOARFA(4,3) = 4
C
C     FACE 4 S456
      NBARFA(4) = 3
      NOARFA(1,4) = 7
      NOARFA(2,4) = 8
      NOARFA(3,4) = 9
C
C     FACE 5 S2365
      NBARFA(5) = 4
      NOARFA(1,5) = 2
      NOARFA(2,5) = 6
      NOARFA(3,5) = 8
      NOARFA(4,5) = 5
      RETURN
C
C     ==================================================================
C     HEXAEDRE
C     ==================================================================
C     FACE 1 S1432
  700 NBARFA(1) = 4
      NOARFA(1,1) = 4
      NOARFA(2,1) = 3
      NOARFA(3,1) = 2
      NOARFA(4,1) = 1
C
C     FACE 2 S1584
      NBARFA(2) = 4
      NOARFA(1,2) = 5
      NOARFA(2,2) = 12
      NOARFA(3,2) = 8
      NOARFA(4,2) = 4
C
C     FACE 3 S1265
      NBARFA(3) = 4
      NOARFA(1,3) = 1
      NOARFA(2,3) = 6
      NOARFA(3,3) = 9
      NOARFA(4,3) = 5
C
C     FACE 4 S5678
      NBARFA(4) = 4
      NOARFA(1,4) = 9
      NOARFA(2,4) = 10
      NOARFA(3,4) = 11
      NOARFA(4,4) = 12
C
C     FACE 5 S2376
      NBARFA(5) = 4
      NOARFA(1,5) = 2
      NOARFA(2,5) = 7
      NOARFA(3,5) = 10
      NOARFA(4,5) = 6
C
C     FACE 6 S4873
      NBARFA(6) = 4
      NOARFA(1,6) = 8
      NOARFA(2,6) = 11
      NOARFA(3,6) = 7
      NOARFA(4,6) = 3
      RETURN
C
C     ==================================================================
C     6-CUBE   REDUIT ICI AU PREMIER 3-CUBE...
C     ==================================================================
C     FACE 1
  800 NBARFA(1) = 4
      NOARFA(1,1) = 1
      NOARFA(2,1) = 4
      NOARFA(3,1) = 3
      NOARFA(4,1) = 2
C
C     FACE 2 S1584
      NBARFA(2) = 4
      NOARFA(1,2) = 5
      NOARFA(2,2) = 12
      NOARFA(3,2) = 8
      NOARFA(4,2) = 4
C
C     FACE 3 S1265
      NBARFA(3) = 4
      NOARFA(1,3) = 1
      NOARFA(2,3) = 6
      NOARFA(3,3) = 9
      NOARFA(4,3) = 5
C
C     FACE 4 S5678
      NBARFA(4) = 4
      NOARFA(1,4) = 9
      NOARFA(2,4) = 10
      NOARFA(3,4) = 11
      NOARFA(4,4) = 12
C
C     FACE 5 S2376
      NBARFA(5) = 4
      NOARFA(1,5) = 2
      NOARFA(2,5) = 7
      NOARFA(3,5) = 10
      NOARFA(4,5) = 6
C
C     FACE 6 S4873
      NBARFA(6) = 4
      NOARFA(1,6) = 8
      NOARFA(2,6) = 11
      NOARFA(3,6) = 7
      NOARFA(4,6) = 3
      RETURN
C
C     ==================================================================
C     PYRAMIDE
C     ==================================================================
C     FACE 1 S1432
  900 NBARFA(1) = 4
      NOARFA(1,1) = 4
      NOARFA(2,1) = 3
      NOARFA(3,1) = 2
      NOARFA(4,1) = 1
C
C     FACE S125  S235  S345  S415
      DO 910 I=2,5
         NBARFA(I) = 3
         NOARFA(1,I) = I-1
         NOARFA(2,I) = I+4
         NOARFA(3,I) = I+3
 910  CONTINUE
      NOARFA(2,5) = 5
      RETURN
      END
