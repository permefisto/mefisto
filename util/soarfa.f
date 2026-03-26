      SUBROUTINE SOARFA( NCOGEL, NOSOAR, NOSOFA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FOURNIR SELON LE CODE GEOMETRIQUE NCOGEL (1:POINT,2:SEGMENT,..)
C ----- LE NO ELEMENTAIRE DES SOMMETS DE CHAQUE ARETE ET FACE
C
C ENTREES:
C --------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:POINT 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE
C          5:TETRAEDRE 6:PENTAEDRE 7:HEXAEDRE  8:6-CUBE 9:PYRAMIDE
C         >9:ERREUR
C
C SORTIES:
C --------
C NOSOAR : NOSOAR(I,J) NO DU I-EME SOMMET DE L ARETE J
C NOSOFA : NOSOFA(I,J) NO DU I-EME SOMMET DE LA FACE J
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C ......................................................................
      include"./incl/gsmenu.inc"
      INTEGER        NOSOAR(2,12), NOSOFA(4,6)
C
C     ******************************************************************
      GOTO ( 1, 200, 300, 400, 500, 600, 700, 800, 900, 1 ), NCOGEL
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'SOARFA:CODE GEOMETRIQUE='//KERR(MXLGER)(1:4)
     %          //' AU LIEU D ETRE COMPRIS ENTRE 2 ET 7'
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     SEGMENT
C     ==================================================================
  200 NOSOAR(1,1) = 1
      NOSOAR(2,1) = 2
      RETURN
C
C     ==================================================================
C     TRIANGLE
C     ==================================================================
  300 NBSOM = 3
  310 DO 350 I=1,NBSOM
         NOSOAR(1,I) = I
         NOSOAR(2,I) = I + 1
         NOSOFA(I,1) = I
  350 CONTINUE
      NOSOAR(2,NBSOM) = 1
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
C
C     FACE 1
  500 NOSOFA(1,1) = 1
      NOSOFA(2,1) = 3
      NOSOFA(3,1) = 2
C
C     FACE 2
      NOSOFA(1,2) = 1
      NOSOFA(2,2) = 4
      NOSOFA(3,2) = 3
C
C     FACE 3
      NOSOFA(1,3) = 1
      NOSOFA(2,3) = 2
      NOSOFA(3,3) = 4
C
C     FACE 4
      NOSOFA(1,4) = 2
      NOSOFA(2,4) = 3
      NOSOFA(3,4) = 4
C
C     LES ARETES
      DO 510 I=1,3
          NOSOAR(1,I+3) = I
          NOSOAR(2,I+3) = 4
  510 CONTINUE
      GOTO 300
C
C     ==================================================================
C     PENTAEDRE
C     ==================================================================
C
C     FACE 1
  600 NOSOFA(1,1) = 1
      NOSOFA(2,1) = 3
      NOSOFA(3,1) = 2
C
C     FACE 2
      NOSOFA(1,2) = 1
      NOSOFA(2,2) = 4
      NOSOFA(3,2) = 6
      NOSOFA(4,2) = 3
C
C     FACE 3
      NOSOFA(1,3) = 1
      NOSOFA(2,3) = 2
      NOSOFA(3,3) = 5
      NOSOFA(4,3) = 4
C
C     FACE 4
      NOSOFA(1,4) = 4
      NOSOFA(2,4) = 5
      NOSOFA(3,4) = 6
C
C     FACE 5
      NOSOFA(1,5) = 2
      NOSOFA(2,5) = 3
      NOSOFA(3,5) = 6
      NOSOFA(4,5) = 5
C
C     LES ARETES
      DO 650 I=1,3
         NOSOAR(1,I+3) = I
         NOSOAR(2,I+3) = I + 3
         NOSOAR(1,I+6) = I + 3
         NOSOAR(2,I+6) = I + 4
  650 CONTINUE
      NOSOAR(2,9) = 4
      GOTO 300
C
C     ==================================================================
C     HEXAEDRE
C     ==================================================================
C
C     FACE 1
  700 NOSOFA(1,1) = 1
      NOSOFA(2,1) = 4
      NOSOFA(3,1) = 3
      NOSOFA(4,1) = 2
C
C     FACE 2
      NOSOFA(1,2) = 1
      NOSOFA(2,2) = 5
      NOSOFA(3,2) = 8
      NOSOFA(4,2) = 4
C
C     FACE 3
      NOSOFA(1,3) = 1
      NOSOFA(2,3) = 2
      NOSOFA(3,3) = 6
      NOSOFA(4,3) = 5
C
C     FACE 4
      NOSOFA(1,4) = 5
      NOSOFA(2,4) = 6
      NOSOFA(3,4) = 7
      NOSOFA(4,4) = 8
C
C     FACE 5
      NOSOFA(1,5) = 2
      NOSOFA(2,5) = 3
      NOSOFA(3,5) = 7
      NOSOFA(4,5) = 6
C
C     FACE 6
      NOSOFA(1,6) = 4
      NOSOFA(2,6) = 8
      NOSOFA(3,6) = 7
      NOSOFA(4,6) = 3
C
C     LES ARETES
      DO 750 I=1,4
         NOSOAR(1,I)  = I
         NOSOAR(2,I)  = I + 1
         I1           = I + 4
         NOSOAR(1,I1) = I
         NOSOAR(2,I1) = I + 4
         I1           = I + 8
         NOSOAR(1,I1) = I + 4
         NOSOAR(2,I1) = I + 5
  750 CONTINUE
C
C     MISE A JOUR
      NOSOAR(2,4) = 1
      NOSOAR(2,12) = 5
      RETURN
C
C     ==================================================================
C     6-CUBE => Projection orthogonale dans R3 du 6-cube
C     ==================================================================
C
C     FACE 1
  800 NOSOFA(1,1) = 1
      NOSOFA(2,1) = 3
      NOSOFA(3,1) = 4
      NOSOFA(4,1) = 2
C
C     FACE 2
      NOSOFA(1,2) = 1
      NOSOFA(2,2) = 5
      NOSOFA(3,2) = 7
      NOSOFA(4,2) = 3
C
C     FACE 3
      NOSOFA(1,3) = 1
      NOSOFA(2,3) = 2
      NOSOFA(3,3) = 6
      NOSOFA(4,3) = 5
C
C     FACE 4
      NOSOFA(1,4) = 5
      NOSOFA(2,4) = 6
      NOSOFA(3,4) = 8
      NOSOFA(4,4) = 7
C
C     FACE 5
      NOSOFA(1,5) = 2
      NOSOFA(2,5) = 3
      NOSOFA(3,5) = 8
      NOSOFA(4,5) = 6
C
C     FACE 6
      NOSOFA(1,6) = 3
      NOSOFA(2,6) = 7
      NOSOFA(3,6) = 8
      NOSOFA(4,6) = 4
C
C     LES ARETES
      NOSOAR(1,1) = 1
      NOSOAR(2,1) = 2
      NOSOAR(1,2) = 2
      NOSOAR(2,2) = 4
      NOSOAR(1,3) = 4
      NOSOAR(2,3) = 3
      NOSOAR(1,4) = 3
      NOSOAR(2,4) = 1
C
      NOSOAR(1,5) = 1
      NOSOAR(2,5) = 5
      NOSOAR(1,6) = 2
      NOSOAR(2,6) = 6
      NOSOAR(1,7) = 4
      NOSOAR(2,7) = 8
      NOSOAR(1,8) = 3
      NOSOAR(2,8) = 7
C
      NOSOAR(1, 9) = 5
      NOSOAR(2, 9) = 6
      NOSOAR(1,10) = 6
      NOSOAR(2,10) = 8
      NOSOAR(1,11) = 8
      NOSOAR(2,11) = 7
      NOSOAR(1,12) = 7
      NOSOAR(2,12) = 5
C
      RETURN
C
C     ==================================================================
C     PYRAMIDE DE BASE QUADRANGULAIRE 1234 et de SOMMET 5
C     ==================================================================
C
C     FACE 1
  900 NOSOFA(1,1) = 1
      NOSOFA(2,1) = 4
      NOSOFA(3,1) = 3
      NOSOFA(4,1) = 2
C
C     FACE 2
      NOSOFA(1,2) = 1
      NOSOFA(2,2) = 2
      NOSOFA(3,2) = 5
C
C     FACE 3
      NOSOFA(1,3) = 2
      NOSOFA(2,3) = 3
      NOSOFA(3,3) = 5
C
C     FACE 4
      NOSOFA(1,4) = 3
      NOSOFA(2,4) = 4
      NOSOFA(3,4) = 5
C
C     FACE 5
      NOSOFA(1,5) = 4
      NOSOFA(2,5) = 1
      NOSOFA(3,5) = 5
C
C     LES ARETES
      DO 910 I=1,4
         NOSOAR(1,I  )  = I
         NOSOAR(2,I  )  = I+1
         NOSOAR(1,I+4)  = I
         NOSOAR(2,I+4)  = 5
 910  CONTINUE
      NOSOAR(2,4) = 1
C
      RETURN
      END
