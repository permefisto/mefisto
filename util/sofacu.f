      SUBROUTINE SOFACU( NCOGEL, NBSOFA, NOSOFA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE CODE GEOMETRIQUE NCOGEL DU "CUBE"
C -----    (5:TETRAEDRE, 6:PENTAEDRE, 7:HEXAEDRE, 8:6-CUBE, 9:PYRAMIDE)
C          LE NUMERO DES SOMMETS DES FACES DU "CUBE"
C          "CUBE" SIGNIFIE TETRAEDRE OU PENTAEDRE OU HEXAEDRE OU 6-CUBE
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
C NBSOFA : NOMBRE DE SOMMETS DES FACES DE L'EF DE TYPE NCOGEL
C NOSOFA : NOSOFA(I,J) NO DU I-EME SOMMET DE LA FACE J DU "CUBE"
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C ......................................................................
      include"./incl/gsmenu.inc"
      INTEGER        NBSOFA(6), NOSOFA(4,6)
C
C     ******************************************************************
      GOTO ( 1, 1, 1, 1, 500, 600, 700, 800, 900, 1 ) , NCOGEL
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'SOFACU: CODE GEOMETRIQUE='//KERR(MXLGER)(1:4)
     %          //' AU LIEU D ETRE COMPRIS ENTRE 5 ET 8'
      CALL LEREUR
      NBSOFA(1) = 0
      RETURN
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
C     LE NOMBRE DE SOMMETS DES FACES
      NBSOFA(1) = 3
      NBSOFA(2) = 3
      NBSOFA(3) = 3
      NBSOFA(4) = 3
      RETURN
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
C     LE NOMBRE DE SOMMETS DES FACES
      NBSOFA(1) = 3
      NBSOFA(2) = 4
      NBSOFA(3) = 4
      NBSOFA(4) = 3
      NBSOFA(5) = 4
      RETURN
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
C     LE NOMBRE DE SOMMETS DES FACES
      NBSOFA(1) = 4
      NBSOFA(2) = 4
      NBSOFA(3) = 4
      NBSOFA(4) = 4
      NBSOFA(5) = 4
      NBSOFA(6) = 4
      RETURN
C
C     ==================================================================
C     6-CUBE  avec ICI PROJECTION EN XYZ (<=> avec OUBLI de UVW)
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
C     LE NOMBRE DE SOMMETS DES FACES
      NBSOFA(1) = 4
      NBSOFA(2) = 4
      NBSOFA(3) = 4
      NBSOFA(4) = 4
      NBSOFA(5) = 4
      NBSOFA(6) = 4
C
      RETURN
C
C     ==================================================================
C     PYRAMIDE DE BASE QUADRANGULAIRE
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
C     LE NOMBRE DE SOMMETS DES FACES
      NBSOFA(1) = 4
      NBSOFA(2) = 3
      NBSOFA(3) = 3
      NBSOFA(4) = 3
      NBSOFA(5) = 3
      RETURN
      END
