      SUBROUTINE TGFACU( NCOGEL, NBTGFA, NOTGFA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE CODE GEOMETRIQUE NCOGEL DU "CUBE"
C -----    (5:TETRAEDRE,6:PENTAEDRE,7:HEXAEDRE, 8:6-CUBE)
C          LE NUMERO DES TANGENTES DES FACES DU "CUBE"
C          "CUBE" SIGNIFIE TETRAEDRE OU PENTAEDRE OU HEXAEDRE
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
C NBTGFA : NOMBRE DE TANGENTES DES FACES DE L'EF DE TYPE NCOGEL
C NOTGFA : NOTGFA(I,J) NO DE LA I-EME TANGENTE DE LA FACE J
C          DU "CUBE" SELON L'ORDRE
C
C    TRIANGLE  :
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S1), NT6(S3S2)
C    QUADRANGLE:
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C    TETRAEDRE :
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S4),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S4),
C                NO TANGENT10(S4S1), NT11(S4S2), NT12(S4S3)
C    PENTAEDRE :
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3)
C    HEXAEDRE  :
C                NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
C
C    6-CUBE    : PAS DE TG
C
C    PYRAMIDE  : PAS DE TG CAR SOMMET 5 AVEC 4 ARETES!
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      INTEGER        NBTGFA(6), NOTGFA(8,6)
C
C     ******************************************************************
      GOTO ( 1 , 1 , 1 , 1 , 500 , 600 , 700 , 800 , 800 , 1 ) , NCOGEL
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'TGFACU: CODE GEOMETRIQUE='//KERR(MXLGER)(1:4)
     %       //' AU LIEU D ETRE COMPRIS ENTRE 5 ET 9'
      CALL LEREUR
      NBTGFA(1) = 0
      RETURN
C
C     ==================================================================
C     TETRAEDRE
C     ==================================================================
C
C     FACE 1 S1 S3 S2
  500 NOTGFA(1,1) = 2
      NOTGFA(2,1) = 1
      NOTGFA(3,1) = 8
      NOTGFA(4,1) = 7
      NOTGFA(5,1) = 5
      NOTGFA(6,1) = 4
C
C     FACE 2 S1 S4 S3
      NOTGFA(1,2) = 3
      NOTGFA(2,2) = 2
      NOTGFA(3,2) = 12
      NOTGFA(4,2) = 10
      NOTGFA(5,2) = 7
      NOTGFA(6,2) = 9
C
C     FACE 3 S1 S2 S4
      NOTGFA(1,3) = 1
      NOTGFA(2,3) = 3
      NOTGFA(3,3) = 6
      NOTGFA(4,3) = 5
      NOTGFA(5,3) = 10
      NOTGFA(6,3) = 11
C
C     FACE 4 S2 S3 S4
      NOTGFA(1,4) = 2
      NOTGFA(2,4) = 3
      NOTGFA(3,4) = 9
      NOTGFA(4,4) = 7
      NOTGFA(5,4) = 10
      NOTGFA(6,4) = 12
C
C     LE NOMBRE DE TANGENTES DES FACES
      NBTGFA(1) = 6
      NBTGFA(2) = 6
      NBTGFA(3) = 6
      NBTGFA(4) = 6
      RETURN
C
C     ==================================================================
C     PENTAEDRE
C     ==================================================================
C
C     FACE 1 S1 S3 S2
  600 NOTGFA(1,1) = 2
      NOTGFA(2,1) = 1
      NOTGFA(3,1) = 8
      NOTGFA(4,1) = 7
      NOTGFA(5,1) = 5
      NOTGFA(6,1) = 4
C
C     FACE 2 S1 S4 S6 S3
      NOTGFA(1,2) = 3
      NOTGFA(2,2) = 2
      NOTGFA(3,2) = 11
      NOTGFA(4,2) = 12
      NOTGFA(5,2) = 18
      NOTGFA(6,2) = 16
      NOTGFA(7,2) = 7
      NOTGFA(8,2) = 9
C
C     FACE 3 S1 S2 S5 S4
      NOTGFA(1,3) = 1
      NOTGFA(2,3) = 3
      NOTGFA(3,3) = 6
      NOTGFA(4,3) = 5
      NOTGFA(5,3) = 14
      NOTGFA(6,3) = 15
      NOTGFA(7,3) = 12
      NOTGFA(8,3) = 10
C
C     FACE 4 S4 S5 S6
      NOTGFA(1,4) = 10
      NOTGFA(2,4) = 11
      NOTGFA(3,4) = 13
      NOTGFA(4,4) = 14
      NOTGFA(5,4) = 16
      NOTGFA(6,4) = 17
C
C     FACE 5 S2 S3 S6 S5
      NOTGFA(1,5) = 4
      NOTGFA(2,5) = 6
      NOTGFA(3,5) = 9
      NOTGFA(4,5) = 8
      NOTGFA(5,5) = 17
      NOTGFA(6,5) = 18
      NOTGFA(7,5) = 15
      NOTGFA(8,5) = 13
C
C     LE NOMBRE DE TANGENTES DES FACES
      NBTGFA(1) = 6
      NBTGFA(2) = 8
      NBTGFA(3) = 8
      NBTGFA(4) = 6
      NBTGFA(5) = 8
      RETURN
C
C     ==================================================================
C     HEXAEDRE
C     ==================================================================
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 1 S1 S4 S3 S2
 700  NOTGFA(1,1) =  1
      NOTGFA(2,1) =  2
      NOTGFA(3,1) =  4
      NOTGFA(4,1) =  5
      NOTGFA(5,1) =  7
      NOTGFA(6,1) =  8
      NOTGFA(7,1) = 10
      NOTGFA(8,1) = 11
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 2 S1 S5 S8 S4
      NOTGFA(1,2) =  2
      NOTGFA(2,2) =  3
      NOTGFA(3,2) = 12
      NOTGFA(4,2) = 10
      NOTGFA(5,2) = 22
      NOTGFA(6,2) = 24
      NOTGFA(7,2) = 15
      NOTGFA(8,2) = 14
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 3 S1 S2 S6 S5
      NOTGFA(1,3) =  1
      NOTGFA(2,3) =  3
      NOTGFA(3,3) =  6
      NOTGFA(4,3) =  5
      NOTGFA(5,3) = 17
      NOTGFA(6,3) = 18
      NOTGFA(7,3) = 15
      NOTGFA(8,3) = 13
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 4 S5 S6 S7 S8
      NOTGFA(1,4) = 13
      NOTGFA(2,4) = 14
      NOTGFA(3,4) = 16
      NOTGFA(4,4) = 17
      NOTGFA(5,4) = 19
      NOTGFA(6,4) = 20
      NOTGFA(7,4) = 22
      NOTGFA(8,4) = 23
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 5 S2 S3 S7 S6
      NOTGFA(1,5) =  4
      NOTGFA(2,5) =  6
      NOTGFA(3,5) =  9
      NOTGFA(4,5) =  8
      NOTGFA(5,5) = 20
      NOTGFA(6,5) = 21
      NOTGFA(7,5) = 18
      NOTGFA(8,5) = 16
C
C     LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 6 S4 S8 S7 S3
      NOTGFA(1,6) = 11
      NOTGFA(2,6) = 12
      NOTGFA(3,6) =  9
      NOTGFA(4,6) =  7
      NOTGFA(5,6) = 15
      NOTGFA(6,6) = 21
      NOTGFA(7,6) = 24
      NOTGFA(8,6) = 23
C
C     LE NOMBRE DE TANGENTES DES FACES
      NBTGFA(1) = 8
      NBTGFA(2) = 8
      NBTGFA(3) = 8
      NBTGFA(4) = 8
      NBTGFA(5) = 8
      NBTGFA(6) = 8
      RETURN
C
C     ==================================================================
C     6-CUBE ou PYRAMIDE
C     ==================================================================
 800  NBTGFA(1) = 0
      RETURN
      END
