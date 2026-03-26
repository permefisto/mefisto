      SUBROUTINE TGAREF( NCOGEL, NOARET, NOTGAR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR LE NUMERO DES 2 TANGENTES DE L'ARETE NOARET
C -----    DE L'EF DE CODE GEOMETRIQUE NCOGEL
C
C ENTREES:
C --------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:POINT      2:SEGMENT    3:TRIANGLE  4:QUADRANGLE
C          5:TETRAEDRE  6:PENTAEDRE  7:HEXAEDRE  8:6-CUBE SI>8:ERREUR
C NOARET : NUMERO DE L'ARETE DE L'EF
C
C SORTIES:
C --------
C NOTGAR : NOTGAR(I) NO DE LA I-EME TANGENTE DE L'ARETE J SELON L'ORDRE
C
C    ARETE     : ARETE : S1S2
C                NO TANGENTE1(S1S2), NT2(S2S1)
C
C    TRIANGLE  : ARETES : S1S2, S2S3, S3S1
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S1), NT6(S3S2)
C
C    QUADRANGLE: ARETES : S1S2 S2S3 S3S4 S4S1
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C
C    TETRAEDRE : ARETES : S1S2, S2S3, S3S1, S1S4, S2S4, S3S4
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S4),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S4),
C                NO TANGENT10(S4S1), NT11(S4S2), NT12(S4S3)
C
C    PENTAEDRE : ARETES : S1S2, S2S3, S3S1, S1S4, S2S5, S3S6, S4S5, S5S6, S6S1
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3)
C
C    HEXAEDRE  : ARETES : S1S2, S2S3, S3S4, S4S1, S1S5, S2S6, S3S7, S4S8,
C                         S5S6, S6S7, S7S8, S8S5
C                NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
C
C    6-CUBE    : IDEM HEXAEDRE => PROJECTION ORTHOGONALE SUR LES XYZ
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C ......................................................................
      include"./incl/gsmenu.inc"
      INTEGER        NOTGAR(2)
C
C     ******************************************************************
      GOTO ( 1, 20, 30, 40, 50, 60, 70, 70, 1 ) , NCOGEL
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'TGAREF:CODE GEOMETRIQUE=' // KERR(MXLGER)(1:4)
     %       // ' INCORRECT'
      CALL LEREUR
      NOTGAR(1) = 0
      NOTGAR(2) = 0
      RETURN
C
C     ARETE
 20   NOTGAR(1) = 1
      NOTGAR(2) = 2
      RETURN
C
C     TRIANGLE
 30   GOTO( 31, 32, 33 ),NOARET
C     ARETE S1S2
 31   NOTGAR(1) = 1
      NOTGAR(2) = 4
      RETURN
C     ARETE S2S3
 32   NOTGAR(1) = 3
      NOTGAR(2) = 6
      RETURN
C     ARETE S3S1
 33   NOTGAR(1) = 5
      NOTGAR(2) = 2
      RETURN
C
C     QUADRANGLE
 40   GOTO( 41, 42, 43, 44 ),NOARET
C     ARETE S1S2
 41   NOTGAR(1) = 1
      NOTGAR(2) = 4
      RETURN
C     ARETE S2S3
 42   NOTGAR(1) = 3
      NOTGAR(2) = 6
      RETURN
C     ARETE S3S4
 43   NOTGAR(1) = 5
      NOTGAR(2) = 8
      RETURN
C     ARETE S4S1
 44   NOTGAR(1) = 7
      NOTGAR(2) = 2
      RETURN
C
C     TETRAEDRE
 50   GOTO( 51, 52, 53, 54, 55, 56 ),NOARET
C     ARETE S1S2
 51   NOTGAR(1) = 1
      NOTGAR(2) = 5
      RETURN
C     ARETE S2S3
 52   NOTGAR(1) = 4
      NOTGAR(2) = 8
      RETURN
C     ARETE S3S1
 53   NOTGAR(1) = 7
      NOTGAR(2) = 2
      RETURN
C     ARETE S1S4
 54   NOTGAR(1) = 3
      NOTGAR(2) = 10
      RETURN
C     ARETE S2S4
 55   NOTGAR(1) = 6
      NOTGAR(2) = 11
      RETURN
C     ARETE S3S4
 56   NOTGAR(1) = 9
      NOTGAR(2) = 12
      RETURN
C
C     PENTAEDRE
 60   GOTO( 61, 62, 63, 64, 65, 66, 67, 68, 69 ),NOARET
C     ARETE S1S2
 61   NOTGAR(1) = 1
      NOTGAR(2) = 5
      RETURN
C     ARETE S2S3
 62   NOTGAR(1) = 4
      NOTGAR(2) = 8
      RETURN
C     ARETE S3S1
 63   NOTGAR(1) = 7
      NOTGAR(2) = 2
      RETURN
C     ARETE S1S4
 64   NOTGAR(1) = 3
      NOTGAR(2) = 12
      RETURN
C     ARETE S2S5
 65   NOTGAR(1) = 6
      NOTGAR(2) = 15
      RETURN
C     ARETE S3S6
 66   NOTGAR(1) = 9
      NOTGAR(2) = 18
      RETURN
C     ARETE S4S5
 67   NOTGAR(1) = 10
      NOTGAR(2) = 14
      RETURN
C     ARETE S5S6
 68   NOTGAR(1) = 13
      NOTGAR(2) = 17
      RETURN
C     ARETE S6S4
 69   NOTGAR(1) = 16
      NOTGAR(2) = 11
      RETURN
C
C     HEXAEDRE
 70   GOTO( 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82 ),NOARET
C     ARETE S1S2
 71   NOTGAR(1) = 1
      NOTGAR(2) = 5
      RETURN
C     ARETE S2S3
 72   NOTGAR(1) = 4
      NOTGAR(2) = 9
      RETURN
C     ARETE S3S4
 73   NOTGAR(1) = 7
      NOTGAR(2) = 11
      RETURN
C     ARETE S4S1
 74   NOTGAR(1) = 10
      NOTGAR(2) = 2
      RETURN
C     ARETE S1S5
 75   NOTGAR(1) = 3
      NOTGAR(2) = 15
      RETURN
C     ARETE S2S6
 76   NOTGAR(1) = 6
      NOTGAR(2) = 18
      RETURN
C     ARETE S3S7
 77   NOTGAR(1) = 8
      NOTGAR(2) = 21
      RETURN
C     ARETE S4S8
 78   NOTGAR(1) = 12
      NOTGAR(2) = 24
      RETURN
C     ARETE S5S6
 79   NOTGAR(1) = 13
      NOTGAR(2) = 17
      RETURN
C     ARETE S6S7
 80   NOTGAR(1) = 16
      NOTGAR(2) = 20
      RETURN
C     ARETE S7S8
 81   NOTGAR(1) = 19
      NOTGAR(2) = 23
      RETURN
C     ARETE S8S5
 82   NOTGAR(1) = 22
      NOTGAR(2) = 14
      RETURN
      END
