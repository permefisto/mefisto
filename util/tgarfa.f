      SUBROUTINE TGARFA( NCOGEL, NOTGAR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE CODE GEOMETRIQUE NCOGEL DE LA FACE
C -----   (3:TRIANGLE,4:QUADRANGLE)
C          LE NUMERO DES TANGENTES DES NCOGEL ARETES DE LA FACE
C
C ENTREES:
C --------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:POINT 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE
C          5:TETRAEDRE 6:PENTAEDRE 7:HEXAEDRE >7:ERREUR
C
C SORTIES:
C --------
C NOTGAR : NOTGAR(I,J) NO DE LA I-EME TANGENTE DE L'ARETE J SELON L'ORDRE
C
C    TRIANGLE  :
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S1), NT6(S3S2)
C    QUADRANGLE:
C                NO TANGENTE1(S1S2), NT2(S1S4),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C ......................................................................
      include"./incl/gsmenu.inc"
      INTEGER        NOTGAR(2,NCOGEL)
C
C     ******************************************************************
      GOTO ( 1 , 1 , 300 , 400 , 1 , 1 , 1 ) , NCOGEL
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'TGARFA:CODE GEOMETRIQUE='//KERR(MXLGER)(1:4)
     %          //' AU LIEU D ETRE COMPRIS ENTRE 3 ET 4'
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     TRIANGLE
C     =================================================================
 300  NOTGAR(1,1) = 1
      NOTGAR(2,1) = 4
      NOTGAR(1,2) = 3
      NOTGAR(2,2) = 6
      NOTGAR(1,3) = 5
      NOTGAR(2,3) = 2
      RETURN
C
C     ==================================================================
C     QUADRANGLE
C     =================================================================
 400  NOTGAR(1,1) = 1
      NOTGAR(2,1) = 4
      NOTGAR(1,2) = 3
      NOTGAR(2,2) = 6
      NOTGAR(1,3) = 8
      NOTGAR(2,3) = 5
      NOTGAR(1,4) = 7
      NOTGAR(2,4) = 2
      RETURN
      END
