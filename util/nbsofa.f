      FUNCTION NBSOFA ( NCOGEL , K )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER LE NOMBRE DE SOMMETS DE LA FACE K
C -----     D UN EF DE CODE GEOMETRIQUE NCOGEL
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
C K      : NO DE LA FACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      INTEGER  NPENTA(5)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DATA     NPENTA / 3 , 4 , 4 , 3 , 4 /
C
      IF( NCOGEL .LE. 0 .OR. NCOGEL .GT. 9 ) GOTO 1000
C
      GOTO ( 1000, 1000, 1000, 1000, 50, 60, 70, 80, 90 ), NCOGEL
C
C     ERREUR
 1000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'CODE ELEMENT NCOGEL='//KERR(MXLGER)(1:4)
      CALL LEREUR
C
C     TETRAEDRE
C     =========
   50 NBSOFA = 3
      GOTO 2000
C
C     PENTAEDRE
C     =========
   60 NBSOFA = NPENTA( K )
      GOTO 2000
C
C     HEXAEDRE
C     ========
   70 NBSOFA = 4
      GOTO 2000
C
C     6-CUBE  ICI UNE FACE EST UN 5-CUBE => 2**5 SOMMETS
C     ======
CCC   70 NBSOFA = 32
   80 NBSOFA = 4
      GOTO 2000
C
C     PYRAMIDE
C     ========
   90 IF( K .EQ. 1 ) THEN
         NBSOFA = 4
      ELSE
         NBSOFA = 3
      ENDIF
      GOTO 2000
C
 2000 RETURN
      END
