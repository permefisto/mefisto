      SUBROUTINE TRARTR( COULARETE, NOSOTR, PTXYZD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACE DES 3 ARETES D'UN TRIANGLE DEFINI PAR LE NO DES 3 SOMMETS
C -----   DANS PTXYZD
C
C ENTREES:
C --------
C COULARETE: NO COULEUR DE TRACE DES ARETES DU TRIANGLE
C NOSOTR   : NO DES 3 SOMMETS DANS LE TABLEAU PTXYZD
C PTXYZD   : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   AVRIL 2008
C2345X7..............................................................012
      include"./incl/mecoit.inc"
C
      INTEGER           COULARETE, NOSOTR(3)
      DOUBLE PRECISION  PTXYZD(4,*)
      REAL              XYZ(3,4)
C
C     LES 4 SOMMETS DU TRIANGLE
      DO 10 I=1,3
         NS   = NOSOTR(I)
         XYZ(1,I) = REAL( PTXYZD(1,NS) )
         XYZ(2,I) = REAL( PTXYZD(2,NS) )
         XYZ(3,I) = REAL( PTXYZD(3,NS) )
 10   CONTINUE
      XYZ(1,4) = XYZ(1,1)
      XYZ(2,4) = XYZ(2,1)
      XYZ(3,4) = XYZ(3,1)
C
      CALL TRAITS3D( COULARETE, 4, XYZ )
C
      RETURN
      END
