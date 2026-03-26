      SUBROUTINE TRTETRAR( COULARETE, NOSOTE, XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      TRACE DES 6 ARETES D'UN TETRAEDRE DEFINI PAR LE NUMERO
C -----      DE SES 4 SOMMETS DANS XYZSOM
C
C ENTREES:
C --------
C COULARETE: NO COULEUR DE TRACE DES ARETES DU TRIANGLE
C NOSOTE   : NO DES 4 SOMMETS DANS LE TABLEAU XYZSOM
C XYZSOM   : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   AVRIL 2008
C2345X7..............................................................012
      include"./incl/mecoit.inc"

      INTEGER           COULARETE, NOSOTE(4)
      REAL              XYZSOM(3,*)

C     LES 3 ARETES DE LA FACE 1 DU TETRAEDRE
      NS1 = NOSOTE(3)
      DO I = 1, 3
         NS2 = NOSOTE(I)
         CALL TRAIT3D( COULARETE, XYZSOM(1,NS1), XYZSOM(1,NS2) )
         NS1 = NS2
      ENDDO

C     LES 3 ARETES NS4-NS1,NS2,NS3 DU TETRAEDRE
      NS1 = NOSOTE(4)
      DO I = 1, 3
         NS2 = NOSOTE(I)
         CALL TRAIT3D( COULARETE, XYZSOM(1,NS1), XYZSOM(1,NS2) )
      ENDDO

      RETURN
      END
