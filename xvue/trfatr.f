      SUBROUTINE TRFATR( COULFACE, COULARETE, NOSOTR, PTXYZD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACE D'UNE FACE TRIANGULAIRE DEFINIE PAR LE NO DE SES 3
C -----   SOMMETS DANS PTXYZD

C ENTREES:
C --------
C COULFACE : NO COULEUR DE TRACE DU TRIANGLE
C COULARETE: NO COULEUR DE TRACE DES ARETES DU TRIANGLE
C NOSOTR : NO DES 3 SOMMETS DANS LE TABLEAU PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE

C MODIFIE:
C --------
C COULFACE : NO COULEUR DE TRACE DU TRIANGLE SI LCRITR=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   AVRIL 2008
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"

      INTEGER           COULFACE, COULARETE, NOSOTR(3)
      DOUBLE PRECISION  PTXYZD(4,*)
      REAL              X(3), Y(3), Z(3), XYZ(3)
      INTRINSIC         REAL
      CHARACTER*8       NMSOMM

C     PROTECTION DE COULFACE POUR APPEL AVEC UNE COULEUR PRIMAIRE
      NCF = COULFACE

C     TRACE DE LA FACE ET DE SES ARETES
      DO I=1,3
         N   = NOSOTR(I)
         X(I) = REAL( PTXYZD(1,N) )
         Y(I) = REAL( PTXYZD(2,N) )
         Z(I) = REAL( PTXYZD(3,N) )
      ENDDO
      CALL FAP13D( NCF, COULARETE, PREDUF, 3, X, Y, Z )

C     TRACE EVENTUEL DU NO DES SOMMETS
      IF( IAVNSO .NE. 0 ) THEN
         DO I=1,3
            XYZ(1) = X(I)
            XYZ(2) = Y(I)
            XYZ(3) = Z(I)
            WRITE( NMSOMM, '(I8)' ) NOSOTR( I )
            CALL SANSBL( NMSOMM, N )
            CALL TEXTE3D( NCNOIR, XYZ, NMSOMM(1:N) )
         ENDDO
      ENDIF

      RETURN
      END
