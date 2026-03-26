      SUBROUTINE ITEMS0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ANNULE LE NOMBRE DES ITEMS DE TOUT TYPE
C -----    1 POINT 2:LIGNE 3:SURFACE 4:VOLUME 5:OBJET
C          6:SOMMET 7:ARETE 8:FACE EF3D 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1990
C ......................................................................
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)

C     AUCUN ITEM VISIBLE A CET INSTANT
      DO I = 1, 9
C        NOMBRE D'ITEMS NUL
         MCN( MNITEM(I) + 2 ) = 0
      ENDDO

      RETURN
      END
