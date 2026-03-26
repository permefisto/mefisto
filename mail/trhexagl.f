      SUBROUTINE TRHEXAGL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DE L'HEXAEDRE ENGLOBANT
C -----    LES XYZ DES 8 SOMMETS SONT DEFINIS A PARTIR de xyzext.inc
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET Saint PIERRE du PERRAY            Novembre 2020
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"

      REAL     XYZ1(3), XYZ2(3)

C     BOUCLE SUR LES FACES DU TABLEAU LFACES
C     ======================================
      CALL XVEPAISSEUR( 0 )

C     ARETE 1-2 // AXZ X
      XYZ2( 1 ) = COOEXT( 1, 2 )
      XYZ2( 2 ) = COOEXT( 2, 1 )
      XYZ2( 3 ) = COOEXT( 3, 1 )
      CALL TRAIT3D( NCROUG, COOEXT(1,1), XYZ2 )

C     ARETE 2-3 // AXE Y
      XYZ1( 1 ) = COOEXT( 1, 2 )
      XYZ1( 2 ) = COOEXT( 2, 2 )
      XYZ1( 3 ) = COOEXT( 3, 1 )
      CALL TRAIT3D( NCVERT, XYZ2, XYZ1 )

C     ARETE 3-4 // AXZ X
      XYZ2( 1 ) = COOEXT( 1, 1 )
      XYZ2( 2 ) = COOEXT( 2, 2 )
      XYZ2( 3 ) = COOEXT( 3, 1 )
      CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )

C     ARETE 4-1 // AXZ Y
      CALL TRAIT3D( NCVERT, XYZ2, COOEXT(1,1) )

C     ARETE 1-5 // AXE Z
      XYZ2( 1 ) = COOEXT( 1, 1 )
      XYZ2( 2 ) = COOEXT( 2, 1 )
      XYZ2( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCBLEU, COOEXT(1,1), XYZ2 )

C     ARETE 2-6 // AXZ Z
      XYZ1( 1 ) = COOEXT( 1, 2 )
      XYZ1( 2 ) = COOEXT( 2, 1 )
      XYZ1( 3 ) = COOEXT( 3, 1 )

      XYZ2( 1 ) = COOEXT( 1, 2 )
      XYZ2( 2 ) = COOEXT( 2, 1 )
      XYZ2( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )

C     ARETE 3-7 // AXZ Z
      XYZ1( 1 ) = COOEXT( 1, 2 )
      XYZ1( 2 ) = COOEXT( 2, 2 )
      XYZ1( 3 ) = COOEXT( 3, 1 )
      CALL TRAIT3D( NCBLEU, XYZ1, COOEXT(1,2) )

C     ARETE 4-8 // AXZ Z
      XYZ1( 1 ) = COOEXT( 1, 1 )
      XYZ1( 2 ) = COOEXT( 2, 2 )
      XYZ1( 3 ) = COOEXT( 3, 1 )

      XYZ2( 1 ) = COOEXT( 1, 1 )
      XYZ2( 2 ) = COOEXT( 2, 2 )
      XYZ2( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )

C     ARETE 5-6 // AXZ X
      XYZ1( 1 ) = COOEXT( 1, 1 )
      XYZ1( 2 ) = COOEXT( 2, 1 )
      XYZ1( 3 ) = COOEXT( 3, 2 )

      XYZ2( 1 ) = COOEXT( 1, 2 )
      XYZ2( 2 ) = COOEXT( 2, 1 )
      XYZ2( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )

C     ARETE 6-7 // AXZ Y
      CALL TRAIT3D( NCVERT, XYZ2, COOEXT(1,2) )

C     ARETE 7-8 // AXZ X
      XYZ2( 1 ) = COOEXT( 1, 1 )
      XYZ2( 2 ) = COOEXT( 2, 2 )
      XYZ2( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCROUG, COOEXT(1,2), XYZ2 )

C     ARETE 8-5 // AXZ Y
      XYZ1( 1 ) = COOEXT( 1, 1 )
      XYZ1( 2 ) = COOEXT( 2, 1 )
      XYZ1( 3 ) = COOEXT( 3, 2 )
      CALL TRAIT3D( NCVERT, XYZ2, XYZ1 )

      RETURN
      END
