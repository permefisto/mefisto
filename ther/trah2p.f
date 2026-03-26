      SUBROUTINE TRAH2P
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DU PROTON (0,-a,0) ET DU PROTON (0,a,0) DE H2+
C -----  TRACE DU REPERE CENTRE DE VECTEURS DE LONGUEUR 1.0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET LABORATOIRE J-L LIONS PARIS ET TAMU  JUIN 2005
C23456---------------------------------------------------------------012
      REAL       a
      PARAMETER (a=4.0)
C     A EST LA DEMI-DISTANCE ENTRE LES NOYAUX REDUITS ICI A 1 PROTON
C
      include"./incl/trvari.inc"
      REAL  ZERO(3), XYZ1(3), XYZ2(3)
C
C     LE REPERE
      CALL XVEPAISSEUR( 5 )
      ZERO(1) = 0
      ZERO(2) = 0
      ZERO(3) = 0
C
C     L'AXE DES X EN ROUGE
      XYZ2(1) = 1.0
      XYZ2(2) = 0
      XYZ2(3) = 0
      CALL TRAIT3D( NCROUG, ZERO, XYZ2 )
C
C     L'AXE DES Y EN VERT
      XYZ2(1) = 0
      XYZ2(2) = 1.0
      XYZ2(3) = 0
      CALL TRAIT3D( NCVERT, ZERO, XYZ2 )
C
C     L'AXE DES Z EN BLEU
      XYZ2(1) = 0
      XYZ2(2) = 0
      XYZ2(3) = 1.0
      CALL TRAIT3D( NCBLEU, ZERO, XYZ2 )
C
C     L'AXE ENTRE LES 2 PROTONS
      CALL XVEPAISSEUR( 1 )
      XYZ1(1) =  0
      XYZ1(2) = -a
      XYZ1(3) =  0
C
      XYZ2(1) =  0
      XYZ2(2) =  a
      XYZ2(3) =  0
      CALL TRAIT3D( NCBLAN, XYZ1, XYZ2 )
C
C     1 PROTON
      CALL SYMBOLE3D( NCNOIR, XYZ1, '+' )
C
C     1 PROTON
      CALL SYMBOLE3D( NCNOIR, XYZ2, '+' )
C
      RETURN
      END
