      SUBROUTINE TRAPEZE3DBORD( XYZ1,  XYZ2,  XYZ3,  XYZ4,
     %                          COUL1, COUL2, COUL3, COUL4,
     %                          NC, NF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE TRAPEZE DE SOMMETS XYZ ET DE COULEURS COUL
C -----    SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C          AINSI QUE LES ARETES DU QUADRANGLE SELON LA COULEUR NC
C          ATTENTION (X,Y,Z) EN COORDONNEES OBJETS 3D
C
C ENTREES:
C --------
C XYZi   : 3 COORDONNEES DU SOMMET i
C COULi  : NUMERO REEL DE LA COULEUR AUX 3 SOMMETS DU TRIANGLE
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C NC     : COULEUR DE TRACE DES 4 ARETES
C          SI NC<0  PAS DE TRACE DU BORD
C NF     : SI NF>=0 TRACE DE LA FACE
C          SI NF<0  PAS DE TRACE DE LA FACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :PERRONNET ALAIN ANALYSE NUMERIQUE UPMC LJLL PARIS OCTOBRE 2003
C2345X7..............................................................012
      REAL           XYZ1(3), XYZ2(3), XYZ3(3), XYZ4(3)
      REAL           COUL1,   COUL2,   COUL3,   COUL4
C
C     TRIANGLE 123
      CALL TRIACOU3DBORD( XYZ1,  XYZ2,  XYZ3,
     %                    COUL1, COUL2, COUL3,
     %                    -2, NF )
C
C     TRIANGLE 234
      CALL TRIACOU3DBORD( XYZ1,  XYZ3,  XYZ4,
     %                    COUL1, COUL3, COUL4,
     %                    -2, NF )
C
C     TRACE SELON LA COULEUR NC DES 4 ARETES DU QUADRANGLE
      IF( NC .GE. 0 ) THEN
         CALL TRAIT3D( NC,  XYZ1, XYZ2 )
         CALL TRAIT3D( NC,  XYZ2, XYZ3 )
         CALL TRAIT3D( NC,  XYZ3, XYZ4 )
         CALL TRAIT3D( NC,  XYZ4, XYZ1 )
      ENDIF
      RETURN
      END
