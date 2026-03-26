      SUBROUTINE QUADCOU3DBORD( XYZ1,  XYZ2,  XYZ3,  XYZ4,
     %                          COUL1, COUL2, COUL3, COUL4,
     %                          NC, NF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE QUADRANGLE DE SOMMETS XYZ ET DE COULEURS COUL
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
      REAL           XYZT(3,3), COULB
C
C     LE QUADRANGLE EST DIVISE EN 4 A PARTIR DU BARYCENTRE B
      DO 10 I=1,3
         XYZT(I,1) = (XYZ1(I) + XYZ2(I) + XYZ3(I) + XYZ4(I) ) * 0.25
 10   CONTINUE
      COULB = ( COUL1 + COUL2 + COUL3 + COUL4 ) * 0.25
C
C     TRIANGLE B12
      CALL TRIACOU3DBORD( XYZT,  XYZ1,  XYZ2,
     %                    COULB, COUL1, COUL2,
     %                    NC, NF )
C
C     TRIANGLE B23
      CALL TRIACOU3DBORD( XYZT,  XYZ2,  XYZ3,
     %                    COULB, COUL2, COUL3,
     %                    NC, NF )
C
C     TRIANGLE B34
      CALL TRIACOU3DBORD( XYZT,  XYZ3,  XYZ4,
     %                    COULB, COUL3, COUL4,
     %                    NC, NF )
C
C     TRIANGLE B41
      CALL TRIACOU3DBORD( XYZT,  XYZ4,  XYZ1,
     %                    COULB, COUL4, COUL1,
     %                    NC, NF )
C
C     SI QUADRANGLE CROISE IL FAUT TRACER B24 ET B13
C     TRIANGLE B24
      CALL TRIACOU3DBORD( XYZT,  XYZ2,  XYZ4,
     %                    COULB, COUL2, COUL4,
     %                    NC, NF )
C
C     TRIANGLE B13
      CALL TRIACOU3DBORD( XYZT,  XYZ1,  XYZ3,
     %                    COULB, COUL1, COUL3,
     %                    NC, NF )
C
      RETURN
      END
