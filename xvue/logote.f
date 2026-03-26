      SUBROUTINE LOGOTE( NCF1, NCF2, NCF3, NCF4, NCA,
     %                   NS1, NS2, NS3, NS4, TETRA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER LES 4 FACES D'UN TETRAEDRE DU LOGO DE MEFISTO
C -----
C
C ENTREES:
C --------
C NCFi   : NUMERO DE LA COULEUR DE LA FACE i=1,...,4
C NCA    : NUMERO DE LA COULEUR DES ARETES
C NSi    : NUMERO DES 4 SOMMETS DU TETRADRE DANS TETRA
C TETRA  : 3 COORDONNEES DES 10 SOMMETS ET MILIEUX DU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    JUILLET 1998
C2345X7..............................................................012
      REAL           XYZ(3,3), TETRA(3,10)
C
      DO 11 K=1,3
         XYZ(K,1) = TETRA(K,NS1)
         XYZ(K,2) = TETRA(K,NS2)
         XYZ(K,3) = TETRA(K,NS3)
 11   CONTINUE
      CALL FACE3D( NCF1, NCA, 3, XYZ )
C
      DO 12 K=1,3
         XYZ(K,1) = TETRA(K,NS1)
         XYZ(K,2) = TETRA(K,NS3)
         XYZ(K,3) = TETRA(K,NS4)
 12   CONTINUE
      CALL FACE3D( NCF2, NCA, 3, XYZ )
C
      DO 13 K=1,3
         XYZ(K,1) = TETRA(K,NS1)
         XYZ(K,2) = TETRA(K,NS2)
         XYZ(K,3) = TETRA(K,NS4)
 13   CONTINUE
      CALL FACE3D( NCF3, NCA, 3, XYZ )
C
      DO 14 K=1,3
         XYZ(K,1) = TETRA(K,NS2)
         XYZ(K,2) = TETRA(K,NS3)
         XYZ(K,3) = TETRA(K,NS4)
 14   CONTINUE
      CALL FACE3D( NCF4, NCA, 3, XYZ )
      RETURN
      END
