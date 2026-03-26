C
C     DONNEES DES SOUS-TETRAEDRES D'UN TETRAEDRE ou PYRAMIDE ou PENTAEDRE
C     ou HEXAEDRE SUBDIVISE REGULIEREMENT POUR FAIRE DES TRACES AVEC DES
C     POLYNOMES D'INTERPOLATION DE DEGRE 2
C
      INTEGER           NOSTAR(2,6)
C
      INTEGER           NOSY(4,2)
      INTEGER           NOSP(4,3)
      INTEGER           NOSH(4,5)
C
      REAL              TETR(3,4)
      REAL              PYRA(3,5)
      REAL              PENT(3,6)
      REAL              HEXA(3,8)
C
C     NUMERO DES SOMMETS DES 6 ARETES DU TETRAEDRE
      DATA     NOSTAR/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /
C
C     NOSY -   NUMEROS DES SOMMETS DE 2 SOUS-TETRAEDRES PRIMAIRES
C              ISSUS DU DECOUPAGE DE LA PYRAMIDE DE REFERENCE
      DATA     NOSY/ 1, 2, 3, 5,
     %               1, 3, 4, 5 /
C
C     NOSP -   NUMEROS DES SOMMETS DE 3 SOUS-TETRAEDRES PRIMAIRES
C              ISSUS DU DECOUPAGE DU PENTAEDRE DE REFERENCE
      DATA     NOSP/ 1, 2, 3, 4 ,
     %               5, 2, 4, 3 ,
     %               4, 6, 5, 3 /
C
C     NOSH -   NUMEROS DES SOMMETS DE 5 SOUS-TETRAEDRES PRIMAIRES
C              ISSUS DU DECOUPAGE DE L'HEXAEDRE DE REFERENCE
      DATA     NOSH/ 4, 5, 2, 7,
     %               1, 2, 4, 5,
     %               3, 4, 2, 7,
     %               8, 7, 5, 4,
     %               6, 5, 7, 2 /

C
C     COORDONNEES DES SOMMETS DU TETRAEDRE, PENTAEDRE, HEXAEDRE DE REFERENCE
      DATA TETR/ 0.0, 0.0, 0.0,
     %           1.0, 0.0, 0.0,
     %           0.0, 1.0, 0.0,
     %           0.0, 0.0, 1.0 /
C
      DATA PYRA/ 0.0, 0.0, 0.0,
     %           1.0, 0.0, 0.0,
     %           1.0, 1.0, 0.0,
     %           0.0, 1.0, 0.0,
     %           0.0, 0.0, 1.0 /
C
      DATA PENT/ 0.0, 0.0, 0.0,
     %           1.0, 0.0, 0.0,
     %           0.0, 1.0, 0.0,
     %           0.0, 0.0, 1.0,
     %           1.0, 0.0, 1.0,
     %           0.0, 1.0, 1.0 /
C
      DATA HEXA/ 0.0, 0.0, 0.0,
     %           1.0, 0.0, 0.0,
     %           1.0, 1.0, 0.0,
     %           0.0, 1.0, 0.0,
     %           0.0, 0.0, 1.0,
     %           1.0, 0.0, 1.0,
     %           1.0, 1.0, 1.0,
     %           0.0, 1.0, 1.0 /
