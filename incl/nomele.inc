C
C     FICHIER ~/incl/nomele.inc
C
C     LES NOMS DE 2 MOTS DE 4 CARACTERES DES TYPES D'ELEMENTS FINIS DISPONIBLES
C     L'ORDRE DEFINIT LE NUMERO DU TYPE DE L'ELEMENT FINI : DE 1 A NBMXEL
C
      PARAMETER             (NBMXEL=34)
      CHARACTER*4   NOMELE(2,NBMXEL)
      DATA          NOMELE /
     % 'TRIA','AP1D',  'TRIA','AP2C',  'QUAD','AQ1C',  'QUAD','AQ2C',
     % 'TRIA','MT10',  'TRIA','MT21',  'QUAD','MQ10',  'QUAD','MQ21',
     % 'TETR','M3T1',  'TETR','M3T2',  'HEXA','M3H1',  'HEXA','M3H2',
     % 'TRIA','2P1D',  'TRIA','2P2D',  'TRIA','2P2C',  'QUAD','2Q1C',
     % 'QUAD','2Q2D',  'QUAD','2Q2C',  'TETR','3P1D',  'TETR','3P2C',
     % 'PENT','3R1C',  'PENT','3R2C',  'HEXA','3Q1C',  'HEXA','3Q2C',
     % 'TRIA','HD06',  'TRIA','EQ06',  'NOEU','SOMM',  'SEGM','1P1D',
     % 'TRIA','2P1C',  '6CUB','6Q1C',  'PYRA','3PY1',  'PYRA','3PY2',
     % 'SEGM','1P2D', 'SEGM','1P3D'  /
