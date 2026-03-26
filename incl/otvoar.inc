C     =======================================================================
C     L'OT LO VOISIN PAR LA FACE NOARET DE L'OT NUOT EST
C     SOIT, -( NUMERO DU FRERE + 100 * NO TYPE DU PERE ) S'IL EXISTE UN FRERE
C     SINON,   NUMERO DE FILS ( 1 A 14 POUR UN PERE OCTAEDRE  )
C                             ( 1 A 5  POUR UN PERE TETRAEDRE )
C              + 100 * NO TYPE DU PERE DE L'OT VOISIN
C                ( NO TYPE DE L'OT = 0 SI OCTAEDRE
C                                    1 SI TETRAEDRE T RONDE (T1)
C                                    2 SI TETRAEDRE T       (T2) )
C              DE L'OT VOISIN PAR LA FACE NOARET DE L'OT NUOT
C
C     OCOC = PERE OCTAEDRE    ET FILS OCTAEDRE
C     T1OC = PERE TETRAEDRE 1 ET FILS OCTAEDRE 
C     T1TE = PERE TETRAEDRE 1 ET FILS TETRAEDRE
C     T2OC = PERE TETRAEDRE 2 ET FILS OCTAEDRE 
C     T2TE = PERE TETRAEDRE 2 ET FILS TETRAEDRE   
C     ======================================================================= 
C
      INTEGER  LOOCOC(1:12,1:6), LOOCTE(1:6,7:14)
      INTEGER  LOT1OC(1:12),     LOT1TE(1:6,2:5)
      INTEGER  LOT2OC(1:12),     LOT2TE(1:6,2:5)
C
      INTEGER  NFOCOC(1:12,1:6), NFOCTE(1:6,7:14)
      INTEGER  NFT1OC(1:12),     NFT1TE(1:6,2:5)
      INTEGER  NFT2OC(1:12),     NFT2TE(1:6,2:5)
C
      DATA LOOCOC /
     &    4,   5,   2,   3,   201, 101, 201, 101,    -2,  -3,  -4,  -5,
     &    6, 201,  -1, 101,     5,  -3,  -5,   3,     1, 101,  -6, 201,
     &  201,   6, 101,  -1,     4,   2,  -4,  -2,   101,   1, 201,  -6,
     &   -1, 101,   6, 201,    -3,   5,   3,  -5,    -6, 201,   1, 101,
     &  101,  -1, 201,   6,    -2,  -4,   2,   4,   201,  -6, 101,   1,
     &   -2,  -3,  -4,  -5,   101, 201, 101, 201,     4,   5,   2,   3 /
C
      DATA NFOCOC /
     &    0,   0,   0,   0,     1,   2,   3,   4,     0,   0,   0,   0,
     &    0,   1,   0,   4,     0,   0,   0,   0,     0,   5,   0,   8,
     &    1,   0,   2,   0,     0,   0,   0,   0,     5,   0,   6,   0,
     &    0,   2,   0,   3,     0,   0,   0,   0,     0,   6,   0,   7,
     &    4,   0,   3,   0,     0,   0,   0,   0,     8,   0,   7,   0,
     &    0,   0,   0,   0,     5,   6,   7,   8,     0,   0,   0,   0 /
C
      DATA LOOCTE /
     &  -11,   -10,   203,   202,    -8,   205,  
     &  -12,    -7,   102,   103,    -9,   104,  
     &  -13,   203,   -10,    -8,   202,   204,  
     &  -14,   102,    -7,    -9,   103,   105,  
     &  103,   104,   105,   -12,   -14,    -7,  
     &  203,   -13,   -11,   204,   205,    -8,  
     &  102,   -14,   -12,   105,   104,    -9,  
     &  202,   205,   204,   -11,   -13,   -10 /
C
      DATA NFOCTE /
     &    0,   0,   1,   1,   0,   1,
     &    0,   0,   2,   2,   0,   2,
     &    0,   3,   0,   0,   3,   3,
     &    0,   4,   0,   0,   4,   4,
     &    5,   5,   5,   0,   0,   0,
     &    6,   0,   0,   6,   6,   0,
     &    7,   0,   0,   7,   7,   0,
     &    8,   8,   8,   0,   0,   0 /
C    
      DATA LOT1OC /
     &    5,   4,   3,   2,     6,   1,   6,   1,     3,   2,   5,   4 /
C    
      DATA NFT1OC /
     &    3,   3,   2,   2,     3,   4,   2,   1,     1,   4,   4,   1 /
C
      DATA LOT1TE /
     &  204,     8,   202,   203,    10,    13,  
     &  205,   202,     8,    10,   203,    11,  
     &    8,   204,   205,    13,    11,   203,  
     &   10,    11,    13,   205,   204,   202 /
C
      DATA NFT1TE /
     &    0,   1,   0,   0,   4,   3,
     &    0,   0,   1,   4,   0,   2,
     &    1,   0,   0,   3,   2,   0,
     &    4,   2,   3,   0,   0,   0 /
C
      DATA LOT2OC /
     &    3,   2,   5,   4,     1,   6,   1,   6,     5,   4,   3,   2 /
C
      DATA NFT2OC /
     &    2,   3,   3,   2,     1,   3,   4,   2,     1,   1,   4,   4 /
C
      DATA LOT2TE /
     &  105,     9,   103,   102,     7,    14,  
     &  104,   103,     9,     7,   102,    12,  
     &    9,   105,   104,    14,    12,   102,  
     &    7,    12,    14,   104,   105,   103 /
C
      DATA NFT2TE /
     &    0,   1,   0,   0,   4,   3,
     &    0,   0,   1,   4,   0,   2,
     &    1,   0,   0,   3,   2,   0,
     &    4,   2,   3,   0,   0,   0 /
