C     =======================================================================
C     L'OT LS OPPOSE PAR LE SOMMET NOSMMT DE L'OT NUOT EST
C     SOIT, -( NUMERO DU FRERE + 100 * NO TYPE DU PERE ) S'IL EXISTE UN FRERE
C     SINON,   NUMERO DE FILS ( 1 A 14 POUR UN PERE OCTAEDRE  )
C                             ( 1 A 5  POUR UN PERE TETRAEDRE )
C              + 100 * NO TYPE DU PERE DE L'OT OPPOSE
C                ( NO TYPE DE L'OT = 0 SI OCTAEDRE
C                                    1 SI TETRAEDRE T RONDE (T1)
C                                    2 SI TETRAEDRE T       (T2) )
C              DE L'OT OPPOSE PAR LE SOMMET NOSMMT DE L'OT NUOT
C         
C     LA.... EST LE NUMERO D'ARETE POUR L'OT PERE DONT LE SOMMET DE L'OT
C     FILS EST LE MILIEU D'ARETE DU PERE
C
C     OCOC = PERE OCTAEDRE    ET FILS OCTAEDRE
C     T1OC = PERE TETRAEDRE 1 ET FILS OCTAEDRE 
C     T1TE = PERE TETRAEDRE 1 ET FILS TETRAEDRE
C     T2OC = PERE TETRAEDRE 2 ET FILS OCTAEDRE 
C     T2TE = PERE TETRAEDRE 2 ET FILS TETRAEDRE   
C     ======================================================================= 
C
      INTEGER  LSOCOC(1:6,1:6), LSOCTE(1:4,7:14)
      INTEGER  LST1OC(1:6),     LST1TE(1:4,2:5)
      INTEGER  LST2OC(1:6),     LST2TE(1:4,2:5) 
C
      INTEGER  LAOCOC(1:6,1:6), LAOCTE(1:4,7:14)
      INTEGER  LAT1OC(1:6),     LATETE(1:4,2:5)
      INTEGER  LAT2OC(1:6)
C
      DATA LSOCOC /
     &    6,    6,    6,    6,   6,  -6,
     &    4,    4,    4,   -4,   4,   4,
     &    5,    5,    5,    5,  -5,   5,
     &    2,   -2,    2,    2,   2,   2,
     &    3,    3,   -3,    3,   3,   3,
     &   -1,    1,    1,    1,   1,   1 /            
C
      DATA LAOCOC /
     &    0,    1,    2,    3,   4,   0,
     &    1,    0,    5,    0,   8,   9,
     &    2,    5,    0,    6,   0,  10,
     &    3,    0,    6,    0,   7,  11,
     &    4,    8,    0,    7,   0,  12,
     &    0,    9,   10,   11,  12,   0 /
C  
      DATA LSOCTE /
     &   13,  -13,   13,   13,
     &   14,  -14,   14,   14,
     &  -11,   11,   11,   11,
     &  -12,   12,   12,   12,
     &    9,    9,    9,   -9,
     &   10,   10,  -10,   10,
     &    7,    7,   -7,    7,
     &    8,    8,    8,   -8  /       
C  
      DATA LAOCTE /
     &    5,    0,    1,    2,
     &    6,    0,    2,    3,
     &    0,    7,    4,    3,
     &    0,    8,    1,    4,
     &   10,    9,    5,    0,
     &   10,   11,    0,    6,
     &   11,   12,    0,    7,
     &    9,   12,    8,    0  /
C    
      DATA LST1OC /
     &  201,  201,  201,  201,  201,  201 / 
C    
      DATA LAT1OC /
     &    6,    3,    4,    5,    2,    1 /
C  
      DATA LST1TE /
     &  205,  205,  205,  205,
     &  204,  204,  204,  204,
     &  202,  202,  202,  202,
     &  203,  203,  203,  203  /       
C  
      DATA LATETE /
     &    0,    1,    3,    4,
     &    1,    0,    2,    5,
     &    3,    2,    0,    6,
     &    4,    5,    6,    0  /
C    
      DATA LST2OC /
     &  101,  101,  101,  101,  101,  101 / 
C    
      DATA LAT2OC /
     &    6,    2,    3,    4,    5,    1 /
C  
      DATA LST2TE /
     &  104,  104,  104,  104,
     &  105,  105,  105,  105,
     &  103,  103,  103,  103,
     &  102,  102,  102,  102  / 

