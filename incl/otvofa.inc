C     =======================================================================
C     L'OT LV VOISIN PAR LA FACE NOFACE DE L'OT NUOT EST  
C     SOIT, -( NUMERO DU FRERE + 100 * NO TYPE DU PERE ) S'IL EXISTE UN FRERE
C     SINON,   NUMERO DE FILS ( 1 A 14 POUR UN PERE OCTAEDRE  )
C                             ( 1 A 5  POUR UN PERE TETRAEDRE )
C              + 100 * NO TYPE DU PERE DE L'OT VOISIN
C                ( NO TYPE DE L'OT = 0 SI OCTAEDRE
C                                    1 SI TETRAEDRE T RONDE (T1)
C                                    2 SI TETRAEDRE T       (T2) )
C              DE L'OT VOISIN PAR LA FACE NOFACE DE L'OT NUOT 
C
C     LF.... NUMERO DE LA FACE DANS L'OT PERE POUR POURSUIVRE LA RECHERCHE
C            LORSQUE LE PERE N'EST PAS DU MEME TYPE QUE LE FILS
C
C     OCOC = PERE OCTAEDRE    ET FILS OCTAEDRE
C     T1OC = PERE TETRAEDRE 1 ET FILS OCTAEDRE 
C     T1TE = PERE TETRAEDRE 1 ET FILS TETRAEDRE
C     T2OC = PERE TETRAEDRE 2 ET FILS OCTAEDRE 
C     T2TE = PERE TETRAEDRE 2 ET FILS TETRAEDRE   
C     ======================================================================= 
C
      INTEGER  LVOCOC(1:8,1:6), LVOCTE(1:4,7:14), LFOCTE(1:4,7:14)
      INTEGER  LVT1OC(1:8),     LVT1TE(1:4,2:5),  LFT1OC(1:8)
      INTEGER  LVT2OC(1:8),     LVT2TE(1:4,2:5),  LFT2OC(1:8) 
C     
      DATA LVOCOC /
     &  205,  104,  204,  105,     -7,   -8,   -9,  -10,
     &  203,   -7,  -10,  102,    104,  -11,  -14,  204,
     &  202,  102,   -8,   -7,    105,  204,  -12,  -11,
     &   -8,  103,  202,   -9,    -12,  205,  105,  -13,
     &  -10,   -9,  203,  103,    -14,  -13,  104,  205,
     &  -11,  -12,  -13,  -14,    103,  203,  102,  202  /
C  
      DATA LVOCTE /
     &   -2,   -1,  201,   -3,
     &   -3,   -1,  101,   -4,
     &   -5,  201,   -1,   -4,
     &   -2,  101,   -1,   -5,
     &  101,   -2,   -3,   -6,
     &   -6,   -4,   -3,  201,
     &   -6,   -5,   -4,  101,
     &  201,   -5,   -2,   -6  /  
C  
      DATA LFOCTE /
     &   0,   0,   1,   0,  
     &   0,   0,   2,   0,  
     &   0,   3,   0,   0,  
     &   0,   4,   0,   0,  
     &   5,   0,   0,   0,  
     &   0,   0,   0,   6,  
     &   0,   0,   0,   7,  
     &   8,   0,   0,   0 /
C    
      DATA LVT1OC /
     &   13, -105,   11, -104,    -102,   10, -103,   8  /   
C    
      DATA LFT1OC /
     &   3,  0,  2,  0,   0,  4,  0,  1 /
C                                              
      DATA LVT1TE /
     &    3, -101,    6,    2,
     &    4,    6, -101,    5,
     &    1,    2,    5, -101,
     & -101,    3,    4,    1  /
C    
      DATA LVT2OC /
     & -204,   14, -205,   12,      9, -202,    7, -203  /   
C    
      DATA LFT2OC /
     &   0,  3,  0,  2,    1,  0, 4,  0 /
C  
      DATA LVT2TE /
     &    4, -201,    6,    3,
     &    5,    6, -201,    2,
     &    1,    3,    2, -201,
     & -201,    4,    5,    1  /
