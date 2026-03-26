C
C     FICHIER ~/incl/kmopbi.inc
C
      PARAMETER   ( N1OPBI=-1, NBOPBI=20 )
      INTEGER       NCOPBI(N1OPBI:NBOPBI)
      CHARACTER*6   NMOPBI(N1OPBI:NBOPBI)
C
      DATA          NMOPBI / 'OR ' , 'AND ',
     %                       'OX ' , 'OU ' , 'ET '  , '<'  , '<=' ,
     %                       '='   , '<>'  , '>='   , '>'  , '+'  ,
     %                       '-'   , '*'   , '/'    , '**' , 'MIN',
     %                       'MAX' , 'MOD' , 'ISIGN', 'NUOBJT' ,
     %                       'DIST2P' /
C
      DATA          NCOPBI / 3 , 4 , 
     %                       3 , 3 , 3 , 1 , 2 , 1 , 2 , 2 , 1 , 1 ,
     %                       1 , 1 , 1 , 2 , 3 , 3 , 3 , 5 , 6 , 6 /
C
