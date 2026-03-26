C
C     FICHIER ~/incl/kmopun.inc
C
      PARAMETER   ( N1OPUN=0, NBOPUN=19 )
      INTEGER       NCOPUN( N1OPUN:NBOPUN )
      CHARACTER*6   NMOPUN( N1OPUN:NBOPUN )
C
      DATA          NMOPUN / 'NOT ', 'NON ', '+','-', 'ABS', 'SQRT',
     %                       'LOG' , 'LOG10','EXP' , 'ASIN', 'ACOS',
     %                       'ATAN', 'SIN'  ,'COS' , 'TAN' , 'SINH',
     %                       'COSH', 'TANH' ,'NINT', 'LLIGNE' /
C
      DATA          NCOPUN / 4 , 4 , 1 , 1 , 3 , 4 , 3 , 5 , 3 , 4 , 4 ,
     %                       4 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 6 /
C
