      FUNCTION VECSCOSTH( VEC1, VEC2 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C GOAL   : CALCULATE COSINE THETA WHICH IS BETWEEN VEC1 AND VEC2
C -----
C
C INPUT PARAMETERS :
C --------
C VEC1, VEC2 : THE 3 COORDINATES OF 2 VECTORS
C
C OUTPUT PARAMETER:
C --------
C PTSVECORI : THE INNER PRODUCT OF VEC1 AND VEC2 DIVIDED BY THE
C             MULTIPLICATION OF THEIR LENGTH
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTHOR : Dai-Ni Hsieh National Taiwan University TAIPEI TAIWAN 01/2010
C23456---------------------------------------------------------------012
      REAL  VEC1(3), VEC2(3)
C
      VECSCOSTH = PROSCR( VEC1, VEC2, 3 ) /
     %              sqrt( PROSCR( VEC1, VEC1, 3 ) *
     %                    PROSCR( VEC2, VEC2, 3 )   )
      RETURN
      END
