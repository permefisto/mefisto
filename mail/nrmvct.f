      SUBROUTINE NRMVCT(V1,V)
C ********************************************************F3BI**********
C     BUT : NORMALISATION DU VECTEUR V1       V=V1/MODUL(V1)
C **********************************************************************
      REAL V1(3),V(3),D
C
      D=SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
      IF(D.EQ.0.)THEN
C      PRINT*,'ERREUR EN NORMALISATION DE VECTEUR'
      RETURN
      ENDIF
      V(1)=V1(1)/D
      V(2)=V1(2)/D
      V(3)=V1(3)/D
      END
