      SUBROUTINE VALMAX( NBV, VAL, XYZ )
      DOUBLE PRECISION VAL(NBV), VMX,V
      REAL             XYZ(3,NBV)
C
      VMX = -1D100
      NB = 0
      DO N=1,NBV
         V = ABS( VAL(N) )
         IF( V .GT. 100D0 ) THEN
            nb = nb + 1
ccc            print *,'Vit(',N,')=',VAL(N),' en XYZ=',(XYZ(K,N),K=1,3)
         ENDIF
         IF( V .GT. VMX ) THEN
            VMX = V
            NMX = N
         ENDIF
      ENDDO
C
      print *,'VALMAX Vit(',NMX,')=',VAL(NMX),' en XYZ=',
     %        (XYZ(K,NMX),K=1,3),'  ',nb,' valeurs > 100'
C
      RETURN
      END
