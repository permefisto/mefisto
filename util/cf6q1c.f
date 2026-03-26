      SUBROUTINE CF6Q1C( POLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 64 COEFFICIENTS DES 64 POLYNOMES DE BASE
C -----    DE L'EF 6CUBE 6Q1C
C SORTIE :
C --------
C POLY   : POLY(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY(I,J,K,L,M,N,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                             U**(L-1) V**(M-1) W**(N-1)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C ......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , NUNITE(30)
      DOUBLE PRECISION  POLY(64,64)
      DOUBLE PRECISION  CI(2), CJ(2), CK(2), CL(2), CM(2), CN(2)
      INTEGER           IC(2), JC(2), KC(2), LC(2), MC(2), NC(2)
C
C     MISE A ZERO DES 64*64 COEFFICIENTS
      CALL AZEROD( 64*64, POLY )
C
      NP = 0
      DO 160 N=1,2
         IF( MOD(N,2) .EQ. 0 ) THEN
            NBNC  = 1
            NC(1) = 2
            CN(1) = 1
         ELSE
            NBNC  = 2
            NC(1) = 1
            CN(1) = 1
            NC(2) = 2
            CN(2) =-1
         ENDIF
         DO 150 M=1,2
            IF( MOD(M,2) .EQ. 0 ) THEN
               NBMC  = 1
               MC(1) = 2
               CM(1) = 1
            ELSE
               NBMC  = 2
               MC(1) = 1
               CM(1) = 1
               MC(2) = 2
               CM(2) =-1
            ENDIF
            DO 140 L=1,2
               IF( MOD(L,2) .EQ. 0 ) THEN
                  NBLC  = 1
                  LC(1) = 2
                  CL(1) = 1
               ELSE
                  NBLC  = 2
                  LC(1) = 1
                  CL(1) = 1
                  LC(2) = 2
                  CL(2) =-1
               ENDIF
               DO 130 K=1,2
                  IF( MOD(K,2) .EQ. 0 ) THEN
                     NBKC  = 1
                     KC(1) = 2
                     CK(1) = 1
                  ELSE
                     NBKC  = 2
                     KC(1) = 1
                     CK(1) = 1
                     KC(2) = 2
                     CK(2) =-1
                  ENDIF
                  DO 120 J=1,2
                     IF( MOD(J,2) .EQ. 0 ) THEN
                        NBJC  = 1
                        JC(1) = 2
                        CJ(1) = 1
                     ELSE
                        NBJC  = 2
                        JC(1) = 1
                        CJ(1) = 1
                        JC(2) = 2
                        CJ(2) =-1
                     ENDIF
                     DO 110 I=1,2
                        IF( MOD(I,2) .EQ. 0 ) THEN
                           NBIC  = 1
                           IC(1) = 2
                           CI(1) = 1
                        ELSE
                           NBIC  = 2
                           IC(1) = 1
                           CI(1) = 1
                           IC(2) = 2
                           CI(2) =-1
                        ENDIF
C
C                       TRAITEMENT DU POLYNOME(I,J,K,L,M,N)
C                       NUMERO DU POLYNOME TRAITE
                        NP = NP + 1
C
C                       SES COEFFICIENTS NON NULS
                        DO 60 NN=1,NBNC
                           DO 50 MM=1,NBMC
                              DO 40 LL=1,NBLC
                                 DO 30 KK=1,NBKC
                                    DO 20 JJ=1,NBJC
                                       DO 10 II=1,NBIC
                                          NUC =       IC(II)
     %                                        +  2 * (JC(JJ)-1)
     %                                        +  4 * (KC(KK)-1)
     %                                        +  8 * (LC(LL)-1)
     %                                        + 16 * (MC(MM)-1)
     %                                        + 32 * (NC(NN)-1)
                                          POLY(NUC,NP) = CI(II) * CJ(JJ)
     %                                                 * CK(KK) * CL(LL)
     %                                                 * CM(MM) * CN(NN)
 10                                    CONTINUE
 20                                 CONTINUE
 30                              CONTINUE
 40                           CONTINUE
 50                        CONTINUE
 60                     CONTINUE
 110                 CONTINUE
 120              CONTINUE
 130           CONTINUE
 140        CONTINUE
 150     CONTINUE
 160  CONTINUE
C
ccc      WRITE(IMPRIM,10500)
ccc10500 FORMAT(/,'LES 64 COEFFICIENTS DES 64 POLYNOMES 6Q1C'/130(1H@)/)
ccc      DO 500 NP=1,64
ccc         WRITE(IMPRIM,10501)  (NUC,NP,POLY(NUC,NP),NUC=1,64)
ccc10501    FORMAT( 8('  6Q1C(',I2,',',I2,')=',G12.3) )
ccc         WRITE(IMPRIM,*)
ccc 500  CONTINUE
      END
