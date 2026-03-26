      SUBROUTINE CF3Q1C( POLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 8 COEFFICIENTS DES 8 POLYNOMES DE BASE
C -----    DE L'EF 3CUBE 3Q1C
C          ATTENTION: PRESENTATION ICI AVEC UNE NUMEROTATION PAR AXES
C             Z
C             5-------7
C            /|      /|
C          6/-------8 |
C           | 1-----|-3Y
C           |/      |/
C          2/-------4
C         X
C SORTIE :
C --------
C POLY   : POLY(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY(I,J,K,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L LIONS UPMC Paris Novembre2006
C ......................................................................
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  POLY(8,8)
      DOUBLE PRECISION  CI(2), CJ(2), CK(2)
      INTEGER           IC(2), JC(2), KC(2)
C
C     MISE A ZERO DES 8*8 COEFFICIENTS
      CALL AZEROD( 8*8, POLY )
C
      NP = 0
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
C              TRAITEMENT DU POLYNOME(I,J,K)
C              NUMERO DU POLYNOME TRAITE
               NP = NP + 1
C
C              SES COEFFICIENTS NON NULS
               DO 30 KK=1,NBKC
                  DO 20 JJ=1,NBJC
                     DO 10 II=1,NBIC
                        NUC = IC(II) + 2 * (JC(JJ)-1) + 4 * (KC(KK)-1)
                        POLY(NUC,NP) = CI(II) * CJ(JJ) * CK(KK)
 10                  CONTINUE
 20               CONTINUE
 30            CONTINUE
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE
C
ccc      WRITE(IMPRIM,10500)
ccc10500 FORMAT(/,'LES 8 COEFFICIENTS DES 8 POLYNOMES 3Q1C'/130(1H@)/)
ccc      DO 500 NP=1,8
ccc         WRITE(IMPRIM,10501)  (NUC,NP,POLY(NUC,NP),NUC=1,8)
ccc10501    FORMAT( 8('  3Q1C(',I2,',',I2,')=',G12.3) )
ccc         WRITE(IMPRIM,*)
ccc 500  CONTINUE
      END
