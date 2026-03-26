      SUBROUTINE PN6DDE( NODERI, N1, P, DP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHE DES COEFFICIENTS DU POLYNOME DERIVE DP DE P
C -----
C
C ENTREES:
C --------
C NODERI : 1 DP EST LE POLYNOME DERIVE PAR RAPPORT A X DU POLYNOME P
C          2 DP EST LE POLYNOME DERIVE PAR RAPPORT A Y DU POLYNOME P
C          3 DP EST LE POLYNOME DERIVE PAR RAPPORT A Z DU POLYNOME P
C          4 DP EST LE POLYNOME DERIVE PAR RAPPORT A U DU POLYNOME P
C          5 DP EST LE POLYNOME DERIVE PAR RAPPORT A V DU POLYNOME P
C          6 DP EST LE POLYNOME DERIVE PAR RAPPORT A W DU POLYNOME P
C N1     : DEGRE+1 DU POLYNOME P EN CHACUNE DES VARIABLES
C P      : TABLEAU A 6 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J,K,L,M,N)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                        U**(L-1) V**(M-1) W**(N-1)
C
C SORTIES:
C --------
C DP     : LE POLYNOME DERIVE DE P PAR RAPPORT A X(NODERI)
C
C ATTENTION: DP DOIT ETRE DIFFERENT DE P
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  P(N1,N1,N1,N1,N1,N1),DP(N1,N1,N1,N1,N1,N1)
      COMMON / UNITES / LECTEU , IMPRIM , NUNITE(30)
C
      IF( N1 .LE. 1 ) GOTO 1000
C
C     CHOIX DE LA DERIVATION PAR RAPPORT A LA NODERI VARIABLE
      GOTO( 100, 200, 300, 400, 500, 600 ), NODERI
C
C     DERIVATION SELON X
 100  DO 160 N=1,N1
         DO 150 M=1,N1
            DO 140 L=1,N1
               DO 130 K=1,N1
                  DO 120 J=1,N1
                     DO 110 I=2,N1
                        I1 = I - 1
                        DP(I1,J,K,L,M,N) = I1 * P(I,J,K,L,M,N)
 110                 CONTINUE
                     DP(N1,J,K,L,M,N) = 0.D0
 120              CONTINUE
 130           CONTINUE
 140        CONTINUE
 150     CONTINUE
 160  CONTINUE
      GOTO 900
C
C     DERIVATION SELON Y
 200  DO 260 N=1,N1
         DO 250 M=1,N1
            DO 240 L=1,N1
               DO 230 K=1,N1
                  DO 220 I=1,N1
                     DO 210 J=2,N1
                        I1 = J - 1
                        DP(I,I1,K,L,M,N) = I1 * P(I,J,K,L,M,N)
 210                 CONTINUE
                     DP(I,N1,K,L,M,N) = 0.D0
 220              CONTINUE
 230           CONTINUE
 240        CONTINUE
 250     CONTINUE
 260  CONTINUE
      GOTO 900
C
C     DERIVATION SELON Z
 300  DO 360 N=1,N1
         DO 350 M=1,N1
            DO 340 L=1,N1
               DO 330 J=1,N1
                  DO 320 I=1,N1
                     DO 310 K=2,N1
                        I1 = K - 1
                        DP(I,J,I1,L,M,N) = I1 * P(I,J,K,L,M,N)
 310                 CONTINUE
                     DP(I,J,N1,L,M,N) = 0.D0
 320              CONTINUE
 330           CONTINUE
 340        CONTINUE
 350     CONTINUE
 360  CONTINUE
      GOTO 900
C
C     DERIVATION SELON U
 400  DO 460 N=1,N1
         DO 450 M=1,N1
            DO 440 K=1,N1
               DO 430 J=1,N1
                  DO 420 I=1,N1
                     DO 410 L=2,N1
                        I1 = L - 1
                        DP(I,J,K,I1,M,N) = I1 * P(I,J,K,L,M,N)
 410                 CONTINUE
                     DP(I,J,K,N1,M,N) = 0.D0
 420              CONTINUE
 430           CONTINUE
 440        CONTINUE
 450     CONTINUE
 460  CONTINUE
      GOTO 900
C
C     DERIVATION SELON V
 500  DO 560 N=1,N1
         DO 550 L=1,N1
            DO 540 K=1,N1
               DO 530 J=1,N1
                  DO 520 I=1,N1
                     DO 510 M=2,N1
                        I1 = M - 1
                        DP(I,J,K,L,I1,N) = I1 * P(I,J,K,L,M,N)
 510                 CONTINUE
                     DP(I,J,K,L,N1,N) = 0.D0
 520              CONTINUE
 530           CONTINUE
 540        CONTINUE
 550     CONTINUE
 560  CONTINUE
      GOTO 900
C
C     DERIVATION SELON W
 600  DO 660 M=1,N1
         DO 650 L=1,N1
            DO 640 K=1,N1
               DO 630 J=1,N1
                  DO 620 I=1,N1
                     DO 610 N=2,N1
                        I1 = N - 1
                        DP(I,J,K,L,M,I1) = I1 * P(I,J,K,L,M,N)
 610                 CONTINUE
                     DP(I,J,K,L,M,N1) = 0.D0
 620              CONTINUE
 630           CONTINUE
 640        CONTINUE
 650     CONTINUE
 660  CONTINUE
      GOTO 900
C
C     ERREURS
 1000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') N1
      KERR(1) = 'DEGRE INCORRECT=' // KERR(MXLGER)(1:4)
      CALL LEREUR
C
  900 RETURN
      END
