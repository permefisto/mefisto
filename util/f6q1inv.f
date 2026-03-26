      SUBROUTINE F6Q1INV( POLY,   DPOLY, XWST6C, XYZUVW,
     %                    abcdef, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER les 6 COORDONNEES abcdef telles que F6(abcdef)=XYZUVW
C -----   ou F6 envoie le 6-cube de reference sur le 6 cube courant
C         de sommets XWST6C(64,6)
C ENTREES:
C --------
C POLY   : POLY(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY(I,J,K,L,M,N,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                             U**(L-1) V**(M-1) W**(N-1)
C DPOLY  : DPOLY(L,Q,K) COEFFICIENT L DU Q-EME POLYNOME DERIVE DANS
C                       la DIRECTION K
C          DPOLY(I,J,K,L,M,N,Q,K)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                                U**(L-1) V**(M-1) W**(N-1)
C          K=1,...,6
C XWST6C : LES 6 COORDONNEES DES 64 SOMMETS DU 6-CUBE COURANT
C XYZUVW : LES 6 COORDONNEES DU POINT DU 6-CUBE COURANT de abcdef A CALCULER
C
C SORTIES:
C --------
C abcdef : LES 6 COORDONNEES DU POINT DU 6-CUBE DE REFERENCE
C          tel que F6(abcdef) = XYZUVW
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L LIONS UPMC Paris Novembre2006
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  XWST6C(64,6), XYZUVW(6)
      DOUBLE PRECISION  abcdef(6),POLY(64,64),DPOLY(64,64,6),
     %                  FE(6), DFE(6,6),
     %                  DFEM1(6,6)
      DOUBLE PRECISION  V, XK, A, B, C, D, E, F, DIF
C
C     CALCULER LES 64 COEFFICIENTS DES 64 POLYNOMES DE BASE
C     DE L'EF 6CUBE 6Q1C
      CALL CF6Q1C( POLY )
C     LES 64 COEFFICIENTS DES 64 POLYNOMES DERIVES DANS LES 6 DIRECTIONS
      DO 3 I=1,64
         DO 2 K=1,6
            CALL PN6DDE( K, N1, POLY(1,I), DPOLY(1,I,K) )
 2       CONTINUE
 3    CONTINUE
C
C      LES 6 XYZUVW DES 64 NOEUDS=SOMMETS
      DO 31 J=1,64
         WRITE(IMPRIM,10031) (J,K,X(J,K),K=1,6)
 31   CONTINUE
10031 FORMAT(6('  XYZUVW(',I2,',',I1,')=',D15.8))
C
C     LES 6 COORDONNEES SUR LE 6-CUBE DE REFERENCE DU POINT DE DEPART
      DO 4 I=1,6
         abcdef(I) = 0.5d0
 4    CONTINUE
C
C     LES ITERATIONS
 10   A = abcdef(1)
      B = abcdef(2)
      C = abcdef(3)
      D = abcdef(4)
      E = abcdef(5)
      F = abcdef(6)
C
C
C        FE,DFEM1,DELTA  AU POINT abcdef
C        ===============================
         DO 30 K=1,6
C
C           LA COORDONNEE K DE FE AU POINT abcdef
            FE(K) = 0.D0
C           LA COLONNE K DE DFE AU POINT abcdef
            DO 5 J=1,6
               DFE(J,K) = 0.D0
 5          CONTINUE
c
            DO 20 I=1,64
C
C              VALEUR DU POLYNOME I AU POINT abcdef
               CALL PN6DVA( 2, POLY(1,I), A, B, C, D, E, F,  V )
C
C              COORDONNEE K du SOMMET I
               XK = XWST6C(I,K)
C
C              FE(abcdef)
               FE(K) = FE(K) + V * XK
C
C              VALEUR en abcdef des 6 DERIVEES PREMIERES DU POLYNOME I
               DO 15 J=1,6
                  CALL PN6DVA( 2, DPOLY(1,I,J), A,B,C,D,E,F,  V )
                  DFE(J,K) = DFE(J,K) + V * XK
 15            CONTINUE
C
 20         CONTINUE
 30      CONTINUE
C
CCC         WRITE(IMPRIM,10029) (I,L,POLY(I,L),I=1,64)
CCC10029    FORMAT(6('  POLY(',I2,',',I2,')=',D15.8))
C
CCC         WRITE(IMPRIM,10032) (L,K,F(L,K),K=1,6)
CCC10032    FORMAT(6('  F(',I2,',',I1,')=',D15.8))
C
CCC         DO 33 J=1,6
CCC            WRITE(IMPRIM,10033) (J,K,DF(J,K),K=1,6)
CCC 33      CONTINUE
CCC10033    FORMAT(6('  DF(',I1,',',I1,')=',D15.8))
CCC         DO 34 K=1,64
CCC            WRITE(IMPRIM,10034) (J,K,L,DPOLY(J,K,L),J=1,6)
CCC 34      CONTINUE
CCC10034    FORMAT(6('  DPOLY(',I1,',',I2,',',I2,')=',D15.8))
C
C     [DF]-1 AU POINT abcdef
      CALL M66INV( DFE,  DFEM1, V, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     abcdef m+1 = abcdef m - DFEM1 ( FE(abcdef m) - XYZUVW )
      DIF = 0d0
      DO 50 K=1,6
         XK = (FE(1) - XYZUVW(1)) * DFEM1(1,K)
     %      + (FE(2) - XYZUVW(2)) * DFEM1(2,K)
     %      + (FE(3) - XYZUVW(3)) * DFEM1(3,K)
     %      + (FE(4) - XYZUVW(4)) * DFEM1(4,K)
     %      + (FE(5) - XYZUVW(5)) * DFEM1(5,K)
     %      + (FE(6) - XYZUVW(6)) * DFEM1(6,K)
         abcdef(K) = abcdef(K) - XK
         DIF = DIF + XK * XK
 50   CONTINUE
      IF( DIF .GT. 1D-8 ) GOTO 10
C
      RETURN
      END
