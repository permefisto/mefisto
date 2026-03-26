      SUBROUTINE F3Q1INV( POLY, DPOLY, XYZS3C, XYZ,
     %                    abc,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER les 3 COORDONNEES abc telles que F3(abc)=XYZ
C -----    ou F3 envoie le 3-cube de reference sur le 3 cube courant
C          de sommets XYZS3C(3,8)
C
C ENTREES:
C --------
C POLY   : POLY(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY(I,J,K,L,M,N,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C DPOLY  : DPOLY(L,Q,K) COEFFICIENT L DU Q-EME POLYNOME DERIVE DANS
C                       la DIRECTION K
C          DPOLY(I,J,K,L,M,N,Q,D)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C          D=1,...,3
C XYZS3C : LES 3 COORDONNEES DES 8 SOMMETS DU 3-CUBE COURANT
C XYZ    : LES 3 COORDONNEES DU POINT DU 3-CUBE de abc A CALCULER
C
C SORTIES:
C --------
C abc    : LES 3 COORDONNEES DU POINT DU 3-CUBE DE REFERENCE
C          tel que F3(abc) = XYZ
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L LIONS UPMC Paris Novembre2006
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  XYZS3C(3,8), XYZ(3)
      DOUBLE PRECISION  abc(3), POLY(8,8), DPOLY(8,8,3),
     %                  FE(3), DFE(3,3), DFEM1(3,3)
      DOUBLE PRECISION  V, XK, A, B, C, DIF
C
cccC     LES 3 XYZ DES 8 NOEUDS=SOMMETS
ccc      DO 31 K=1,8
ccc         WRITE(IMPRIM,10031) (J,K,XYZS3C(J,K),J=1,3)
ccc 31   CONTINUE
ccc10031 FORMAT(3('  XYZS3C(',I2,',',I1,')=',D15.8))
ccc      print *,'FE-1(',XYZ,')'
C
      IERR = 0
C
C     LES 3 COORDONNEES SUR LE 3-CUBE DE REFERENCE DU POINT DE DEPART
      DO 4 I=1,3
         abc(I) = 0.5d0
 4    CONTINUE
C
C     LES ITERATIONS
 10   A = abc(1)
      B = abc(2)
      C = abc(3)
C
C     FE,DFEM1,DELTA  AU POINT abc
C     ============================
      DO 30 K=1,3
C
C        LA COORDONNEE K DE FE AU POINT abc
         FE(K) = 0.D0
C        LA COLONNE K DE DFE AU POINT abc
         DO 5 J=1,3
            DFE(J,K) = 0.D0
 5       CONTINUE
c
         DO 20 I=1,8
C
C           VALEUR DU POLYNOME I AU POINT abc
            CALL PN3DVA( 2, POLY(1,I), A, B, C,  V )
C
C           COORDONNEE K du SOMMET I
            XK = XYZS3C(K,I)
C
C           FE(abc)
            FE(K) = FE(K) + V * XK
C
C           VALEUR en abc des 3 DERIVEES PREMIERES DU POLYNOME I
            DO 15 J=1,3
               CALL PN3DVA( 2, DPOLY(1,I,J), A,B,C,  V )
               DFE(J,K) = DFE(J,K) + V * XK
 15         CONTINUE
C
 20      CONTINUE
C
 30   CONTINUE
C
C     [DF]-1 AU POINT abc
      CALL M33INV( DFE,  V, DFEM1 )
      IF( IERR .NE. 0 ) RETURN
C
C     Methode de NEWTON sur FE(abc) - XYZ = 0
C     abc m+1 = abc m - DFEM1 ( FE(abc m) - XYZ )
      DIF = 0d0
      DO 50 K=1,3
         XK = (FE(1) - XYZ(1)) * DFEM1(1,K)
     %      + (FE(2) - XYZ(2)) * DFEM1(2,K)
     %      + (FE(3) - XYZ(3)) * DFEM1(3,K)
         abc(K) = abc(K) - XK
         DIF = DIF + XK * XK
 50   CONTINUE
      IF( DIF .GT. 1D-8 ) GOTO 10
C
ccc      print *,'FE(',abc,')=',XYZ
C
      RETURN
      END
