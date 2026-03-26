         SUBROUTINE  PEN9AR( XYZ, M, NZ, XYSCA1, COSO,  XYZA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES COORDONNEES DES POINTS SUR LES ARETES
C -----    DU PENTAEDRE RECTANGLE UNITE DE LA FORMULE DE
C          L'INTERPOLATION TRANSFINIE DEGRE 1 (cf CRAS A.PERRONNET)
C
C ENTREES:
C --------
C XYZ    : LES 3 COORDONNEES DU POINT DANS LE PENTAEDRE UNITE
C M      : NOMBRE DE SOMMETS PAR ARETE DES FACES TRIANGULAIRES
C NZ     : NOMBRE DE SOMMETS PAR ARETE EN Z DES FACES QUADRILATERES
C XYSCA1 : LES 2 COORDONNEES DES POINTS SUR LE CARRE UNITE
C          DES FACES 2 3 4 DU PENTAEDRE
C COSO   : LES COORDONNEES DES SOMMETS SUR LES FACES DU PENTAEDRE COURBE
C
C SORTIES:
C --------
C XYZA   : LES 3 COORDONNEES DES SOMMETS DES ARETES DU PENTAEDRE COURBE
C          NECESSAIRE DANS LA FORMULE DE L'INTERPOLATION TRANSFINIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL     XYZ(3)
      REAL     XYSCA1(2,M,NZ,2:4)
      REAL     COSO(3,*)
      REAL     XYZA(3,15)
C
C     LA FONCTION FORMULE DES NUMEROS DES SOMMETS DU PENTAEDRE
      NUSOPE(I,J,K) = J + ( I * I - I + (K-1) * ( M * M + M ) ) / 2
C
C     LES COORDONNEES BARYCENTRIQUES DU TRIANGLE DE SECTION
      CB2 = XYZ(1)
      CB3 = XYZ(2)
      CB1 = 1 - CB2 - CB3
      Z   = XYZ(3)
C
C     SUR L'ARETE 1: S1 S2 DU PENTAEDRE = ARETE 1 DE LA FACE CARREE 2
C     ===============================================================
C
C     S(1-CB2,CB2,0,0) = A1F2(CB2)
C     ----------------------------
      R2 = 0
      R  = CB2
      DO 10 I=2,M
         R2 = XYSCA1(1,I,1,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 11
 10   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU PENTAEDRE DE SOMMETS L1, L2
 11   L1 = NUSOPE(I-1,1,1)
      L2 = NUSOPE(I  ,1,1)
      R1 = XYSCA1(1,I-1,1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,1) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,1) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,1) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(CB1,1-CB1,0,0) = A1F2(1-CB1)
C     ------------------------------
      R = 1-CB1
      DO 15 I=2,M
         R2 = XYSCA1(1,I,1,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 16
 15   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU PENTAEDRE DE SOMMETS L1, L2
 16   L1 = NUSOPE(I-1,1,1)
      L2 = NUSOPE(I  ,1,1)
      R1 = XYSCA1(1,I-1,1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,3) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,3) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,3) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 4: S4 S5 DU PENTAEDRE = ARETE 3 DE LA FACE CARREE 2
C     ===============================================================
C
C     S(1-CB2,CB2,0,1) = A3F2(CB2)
C     ----------------------------
      R = CB2
      DO 40 I=2,M
         R2 = XYSCA1(1,I,NZ,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 41
 40   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU PENTAEDRE DE SOMMETS L1, L2
 41   L1 = NUSOPE(I-1,1,NZ)
      L2 = NUSOPE(I  ,1,NZ)
      R1 = XYSCA1(1,I-1,NZ,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,2) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,2) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,2) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(CB1,1-CB1,0,1) = A3F2(1-CB1)
C     ------------------------------
      R = 1-CB1
      DO 45 I=2,M
         R2 = XYSCA1(1,I,NZ,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 46
 45   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU PENTAEDRE DE SOMMETS L1, L2
 46   L1 = NUSOPE(I-1,1,NZ)
      L2 = NUSOPE(I  ,1,NZ)
      R1 = XYSCA1(1,I-1,NZ,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,4) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,4) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,4) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 2: S2 S3 DU PENTAEDRE = ARETE 1 DE LA FACE CARREE 3
C     ===============================================================
C
C     S(0,1-CB3,CB3,0) = A1F3(CB3)
C     ----------------------------
      R = CB3
      DO 20 I=2,M
         R2 = XYSCA1(1,I,1,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 21
 20   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 2 DU PENTAEDRE DE SOMMETS L1, L2
 21   L1 = NUSOPE(M,I-1,1)
      L2 = NUSOPE(M,I  ,1)
      R1 = XYSCA1(1,I-1,1,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,5) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,5) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,5) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(0,CB2,1-CB2,0) = A1F3(1-CB2)
C     ------------------------------
      R = 1-CB2
      DO 25 I=2,M
         R2 = XYSCA1(1,I,1,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 26
 25   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 2 DU PENTAEDRE DE SOMMETS L1, L2
 26   L1 = NUSOPE(M,I-1,1)
      L2 = NUSOPE(M,I  ,1)
      R1 = XYSCA1(1,I-1,1,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,7) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,7) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,7) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 5: S5 S6 DU PENTAEDRE = ARETE 3 DE LA FACE CARREE 3
C     ===============================================================
C
C     S(0,1-CB3,CB3,1) = A3F3(CB3)
C     ----------------------------
      R = CB3
      DO 50 I=2,M
         R2 = XYSCA1(1,I,NZ,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 51
 50   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU PENTAEDRE DE SOMMETS L1, L2
 51   L1 = NUSOPE(M,I-1,NZ)
      L2 = NUSOPE(M,I  ,NZ)
      R1 = XYSCA1(1,I-1,NZ,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,6) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,6) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,6) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(0,CB2,1-CB2,1) = A3F3(1-CB2)
C     ------------------------------
      R = 1-CB2
      DO 55 I=2,M
         R2 = XYSCA1(1,I,NZ,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 56
 55   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU PENTAEDRE DE SOMMETS L1, L2
 56   L1 = NUSOPE(M,I-1,NZ)
      L2 = NUSOPE(M,I  ,NZ)
      R1 = XYSCA1(1,I-1,NZ,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,8) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,8) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,8) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 3: S3 S1 DU PENTAEDRE = ARETE 1 DE LA FACE CARREE 4
C     ===============================================================
C
C     S(CB1,0,1-CB1,0) = A1F4(1-CB1) MAIS F4 AVEC Y EN SENS INVERSE => CB1
C     ------------------------------
      R = CB1
      DO 30 I=2,M
         R2 = XYSCA1(1,I,1,4)
         IF( R .LE. R2 + EPSXYZ ) GOTO 31
 30   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 3 DU PENTAEDRE DE SOMMETS L1, L2
C     MAIS I DANS LA FACE EST M+1-I DANS LE PENTAEDRE
 31   L1 = NUSOPE(M+2-I,M+2-I,1)
      L2 = NUSOPE(M+1-I,M+1-I,1)
      R1 = XYSCA1(1,I-1,1,4)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,9) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,9) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,9) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(1-CB3,0,CB3,0) = A1F4(CB3) MAIS F4 AVEC Y EN SENS INVERSE => 1-CB3
C     ----------------------------
      R = 1-CB3
      DO 35 I=2,M
         R2 = XYSCA1(1,I,1,4)
         IF( R .LE. R2 + EPSXYZ ) GOTO 36
 35   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 3 DU PENTAEDRE DE SOMMETS L1, L2
C     MAIS I DANS LA FACE EST M+1-I DANS LE PENTAEDRE
 36   L1 = NUSOPE(M+2-I,M+2-I,1)
      L2 = NUSOPE(M+1-I,M+1-I,1)
      R1 = XYSCA1(1,I-1,1,4)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,11) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,11) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,11) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 6: S6 S4 DU PENTAEDRE = ARETE 3 DE LA FACE CARREE 4
C     ===============================================================
C
C     S(CB1,0,1-CB1,1) = A3F4(1-CB1) MAIS F4 AVEC Y EN SENS INVERSE => CB1
C     ------------------------------
      R = CB1
      DO 60 I=2,M
         R2 = XYSCA1(1,I,NZ,4)
         IF( R .LE. R2 + EPSXYZ ) GOTO 61
 60   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU PENTAEDRE DE SOMMETS L1, L2
C     MAIS I DANS LA FACE EST M+1-I DANS LE PENTAEDRE
 61   L1 = NUSOPE(M+2-I,M+2-I,NZ)
      L2 = NUSOPE(M+1-I,M+1-I,NZ)
      R1 = XYSCA1(1,I-1,NZ,4)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,10) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,10) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,10) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(1-CB3,0,CB3,1) = A3F4(CB3) MAIS F4 AVEC Y EN SENS INVERSE => 1-CB3
C     ----------------------------
      R = 1-CB3
      DO 65 I=2,M
         R2 = XYSCA1(1,I,NZ,4)
         IF( R .LE. R2 + EPSXYZ ) GOTO 66
 65   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU PENTAEDRE DE SOMMETS L1, L2
C     MAIS I DANS LA FACE EST M+1-I DANS LE PENTAEDRE
 66   L1 = NUSOPE(M+2-I,M+2-I,NZ)
      L2 = NUSOPE(M+1-I,M+1-I,NZ)
      R1 = XYSCA1(1,I-1,NZ,4)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,12) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,12) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,12) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 7: S1 S4 DU PENTAEDRE = ARETE 4 DE LA FACE CARREE 2
C     ===============================================================
C
C     S(1,0,0,Z) = A4F2(Z)
C     --------------------
      R = Z
      DO 70 I=2,NZ
         R2 = XYSCA1(2,1,I,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 71
 70   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 7 DU PENTAEDRE DE SOMMETS L1, L2
 71   L1 = NUSOPE(1,1,I-1)
      L2 = NUSOPE(1,1,I  )
      R1 = XYSCA1(2,1,I-1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,13) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,13) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,13) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 8: S2 S5 DU PENTAEDRE = ARETE 4 DE LA FACE CARREE 3
C     ===============================================================
C
C     S(0,1,0,Z) = A4F3(Z)
C     --------------------
      DO 80 I=2,NZ
         R2 = XYSCA1(2,1,I,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 81
 80   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 7 DU PENTAEDRE DE SOMMETS L1, L2
 81   L1 = NUSOPE(M,1,I-1)
      L2 = NUSOPE(M,1,I  )
      R1 = XYSCA1(2,1,I-1,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,14) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,14) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,14) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 9: S3 S6 DU PENTAEDRE = ARETE 4 DE LA FACE CARREE 4
C     ===============================================================
C
C     S(0,0,1,Z) = A4F4(Z)
C     --------------------
      DO 90 I=2,NZ
         R2 = XYSCA1(2,1,I,4)
         IF( R .LE. R2 + EPSXYZ ) GOTO 91
 90   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 7 DU PENTAEDRE DE SOMMETS L1, L2
 91   L1 = NUSOPE(M,M,I-1)
      L2 = NUSOPE(M,M,I  )
      R1 = XYSCA1(2,1,I-1,4)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      XYZA(1,15) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,15) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,15) = A1 * COSO(3,L1) + A * COSO(3,L2)
      END
