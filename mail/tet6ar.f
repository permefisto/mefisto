      SUBROUTINE TET6AR( XYZ, M, NBSTTR, XYSFTR, COSO,
     %                   XYZA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES COORDONNEES DES POINTS SUR LES ARETES
C -----    DU TETRAEDRE RECTANGLE UNITE DE LA FORMULE DE
C          L'INTERPOLATION TRANSFINIE DEGRE 1
C          (cf CRAS A.PERRONNET)
C
C ENTREES:
C --------
C XYZ     : LES 3 COORDONNEES DU POINT DANS LE TETRAEDRE UNITE
C M       : LE NOMBRE DE SOMMETS PAR ARETE DES FACES TRIANGULAIRES
C NBSTTR  : LE NOMBRE DE SOMMETS PAR TRIANGLE NBSTTR=M*(M+1)/2
C XYSFTR  : LES 2 COORDONNEES DES POINTS SUR LE CARRE UNITE
C           DES FACES 2 3 4 DU TETRAEDRE
C COSO    : LES COORDONNEES DES SOMMETS SUR LES FACES DU TETRAEDRE COURBE
C
C SORTIES:
C --------
C XYZA   : LES 3 COORDONNEES DES SOMMETS DES ARETES DU TETRAEDRE COURBE
C          NECESSAIRE DANS LA FORMULE DE L'INTERPOLATION TRANSFINIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL     XYZ(3)
      REAL     XYSFTR(2,NBSTTR,1:4)
      REAL     COSO(3,*)
      REAL     XYZA(3,12)
C
C     LA FONCTION FORMULE DU NUMERO DES SOMMETS DU TRIANGLE
      NUSOTR(I,J)   = (I*I-I)/2 + J
C     LA FONCTION FORMULE DU NUMERO DES SOMMETS DU TETRAEDRE
      NUSOTE(I,J,K) = (I-1)*I*(I+1)/6 + (J+K-2)*(J+K-1)/2 + K
C
C     LES COORDONNEES BARYCENTRIQUES DU TRIANGLE DE SECTION
      CB2 = XYZ(1)
      CB3 = XYZ(2)
      CB4 = XYZ(3)
      CB1 = 1 - CB2 - CB3 - CB4
C
C     SUR L'ARETE 1: S1 S2 DU TETRAEDRE = ARETE 1 DE LA FACE 1
C     ========================================================
C
C     S(1-CB2,CB2,0,0) = A1F1(CB2)
C     ----------------------------
      R2 = 0
      R  = CB2
      DO 10 I=2,M
         L2 = NUSOTR(I,1)
         R2 = XYSFTR(1,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 11
 10   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU TETRAEDRE DE SOMMETS L1, L2
 11   L1 = NUSOTR(I-1,1)
      R1 = XYSFTR(1,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,1,1)
      L2 = NUSOTE(I  ,1,1)
      XYZA(1,1) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,1) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,1) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(CB1,1-CB1,0,0) = A1F1(1-CB1)
C     ------------------------------
      R = 1-CB1
      DO 15 I=2,M
         L2 = NUSOTR(I,1)
         R2 = XYSFTR(1,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 16
 15   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 1 DU TETRAEDRE DE SOMMETS L1, L2
 16   L1 = NUSOTR(I-1,1)
      R1 = XYSFTR(1,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,1,1)
      L2 = NUSOTE(I  ,1,1)
      XYZA(1,2) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,2) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,2) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 2: S2 S3 DU TETRAEDRE = ARETE 1 DE LA FACE 1
C     ========================================================
C
C     S(0,1-CB3,CB3,0) = A2F1(CB3)
C     ----------------------------
      R = CB3
      DO 20 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 21
 20   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 2 DU TETRAEDRE DE SOMMETS L1, L2
 21   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,I-1,1)
      L2 = NUSOTE(M,I  ,1)
      XYZA(1,3) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,3) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,3) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(0,CB2,1-CB2,0) = A2F1(1-CB2)
C     ------------------------------
      R = 1-CB2
      DO 25 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 26
 25   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 2 DU TETRAEDRE DE SOMMETS L1, L2
 26   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,I-1,1)
      L2 = NUSOTE(M,I  ,1)
      XYZA(1,4) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,4) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,4) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 3: S3 S1 DU TETRAEDRE = ARETE 3 DE LA FACE 1
C     ========================================================
C
C     S(CB1,0,1-CB1,0) = A3F1(1-CB1)
C     ------------------------------
      R = 1-CB1
      DO 30 I=2,M
         L2 = NUSOTR(I,I)
         R2 = XYSFTR(2,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 31
 30   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 3 DU TETRAEDRE DE SOMMETS L1, L2
 31   L1 = NUSOTR(I-1,I-1)
      R1 = XYSFTR(2,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,I-1,1)
      L2 = NUSOTE(I  ,I  ,1)
      XYZA(1,5) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,5) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,5) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(1-CB3,0,CB3,0) = A3F1(CB3)
C     ----------------------------
      R = CB3
      DO 35 I=2,M
         L2 = NUSOTR(I,I)
         R2 = XYSFTR(2,L2,1)
         IF( R .LE. R2 + EPSXYZ ) GOTO 36
 35   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 3 DU TETRAEDRE DE SOMMETS L1, L2
 36   L1 = NUSOTR(I-1,I-1)
      R1 = XYSFTR(2,L1,1)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,I-1,1)
      L2 = NUSOTE(I  ,I  ,1)
      XYZA(1,6) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,6) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,6) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 4: S1 S4 DU TETRAEDRE = ARETE 3 DE LA FACE 2
C     ========================================================
C
C     S(1-CB4,0,0,CB4) = A3F2(CB4)
C     ----------------------------
      R = CB4
      DO 40 I=2,M
         L2 = NUSOTR(I,I)
         R2 = XYSFTR(2,L2,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 41
 40   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 4 DU TETRAEDRE DE SOMMETS L1, L2
 41   L1 = NUSOTR(I-1,I-1)
      R1 = XYSFTR(2,L1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,1,I-1)
      L2 = NUSOTE(I  ,1,I  )
      XYZA(1,7) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,7) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,7) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(CB1,0,0,1-CB1) = A3F2(1-CB1)
C     ------------------------------
      R = 1-CB1
      DO 45 I=2,M
         L2 = NUSOTR(I,I)
         R2 = XYSFTR(2,L2,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 46
 45   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 4 DU TETRAEDRE DE SOMMETS L1, L2
 46   L1 = NUSOTR(I-1,I-1)
      R1 = XYSFTR(2,L1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(I-1,1,I-1)
      L2 = NUSOTE(I  ,1,I  )
      XYZA(1,8) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,8) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,8) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 5: S2 S4 DU TETRAEDRE = ARETE 2 DE LA FACE 3
C     ========================================================
C
C     S(0,1-CB4,0,CB4) = A2F3(CB4) => 1-CB4 CAR ARETE S2 VERS S4
C     ----------------------------
      R = 1-CB4
      DO 50 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 51
 50   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU TETRAEDRE DE SOMMETS L1, L2
 51   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,1,M+2-I)
      L2 = NUSOTE(M,1,M+1-I)
      XYZA(1,9) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,9) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,9) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(0,CB2,0,1-CB2) = A2F3(1-CB2) => CB2 CAR ARETE S2 VERS S4
C     ------------------------------
      R = CB2
      DO 55 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,3)
         IF( R .LE. R2 + EPSXYZ ) GOTO 56
 55   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 5 DU TETRAEDRE DE SOMMETS L1, L2
 56   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,3)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,1,M+2-I)
      L2 = NUSOTE(M,1,M+1-I)
      XYZA(1,10) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,10) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,10) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     SUR L'ARETE 6: S3 S4 DU TETRAEDRE = ARETE 2 DE LA FACE 2
C     ========================================================
C
C     S(0,0,1-CB4,CB4) = A2F2(CB4)
C     ----------------------------
      R = CB4
      DO 60 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 61
 60   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 6 DU TETRAEDRE DE SOMMETS L1, L2
 61   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,M+2-I,I-1)
      L2 = NUSOTE(M,M+1-I,I  )
      XYZA(1,11) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,11) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,11) = A1 * COSO(3,L1) + A * COSO(3,L2)
C
C     S(0,0,CB3,1-CB3) = A2F2(1-CB3)
C     ------------------------------
      R = 1-CB3
      DO 65 I=2,M
         L2 = NUSOTR(M,I)
         R2 = XYSFTR(2,L2,2)
         IF( R .LE. R2 + EPSXYZ ) GOTO 66
 65   CONTINUE
C     R EST DANS LE SEGMENT [I-1,I] DE L'ARETE 6 DU TETRAEDRE DE SOMMETS L1, L2
 66   L1 = NUSOTR(M,I-1)
      R1 = XYSFTR(2,L1,2)
      A  = ( R - R1 ) / ( R2 - R1 )
      A1 = 1 - A
      L1 = NUSOTE(M,M+2-I,I-1)
      L2 = NUSOTE(M,M+1-I,I  )
      XYZA(1,12) = A1 * COSO(1,L1) + A * COSO(1,L2)
      XYZA(2,12) = A1 * COSO(2,L1) + A * COSO(2,L2)
      XYZA(3,12) = A1 * COSO(3,L1) + A * COSO(3,L2)
      END
