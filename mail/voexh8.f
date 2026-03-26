         SUBROUTINE VOEXH8( XYZ, CUXYZ, HEXYZ, NBX, NBY, NBZ, A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :         CALCULER LES COORDONNEES DE L'IMAGE D'UN SOMMET
C -----         INTERNE DU CUBE UNITE SUR CHAQUE ARETE DE L'HEXAEDRE
C
C ENTREES :
C ---------
C XYZ         : COORDONNEES DU POINT A TRAITER DANS LE CUBE
C CUXYZ       : COORDONNEES DES SOMMETS DES FACES DU CUBE
C HEXYZ       : COORDONNEES DES SOMMETS DES FACES DE L'HEXAEDRE
C NBX,NBY,NBZ : NOMBRE DE POINTS DANS CHAQUE DIRECTION
C
C SORTIES :
C ---------
C A           : COORDONNEES DES SOMMETS SUR LES 12 ARETES DE L'HEXAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL     XYZ(3)
      REAL     HEXYZ(3,6,*)
      REAL     CUXYZ(3,6,*)
      REAL     A(3,12)
C
C     ------------------------------------------------------------------
C     LES 4 ARETES DANS LA DIRECTION X SUR LE CUBE UNITE
C     ------------------------------------------------------------------
C     LE POINT SUR L'ARETE 1 (ARETE 1 DE LA FACE 1 DU CUBE) DE L'HEXAEDRE
      NF = 1
      R  = XYZ(1)
      DO 15 K=2,NBX
         IF( R .LE. CUXYZ(1,NF,K)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [Xk-1,Xk]
            ALFA  = (R-CUXYZ(1,NF,K-1))/(CUXYZ(1,NF,K)-CUXYZ(1,NF,K-1))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(Xk-1) et HEXYZ(k)
            DO 12 J=1,3
               A(J,1) = ALFA1 * HEXYZ(J,NF,K-1) + ALFA * HEXYZ(J,NF,K)
 12         CONTINUE
            GOTO 30
         ENDIF
 15   CONTINUE
C
C     LE POINT SUR L'ARETE 3 (ARETE 3 DE LA FACE 1 DU CUBE) DE L'HEXAEDRE
C     DECALAGE POUR ARRIVER AU PREMIER POINT DE L'ARETE 2 DANS LA FACE 1
 30   N0 = NBX * NBY - NBX
      DO 35 K=N0+2,N0+NBX
         IF( R .LE. CUXYZ(1,NF,K)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [Xk-1,Xk]
            ALFA  = (R-CUXYZ(1,NF,K-1))/(CUXYZ(1,NF,K)-CUXYZ(1,NF,K-1))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(Xk-1) et HEXYZ(k)
            DO 32 J=1,3
               A(J,3) = ALFA1 * HEXYZ(J,NF,K-1) + ALFA * HEXYZ(J,NF,K)
 32         CONTINUE
            GOTO 50
         ENDIF
 35   CONTINUE
C
C     LE POINT SUR L'ARETE 5 (ARETE 1 DE LA FACE 4 DU CUBE) DE L'HEXAEDRE
 50   NF = 4
      DO 55 K=2,NBX
         IF( R .LE. CUXYZ(1,NF,K)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [Xk-1,Xk]
            ALFA  = (R-CUXYZ(1,NF,K-1))/(CUXYZ(1,NF,K)-CUXYZ(1,NF,K-1))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(Xk-1) et HEXYZ(k)
            DO 52 J=1,3
               A(J,5) = ALFA1 * HEXYZ(J,NF,K-1) + ALFA * HEXYZ(J,NF,K)
 52         CONTINUE
            GOTO 70
         ENDIF
 55   CONTINUE
C
C     LE POINT SUR L'ARETE 7 (ARETE 3 DE LA FACE 4 DU CUBE) DE L'HEXAEDRE
C     DECALAGE POUR ARRIVER AU PREMIER POINT DE L'ARETE 3 DANS LA FACE 4
 70   DO 75 K=N0+2,N0+NBX
         IF( R .LE. CUXYZ(1,NF,K)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [Xk-1,Xk]
            ALFA  = (R-CUXYZ(1,NF,K-1))/(CUXYZ(1,NF,K)-CUXYZ(1,NF,K-1))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(Xk-1) et HEXYZ(k)
            DO 72 J=1,3
               A(J,7) = ALFA1 * HEXYZ(J,NF,K-1) + ALFA * HEXYZ(J,NF,K)
 72         CONTINUE
            GOTO 20
         ENDIF
 75   CONTINUE
C
C     ------------------------------------------------------------------
C     LES 4 ARETES DANS LA DIRECTION Y SUR LE CUBE UNITE
C     ------------------------------------------------------------------
C     LE POINT SUR L'ARETE 2 (ARETE 2 DE LA FACE 1 DU CUBE) DE L'HEXAEDRE
 20   NF = 1
      N1 = NBX
      R  = XYZ(2)
      DO 25 K=2,NBY
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(2,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [YN0,YN1]
            ALFA  = (R-CUXYZ(2,NF,N0))/(CUXYZ(2,NF,N1)-CUXYZ(2,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(YN0) et HEXYZ(YN1)
            DO 22 J=1,3
               A(J,2) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 22         CONTINUE
            GOTO 40
         ENDIF
 25   CONTINUE
C
C     LE POINT SUR L'ARETE 4 (ARETE 2 DE LA FACE 1 DU CUBE) DE L'HEXAEDRE
 40   N1 = 1
      DO 45 K=2,NBY
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(2,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [YN0,YN1]
            ALFA  = (R-CUXYZ(2,NF,N0))/(CUXYZ(2,NF,N1)-CUXYZ(2,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(YN0) et HEXYZ(YN1)
            DO 42 J=1,3
               A(J,4) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 42         CONTINUE
            GOTO 60
         ENDIF
 45   CONTINUE
C
C     LE POINT SUR L'ARETE 6 (ARETE 2 DE LA FACE 4 DU CUBE) DE L'HEXAEDRE
 60   NF = 4
      N1 = NBX
      DO 65 K=2,NBY
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(2,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [YN0,YN1]
            ALFA  = (R-CUXYZ(2,NF,N0))/(CUXYZ(2,NF,N1)-CUXYZ(2,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(YN0) et HEXYZ(YN1)
            DO 62 J=1,3
               A(J,6) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 62         CONTINUE
            GOTO 80
         ENDIF
 65   CONTINUE
C
C     LE POINT SUR L'ARETE 8 (ARETE 4 DE LA FACE 4 DU CUBE) DE L'HEXAEDRE
 80   N1 = 1
      DO 85 K=2,NBY
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(2,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [YN0,YN1]
            ALFA  = (R-CUXYZ(2,NF,N0))/(CUXYZ(2,NF,N1)-CUXYZ(2,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(YN0) et HEXYZ(YN1)
            DO 82 J=1,3
               A(J,8) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 82         CONTINUE
            GOTO 90
         ENDIF
 85   CONTINUE
C
C     ------------------------------------------------------------------
C     LES 4 ARETES DANS LA DIRECTION Z SUR LE CUBE UNITE
C     ------------------------------------------------------------------
C     LE POINT SUR L'ARETE 9 (ARETE 4 DE LA FACE 2 DU CUBE) DE L'HEXAEDRE
 90   NF = 2
      N1 = 1
      R  = XYZ(3)
      DO 95 K=2,NBZ
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(3,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [ZN0,ZN1]
            ALFA  = (R-CUXYZ(3,NF,N0))/(CUXYZ(3,NF,N1)-CUXYZ(3,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(ZN0) et HEXYZ(ZN1)
            DO 92 J=1,3
               A(J,9) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 92         CONTINUE
            GOTO 100
         ENDIF
 95   CONTINUE
C
C     LE POINT SUR L'ARETE 10 (ARETE 3 DE LA FACE 2 DU CUBE) DE L'HEXAEDRE
 100  N1 = NBX
      DO 105 K=2,NBZ
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(3,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [ZN0,ZN1]
            ALFA  = (R-CUXYZ(3,NF,N0))/(CUXYZ(3,NF,N1)-CUXYZ(3,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(ZN0) et HEXYZ(ZN1)
            DO 102 J=1,3
               A(J,10) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 102        CONTINUE
            GOTO 110
         ENDIF
 105  CONTINUE
C
C     LE POINT SUR L'ARETE 11 (ARETE 1 DE LA FACE 5 DU CUBE) DE L'HEXAEDRE
 110  NF = 5
      N1 = NBX
      DO 115 K=2,NBZ
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(3,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [ZN0,ZN1]
            ALFA  = (R-CUXYZ(3,NF,N0))/(CUXYZ(3,NF,N1)-CUXYZ(3,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(ZN0) et HEXYZ(ZN1)
            DO 112 J=1,3
               A(J,11) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 112        CONTINUE
            GOTO 120
         ENDIF
 115  CONTINUE
C
C     LE POINT SUR L'ARETE 12 (ARETE 3 DE LA FACE 5 DU CUBE) DE L'HEXAEDRE
 120  N1 = 1
      DO 125 K=2,NBZ
         N0 = N1
         N1 = N1 + NBX
         IF( R .LE. CUXYZ(3,NF,N1)+EPSXYZ ) THEN
C           LE POINT R EST DANS L'INTERVALLE [ZN0,ZN1]
            ALFA  = (R-CUXYZ(3,NF,N0))/(CUXYZ(3,NF,N1)-CUXYZ(3,NF,N0))
            ALFA1 = 1 - ALFA
C           INTERPOLATION LINEAIRE SUR LA CORDE ENTRE HEXYZ(ZN0) et HEXYZ(ZN1)
            DO 122 J=1,3
               A(J,12) = ALFA1 * HEXYZ(J,NF,N0) + ALFA * HEXYZ(J,NF,N1)
 122        CONTINUE
            GOTO 300
         ENDIF
 125  CONTINUE
C
 300  RETURN
      END
