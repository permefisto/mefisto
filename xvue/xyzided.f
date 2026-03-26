      SUBROUTINE XYZIDED( XYZ1, XYZ2, IDENTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NON 2 POINTS DEFINIS PAR 3 COORDONNEES
C -----    REELLES DOUBLE PRECISION
C
C ENTREES :
C ---------
C XYZ1   : LE PREMIER POINT
C XYZ2   : LE SECOND  POINT
C
C SORTIE :
C --------
C IDENTQ : 1 SI LES 2 POINTS SONT JUGES IDENTIQUES
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      DOUBLE PRECISION  XYZ1(3), XYZ2(3), D1, D2, DEPZERO, DEPSXYZ
      INTRINSIC         ABS, DBLE
C
C     SEUILS EN DOUBLE PRECISION
      DEPZERO = DBLE( EPZERO )
      DEPSXYZ = DBLE( EPSXYZ )
C
C     COMPARAISON DES 3 COORDONNEES
      DO K=1,3
         D1 = ABS( XYZ1(K) )
         D2 = ABS( XYZ2(K) )
         IF     ( D1 .LE. DEPZERO ) THEN
            IF  ( D2 .GT. DEPZERO ) GOTO 30
         ELSE IF( D2 .LE. DEPZERO ) THEN
            GOTO 30
         ELSE IF( ABS( XYZ1(K)-XYZ2(K) ) .GT. (D1+D2) * DEPSXYZ ) THEN
            GOTO 30
         ENDIF
      ENDDO
C
C     LES 2 POINTS XYZ1 ET XYZ2 SONT IDENTIQUES
      IDENTQ = 1
      GOTO 9999
C
C     LES POINTS SONT DIFFERENTS
 30   IDENTQ = 0
C
 9999 RETURN
      END
