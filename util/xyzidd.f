      SUBROUTINE XYZIDD( XYZ1 , XYZ2 , IDENTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NOM 2 POINTS DEFINIS PAR 3 COORDONNEES
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
      DOUBLE PRECISION  XYZ1(3), XYZ2(3), R1, R2
      INTRINSIC         ABS
C
C     COMPARAISON DES 3 COORDONNEES
      DO 10 K=1,3
         R1 = ABS( XYZ1(K) )
         R2 = ABS( XYZ2(K) )
         IF     ( R1 .LE. EPZERO ) THEN
            IF  ( R2 .GT. EPZERO ) GOTO 30
         ELSE IF( R2 .LE. EPZERO ) THEN
            GOTO 30
         ELSE IF( ABS( XYZ1(K)-XYZ2(K) ) .GT. R1 * EPSXYZ ) THEN
            GOTO 30
         ENDIF
 10   CONTINUE
C
C     LES 2 POINTS XYZ1 ET XYZ2 SONT IDENTIQUES
      IDENTQ = 1
      RETURN
C
C     LES POINTS SONT DIFFERENTS
  30  IDENTQ = 0
      RETURN
      END
