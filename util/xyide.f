      SUBROUTINE XYIDE( XY1, XY2, IDENTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NOM 2 POINTS DEFINIS PAR 3 COORDONNEES
C -----

C ENTREES :
C ---------
C XY1    : LE PREMIER POINT
C XY2    : LE SECOND  POINT

C SORTIE :
C --------
C IDENTQ : =1 SI LES 2 POINTS SONT JUGES IDENTIQUES
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL              XY1(2), XY2(2)

C     COMPARAISON DES 2 COORDONNEES
      DO 10 K=1,2
         R1 = ABS( XY1(K) )
         R2 = ABS( XY2(K) )
         IF     ( R1 .LE. EPZERO ) THEN
            IF  ( R2 .GT. EPZERO ) GOTO 30
         ELSE IF( R2 .LE. EPZERO ) THEN
            GOTO 30
         ELSE IF( ABS( XY1(K)-XY2(K) ) .GT. R1 * EPSXYZ ) THEN
            GOTO 30
         ENDIF
 10   CONTINUE

C     LES 2 POINTS XY1 ET XY2 SONT IDENTIQUES
      IDENTQ = 1
      RETURN

C     LES POINTS SONT DIFFERENTS
  30  IDENTQ = 0

      RETURN
      END
