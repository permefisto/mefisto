      SUBROUTINE XYZIDE( XYZ1, XYZ2, IDENTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NON 2 POINTS DEFINIS PAR 3 COORDONNEES
C -----

C ENTREES :
C ---------
C XYZ1   : LE PREMIER POINT
C XYZ2   : LE SECOND  POINT

C SORTIE :
C --------
C IDENTQ : =1 SI LES 2 POINTS SONT JUGES IDENTIQUES
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO,  EPSXYZ
      REAL              XYZ1(3), XYZ2(3)

C     COMPARAISON DES 3 COORDONNEES
      DO K=1,3
         R1 = ABS( XYZ1(K) )
         R2 = ABS( XYZ2(K) )
         IF     ( R1 .LE. EPZERO ) THEN
            IF  ( R2 .GT. EPZERO ) GOTO 30
         ELSE IF( R2 .LE. EPZERO ) THEN
            GOTO 30
         ELSE IF( ABS( XYZ1(K)-XYZ2(K) ) .GT. (R1+R2) * EPSXYZ ) THEN
            GOTO 30
         ENDIF
      ENDDO

C     LES 2 POINTS XYZ1 ET XYZ2 SONT IDENTIQUES
      IDENTQ = 1
      GOTO 9999

C     LES POINTS SONT DIFFERENTS
 30   IDENTQ = 0

 9999 RETURN
      END
