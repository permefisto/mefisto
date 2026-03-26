      SUBROUTINE PTDATE( POINT, XYZSOM, NOSOTE, NSIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE POINT EST IL DANS LE TETRAEDRE NOSOTE
C -----    Cf ptdste.f pour une autre version

C ENTREES:
C --------
C POINT  : LES 3 COORDONNEES DU POINT
C NOSOTE : LE NUMERO DES 4 SOMMETS DU TETRAEDRE

C SORTIES:
C --------
C NSIGNE : 1 SI LE POINT EST DANS LE TETRAEDRE
C          0 SI LE TETRAEDRE EST DEGENERE OU NE CONTIENT PAS LE POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1991
C....................................................................012
      INTEGER           NOSOTE(4)
      REAL              POINT(3), XYZSOM(3,*)
      DOUBLE PRECISION  COORTE(4,5)

C     CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT DANS NT
      DO 20 J=1,4
         DO 10 I=1,3
            COORTE(I,J) = XYZSOM(I,NOSOTE(J))
  10     CONTINUE
C        LA DERNIERE LIGNE DE LA MATRICE
         COORTE(4,J) = 1D0
C        LE SECOND MEMBRE
         IF( J .NE. 4 ) COORTE(J,5) = POINT(J)
  20  CONTINUE
      COORTE(4,5) = 1D0

C     INVERSION DU SYSTEME LINEAIRE
      CALL GAUSPT( 4 , 1 , COORTE , I )
      IF( I .NE. 0 ) THEN
         NSIGNE = 0
         RETURN
      ENDIF

C     TOUTES SES COORDONNEES BARYCENTRIQUES SONT COMPRISES ENTRE 0 ET 1
C     => LE POINT EST DANS LE TETRAEDRE
      DO 30 I=1,4
         IF( COORTE(I,5) .LT. -1D-15 .OR.
     %       COORTE(I,5) .GT.  1.D0 )
     %       GOTO 9000
 30   CONTINUE
      NSIGNE = 1
      RETURN

C     POINT NON DANS LE TETRAEDRE
 9000 NSIGNE = 0
      RETURN
      END
