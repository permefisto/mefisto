      SUBROUTINE MIMXPP( NBPOIN , NPPOIN , COINLI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA BORNE INFERIEURE ET SUPERIEURE DES 3 COORDONNEES
C ----- DES POINTS POOL DE NUMERO DANS LE TABLEAU NPPOIN
C
C ENTREES:
C --------
C NBPOIN : NOMBRE DE POINTS POOL
C NPPOIN : NUMERO DES NBPOIN POINTS POOL
C
C SORTIE :
C --------
C COINLI : COINLI(.,1)=MIN  COINLI(.,2)=MAX
C          LES MINIMA ET MAXIMA DES COORDONNEES DES POINTS
C...............................................................................
C PROGRAMMEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS MARS 1987
C...............................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NPPOIN(NBPOIN)
      REAL              COINLI(6,2),XYZ(3)
C
C     INITIALISATION DE COINLI AU MIN ET MAX INVERSES
      CALL MIMXCO( COINLI )
C
C     LA BOUCLE SUR LES POINTS
      DO 20 I=1,NBPOIN
C        LE NUMERO POOL DU POINT
         NPP = NPPOIN( I )
         IF( NPP .GT. 0 ) THEN
C           LES COORDONNEES DE CE POINT
            CALL COORLE( XYZ , J , NUPOIN , NPP )
            DO 10 J=1,3
               C = XYZ(J)
               IF( C .LT. COINLI(J,1) ) THEN
                   COINLI(J,1) = C
               ENDIF
               IF( C .GT. COINLI(J,2) ) THEN
                   COINLI(J,2) = C
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
      END
