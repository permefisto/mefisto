      SUBROUTINE XYZIDEDS( XYZ, NBP, POINTS, NUMPT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NON 1 POINT XYZ PARMI LES NBP POINTS
C -----
C
C ENTREES:
C --------
C XYZ    : LES 3 COORDONNEES DU POINT A COMPARER
C NBP    : LE NOMBRE DE POINTS DU TABLEAU POINTS
C POINTS : LES 3 COORDONNEES DES POINTS
C
C SORTIE :
C --------
C NUMPT  : >0 NUMERO DANS POINTS DU POINT IDENTIFIE A XYZ
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Octobre 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION XYZ(3), POINTS(3,NBP)
C
      DO NUMPT = 1, NBP
         CALL XYZIDED( XYZ, POINTS(1,NUMPT), IDENTQ )
C        IDENTQ=1 SI LES 2 POINTS SONT JUGES IDENTIQUES, 0 SINON
         IF( IDENTQ .EQ. 1 ) GOTO 9999
      ENDDO
      NUMPT = 0
C
 9999 RETURN
      END
