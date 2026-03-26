      SUBROUTINE NUAR2STR( NS1, NS2, NOSOTR, NA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: RETROUVER LE NUMERO D'ARETE DE SOMMETS NS1 NS2 DU TRIANGLE NOSOTR
C ----

C ENTREES:
C --------
C NS1 NS2: NUMERO DES 2 SOMMETS DE L'ARETE
C NOSOTR : NUMERO DES 3 SOMMETS DU TRIANGLE

C SORTIE :
C --------
C NA     : NUMERO DE 1 A 3 DE L'ARETE NS1-NS2 DU TRIANGLE NOSOTR
C          =0 SI PAS D'ARETE NS1-NS2 DANS NOSOTR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Veulettes sur mer               Fevrier 2020
C23456...............................................................012
      INTEGER  NOSOTR(3)

      DO NA = 1, 3

C        LE NO DES 2 SOMMETS DE L'ARETE NA DE NOSOTR
         NSA1 = NOSOTR( NA )
         IF( NA .EQ. 3 ) THEN
            NA1 = 1
         ELSE
            NA1 = NA + 1
         ENDIF
         NSA2 = NOSOTR( NA1 )

         IF( ( NSA1.EQ.NS1 .AND. NSA2.EQ.NS2 ) .OR.
     %       ( NSA1.EQ.NS2 .AND. NSA2.EQ.NS1 ) ) GOTO 9999

      ENDDO

C     NS1-NS2 N'EST PAS UNE ARETE DE NOSOTR
      NA = 0

 9999 RETURN
      END
