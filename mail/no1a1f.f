      SUBROUTINE  NO1A1F( NS1, NS2,  NOSOFA, NA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO LOCAL DE L'ARETE DE SOMMETS NOSOFA
C -----    PARMI LES 3 SOMMETS NOSOFA D'UN TRIANGLE
C
C ENTREES:
C --------
C NS1,NS2 : NUMERO DES 2 SOMMETS DE L'ARETE
C NOSOFA  : NUMERO DES 3 SOMMETS DU TRIANGLE
C
C SORTIE :
C --------
C NA : NUMERO DE 1 A 3 DE L'ARETE RETROUVEE
C      0 SI ARETE NON RETROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  SEPTEMBRE 1991
C2345X7..............................................................012
       INTEGER NOSOFA(1:3)
C
C      BOUCLE SUR LES 3 ARETES DU TRIANGLE
       DO 10 NA=1,3
C         ARETE I DE SOMMET NOSOFA(NA) ET NOSOFA(NA+1)
          IF( NA .NE. 3 ) THEN
             NA1 = NA + 1
          ELSE
             NA1 = 1
          ENDIF
          IF( (NS1 .EQ. NOSOFA(NA) .AND. NS2 .EQ. NOSOFA(NA1)) .OR.
     %        (NS2 .EQ. NOSOFA(NA) .AND. NS1 .EQ. NOSOFA(NA1)) ) RETURN
 10   CONTINUE
C
C     ARETE NON RETROUVEE
      NA = 0
      END
