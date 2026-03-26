      SUBROUTINE ARCO2TET( NSTET1, NSTET2, NAR1, NAR2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TROUVER DANS CHACUN DES TETRAEDRES LE NO LOCAL DE L'ARETE
C ----   COMMUNE AUX 2 TETRAEDRES

C ENTREES:
C --------
C NSTET1 : NO DES 4 SOMMETS DU TETRAEDRE 1
C NSTET2 : NO DES 4 SOMMETS DU TETRAEDRE 2

C SORTIES:
C --------
C NAR1   : NO DE 1 A 6 DE L'ARETE COMMUNE DANS LE TETRAEDRE 1
C          0 SI PAS D'ARETE COMMUNE
C NAR2   : NO DE 1 A 6 DE L'ARETE COMMUNE DANS LE TETRAEDRE 2
C          0 SI PAS D'ARETE COMMUNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY           Janvier 2018
C2345X7..............................................................012
      INTEGER  NSTET1(4), NSTET2(4)

C     NO DES 2 SOMMETS DES 6 ARETES D'UN TETRAEDRE
      INTEGER  NOSOARTE(2,6)
      DATA     NOSOARTE / 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      DO NAR1=1,6

C        NO DES 2 SOMMETS DE L'ARETE NAR1 DE NTET1
         NS1T1 = NSTET1( NOSOARTE(1,NAR1) )
         NS2T1 = NSTET1( NOSOARTE(2,NAR1) )
         DO NAR2=1,6
            NS1T2 = NSTET2( NOSOARTE(1,NAR2) )
            NS2T2 = NSTET2( NOSOARTE(2,NAR2) )
            IF( ( NS1T1 .EQ. NS2T2 .AND. NS2T1 .EQ. NS1T2) .OR.
     %          ( NS1T1 .EQ. NS1T2 .AND. NS2T1 .EQ. NS2T2) ) THEN
C              ARETE COMMUNE RETROUVEE
               RETURN
            ENDIF
         ENDDO
      ENDDO

      NAR1 = 0
      NAR2 = 0

      RETURN
      END
