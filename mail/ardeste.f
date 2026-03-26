      SUBROUTINE ARDESTE( NS1, NS2, NBTECF, NOTECF, NOTETR, NTE, NA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER UNE ARETE NS1-NS2 PARMI LES ARETES DES NBTECF
C -----    TETRAEDRES DE NOTETR

C ENTREES:
C --------
C NS1 NS2: NO DES 2 SOMMETS DE L'ARETE A RECHERCHER
C NBTECF : NOMBRE DE TETRAEDRES
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES:
C --------
C NTE    : NUMERO NOTETR D'UN TETRAEDRE D'ARETE NS1-NS2
C          =0  SI AUCUNE ARETE DES TETRAEDRES N'EST L'ARETE NS1-NS2
C NA     : NUMERO DE 1 A 6 DE L'ARETE NS1-NS2 DE NTE
C          =0 SI NON RETROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  VEULETTES sur MER               Fevrier 2018
C23456...............................................................012
      INTEGER  NOTECF(NBTECF), NOTETR(8,*), NOSOARTE(2,6)
      DATA     NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /

      DO N = 1, NBTECF

         NTE = NOTECF( N )

         IF( NTE.GT.0 .AND. NOTETR(1,NTE).GT.0 ) THEN

            DO NA = 1, 6

               IF( NA .EQ. 3 ) THEN
                  NN = 1
               ELSE
                  NN = NA+1
               ENDIF
               NSA1 = NOTETR( NA, NTE )
               NSA2 = NOTETR( NN, NTE )

               IF( ( NSA1 .EQ. NS1 .AND. NSA2 .EQ. NS2 ) .OR.
     %             ( NSA1 .EQ. NS2 .AND. NSA2 .EQ. NS1 ) ) THEN
C                 ARETE NS1-NS2 EST L'ARETE NA DU TETRAEDRE NTE
                  GOTO 9999
               ENDIF

            ENDDO

         ENDIF

      ENDDO

C     ARETE NS1-NS2 NON RETROUVEE DANS NOTECF
      NTE = 0
      NA  = 0

 9999 RETURN
      END
