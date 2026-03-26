      SUBROUTINE TE3TOP( NTE, NOTETR, NTEOP1, NTEOP4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     RECHERCHER SI PARMI LES 4 TETRAEDRES OPPOSES AUX FACES
C ----     IL EXISTE UN TETRAEDRE TRIPLE NTEOP1

C ENTREES:
C --------
C NTE    : NUMERO NOTETR DU TETRAEDRE A TRAITER
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES:
C --------
C NTEOP1 : >0   NUMERO NOTETR DU TETRAEDRE OPPOSE TRIPLE
C           =-1 PAS DE TETRAEDRE OPPOSE TRIPLE
C NTEOP4 : SI NTEOP1>0 TETRAEDRE OPPOSE TRIPLE EXISTE NTEOP4 EST
C             LE 4-EME TETRAEDRE OPPOSE A UNE FACE DU TETRAEDRE NTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  Saint PIERRE DU PERRAY              Mars 2019
C....................................................................012
      INTEGER  NOTETR(8,*)

C     BOUCLE SUR LES 4 TETRAEDRES OPPOSES DE NTE
      DO K=5,6
         NTEOP1 = NOTETR(K,NTE)
         IF( NTEOP1 .GT. 0 ) THEN
            DO L=K+1,8
               NTEOP2 = NOTETR(L,NTE)
               IF( NTEOP2 .EQ. NTEOP1 ) THEN
                  DO M=L+1,8
                     NTEOP3 = NOTETR(M,NTE)
                     IF( NTEOP3 .EQ. NTEOP1 ) THEN
C                        NTEOP1 TETRAEDRE OPPOSE A 3 FACES DE NTE

C                        NTEOP4 TETRAEDRE OPPOSE A LA FACE RESTANTE
                         DO I=5,8
                            NTEOP4 = NOTETR(I,NTE)
                            IF( NTEOP4 .NE. NTEOP1 ) THEN
                               GOTO 9999
                            ENDIF
                         ENDDO

                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
       ENDDO

       NTEOP1 = -1
       NTEOP4 = -1

 9999  RETURN
       END
