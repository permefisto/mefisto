      SUBROUTINE TETOPDOU( NTE, NOTETR, NF1, NF2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     UN TETRAEDRE OPPOSE AU TETRAEDRE NTE EST IL DOUBLE OPPOSE?
C -----

C ENTREES :
C ---------
C NTE     : NUMERO DANS NOTETR DU TETRAEDRE A EXAMINER
C NOTETR  : LISTE DES TETRAEDRES
C           SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C           TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C           OPPOSE est de l'AUTRE COTE DE LA FACE
C           1: 123      2: 234      3: 341      4: 412

C SORTIE  :
C ---------
C NF1     : =0 SI PAS DE TETRAEDRE OPPOSE DOUBLE
C           >0 NO DE 1     A 4 DE LA FACE DE NTE
C              DE TETRAEDRE OPPOSE DOUBLE
C NF2     : =0 SI PAS DE TETRAEDRE OPPOSE DOUBLE
C           >0 NO DE NF1+1 A 4 DE L'AUTRE FACE DE NTE
C              DE TETRAEDRE OPPOSE DOUBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Veulettes sur mer              Septembre 2019
C2345X7..............................................................012
      INTEGER NOTETR(8,*)

      IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

         DO NF1=1,4

            NTOP1 = NOTETR( 4+NF1, NTE )

            IF( NTOP1 .GT. 0 ) THEN

               DO NF2=NF1+1,4

                  NTOP2 = NOTETR( 4+NF2, NTE )
                  IF( NTOP1 .EQ. NTOP2 ) THEN

C                    NTEOP1 est un TETRAEDRE OPPOSE DOUBLE POUR NTE
                     NONOUI = 1
                     PRINT*,'tetopdou: le TETRAEDRE',NTE,
     %                      ' a 2 TETRAEDRES OPPOSES IDENTIQUES'
                     PRINT*,'tetopdou: NOTETR(',NTE,')=',
     %                      (NOTETR(kk,NTE),kk=1,8)
                     GOTO 9999

                  ENDIF

               ENDDO
            ENDIF

         ENDDO

      ENDIF

C     PAS DE TETRAEDRE OPPOSE DOUBLE
      NF1 = 0
      NF2 = 0

 9999 RETURN
      END


   
