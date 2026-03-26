      SUBROUTINE NOSTDA1S( NST,    NBTE1S, NOTE1S, NOTETR,
     %                     MXSTAR, NBSTAR, NOSTAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DU TABLEAU DES NO DES SOMMETS (NON NST) DES
C -----    ARETES DES TETRAEDRES DE SOMMET NST

C ENTREES:
C --------
C NST    : NO PTXYZD DU SOMMET DES TETRAEDRES
C NBTE1S : NOMBRE DE TETRAEDRES DE SOMMET NST
C NOTE1S : NO NOTETR DES NBTE1S TETRAEDRES DE SOMMET NST
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C MXSTAR : NOMBRE MAXIMAL DE NO A STOCKER DANS NOSTAR

C SORTIES:
C --------
C NBSTAR : NOMBRE DE SOMMETS DES ARETES ISSUES DE NST DES TETRAEDRES
C NOSTAR : NO NOTETR DES SOMMETS DES ARETES ISSUES DE NST DES TETRAEDRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Janvier 2019
C....................................................................012
      INTEGER  NOTETR(8,*), NOTE1S(NBTE1S), NOSTAR(MXSTAR)

      NBSTAR = 0

      DO 20 N = 1, NBTE1S

C        NO NOTETR DU TETRAEDRE N
         NTE = NOTE1S( N )
         IF( NOTETR(1,NTE) .LE. 0 ) GOTO 20

         DO 10 K = 1, 4

C           NO PTXYZD DU SOMMET K DU TETRAEDRE NTE
            NSK = NOTETR( K, NTE )

            IF( NSK .NE. NST ) THEN

C              AJOUT DE NSK SI CE N'EST PAS ENCORE FAIT
               DO L = 1, NBSTAR
                  IF( NSK .EQ. NOSTAR(L) ) GOTO 10
               ENDDO

               IF( NBSTAR .GE. MXSTAR ) THEN
                PRINT*,'nostda1s: SATURATION du TABLEAU NOSTAR MXSTAR=',
     %                  MXSTAR
                  GOTO 9999
               ENDIF
               NBSTAR = NBSTAR + 1
               NOSTAR( NBSTAR ) = NSK

            ENDIF

 10      ENDDO

 20   ENDDO

 9999 RETURN
      END
