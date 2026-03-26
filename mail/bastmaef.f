      SUBROUTINE BASTMAEF( NBSTEF, XYZSOM, NOSTAB, N1STEF, NBEF, NOSTEF)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    BARYCENTRAGE/2 DES SOMMETS D'UN MAILLAGE D'ELEMENTS FINIS
C -----    EXCEPTES LES SOMMETS DE NUMERO NOSTAB NON NUL

C ENTREES:
C --------
C NBSTEF : NOMBRE DE SOMMETS DES ELEMENTS FINIS DU MAILLAGE
C NOSTAB : =0 SI SOMMET A BARYCENTRER
C          >0 SINON
C N1STEF : MAXIMUM DU PREMIER INDICE DU TABLEAU NOSTEF
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NOSTEF : NUMERO DE 1 A NBSTEF DES N1STEF SOMMETS DU MAILLAGE
C           ou =0 SI PAS DE SOMMET

C MODIFIE:
C --------
C XYZSOM : EN ENTREE XYZ DES NBSTEF SOMMETS DU MAILLAGE
C          EN SORTIE XYZ DES NBSTEF SOMMETS BARYCENTRES D'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY              Mai 2019
C23456...............................................................012
      INTEGER  NOSTAB(NBSTEF),   NOSTEF(N1STEF,NBEF)
      REAL     XYZSOM(3,NBSTEF), XYZBAR(3)

      NBBAR = 0
      DO 20 NS = 1, NBSTEF

         IF( NOSTAB( NS ) .NE. 0 ) GOTO 20

C        NS EST UN SOMMET INTERNE DEPLACE AU BARYCENTRE DE SES VOISINS
         XYZBAR( 1 ) = 0.
         XYZBAR( 2 ) = 0.
         XYZBAR( 3 ) = 0.
         NBSTBA = 0

C        RECHERCHE EXHAUSTIVE DES TRIANGLES DE SOMMET NS
         DO 10 NT = 1, NBEF

            DO I = 1, N1STEF

C              NUMERO DU SOMMET I DU TRIANGLE NT
               NSTI = NOSTEF( I, NT )
               IF( NSTI .GT. 0 ) THEN

                  IF( NSTI .EQ. NS ) THEN
C                    AJOUT AU BARYCENTRE DES COORDONNEES
C                    DES AUTRES SOMMETS
                     DO J = 1, N1STEF
                        IF( J .NE. I ) THEN
C                          NUMERO DU SOMMET J DU TRIANGLE NT
                           NSTJ = NOSTEF( J, NT )
                           IF( NSTJ .GT. 0 ) THEN
                              DO K=1,3
                                 XYZBAR(K) = XYZBAR(K) + XYZSOM(K,NSTJ)
                              ENDDO
                              NBSTBA = NBSTBA + 1
                           ENDIF
                        ENDIF
                     ENDDO
                     GOTO 10
                  ENDIF

               ENDIF
            ENDDO

C           PASSAGE AU TRIANGLE SUIVANT
 10      ENDDO

C        NS EST PLACE AU BARYCENTRE DE SES VOISINS ET DE LUI MEME
         DO K=1,3
            XYZSOM( K, NS )=( XYZSOM( K, NS ) + XYZBAR( K )/NBSTBA ) / 2
         ENDDO
         NBBAR = NBBAR + 1

 20   ENDDO

      PRINT*,'bastmaef: BARYCENTRAGE de',NBBAR,
     %       ' SOMMETS du MAILLAGE'
      RETURN
      END
