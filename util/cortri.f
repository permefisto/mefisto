      SUBROUTINE CORTRI( NT, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    VERIFIER LA COHERENCE DU TRIANGLE NT DE NOTRIA
C -----
C ENTREES :
C ---------
C NT     : NUMERO DANS NOTRIA DES NBT TRIANGLES
C NOTRIA : NUMERO DES 3 SOMMETS ET TRIANGLES OPPOSES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS         MAI 1995
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOTRIA(1:6,*)
C
      IF( NT .GT. 0 ) THEN
C
C        COHERENCE DES SOMMETS
         IF( NOTRIA(1,NT) .NE. 0 ) THEN
            IF( NOTRIA(1,NT) .EQ. NOTRIA(2,NT)  .OR.
     %          NOTRIA(1,NT) .EQ. NOTRIA(3,NT)  .OR.
     %          NOTRIA(2,NT) .EQ. NOTRIA(3,NT)  ) THEN
               WRITE(IMPRIM,*) 'PB AVEC TRIANGLE ',NT
               WRITE(IMPRIM,*) 'DE SOMMETS ',(NOTRIA(K,NT),K=1,6)
               CALL XVPAUSE
            ENDIF
C
C           COHERENCE DES TRIANGLES OPPOSES
            DO 8 J=4,5
               IF( NOTRIA(J,NT) .EQ. 0 ) GOTO 8
               DO 6 I=J+1,6
                  IF( NOTRIA(I,NT) .EQ. 0 ) GOTO 6
                  IF( NOTRIA(I,NT) .EQ. NOTRIA(J,NT) ) THEN
                     WRITE(IMPRIM,*) 'PB AVEC TRIANGLE ',NT
                     WRITE(IMPRIM,*) 'DE TRIANGLES OPPOSES ',
     %                               (NOTRIA(K,NT),K=1,6)
                     CALL XVPAUSE
                  ENDIF
 6             CONTINUE
 8          CONTINUE
         ENDIF
      ENDIF
      END
