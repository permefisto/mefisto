      SUBROUTINE UNITABL( LETABL, NBVAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER LES VARIABLES NEGATIVES ou NULLES et les VARIABLES
C -----    STRICTEMENT POSITIVES MULTIPLES du TABLEAU LETABL POUR
C          OBTENIR UNE SEULE PRESENCE DANS LE TABLEAU LETABL

C ENTREE :
C --------
C NBVAR  : NOMBRE DE VARIABLES DU TABLEAU EN ENTREE

C MODIFIE:
C --------
C LETABL : LE TABLEAU DES VARIABLES AVANT ET APRES AVOIR OBTENU L'UNICITE

C SORTIE :
C --------
C NBVAR  : NOMBRE DE VARIABLES UNIQUES DU TABLEAU EN SORTIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY              Mai 2020
C2345X7..............................................................012
      INTEGER  LETABL( NBVAR )

      NBVAR0 = NBVAR

      IF( NBVAR0 .GT. 32 ) THEN

C        TRI CROISSANT DES NBVAR0 ENTIERS LETABL
         CALL TRITAENT( NBVAR0, LETABL )

C        SUPPRESSION DES MEMES VALEURS >0
         NBVAR = 0
         N1 = 0
 10      IF( N1 .GE. NBVAR0 ) GOTO 9999
         N1 = N1 + 1
         NV = LETABL( N1 )
         IF( NV .LE. 0 ) GOTO 10

         NBVAR = 1
         LETABL( 1 ) = NV

         DO 20 I = N1+1, NBVAR0

            NV = LETABL( I )
            IF( NV .LE. 0 ) GOTO 20

            IF( LETABL(NBVAR) .NE. NV ) THEN
               NBVAR = NBVAR + 1
               LETABL(NBVAR) = NV
            ENDIF

 20      ENDDO

      ELSE

C        UNICITE DES NBVAR0 VARIABLES DU TABLEAU LETABL
         NBVAR  = 0
         DO 50 I = 1, NBVAR0

            NV = LETABL( I )
            IF( NV .LE. 0 ) GOTO 50

            DO J = 1, NBVAR
               IF( NV .EQ. LETABL(J) ) GOTO 50
            ENDDO

            NBVAR = NBVAR + 1
            LETABL( NBVAR ) = NV

 50      ENDDO

      ENDIF

 9999 RETURN
      END
