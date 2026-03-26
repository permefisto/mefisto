      SUBROUTINE COMPENTP( NBENT0, ENTIERS,  NBENTP1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    COMPRESSION DU TABLEAU ENTIERS DES VALEURS NEGATIVES OU NULLES
C -----

C ENTREES:
C --------
C NBENT0 : NOMBRE DE VALEURS ENTIERES DU TABLEAU ENTIERS

C MODIFIE:
C --------
C ENTIERS : TABLEAU D'ENTIERS COMPRESSE EN OUBLIANT LES ENTIERS <=0

C SORTIE :
C --------
C NBENTP1: NOMBRE D'ENTIERS STRICTEMENT POSITIFS DU TABLEAU ENTIERS
C+++++++++++++++++++++++++++ +++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY      Mai 2015
C23456...............................................................012
      INTEGER  ENTIERS(NBENT0)

C     COMPRESSION DU TABLEAU D'ENTIERS
      NBENTP1 = 0
      DO K = 1, NBENT0
         N = ENTIERS( K )
         IF( N .GT. 0 ) THEN
            NBENTP1 = NBENTP1 + 1
            ENTIERS( NBENTP1 ) = N
         ENDIF
      ENDDO

      RETURN
      END
