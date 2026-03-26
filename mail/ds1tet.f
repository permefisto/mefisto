      SUBROUTINE DS1TET( NT, NBSOTE, N1TEVI, NSTETR,
     %                       NO1TSO, N1TESO, NOTESO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LE TETRAEDRE NT DU TABLEAU NSTETR ET NOTESO
C -----
C
C ENTREES:
C --------
C NT     : NUMERO DANS NSTETR DU TETRAEDRE A SUPPRIMER
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C N1TEVI : NUMERO 1-ER TETRAEDRE VIDE DE NSTETR
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE si OCCUPE
C          ou NSTETR(2,NT) TETRAEDRE VIDE SUIVANT
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C N1TESO : NUMERO 1-ER TETRAEDRE VIDE DANS NOTESO
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      VIDE OU OCCUPE
C                      0 SI C'EST LE DERNIER
C IERR   : =0 SI PAS D'ERREUR
C          >0 SINON ( SATURATION DES FACES DE L'ETOILE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1991
C....................................................................012
      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*)
C
C     LE CHAINAGE NOTESO DE NT EST SUPPRIME DES TETRAEDRES DES 4 SOMMETS
      DO I=1,4
C
C        LE I-EME SOMMET DU TETRAEDRE NT
         NS = NSTETR(I,NT)
         IF( NS .LE. 0 ) GOTO 9999
C
C        RECHERCHE AVEC PRECEDENT DU TETRAEDRE NT
         N0 = 0
         N1 = NO1TSO( NS )
 10      IF( N1 .GT. 0 ) THEN
            NT1 = NOTESO( 1, N1 )
            IF( NT1 .NE. NT ) THEN
C              TETRAEDRE NON RETROUVE : PASSAGE AU SUIVANT
               N0 = N1
               N1 = NOTESO( 2, N1 )
               GOTO 10
            ELSE
C              TETRAEDRE RETROUVE : CHAINAGE SUPPRIME
               IF( N0 .EQ. 0 ) THEN
C                 C'EST LE PREMIER DU CHAINAGE
                  NO1TSO( NS ) = NOTESO( 2, N1 )
               ELSE
C                 IL EXISTE UN PRECEDENT
                  NOTESO( 2, N0 ) = NOTESO( 2, N1 )
               ENDIF
C              N1 REDEVIENT LIBRE DANS NOTESO
               NOTESO( 2, N1 ) = N1TESO
               N1TESO = N1
            ENDIF
         ENDIF
      ENDDO
C
C     NT EST AJOUTE AUX TETRAEDRES VIDES DE NSTETR
      NSTETR( 1, NT ) = 0
      NSTETR( 2, NT ) = N1TEVI
      N1TEVI = NT
C
 9999 RETURN
      END
