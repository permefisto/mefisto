      SUBROUTINE TEOPFA( NBTETR, NSTETR, NTE0, NTEOPPO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RECHERCHE DU TETRAEDRE OPPOSE AUX 4 FACES DU TETRAEDRE NTE0
C -----
C ENTREES:
C --------
C NBTETR : NOMBRE DE TETRAEDRES
C NSTETR : NO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NTE0   : NO NSTETR DU TETRAEDRE INITIAL

C SORTIE :
C --------
C NTEOPPO: NO DES 4 TETRAEDRES OPPOSES AUX FACES DU TETRAEDRE NTE0
C          0 SI PAS DE TETRAEDRE OPPOSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         Novembre 2021
C2345X7..............................................................012
      INTEGER  NSTETR(8,NBTETR), NTEOPPO(4)

      INTEGER  NOSOFA0(3), NOSOFA1(3)
      INTEGER  NOSOFATE(3,4)
      DATA     NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /

      DO K0=1,4
         NTEOPPO( K0 ) = 0
      ENDDO

C     PARCOURS DES 4 FACES DU TETRAEDRE NTE0
      DO 10 K0=1,4

C        LES 3 NO DES SOMMETS DE LA FACE K DE NTE0
         DO M0=1,3
            NOSOFA0(M0) = NSTETR( NOSOFATE(M0,K0), NTE0 )
         ENDDO
C        TRI CROISSANT
         CALL TRI3NO( NOSOFA0, NOSOFA0 )

C        PARCOURS BESTIAL PAR MANQUE DE DONNEES
         DO 5 NTE1=1,NBTETR

            IF( NTE1 .EQ. NTE0 ) GOTO 5

C           PARCOURS DES 4 FACES DU TETRAEDRE NTE1
            DO K1=1,4

C              LES 3 NO DES SOMMETS DE LA FACE K1 DE NTE1
               DO M1=1,3
                  NOSOFA1(M1) = NSTETR( NOSOFATE(M1,K1), NTE1 )
               ENDDO
C              TRI CROISSANT
               CALL TRI3NO( NOSOFA1, NOSOFA1)

C              COMPARAISON DES 3 NO DE SOMMETS DES 2 FACES
               IF( NOSOFA0(1) .EQ. NOSOFA1(1) ) THEN
                  IF( NOSOFA0(2) .EQ. NOSOFA1(2) ) THEN
                     IF( NOSOFA0(3) .EQ. NOSOFA1(3) ) THEN
C                       FACE RETROUVEE
                        NTEOPPO( K0 ) = NTE1
                        GOTO 10
                     ENDIF
                  ENDIF
               ENDIF

            ENDDO
 5       ENDDO

 10   ENDDO

      RETURN
      END


