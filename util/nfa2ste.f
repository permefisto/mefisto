      SUBROUTINE NFA2STE( NS1, NS2, NF1, NF2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DANS UN TETRAEDRE, RETROUVER LE NUMERO NF1<NF2 ( 1 A 4 )
C -----   DES 2 FACES D'ARETE COMMUNE NS1 (1 A 3 ) < NS2 ( 2 A 4 )
C         POUR LES NUMEROTATIONS LOCALES AU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY      Mai 2018
C2345X7..............................................................012
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

C     RECHERCHE DE LA PREMIERE FACE D'ARETE NS1-NS2
      DO NF = 1, 4
         DO K1 = 1, 3
            IF( NOSOFATE(K1,NF) .EQ. NS1 ) THEN
               DO K2 = 1, 3
                  IF( K2 .NE. K1 ) THEN
                     IF( NOSOFATE(K2,NF) .EQ. NS2 ) THEN
                        NF1 = NF
                        GOTO 10
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

C     RECHERCHE DE LA SECONDE FACE D'ARETE NS1-NS2
 10   DO NF = 1, 4
         IF( NF .NE. NF1 ) THEN
            DO K1 = 1, 3
               IF( NOSOFATE(K1,NF) .EQ. NS1 ) THEN
                  DO K2 = 1, 3
                     IF( K2 .NE. K1 ) THEN
                        IF( NOSOFATE(K2,NF) .EQ. NS2 ) THEN
                           NF2 = NF
                           GOTO 9999
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

 9999 RETURN
      END
