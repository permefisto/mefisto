      SUBROUTINE TR2ADJ1AR( NT0, NOTRIA, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     2 TRIANGLES ADJACENTS A NT0 ONT ILS  UNE ARETE COMMUNE?
C ----

C ENTREES:
C --------
C NT0    : NUMERO NOTRIA DU TRIANGLE D'ANGLE TROP GRAND
C NOTRIA : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES

C SORTIE :
C --------
C NONOUI : =1 S'IL EXISTE UNE ARETE COMMUNE A 2 TRIANGLES ADJACENTS A NT0
C          =0 S'IL N'EXISTE PAS UNE TELLE ARETE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET   Saint Pierre du Perray            Novembre 2019
C2345X7..............................................................012
      INTEGER  NOTRIA(6,*)

C     RECHERCHE DE 2 TRIANGLES OPPOSES A NT0 AYANT UNE ARETE COMMUNE
      NTOP1 = NOTRIA( 6, NT0 )
      DO J=1,3
            NTOP2 = NOTRIA( 3+J, NT0 )
            IF( NTOP1 .GT. 0 .AND. NTOP2 .GT. 0 ) THEN
C              RECHERCHE D'UNE MEME ARETE COMMUNE
               DO J1=1,3
                  NST1 = NOTRIA( J1, NTOP1 )
                  DO J2=1,3
                     NST12 = NOTRIA( J2, NTOP2 )
                     IF( NST1 .EQ. NST12 ) THEN
C                       NST1 SOMMET COMMUN A NTOP1 et NTOP2

                        DO K1=1,3
                           IF( K1 .NE. J1 ) THEN
                              NST2 = NOTRIA( K1, NTOP1 )
                              DO K2=1,3
                                 IF( K2 .NE. J2 ) THEN
                                    NST22 = NOTRIA( K2, NTOP2 )
                                    IF( NST2 .EQ. NST22 ) THEN
C                                      NST1 et NST2 SOMMETS COMMUNS
C                                      A NTOP1 et NTOP2
C                                      => PAS D'IDENTIFICATION DE NST1 et NST2
C                                      DE PLUS, COMME L'EPAISSEUR DU MAILLAGE
C                                      EST FAIBLE NT0 N'EST PAS MODIFIE
                                       NONOUI = 1
                                       GOTO 9999
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
            NTOP1 = NTOP2
      ENDDO

C     PAS DE TELLE ARETE
      NONOUI = 0

 9999 RETURN
      END
