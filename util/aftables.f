      SUBROUTINE AFTABLES( NOMTABL, NBLIGN, NBCOLO, NBTABL, TABLES,
     %                     NL1, NL2, NC1, NC2, NM1, NM2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES COEFFICIENTS
C -----    DES LIGNES NL1 A NL2
C          DES COLONNES NC1 A NC2
C          DES TABLES NM1 A NM2
C          DES TABLES(NBLIGN,NBCOLO,NBTABL) DE NOM NOMTABL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A&M University at QATAR       Mars 2012
C23456---------------------------------------------------------------012
      CHARACTER*(*)     NOMTABL
      DOUBLE PRECISION  TABLES( NBLIGN, NBCOLO, NBTABL )
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      NBC = NUDCNB( NOMTABL )
      DO NM=NM1, NM2
         DO NC=NC1,NC2
            DO NL=NL1,NL2
               WRITE(IMPRIM,10000)
     %               NOMTABL(1:NBC), NL, NC, NM, TABLES(NL,NC,NM)
            ENDDO
         ENDDO
      ENDDO
C
10000 FORMAT(A,I9,',',I9,',',I9,')=',G25.17)
C
      RETURN
      END
