      SUBROUTINE AFTABLE( NOMTABL, NBLIGN, NBCOLO, TABLE,
     %                    NL1, NL2, NC1, NC2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES COEFFICIENTS
C -----    DES LIGNES NL1 A NL2
C          DES COLONNES NC1 A NC2
C          DE LA TABLE(NBLIGN,NBCOLO) DE NOM NOMTABL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A&M University at QATAR       Mars 2012
C23456---------------------------------------------------------------012
      CHARACTER*(*)     NOMTABL
      DOUBLE PRECISION  TABLE( NBLIGN, NBCOLO )
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      NBC = NUDCNB( NOMTABL )
      DO NL=NL1,NL2
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,10000)
     %        (NOMTABL(1:NBC), NL, NC, TABLE(NL,NC),NC=NC1,NC2)
      ENDDO
C
10000 FORMAT(A,'(',I9,',',I9,')=',G25.17)
C
      RETURN
      END
