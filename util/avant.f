      LOGICAL FUNCTION AVANT( NTAB1 , NTAB2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LA DATE DE NTAB1 EST ELLE ANTERIEURE A CELLE DE NTAB2 ?
C -----
C
C ENTREES :
C ---------
C NTAB1 , NTAB2 : LES 2 TABLEAUX CONTENANT LA DATE DANS LES 2-ERS MOTS
C
C SORTIE :
C --------
C AVANT  : .TRUE. SI DATE DE NTAB1 =< DATE NTAB2    .FALSE. SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS          MARS 1996
C.....................................................................
      INTEGER           NTAB1(1:2), NTAB2(1:2)
C
      DOUBLE PRECISION   DATE1,       DATE2
      INTEGER           NDATE1(1:2), NDATE2(1:2)
      EQUIVALENCE      (DATE1,NDATE1(1)), (DATE2,NDATE2(1))
C
C     LES DATES : 2 ENTIERS => 1 DOUBLE PRECISION
      NDATE1(1) = NTAB1(1)
      NDATE1(2) = NTAB1(2)
C
      NDATE2(1) = NTAB2(1)
      NDATE2(2) = NTAB2(2)
C
C     LA COMPARAISON ENTRE 2 DOUBLE PRECISION
      AVANT     = DATE1 .LE. DATE2
CCC      PRINT *,'AVANT: DATES=',DATE1,' ',DATE2,' RESULTAT=',AVANT
      END
