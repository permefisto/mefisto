      FUNCTION CHARX(ICHARX)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REALISER LA CONVERSION : INTEGER AX -> CHARACTER*X
C -----
C **********************************************************************
C ATTENTION : CETTE FONCTION EST DEPENDANTE MACHINE  ICI VERSION APOLLO
C             NBCAMO = NOMBRE DE CARACTERES PAR MOT = 1 ENTIER
C **********************************************************************
C ENTREE :
C --------
C ICHARX : ENTIER CONTENANT NBCAMO CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1984
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/nbcamo.inc"
      CHARACTER*(NBCAMO) CHARX, CHAR4
      INTEGER            ICHARX, I
C
      EQUIVALENCE (CHAR4,I)
C
      I      = ICHARX
      CHARX  = CHAR4
C
C     CONVERSION DE L'ENTIER EN NBCAMO CARACTERES
C     ATTENTION (A4) DOIT ETRE REMPLACE CI DESSOUS PAR (A'NBCAMO')
C     ========= SOUS PEINE D'ERREUR
CCC      WRITE (CHARX,'(A4)') ICHARX
      END
