      INTEGER FUNCTION ICHARX(CHARX)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   REALISER LA CONVERSION : CHARACTER*X -> INTEGER AX
C **********************************************************************
C ATTENTION : CETTE FONCTION EST DEPENDANTE MACHINE ICI VERSION IBM
C **********************************************************************
C ENTREE :
C --------
C CHARX  : CHAINE DE X CARACTERES A CONVERTIR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS    OCTOBRE 1984
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/nbcamo.inc"
      CHARACTER*(*)      CHARX
      CHARACTER*(NBCAMO) BUFFER
      INTEGER            I
      EQUIVALENCE       (I,BUFFER)
C
C     LE BUFFER RECOIT CHARX ET PORTE SA LONGUEUR A NBCAMO CARACTERES
      BUFFER = CHARX
      ICHARX = I
C
C     CONVERSION DE BUFFER NBCAMO CARACTERES EN UN ENTIER
C     POUR VERSION /= APOLLO REVOIR (A4) CI DESSOUS A REMPLACER PAR
C     (A NBCAMO)
CCC      READ (BUFFER,'(A4)') ICHARX
      END
