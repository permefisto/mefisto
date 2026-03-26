       INTEGER FUNCTION NTYSOB( KNMSOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RETOURNER LE NUMERO (OU DIMENSION) DU SOUS-OBJET
C -----   DE NOM KNMSOB
C
C ENTREES :
C ---------
C KNMSOB : NOM DU SOUS-OBJET 'SOMMET' 'ARETE' 'FACE' 'CUBE' 'ELEMENT'
C
C SORTIE :
C --------
C NTYSOB : NUMERO COORRESPONDANT   1 OU 2 OU 3 OU 4 OU 5
C          0 SI KNMSOB N'EST PAS UN SOUS-OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS   AVRIL 1987
C.....................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     KNMSOB
      CHARACTER*7       NMSOBJ(1:5)
      DATA              NMSOBJ/ 'SOMMET' , 'ARETE' , 'FACE' ,
     %                          'CUBE'  , 'ELEMENT' /
C
C     LE NUMERO DANS LA LISTE NMSOBJ
      DO 10 NTYSOB = 1 , 5
         IF( NMSOBJ(NTYSOB)(1:4) .EQ. KNMSOB(1:4) ) RETURN
 10   CONTINUE
C
      NBLGRC(NRERR) = 1
      KERR(1) = 'TYPE INCORRECT DE SOUS-OBJET '// KNMSOB
      CALL LEREUR
      NTYSOB = 0
      END
