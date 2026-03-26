      SUBROUTINE GIBB5( NDEG, X1, X2, STK2,
     %                  STK1              )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT DU SP:
C ----------
C RANGER LES X2 ELEMENTS DE STK2 SUIVANT L'ORDRE DECROISSANT DES DEGRES
C DANS STK1 A PARTIR DE L'ADRESSE X1.
C
C ENTREES :
C ---------
C NDEG : TABLEAU DES DEGRES DE CHAQUE NOEUD
C X1   : ADRESSE DANS STK1 A PARTIR DE LAQUELLE ON MET LES ELTS DE STK2
C X2   : NOMBRE D'ELEMENTS DE STK2 A RANGER
C STK2 : TABLEAU DES ELEMENTS A RANGER
C
C MODIFIEES :
C -----------
C STK1 : TABLEAU DANS LEQUEL ON RANGE LES ELEMENTS DE STK2
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR   : BARGACH MOHAMED LAN189 PARIS    OCTOBRE  1980
C MODIFICATIONS : DEFAIX THIERRY                  DECEMBRE 1989
C23456---------------------------------------------------------------012
      INTEGER   STK1, STK2, X1, X2, TEMP
      DIMENSION STK1(*), STK2(*)
      DIMENSION NDEG(*)
C
C     CLASSEMENT DES ELEMENTS DE STK2 PAR TRI A BULLES
C     ================================================
      IND=X2
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) GO TO 30
      DO 20 I=1,IND
         J=I+1
         ISTK2=STK2(I)
         JSTK2=STK2(J)
         IF(NDEG(ISTK2).GT.NDEG(JSTK2)) THEN
            ITEST=1
            TEMP=STK2(I)
            STK2(I)=STK2(J)
            STK2(J)=TEMP
         END IF
   20 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
   30 CONTINUE
C
C     AFFECTATION DANS STK1 DES ELEMENTS DE STK2 RANGES
C     =================================================
      DO 40 I=1,X2
         X1=X1+1
         STK1(X1)=STK2(I)
   40 CONTINUE
      RETURN
      END
