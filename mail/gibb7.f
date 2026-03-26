      SUBROUTINE GIBB7( XC,
     &                  SIZE,STPT,
     &                  TT        )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT DU SP:
C ----------
C  REORDONNE LES 2 TABLEAUX SIZE ET STPT SUIVANT L'ORDRE CROISSANT DES
C  ELEMENTS DE SIZE.
C
C ENTREES :
C ---------
C XC   : NOMBRE DE COMPOSANTES CONNEXES
C
C MODIFIEES :
C -----------
C SIZE : TABLEAU DES TAILLES DES COMPOSANTES CONNEXES
C STPT : TABLEAU DES NUMEROS DU PREMIER NOEUD DE CHAQUE COMPOSANTE
C
C SORTIE :
C --------
C TT   : EGAL A ZERO SI AUCUNE COMPOSANTE CONNEXE, EGAL A UN SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR   : BARGACH MOHAMED LAN189 PARIS    OCTOBRE  1980
C MODIFICATIONS : DEFAIX THIERRY                DECEMBRE 1989
C23456---------------------------------------------------------------012
      INTEGER TEMP,SIZE,XC,STPT,TT
      DIMENSION SIZE(1),STPT(1)
C
      TT=0
      IF(XC.EQ.0) RETURN
      TT=1
      IND=XC
C
C     CLASSEMENT PAR TRI A BULLES
C     ===========================
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) RETURN
      DO 20 I=1,IND
         J=I+1
         IF(SIZE(I).LE.SIZE(J)) THEN
            ITEST=1
            TEMP=SIZE(I)
            SIZE(I)=SIZE(J)
            SIZE(J)=TEMP
            TEMP=STPT(I)
            STPT(I)=STPT(J)
            STPT(J)=TEMP
         END IF
   20 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
      RETURN
      END
