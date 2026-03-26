      SUBROUTINE TRIA( TAB,NTAB,TAB1,TAB2,LTAB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RANGER LES COEFFICIENTS DU TABLEAU TAB DANS L'ORDRE CROISSANT
C ---    ET DE MANIERE COHERENTE CEUX DU TABLEAU NTAB
C
C ENTREES:
C --------
C TAB  : TABLEAU A TRIER , DE LONGUEUR LTAB
C NTAB : TABLEAU AUXILIAIRE, DE LONGUEUR LTAB
C TAB1 : TABLEAU AUXILIAIRE, DE LONGUEUR LTAB
C TAB2 : TABLEAU AUXILIAIRE, DE LONGUEUR LTAB
C
C SORTIES:
C --------
C TAB  : TABLEAU APRES LE TRI
C NTAB : TABLEAU APRES LE TRI
C TAB1 : TABLEAU APRES LE TRI
C TAB2 : TABLEAU APRES LE TRI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION TAB(LTAB),NTAB(LTAB),TAB1(LTAB),TAB2(LTAB)
      COMMON / EPSSSS / EPZERO, EPSXYZ
C
      DO 1 K=LTAB,2,-1
C
         INDIC=0
C
         DO 2 I=1,K-1
            IP1=I+1
            IF( TAB(I) .GT. TAB(IP1)+EPZERO ) THEN
C
C              ECHANGE
               INDIC=1
C
               TABI=TAB(I)
               TAB(I)=TAB(IP1)
               TAB(IP1)=TABI
C
               NTABI=NTAB(I)
               NTAB(I)=NTAB(IP1)
               NTAB(IP1)=NTABI
C
               TABI1=TAB1(I)
               TAB1(I)=TAB1(IP1)
               TAB1(IP1)=TABI1
C
               TABI2=TAB2(I)
               TAB2(I)=TAB2(IP1)
               TAB2(IP1)=TABI2
C
            END IF
2        CONTINUE
C
         IF( INDIC .EQ. 0 ) RETURN
1     CONTINUE
C
      RETURN
      END
