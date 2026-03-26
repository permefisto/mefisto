      SUBROUTINE GIBB8( LVL1,NACU,CCST,SIZE,STPT,IDFLT,IDPTH,XC,
     &                  LVL2,NHIG,NLOW,
     &                  ISDIR      )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT DU SP:
C ----------
C  TROUVER LA STRUCTURE DE NIVEAU FINALE NOTEE 'LVL2' QUI VA SERVIR POUR
C  RENUMEROTAGE.
C ENTREES :
C ---------
C LVL1  : TABLEAU DE LA STRUCTURE CONSTRUITE A PARTIR DE SND1
C NACU  : TABLEAU DONNANT LE NOMBRE D'ELEMENTS DANS L'INTERSECTION DE
C         LA STRUCTURE LVL1 ET DE LA STRUCTURE LVL2 RENVERSEE PAR NIVEAU
C CCST  : TABLEAU DES ELEMENTS DE CHAQUE COMPOSANTE CONNEXE
C SIZE  : TABLEAU DES TAILLES DE CHAQUE COMPOSANTE CONNEXE
C STPT  : TABLEAU DES ADRESSES DANS CCST DE CHAQUE COMPOSANTE CONNEXE
C IDFLT : NUMERO DE LA MEILLEURE STRUCTURE ENTRE LVL1 ET LVL2
C IDPTH : HAUTEUR DE LA STRUCTURE DE NIVEAUX
C XC    : NOMBRE DE COMPOSANTES CONNEXES
C
C MODIFIEES :
C -----------
C LVL2  : EN ENTREE? CONTIENT LA STRUCTURE CONSTRUITE A PARTIR DE SND2
C         ET RENVERSEE; EN SORTIE CONTIENT LA STRUCTURE FINALE.
C NHIG  : TABLEAU DES ELEMENTS hi DU PAS 3 DE L'ALGORITHME II
C NLOW  : TABLEAU DES ELEMENTS li
C
C SORTIE :
C --------
C ISDIR : DEFINIT LE SENS DE PARCOURS DE LVL2 DANS LA RENUMEROTATION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR :  BARGACH MOHAMED LAN189 PARIS    OCTOBRE  1980
C MODIFICATION : DEFAIX THIERRY                 DECEMBRE 1989
C23456---------------------------------------------------------------012
      INTEGER SIZE,END,STPT,XC,CCST
      DIMENSION LVL1(1),LVL2(1),NACU(1),CCST(1),SIZE(1),STPT(1),NHIG(1),
     &NLOW(1)
C
C     BOUCLES SUR LES COMPOSANTES CONNEXES
C     ====================================
      DO 80 I=1,XC
         J=STPT(I)
         END=SIZE(I)+J-1
C
C        INITIALISATION DE NLOW ET NHIG
C        ------------------------------
         DO 10 K=1,IDPTH
            NHIG(K)=NACU(K)
            NLOW(K)=NACU(K)
   10    CONTINUE
C
C        BOUCLE SUR LES NOEUDS DE LA COMP./CONSTRUCTION DE NHIG ET NLOW
C        --------------------------------------------------------------
         DO 20 K=J,END
            INODE=CCST(K)
            LVLNH=LVL1(INODE)
            NHIG(LVLNH)=NHIG(LVLNH)+1
            LVLNL=LVL2(INODE)
            NLOW(LVLNL)=NLOW(LVLNL)+1
   20    CONTINUE
C
C        CALCUL DE MAX1=h0 et MAX2=l0
C        ----------------------------
         MAX1=0
         MAX2=0
         DO 30 K=1,IDPTH
            IF (2*NACU(K).NE.NLOW(K)+NHIG(K)) THEN
               IF (NHIG(K).GT.MAX1) MAX1=NHIG(K)
               IF (NLOW(K).GT.MAX2) MAX2=NLOW(K)
            END IF
   30    CONTINUE
         IT=1
         IF(MAX1.GT.MAX2) IT=2
         IF(MAX1.EQ.MAX2)IT=IDFLT
         IF(IT.EQ.2) THEN
C        ---LVL2 contient la bonne structure, mise a jour de NACU
            DO 70 K=1,IDPTH
               NACU(K)=NLOW(K)
   70       CONTINUE
          ELSE
C         ---LVL1 contient la bonne structure, on la copie dans LVL2
C         ---Mise a jour de NACU
            IF(I.EQ.1) ISDIR=-1
            DO 40 K=J,END
               INODE=CCST(K)
               LVL2(INODE)=LVL1(INODE)
   40       CONTINUE
            DO 50 K=1,IDPTH
               NACU(K)=NHIG(K)
   50       CONTINUE
         END IF
   80 CONTINUE
      RETURN
      END
