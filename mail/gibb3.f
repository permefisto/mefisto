      SUBROUTINE GIBB3( NOE,   LISTVOI, LPVOIS, IANDEG,
     &                  IALVL, IANPER,  IANDLS, IAIWK,
     &                  SND1,  SND2,    IALVL1, IALVL2, IDFLT, IDPTH )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT DU SP:
C ----------
C TROUVER LES 2 EXTREMITES SND1 ET SND2 DU PSEUDO-DIAMETRE AVEC LEURS
C STRUCTURES DE NIVEAUX RESPECTIVES LVL1 ET LVL2 DE HAUTEUR IDPTH,LE
C NUMERO DE LA STRUCTURE QUI DONNE LA LARGEUR DE BANDE LA PLUS FAIBLE
C EST CONTENU DANS IDFLT.
C
C ENTREES :
C ---------
C NOE    : NOMBRE TOTAL DE NOEUDS
C LISTVOI: TABLEAU QUI DONNE LA LISTE DES VOISINS
C LPVOIS : POINTEUR DERNIER NOEUD VOISIN DE CHAQUE NOEUD
C IANDEG : ADRESSE MCN DU TABLEAU DES DEGRES DES NOEUDS.
C
C MODIFIEES :
C -----------
C IALVL  : ADRESSE MCN D'UN TABLEAU DE STRUCTURES DE NIVEAUX INTERMEDIAIRES
C IANPER : ADRESSE MCN DU TABLEAU QUI CONTIENT LES NOEUDS DU DERNIER NIVEAU
C                         DE LA STRUCTURE D'EXTREMITE SND1
C IANDLS : ADRESSE MCN DU TABLEAU QUI CONTIENT LES MEMES NOEUDS QUE LE
C                         PRECEDENT MAIS ORDONNES DANS L'ORDRE DES DEGRES
C                         CROISSANT.
C IAIWK  : ADRESSE MCN DU TABLEAU INVERSE DE LVL
C
C SORTIES :
C ---------
C SND1,SND2     : EXTREMITES DU PSEUDO-DIAMETRE
C IALVL1,IALVL2 : ADRESSE MCN DES TABLEAUX DES STRUCTURES DE NIVEAUX
C                            ASSOCIEES A CES EXTREMITES
C IDFLT         : NUMERO DE LA MEILLEURE DES DEUX STRUCTURES
C IDPTH         : HAUTEUR DE CES STRUCTURES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: BARGACH MOHAMED LAN189 PARIS                    OCTOBRE  1980
C MODIFS : DEFAIX THIERRY                                  DECEMBRE 1989
C MODIFS : ALAIN PERRONNET  Saint PIERRE du PERRAY            Mars  2021
C23456---------------------------------------------------------------012
      INTEGER SND, SND1, TT, SND2, LISTVOI(*)
      include"./incl/pp.inc"
      COMMON  MCN(MOTMCN)
      INTEGER LPVOIS(0:NOE)
C
C     INITIALISATION
C     ==============
      NDXN = 0
      MTW1 = 0
      FLAG = 0
      MTW2 = NOE
      SND  = SND1
C
C     CONSTRUCTION DE LA STRUCTURE D'EXTREMITE SDN
C     ============================================
   10 DO 20 I=1,NOE
         MCN(IALVL-1+I)=0
   20 CONTINUE
      LVLN=1
      CALL GIBB4( SND, LISTVOI, LPVOIS,
     &            MCN(IALVL),MCN(IAIWK),LVLWTH,LVLBOT,LVLN,MAXLW,MTW2 )
      IF(FLAG.GE.1) GO TO 50
      FLAG=1
   30 CONTINUE
C
C     LVL1 RECOIT LA STRUCTURE DE NIVEAU LVL
C     --------------------------------------
      IDPTH=LVLN-1
      MTW1=MAXLW
      DO 40 I=1,NOE
         MCN(IALVL1-1+I)=MCN(IALVL-1+I)
   40 CONTINUE
C
C     ON ISOLE LES NOEUDS DU DERNIER NIVEAU
C     -------------------------------------
      NDXN=1
      NDXL=0
      MTW2=NOE
      TT=LVLBOT+LVLWTH-1
      DO 100 I=LVLBOT,TT
         J=I-LVLBOT+1
         MCN(IANPER-1+J)=MCN(IAIWK-1+I)
  100 CONTINUE
      CALL GIBB5( MCN(IANDEG), NDXL, LVLWTH, MCN(IANPER),
     &            MCN(IANDLS))
      SND=MCN(IANDLS-1+1)
      GO TO 10
   50 CONTINUE
C
C     TRAITEMENT DE LA DERNIERE STRUCTURE CALCULEE
C     --------------------------------------------
      IF(IDPTH.LT.LVLN-1) THEN
C     ---La Derniere Structure construite est de Hauteur superieure
C     ---donc elle est meilleure.
         SND1=SND
         GO TO 30
      END IF
      IF(MAXLW.LT.MTW2) THEN
C     ---La Structure construite a partir de SDN est meilleur en
C     ---epaisseur de bande.
         MTW2=MAXLW
         SND2=SND
         DO 70 I=1,NOE
            MCN(IALVL2-1+I)=MCN(IALVL-1+I)
   70    CONTINUE
      ENDIF
      IF(NDXN.NE.NDXL) THEN
C     ---Il reste des noeuds a explorer
         NDXN=NDXN+1
         SND=MCN(IANDLS-1+NDXN)
         GO TO 10
      END IF
C
C     FIN DE LA RECHERCHE
C     ===================
      IDFLT=1
      IF(MTW2.LE.MTW1) IDFLT=2
C
      RETURN
      END
