      SUBROUTINE GIBB4(IROOT,NVOI,MUNE,
     &                 LVL,IWK,LVLWTH,LVLBOT,LVLN,MAXLW,IBORT)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT DU SP:
C ----------
C  GENERER UNE STRUCTURE DE NIVEAU A PARTIR DE LA RACINE  'IROOT'.
C
C ENTREES :
C ---------
C IROOT : RACINE DONT ON VEUT EXTRAIRE LA STRUCTURE DE NIVEAU.
C NVOI  : TABLEAU QUI DONNE LA LISTE DES VOISINS
C MUNE  : TABLEAU DES POINTEURS ASSOCIES AU TABLEAU PRECEDENT
C
C SORTIES :
C ---------
C LVL    : TABLEAU QUI CONTIENT LA STRUCTURE DE NIVEAUX
C IWK    : TABLEAU INVERSE DU PRECEDENT
C LVLWTH : LONGUEUR DE DERNIER NIVEAU DE LVL
C LVLBOT : ADRESSE DANS IWK DU PREMIER NOEUD DU DERNIER NIVEAU
C LVLN   : HAUTEUR DE LA STRUCTURE LVL
C MAXLW  : PLUS GRANDE LONGUEUR DE NIVEAU RENCONTREE
C IBORT  : PARAMETRE PERMETTANT LE RETOUR SI MAXLW DEVIENT TROP GRAND
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR   : BARGACH MOHAMED LAN189 PARIS    OCTOBRE  1980
C MODIFICATIONS : DEFAIX THIERRY                  DECEMBRE 1989
C23456---------------------------------------------------------------012
      DIMENSION IWK(*),LVL(*),MUNE(*),NVOI(*)
C
C     INITIALISATION ET CONSTRUCTION DU PREMIER NIVEAU
C     ================================================
      MAXLW=0
      ITOP=LVLN
      INOW=LVLN
      LVLBOT=LVLN
      LVLTOP=LVLN+1
      LVLN=1
      LVL(IROOT)=1
      IWK(ITOP)=IROOT
C
C     CONSTRUCTION D'UN NIVEAU A PARTIR DU PRECEDENT
C     ==============================================
   10 LVLN=LVLN+1
C
   20 CONTINUE
C     Boucle sur les noeuds du niveau precedent
      IWKNOW=IWK(INOW)
      MU1=MUNE(IWKNOW)+1
      MU2=MUNE(IWKNOW+1)
      DO 30 J=MU1,MU2
         ITEST=NVOI(J)
         IF(LVL(ITEST).EQ.0) THEN
C        ---Noeud non encore numerote
            LVL(ITEST)=LVLN
            ITOP=ITOP+1
            IWK(ITOP)=ITEST
         END IF
   30 CONTINUE
      INOW=INOW+1
      IF(INOW.LT.LVLTOP) GO TO 20
C
C     TRAITEMENT DU NIVEAU COURANT ET FIN EVENTUELLE
C     ==============================================
      LVLWTH=LVLTOP-LVLBOT
      IF(MAXLW.LT.LVLWTH) MAXLW=LVLWTH
      IF(MAXLW.GE.IBORT) RETURN
      IF(ITOP.GE.LVLTOP) THEN
C     ---Le niveau courant contient au moins un noeud
C     ---Donc la numerotation peut ne pas etre terminee
C     ---On passe au niveau superieur
         LVLBOT=INOW
         LVLTOP=ITOP+1
         GO TO10
      END IF
      RETURN
      END
