      SUBROUTINE JOTHE7(NUOBSD,NBJOIN,MXNOJO,
     S                  NBOBJE,NBNUJO,NUNOJO,NUSD,NDCAL,
     S                  DG,JOIM,NBJOIM,DIAG,LPJOIN,NPJOIN)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SOMMATION DES COEFFICIENTS DIAGONAUX DES NOEUDS SUR LE JOINT
C -----
C
C ENTREES :
C ---------
C NUOBSD : LE NUMERO D'OBJET DU SOUS-DOMAINE TRAITE
C NBJOIN : LE NOMBRE DE JOINTS
C MXNOJO : LE NOMBRE MAXIMUM DE NOEUDS D'UN JOINT
C NBOBJE : LES TYPES ET NUMEROS DES OBJETS ASSOCIES AUX JOINTS
C NBNUJO : LE NOMBRE DE NOEUDS ASSOCIES AUX JOINTS
C NUNOJO : LES NUMEROS DES NOEUDS ASSOCIES AUX JOINTS
C NUSD   : LA LISTE DES NUMEROS DES SOUS-DOMAINES ASSOCIES AUX JOINTS
C NDCAL  : DECALAGE POUR LA NUMEROTATION GLOBALE
C DG     : LA DIAGONALE GLOBALE
C JOIM   : LISTE DES NBJOIM JOINTS MAITRES
C JOIE   : LISTE DES NBJOIE JOINTS ESCLAVES
C
C SORTIE :
C --------
C DIAG   : LA DIAGONALE GLOBALE DES NOEUDS DES JOINTS
C NPJOIN : LE NOMBRE TOTAL DE POINTS SUR LES JOINTS
C LPJOIN : LA LISTE DE LEUR NUMERO GLOBAL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1994
C23456---------------------------------------------------------------012
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSIMTA / NOIMPR
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))
      DIMENSION JOIM(NBJOIM),NUSD(1),NDCAL(1),LPJOIN(1)
      DIMENSION NBOBJE(NBJOIN,4),NBNUJO(NBJOIN,2),
     S          NUNOJO(NBJOIN,2,MXNOJO)
      DOUBLE PRECISION DIAG(1),DG(1)
C
      IM = 0
      IE = 0
      DO 1 J = 1 , NBJOIM
C
C        COTE ESCLAVE / COTE MAITRE
C        --------------------------
C        LE NUMERO DU JOINT MAITRE
         NJ = JOIM(J)
         NJ10 = (NJ-1)*10
C        LES NUMEROS DES OBJETS ASSOCIES AU JOINT
         NOBJE1 = NBOBJE(NJ,2)
         NOBJE2 = NBOBJE(NJ,4)
         IF (NUOBSD.EQ.NOBJE1) THEN
            NBNOJM=NBNUJO(NJ,1)
            NBNOJE=NBNUJO(NJ,2)
            NUSDM=NUSD(NJ*2-1)
            NUSDE=NUSD(NJ*2)
            IM=1
            IE=2
         ELSE IF (NUOBSD.EQ.NOBJE2) THEN
            NBNOJE=NBNUJO(NJ,1)
            NBNOJM=NBNUJO(NJ,2)
            NUSDE=NUSD(NJ*2-1)
            NUSDM=NUSD(NJ*2)
            IE=1
            IM=2
         ENDIF
C=========================================================================
         IF (NOIMPR.GT.10)
     S   WRITE(IMPRIM,1000) NJ,NBNOJM,NBNOJE,NUSDM,NUSDE
C=========================================================================
C
C        PREMIER NOEUD
C        LE NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,1)
C        LE NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
         NPJOIN = NPJOIN + 1
         LPJOIN(NPJOIN) = NUGLMA
         DIAG(NPJOIN) = DG(NUGLMA)
C=========================================================================
         IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,2000) NPJOIN,NUGLMA,DIAG(NPJOIN)
C=========================================================================
         DO 2 K = 2 , NBNOJM - 1
C           LE NUMERO LOCAL COTE MAITRE
            NULOMA = NUNOJO(NJ,IM,K)
C           LE NUMERO GLOBAL COTE MAITRE
            NUGLMA = NULOMA + NDCAL(NUSDM)
C           LE NUMERO LOCAL COTE ESCLAVE
            NULOES = NUNOJO(NJ,IE,K)
C           LE NUMERO GLOBAL COTE ESCLAVE
            NUGLES = NULOES + NDCAL(NUSDE)
            NPJOIN = NPJOIN + 1
            LPJOIN(NPJOIN) = NUGLMA
            DIAG(NPJOIN) =  DG(NUGLMA) + DG(NUGLES)
C=========================================================================
         IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,2000) NPJOIN,NUGLMA,DIAG(NPJOIN)
C=========================================================================
2        CONTINUE
C        DERNIER NOEUD
C        LE NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,NBNOJM)
C        LE NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
         NPJOIN = NPJOIN + 1
         LPJOIN(NPJOIN) = NUGLMA
         DIAG(NPJOIN) = DG(NUGLMA)
C=========================================================================
         IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,2000) NPJOIN,NUGLMA,DIAG(NPJOIN)
C=========================================================================
C
 1    CONTINUE
C
      RETURN
1000  FORMAT(1X/,80(1H-)/
     S       ' THE7 : JOINT',I5,' NBN MAITRES ',I5,
     S       ' NBN ESCLAVES',I5,' SD MA/ES',2I5)
2000  FORMAT(I5,' NOEUD GLOBAL ',I5,' VALEUR ',D15.6)
      END
