      SUBROUTINE JOTHE8(NUOBSD,NBJOIN,MXNOJO,NTDL,
     S                  NBOBJE,NBNUJO,NUNOJO,NUSD,NDCAL,
     S                  VG,U,JOIM,NBJOIM,JOIE,NBJOIE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EXTRACTION DES VALEURS D'UN VECTEUR SUR LE JOINT
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
C VG     : LE VECTEUR LOCAL DU SOUS-DOMAINE
C JOIM   : LISTE DES NBJOIM JOINTS MAITRES
C JOIE   : LISTE DES NBJOIE JOINTS ESCLAVES
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1994
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DIMENSION JOIM(NBJOIM),JOIE(NBJOIE),NUSD(1),NDCAL(1)
      DIMENSION NBOBJE(NBJOIN,4),NBNUJO(NBJOIN,2),
     S          NUNOJO(NBJOIN,2,MXNOJO)
      DOUBLE PRECISION VG(1),U(NTDL)
C
      IM = 0
      IE = 0
      DO 1 J = 1 , NBJOIM
C
C        COTE ESCLAVE / COTE MAITRE
C        --------------------------
C        LE NUMERO DU JOINT MAITRE
         NJ = JOIM(J)
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
         WRITE(IMPRIM,1000) NJ,NBNOJM,NBNOJE,NUSDM,NUSDE
C=========================================================================
         DO 2 K = 1 , NBNOJM
C           LE NUMERO LOCAL COTE MAITRE
            NULOMA = NUNOJO(NJ,IM,K)
C           LE NUMERO GLOBAL COTE MAITRE
            NUGLMA = NULOMA + NDCAL(NUSDM)
C=========================================================================
            WRITE(IMPRIM,3000) K,NULOMA,VG(NUGLMA),U(NULOMA)
C=========================================================================
 2       CONTINUE
 1    CONTINUE
C
      DO 3 J = 1 , NBJOIE
C
C        COTE ESCLAVE / COTE MAITRE
C        --------------------------
C        LE NUMERO DU JOINT MAITRE
         NJ = JOIE(J)
C        LES NUMEROS DES OBJETS ASSOCIES AU JOINT
         NOBJE1 = NBOBJE(NJ,2)
         NOBJE2 = NBOBJE(NJ,4)
         IF (NUOBSD.EQ.NOBJE1) THEN
            NBNOJE=NBNUJO(NJ,1)
            NBNOJM=NBNUJO(NJ,2)
            NUSDE=NUSD(NJ*2-1)
            NUSDM=NUSD(NJ*2)
            IE=1
            IM=2
         ELSE IF (NUOBSD.EQ.NOBJE2) THEN
            NBNOJM=NBNUJO(NJ,1)
            NBNOJE=NBNUJO(NJ,2)
            NUSDM=NUSD(NJ*2-1)
            NUSDE=NUSD(NJ*2)
            IM=1
            IE=2
         ENDIF
C=========================================================================
         WRITE(IMPRIM,4000) NJ,NBNOJE,NBNOJM,NUSDE,NUSDM
C=========================================================================
         DO 4 K = 1 , NBNOJE
C           NUMERO LOCAL COTE ESCLAVE
            NULOES = NUNOJO(NJ,IE,K)
C           LE NUMERO GLOBAL COTE ESCLAVE
            NUGLES = NULOES + NDCAL(NUSDE)
C=========================================================================
            WRITE(IMPRIM,2000) K,NULOES,VG(NUGLES),U(NULOES)
C=========================================================================
 4       CONTINUE
 3    CONTINUE
C
      WRITE(IMPRIM,5000) (I,U(I),I=1,NTDL)
C
      RETURN
1000  FORMAT(1X/,80(1H-)/
     S       ' THE8 : JOINT',I5,' NBN MAITRES ',I5,
     S       ' NBN ESCLAVES ',I5,' SD MA/ES',2I5)
2000  FORMAT(' NOEUD GLOBAL ESCLAVE  ',I5,'(',I5,') VALEUR',2D15.6)
3000  FORMAT(' NOEUD GLOBAL MAITRE   ',I5,'(',I5,') VALEUR',2D15.6)
4000  FORMAT(1X/,80(1H-)/
     S       ' THE8 : JOINT',I5,' NBN ESCLAVES ',I5,
     S       ' NBN MAITRES ',I5,' SD ES/MA',2I5)
5000  FORMAT(1X/,80(1H=)/,100(4(I5,D15.6)/)/)
      END
