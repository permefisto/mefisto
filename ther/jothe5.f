      SUBROUTINE JOTHE5( NUOBSD, NBJOIN, MXNOJO, 
     S                   NBOBJE, NBNUJO, NUNOJO, NUSD,   NDCAL, 
     S                   VG,     VL,     JOIE,   NBJOIE, MNMAJO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CONTRIBUTION D'UN VECTEUR GLOBAL AU VECTEUR LOCAL
C -----  ( MULTIPLICATION PAR LA MATRICE Q )
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
C VG     : LE VECTEUR GLOBAL
C JOIM   : LISTE DES NBJOIM JOINTS MAITRES
C JOIE   : LISTE DES NBJOIE JOINTS ESCLAVES
C MNMAJO : ADRESSES DES MATRICES ASSOCIEES AUX JOINTS
C
C SORTIE :
C --------
C VL     : LE VECTEUR LOCAL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1994
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSIMTA / NOIMPR
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))
      DIMENSION         JOIE(NBJOIE), NUSD(1), NDCAL(1)
      DIMENSION         NBOBJE(NBJOIN,4), NBNUJO(NBJOIN,2),
     S                  NUNOJO(NBJOIN,2,MXNOJO)
      DOUBLE PRECISION  VG(1), VL(1)
C
      IM = 0
      IE = 0
      DO 1 J = 1 , NBJOIE
C
C        COTE ESCLAVE / COTE MAITRE
C        --------------------------
C        LE NUMERO DU JOINT ESCLAVE
         NJ = JOIE(J)
         NJ10 = (NJ-1)*10
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
            NBNOJE=NBNUJO(NJ,2)
            NBNOJM=NBNUJO(NJ,1)
            NUSDE=NUSD(NJ*2)
            NUSDM=NUSD(NJ*2-1)
            IE=2
            IM=1
         ENDIF
C=========================================================================
         IF (NOIMPR.GT.10)
     S   WRITE(IMPRIM,1000) NJ,NBNOJE,NBNOJM,NUSDE,NUSDM
C=========================================================================
C
C        DECLARATION DES TABLEAUX
C        ------------------------
         LONV = NBNOJM
         IF (NBNOJE.GT.LONV) LONV = NBNOJE
         MNVSS = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVSS )
         CALL AZEROD( LONV , MCN(MNVSS) )
         IU = ( MNVSS - 1 ) / 2
         MNVMA = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVMA )
         CALL AZEROD( LONV , MCN(MNVMA) )
         IV = ( MNVMA - 1 ) / 2
         MNVES = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVES )
         CALL AZEROD( LONV , MCN(MNVES) )
         IW = ( MNVES - 1 ) / 2
         MNVSM = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVSM )
         CALL AZEROD( LONV , MCN(MNVSM) )
         IZ = ( MNVSM - 1 ) / 2
C
C        INTERPOLATION   MAITRE > ESCLAVE
C        --------------------------------
         DO 3 K = 1 , NBNOJM
C           LE NUMERO LOCAL COTE MAITRE
            NULOMA = NUNOJO(NJ,IM,K)
C           LE NUMERO GLOBAL COTE MAITRE
            NUGLMA = NULOMA + NDCAL(NUSDM)
            DMCN(IV+K) = VG(NUGLMA)
C=========================================================================
            IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,3000) NULOMA,NUGLMA,VG(NUGLMA)
C=========================================================================
 3       CONTINUE
         MNQQ = MCN(MNMAJO+NJ10+8)
         MNNU = MCN(MNMAJO+NJ10+9)
         CALL SDRES13(MCN(MNQQ),MCN(MNVMA),MCN(MNNU),MCN(MNVES),NBNOJE)
C=========================================================================
         IF (NOIMPR.GT.10) THEN
         DO 5001 K=1,NBNOJM
            WRITE(IMPRIM,5000) K,DMCN(IV+K)
 5001    CONTINUE
         DO 5002 K=1,NBNOJE
            WRITE(IMPRIM,6000) K,DMCN(IW+K)
 5002    CONTINUE
         ENDIF
C=========================================================================
C
C        LE SECOND MEMBRE
C        ----------------
C        PREMIER NOEUD : NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,1)
C        NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
C        LE NUMERO LOCAL COTE ESCLAVE
         NULOES = NUNOJO(NJ,IE,1)
C        LE NUMERO GLOBAL COTE ESCLAVE
         NUGLES = NULOES + NDCAL(NUSDE)
         DMCN(IU+1) = VG(NUGLMA) - VG(NUGLES)
C=========================================================================
         IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,2000) NULOES,VG(NUGLES),NULOMA,VG(NUGLMA)
C=========================================================================
         DO 2 K = 2 , NBNOJE - 1
            DMCN(IU+K) = DMCN(IW+K)
 2       CONTINUE
C        DERNIER NOEUD : NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,NBNOJM)
C        NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
C        LE NUMERO LOCAL COTE ESCLAVE
         NULOES = NUNOJO(NJ,IE,NBNOJE)
C        LE NUMERO GLOBAL COTE ESCLAVE
         NUGLES = NULOES + NDCAL(NUSDE)
         DMCN(IU+NBNOJE) = VG(NUGLMA) - VG(NUGLES)
C=========================================================================
         IF (NOIMPR.GT.10)
     S       WRITE(IMPRIM,2000) NULOES,VG(NUGLES),NULOMA,VG(NUGLMA)
C=========================================================================
C
C        CALCUL DE LA VALEUR DU COTE ESCLAVE
C        -----------------------------------
         MNQ = MCN(MNMAJO+NJ10+5)
         CALL SDRES12(MCN(MNQ),MCN(MNVSS),MCN(MNVSM),NBNOJE)
         MOREE2 = MOTVAR(6)
         MNVMA2 = MNVMA + MOREE2
         MNVSM2 = MNVSM + MOREE2
         MNLU = MCN(MNMAJO+NJ10+6)
         CALL AZEROD( LONV , MCN(MNVMA) )
         CALL SDRES11(MCN(MNLU),MCN(MNVMA2),MCN(MNVSM2),NBNOJE-2)
C=========================================================================
         IF (NOIMPR.GT.10) THEN
         DO 5003 K=1,NBNOJE
            NULOES = NUNOJO(NJ,IE,K)
            WRITE(IMPRIM,7000) K,VL(NULOES)
 5003    CONTINUE
         ENDIF
C=========================================================================
C
C        TRANSFERT DES VALEURS
C        ---------------------
         DO 4 K = 2 , NBNOJE - 1
C           LE NUMERO LOCAL COTE ESCLAVE
            NULOES = NUNOJO(NJ,IE,K)
C           LE NUMERO GLOBAL COTE ESCLAVE
            NUGLES = NULOES + NDCAL(NUSDE)
            VL(NULOES) = DMCN(IV+K)
C=========================================================================
            IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,4000) NULOES,NUGLES,VL(NULOES)
C=========================================================================
 4       CONTINUE
C
C        DESTRUCTION DES TABLEAUX TEMPORAIRES
C        ------------------------------------
         CALL TNMCDS( 'REEL2' , LONV , MNVSM )
         CALL TNMCDS( 'REEL2' , LONV , MNVMA )
         CALL TNMCDS( 'REEL2' , LONV , MNVES )
         CALL TNMCDS( 'REEL2' , LONV , MNVSS )
C
 1    CONTINUE
C
      RETURN
1000  FORMAT(1X/,80(1H-)/
     S       ' THE5 : JOINT',I5,' NBN ESCLAVES',I5,
     S       ' NBN MAITRES ',I5,' SD ES/MA',2I5)
2000  FORMAT(' CORRECTION AU BORD         ',I5,D15.6,I7,D15.6)
3000  FORMAT(' NOEUD GLOBAL MAITRE        ',I5,'(',I5,') VALEUR',D15.6)
4000  FORMAT(' NOEUD LOCAL ESCLAVE MODIFIE',I5,'(',I5,') VALEUR',D15.6)
5000  FORMAT(' VALEUR MAITRE AVANT INTERPOLATION ',I4,D15.6)
6000  FORMAT(' VALEUR MAITRE APRES INTERPOLATION ',I4,D15.6)
7000  FORMAT(' VALEUR ESCLAVE AVANT CORRECTION   ',I4,D15.6)
      END
