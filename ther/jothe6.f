      SUBROUTINE JOTHE6( NUOBSD, NBJOIN, MXNOJO,
     S                   NBOBJE, NBNUJO, NUNOJO, NUSD, NDCAL,
     S                   VG,     JOIM,   NBJOIM, MNMAJO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CONTRIBUTION D'UN VECTEUR LOCAL AU VECTEUR GLOBAL
C ----- ( MULTIPLICATION PAR LA MATRICE Q TRANSPOSEE )
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
C JOIM   : LISTE DES NBJOIM JOINTS MAITRES
C MNMAJO : ADRESSES DES MATRICES ASSOCIEES AUX JOINTS
C
C SORTIE :
C --------
C VG     : LE VECTEUR GLOBAL
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
      DIMENSION JOIM(NBJOIM),NUSD(1),NDCAL(1)
      DIMENSION NBOBJE(NBJOIN,4),NBNUJO(NBJOIN,2),
     S          NUNOJO(NBJOIN,2,MXNOJO)
      DOUBLE PRECISION VG(1)
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
C        DECLARATION DES TABLEAUX
C        ------------------------
         LONV = NBNOJM
         IF (NBNOJE.GT.LONV) LONV = NBNOJE
         MNVSS = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVSS )
         CALL AZEROD( LONV , MCN(MNVSS) )
         IU = ( MNVSS - 1 ) / 2
         MNVES = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVES )
         CALL AZEROD( LONV , MCN(MNVES) )
         IV = ( MNVES - 1 ) / 2
         MNVMA = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVMA )
         CALL AZEROD( LONV , MCN(MNVMA) )
         IW = ( MNVMA - 1 ) / 2
         MNVSM = 0
         CALL TNMCDC( 'REEL2' , LONV , MNVSM )
         CALL AZEROD( LONV , MCN(MNVSM) )
         IZ = ( MNVSM - 1 ) / 2
C
C        INTERPOLATION  ESCLAVE > MAITRE
C        -------------------------------
         DO 3 K = 1 , NBNOJE
C           LE NUMERO LOCAL COTE ESCLAVE
            NULOES = NUNOJO(NJ,IE,K)
C           LE NUMERO GLOBAL COTE ESCLAVE
            NUGLES = NULOES + NDCAL(NUSDE)
            DMCN(IV+K) =  VG(NUGLES)
C=========================================================================
         IF (NOIMPR.GT.10)
     S       WRITE(IMPRIM,3000) NUGLES,NULOES,VG(NUGLES)
C=========================================================================
 3       CONTINUE
         MNQQ = MCN(MNMAJO+NJ10+3)
         MNNU = MCN(MNMAJO+NJ10+4)
         CALL SDRES13(MCN(MNQQ),MCN(MNVES),MCN(MNNU),MCN(MNVMA),NBNOJM)
C=========================================================================
         IF (NOIMPR.GT.10) THEN
         DO 5001 K=1,NBNOJE
            WRITE(IMPRIM,5000) K,DMCN(IV+K)
 5001    CONTINUE
         DO 5002 K=1,NBNOJM
            WRITE(IMPRIM,6000) K,DMCN(IW+K)
 5002    CONTINUE
         ENDIF
C=========================================================================
C
C        LE SECOND MEMBRE
C        ----------------
C        PREMIER NOEUD : NUMERO LOCAL COTE ESCLAVE
         NULOES = NUNOJO(NJ,IE,1)
C        LE NUMERO GLOBAL COTE ESCLAVE
         NUGLES = NULOES + NDCAL(NUSDE)
C        LE NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,1)
C        LE NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
         DMCN(IU+1) = VG(NUGLMA) - VG(NUGLES)
C=========================================================================
         IF (NOIMPR.GT.10)
     S       WRITE(IMPRIM,2000) NULOES,VG(NUGLES),NULOMA,VG(NUGLMA)
C=========================================================================
         DO 2 K = 2 , NBNOJM - 1
            DMCN(IU+K) = DMCN(IW+K)
 2       CONTINUE
C        DERNIER NOEUD : NUMERO LOCAL COTE ESCLAVE
         NULOES = NUNOJO(NJ,IE,NBNOJE)
C        LE NUMERO GLOBAL COTE ESCLAVE
         NUGLES = NULOES + NDCAL(NUSDE)
C        LE NUMERO LOCAL COTE MAITRE
         NULOMA = NUNOJO(NJ,IM,NBNOJM)
C        LE NUMERO GLOBAL COTE MAITRE
         NUGLMA = NULOMA + NDCAL(NUSDM)
C        LA VALEUR DU COTE ESCLAVE
         DMCN(IU+NBNOJM) = VG(NUGLES) - VG(NUGLMA)
C=========================================================================
         IF (NOIMPR.GT.10)
     S      WRITE(IMPRIM,2000) NULOES,VG(NUGLES),NULOMA,VG(NUGLMA)
C=========================================================================
C
C        CALCUL DE LA VALEUR DU COTE MAITRE
C        ----------------------------------
         MNQ = MCN(MNMAJO+NJ10)
         CALL SDRES12(MCN(MNQ),MCN(MNVSS),MCN(MNVSM),NBNOJM)
         MOREE2 = MOTVAR(6)
         MNVES2 = MNVES + MOREE2
         MNVSM2 = MNVSM + MOREE2
         MNLU = MCN(MNMAJO+NJ10+1)
         CALL AZEROD( LONV , MCN(MNVES) )
         CALL SDRES11(MCN(MNLU),MCN(MNVES2),MCN(MNVSM2),NBNOJM-2)
C=========================================================================
         IF (NOIMPR.GT.10) THEN
         DO 5003 K=1,NBNOJM
            NUGLMA = NUNOJO(NJ,IM,K) + NDCAL(NUSDM)
            WRITE(IMPRIM,7000) K,VG(NUGLMA)
 5003    CONTINUE
         ENDIF
C=========================================================================
C
C        TRANSFERT DES VALEURS
C        ---------------------
         DO 4 K = 2 , NBNOJM - 1
C           LE NUMERO LOCAL COTE MAITRE
            NULOMA = NUNOJO(NJ,IM,K)
C           LE NUMERO GLOBAL COTE MAITRE
            NUGLMA = NULOMA + NDCAL(NUSDM)
C           1) CUMUL DES VALEURS : COTE MAITRE + COTE ESCLAVE
C           POUR LA CONTRIBUTION TOTALE DU JOINT AU SYSTEME GLOBAL
C           LES INCONNUES DU JOINT SONT PORTEES PAR LE COTE MAITRE
            VG(NUGLMA) =  VG(NUGLMA) + DMCN(IV+K)
C=========================================================================
            IF (NOIMPR.GT.10)
     S         WRITE(IMPRIM,4000) NUGLMA,NULOMA,VG(NUGLMA)
C=========================================================================
 4       CONTINUE
C
C        CORRECTION GLOBALE
C        ------------------
         DO 5 K = 2 , NBNOJE - 1
C           LE NUMERO GLOBAL COTE ESCLAVE
            NUGLES = NUNOJO(NJ,IE,K) + NDCAL(NUSDE)
C           2) SUPPRESSION DES VALEURS ESCLAVES
C           LES INCONNUES COTE ESCLAVE N'APPARAISSENT PAS DANS
C           LE SYSTEME GLOBAL ( VOIR 1 )
            VG(NUGLES) = 0.D0
 5       CONTINUE
C
C        DESTRUCTION DES TABLEAUX TEMPORAIRES
C        ------------------------------------
         CALL TNMCDS( 'REEL2'  , LONV , MNVSM )
         CALL TNMCDS( 'REEL2'  , LONV , MNVES )
         CALL TNMCDS( 'REEL2'  , LONV , MNVMA )
         CALL TNMCDS( 'REEL2'  , LONV , MNVSS )
C
 1    CONTINUE
C
      RETURN
1000  FORMAT(1X/,80(1H-)/
     S       ' THE6 : JOINT',I5,' NBN MAITRES ',I5,
     S       ' NBN ESCLAVES',I5,' SD MA/ES',2I5)
2000  FORMAT(' CORRECTION AU BORD         ',I5,D15.6,I7,D15.6)
3000  FORMAT(' NOEUD GLOBAL ESCLAVE       ',I5,'(',I5,') VALEUR',D15.6)
4000  FORMAT(' NOEUD GLOBAL MAITRE MODIFIE',I5,'(',I5,') VALEUR',D15.6)
5000  FORMAT(' VALEUR ESCLAVE AVANT INTERPOLATION ',I4,D15.6)
6000  FORMAT(' VALEUR ESCLAVE APRES INTERPOLATION ',I4,D15.6)
7000  FORMAT(' VALEUR MAITRE  AVANT CORRECTION    ',I4,D15.6)
      END
