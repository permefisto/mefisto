      SUBROUTINE RESBIJ(NBS1,NBS2,FCI,FCJ,NBELTR,NBPEXT,ITRAV)
C***********************************************************************
C BUT :    RESOLUTION DE LA BIJECTIVITE
C***********************************************************************
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           NBELTR : NOMBRE DE MAILLES DEGENEREES
C           NBPEXT : NOMBRE DE POINTS EXTERIEURS AU MAILLAGE
C           ITRAV  : VECTEUR DE TRAVAIL ENTIERS CONTENAT LES NUMEROS
C                    DES MAILLES DEGENREES ET DES POINTS EXTERIEURS
C
C SORTIE:
C           FCI    : VECTEUR DES FONCTIONS DE CONTROLE EN I
C           FCJ    : VECTEUR DES FONCTIONS DE CONTROLE EN J
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DIMENSION FCI(NBS1*NBS2),FCJ(NBS1*NBS2)
      DIMENSION ITRAV(NBS1*NBS2)
      INTEGER CHAINE(1000),CHAINL(24)
      WRITE (IMPRIM,*) ' '
      WRITE (IMPRIM,1010)
 1010 FORMAT('CORRECTION DE LA NON BIJECTIVITE DU MAILLAGE')
CCC      WRITE (IMPRIM,1000)
CCC      WRITE (IMPRIM,1010)
CCC      WRITE (IMPRIM,1020)
CCC 1000 FORMAT('**************************************')
CCC 1010 FORMAT('*  CORRECTION DE LA NON BIJECTIVITE  *')
CCC 1020 FORMAT('*                                    *')
CCC 1030 FORMAT('*   A PARTIR DES NOEUDS EXTERIEURS   *')
CCC 1040 FORMAT('*  A PARTIR DES MAILLES DEGENEREES   *')
C***********************************************************************
C LES POINTS ET LES MAILLES DEGENERES TROUVES, ON VA TRAITER LES ZONES
C***********************************************************************
      KTEST = 0
      IF (NBPEXT.NE.0) THEN
CCC        WRITE (IMPRIM,1030)
CCC        WRITE (IMPRIM,1000)
C***********************************************************************
C     ON EST DANS LE CAS OU IL Y A DES POINTS EXTERIEURS AU MAILLAGE
C***********************************************************************
C          RECHERCHE DE ZONES DE POINTS EXTERIEURS VOISINS
C***********************************************************************
C                 CHOIX D UN PREMIER POINT EXTERIEUR
C             ET INITIALISATION DE LA ZONE DE RECHERCHE
C-----------------------------------------------------------------------
        NDEB  = NBELTR+1
        NBRE  = 1
  200   CONTINUE
        NEU  = ITRAV(NDEB)
        J    = 1+(NEU-1)/NBS1
        I    = NEU-NBS1*(J-1)
        IINF = I-1
        ISUP = I+1
        JINF = J-1
        JSUP = J+1
  210   CONTINUE
C-----------------------------------------------------------------------
C      RECHERCHE DE LA ZONE A LAQUELLE APPARTIENT LE PREMIER POINT
C-----------------------------------------------------------------------
        ITEST = 0
        DO 220 KK=NDEB+1,NBELTR+NBPEXT
          NEU = ITRAV(KK)
          IF (NEU.NE.0) THEN
            J = 1+(NEU-1)/NBS1
            I = NEU-NBS1*(J-1)
            IF ((I.GE.IINF).AND.(I.LE.ISUP).AND.
     S          (J.GE.JINF).AND.(J.LE.JSUP)) THEN
              IINF      = MIN(IINF,I-1)
              ISUP      = MAX(ISUP,I+1)
              JINF      = MIN(JINF,J-1)
              JSUP      = MAX(JSUP,J+1)
              ITRAV(KK) = 0
              ITEST     = 1
              NBRE      = NBRE+1
            ENDIF
          ENDIF
  220   CONTINUE
        IF ((ITEST.EQ.1).AND.(NBRE.LT.NBPEXT)) THEN
          GOTO 210
        ENDIF
C***********************************************************************
C LA ZONE TROUVEE ON CHERCHE OU LES FONCT. DE CONT. DOIVENT ETRE PLACEES
C***********************************************************************
        IDINF = IINF-1
        IDSUP = NBS1-ISUP
        JDINF = JINF-1
        JDSUP = NBS2-JSUP
        ID = MIN(IDINF,IDSUP)
        JD = MIN(JDINF,JDSUP)
        ID = MIN(ID,JD)
        IF (ID.EQ.0) THEN
C-----------------------------------------------------------------------
C             LA ZONE A TRAITER TOUCHE AU MOINS UN BORD
C        RECHERCHE DE DES COORDONNEES TOPOLOGIQUES DE CE BORD
C-----------------------------------------------------------------------
          KTEST = 1
          IF (IDINF.EQ.ID) THEN
            I0 = IINF
          ELSE
            IF (IDSUP.EQ.ID) THEN
              I0 = ISUP
            ELSE
              I0 = 0
            ENDIF
          ENDIF
          IF (JDINF.EQ.ID) THEN
            J0 = JINF
          ELSE
            IF (JDSUP.EQ.ID) THEN
              J0 = JSUP
            ELSE
              J0 = 0
            ENDIF
          ENDIF
          IF ((I0*J0).EQ.0) THEN
C-----------------------------------------------------------------------
C     LA ZONE A TRAITER TOUCHE UN BORD MAIS NE TOUCHE PAS DE COIN
C
C        |   |   |   |   |   |
C       -+---+---X---+---+---+--               * NOEUD DU BORD
C        |   |   |   |   |   |
C       -+---+---X---X---X---+--              + NOEUD INTERIEUR
C        |   |   |   |   |   |
C       -*---*---*---*---*---*--> I          X NOEUD EXTERIEUR
C                ^   ^   ^
C                |   |   |
C               FCJ FCJ FCJ
C-----------------------------------------------------------------------
C  CALCUL DES FONCTIONS DE CONTROLE EN FONCTION DE L ETENDUE DE LA ZONE
C-----------------------------------------------------------------------
            IF (I0.NE.0) THEN
              LZONE = ISUP-IINF-1
              BI  = LOG(10.)/(LZONE+0.5)
              AIB = 0.3*(LZONE+3.)/4.
              AI  = AIB/EXP(-BI)
              DO 230 JJ=JINF+1,JSUP-1
                ID = 0
                CALL CALFCI(ID,I0,JJ,AI,BI,NBS1,NBS2,FCI)
  230         CONTINUE
            ELSE
              LZONE = JSUP-JINF-1
              BI  = LOG(10.)/(LZONE+0.5)
              AIB = 0.3*(LZONE+3.)/4.
              AI  = AIB/EXP(-BI)
              DO 240 II=IINF+1,ISUP-1
                ID = 0
                CALL CALFCJ(ID,II,J0,AI,BI,NBS1,NBS2,FCJ)
  240         CONTINUE
            ENDIF
          ELSE
C-----------------------------------------------------------------------
C                  LA ZONE A TRAITER TOUCHE UN COIN
C
C            J
C            ^   |   |   |
C  FCI -->   *---X---+---+--               * NOEUD DU BORD
C            |   |   |   |
C  FCI -->   *---X---X---X--               + NOEUD INTERIEUR
C            |   |   |   |
C  FCI -->   *---*---*---*---> I           X NOEUD EXTERIEUR
C            ^   ^   ^   ^
C            |   |   |   |
C           FCJ FCJ FCJ FCJ
C-----------------------------------------------------------------------
C  CALCUL DES FONCTIONS DE CONTROLE EN FONCTION DE L ETENDUE DE LA ZONE
C-----------------------------------------------------------------------
            N = INT( I-NBS1/2. )
            IF (N.GT.0.) THEN
              IDIR = 1
            ELSE
              IDIR = -1
            ENDIF
            N = INT( J-NBS2/2. )
            IF (N.GT.0.) THEN
              JDIR = 1
            ELSE
              JDIR = -1
            ENDIF
            ILONG = ISUP-IINF-1
            JLONG = JSUP-JINF-1
            LZONE = MAX(ILONG,JLONG)
            AIB   = 0.5*(LZONE+3.)/4.
            BI    = LOG(10.)/LZONE
            AI    = AIB/EXP(-BI)
            ID = 0
            CALL CALFCI(ID,I0,J0,AI,BI,NBS1,NBS2,FCI)
            CALL CALFCJ(ID,I0,J0,AI,BI,NBS1,NBS2,FCJ)
            DO 250 KK=1,LZONE
              BI = LOG(10.)/(KK+0.5)
              AI = AIB/2./EXP(-BI)/(KK+1)
              I1 = I0-KK*IDIR
              CALL CALFCJ(ID,I1,J0,AI,BI,NBS1,NBS2,FCJ)
              J1 = J0-KK*JDIR
              CALL CALFCI(ID,I0,J1,AI,BI,NBS1,NBS2,FCI)
  250       CONTINUE
          ENDIF
        ENDIF
C***********************************************************************
C                 LA ZONE PRECEDENTE AYANT ETE TRAITEE
C             RECHERCHE SI NECESSAIRE D UNE NOUVELLE ZONE
C***********************************************************************
        IF (NBRE.NE.NBPEXT) THEN
          DO 260 II=NDEB+1,NBELTR+NBPEXT
            IF (ITRAV(II).NE.0) THEN
              NDEB = II
              NBRE = NBRE+1
              GOTO 200
            ENDIF
  260     CONTINUE
        ENDIF
        IF (KTEST.EQ.1) THEN
          RETURN
        ENDIF
      ENDIF
C***********************************************************************
C        ON EST DANS LE CAS OU IL N Y A PAS DE POINTS EXTERIEURS
C             ON VA DONC TRAITER LES MAILLES DEGENEREES
C***********************************************************************
C              RECHERCHE DE ZONES DE MAILLES RETOURNEES
C***********************************************************************
C                     CHOIX D UNE PREMIERE MAILLE
C              ET INITIALISATION DE LA ZONE DE RECHERCHE
C-----------------------------------------------------------------------
CCC      WRITE (IMPRIM,1040)
CCC      WRITE (IMPRIM,1000)
      NELDEB = 1
      NBRE   = 0
   70 CONTINUE
      NUMELT = ITRAV(NELDEB)
      IF (NUMELT.NE.0) THEN
        J0 = 1+(NUMELT-1)/(NBS1-1)
        I0 = NUMELT-(NBS1-1)*(J0-1)
        IINF = I0-1
        ISUP = I0+2
        JINF = J0-1
        JSUP = J0+2
        I1 = MAX(1,IINF)
        I2 = MIN(NBS1,ISUP)
        J1 = MAX(1,JINF)
        J2 = MIN(NBS2,JSUP)
        NUM  = 0
        DO 71 I=I1,I2
          NUM         = NUM+1
          NEU         = (J1-1)*NBS1+I
          CHAINE(NUM) = NEU
   71   CONTINUE
        DO 72 J=J1+1,J2
          NUM         = NUM+1
          NEU         = (J-1)*NBS1+I2
          CHAINE(NUM) = NEU
   72   CONTINUE
        DO 73 I=I2-1,I1,-1
          NUM         = NUM+1
          NEU         = (J2-1)*NBS1+I
          CHAINE(NUM) = NEU
   73   CONTINUE
        DO 74 J=J2-1,J1+1,-1
          NUM         = NUM+1
          NEU         = (J-1)*NBS1+I1
          CHAINE(NUM) = NEU
   74   CONTINUE
        LONGCH = NUM
        NBRE   = NBRE+1
      ELSE
        NELDEB = NELDEB+1
        IF (NBRE.LE.NBELTR) THEN
          GOTO 70
        ELSE
          GOTO 150
        ENDIF
      ENDIF
      ITRAV(NELDEB) = -ITRAV(NELDEB)
   80 CONTINUE
C-----------------------------------------------------------------------
C    RECHERCHE DE LA ZONE A LAQUELLE APPARTIENT LA PREMIERE MAILLE
C-----------------------------------------------------------------------
      ITEST = 0
      DO 90 KK=NELDEB+1,NBELTR
        NUMELT = ITRAV(KK)
        IF (NUMELT.GT.0) THEN
          J0 = 1+(NUMELT-1)/(NBS1-1)
          I0 = NUMELT-(NBS1-1)*(J0-1)
          KTEST = 0
          DO 500 I=NELDEB,NBELTR
            IF ((ABS(NUMELT+ITRAV(I)).EQ.1).OR.
     S          ((ABS(ABS(NUMELT+ITRAV(I))-(NBS1-1))).LE.1)) THEN
              KTEST = 1
              ITRAV(KK) = -ITRAV(KK)
              GOTO 501
            ENDIF
  500     CONTINUE
  501     CONTINUE
          IF (KTEST.EQ.1) THEN
            NUM = 0
            I1 = MAX(1,I0-1)
            I2 = MIN(NBS1,I0+2)
            J1 = MAX(1,J0-1)
            J2 = MIN(NBS2,J0+2)
            DO 81 I=I1,I2
              NUM = NUM+1
              NEU = (J1-1)*NBS1+I
              CHAINL(NUM) = NEU
   81       CONTINUE
            DO 82 J=J1+1,J2
              NUM = NUM+1
              NEU = (J-1)*NBS1+I2
              CHAINL(NUM) = NEU
   82       CONTINUE
            DO 83 I=I2-1,I1,-1
              NUM = NUM+1
              NEU = (J2-1)*NBS1+I
              CHAINL(NUM) = NEU
   83       CONTINUE
            DO 84 J=J2-1,J1+1,-1
              NUM = NUM+1
              NEU = (J-1)*NBS1+I1
              CHAINL(NUM) = NEU
   84       CONTINUE
            LONGLO = NUM
C
            NDEBLO = 0
            NDEBCH = 0
            NFINLO = 0
            NFINCH = 0
            DO 86 N=1,LONGCH
              NEUREF = CHAINE(N)
              DO 85 I=1,LONGLO
                IF (CHAINL(I).EQ.NEUREF) THEN
                  IF (NFINCH.EQ.0) THEN
                    NDEBLO = I
                    NDEBCH = N
                  ELSE
                    NFINLO = I
                    NFINCH = N
                    GOTO 86
                  ENDIF
                  GOTO 86
                ENDIF
   85         CONTINUE
              IF (NDEBCH.NE.0) THEN
                NFINCH = N
              ENDIF
   86       CONTINUE
            IF (NDEBLO.GT.NFINLO) THEN
              NFINLO = NFINLO+LONGLO
              DO 88 I=1,LONGLO
                CHAINL(I+LONGLO) = CHAINL(I)
   88         CONTINUE
            ENDIF
            NRAJOU = (NFINLO-NDEBLO)-(NFINCH-NDEBCH)
            LONGCH = LONGCH+NRAJOU
            IF (NRAJOU.GT.0) THEN
              DO 91 I=LONGCH,NFINCH+NRAJOU,-1
                CHAINE(I) = CHAINE(I-NRAJOU)
   91         CONTINUE
            ENDIF
            IF (NRAJOU.LT.0) THEN
              DO 92 I=NFINCH+NRAJOU,LONGCH
                CHAINE(I) = CHAINE(I-NRAJOU)
   92         CONTINUE
            ENDIF
            DO 93 I=1,NFINLO-NDEBLO-1
              CHAINE(NDEBCH+I) = CHAINL(NDEBLO+I)
   93       CONTINUE
            ITEST     = 1
            NBRE      = NBRE+1
CCC            WRITE(IMPRIM,*) ' NBELTR NBRE ',NBELTR,NBRE
          ENDIF
        ENDIF
   90 CONTINUE
      IF ((ITEST.EQ.1).AND.(NBRE.NE.NBELTR)) THEN
        GOTO 80
      ENDIF
C***********************************************************************
C LA ZONE TROUVEE ON CHERCHE OU LES FONCT. DE CONT. DOIVENT ETRE PLACEES
C***********************************************************************
      DO 502 I=NELDEB,NBELTR
        IF (ITRAV(I).LT.0) THEN
          ITRAV(I) = 0
        ENDIF
  502 CONTINUE
C-----------------------------------------------------------------------
C  CALCUL DES FONCTIONS DE CONTROLE EN FONCTION DE L ETENDUE DE LA ZONE
C-----------------------------------------------------------------------
      LZONEI = 5
      BI  = LOG(10.)/(LZONEI+0.5)
      AIB = 0.1*(LZONEI+3.)/4.
      AI  = -AIB/EXP(-BI)/8
      CHAINE(LONGCH+1) = CHAINE(1)
      DO 100 N=1,LONGCH
        NEU1 = CHAINE(N)
        J0 = 1+(NEU1-1)/NBS1
        I0 = NEU1-NBS1*(J0-1)
        NEU2 = CHAINE(N+1)
        NDIF = NEU2-NEU1
        IF (NDIF.EQ.1) THEN
          ID = 1
          J1 = MAX(1,J0-3)
          IF (J1.NE.1) THEN
            CALL CALFCJ(ID,I0,J1,AI,BI,NBS1,NBS2,FCJ)
            CALL CALFCJ(ID,I0+1,J1,AI,BI,NBS1,NBS2,FCJ)
          ENDIF
        ENDIF
        IF (NDIF.EQ.-1) THEN
          ID = -1
          J1 = MIN(NBS2,J0+3)
          IF (J1.NE.NBS2) THEN
            CALL CALFCJ(ID,I0,J1,AI,BI,NBS1,NBS2,FCJ)
            CALL CALFCJ(ID,I0-1,J1,AI,BI,NBS1,NBS2,FCJ)
          ENDIF
        ENDIF
        IF (NDIF.EQ.NBS1) THEN
          JD = -1
          I1 = MIN(NBS1,I0+3)
          IF (I1.NE.NBS1) THEN
            CALL CALFCI(JD,I1,J0,AI,BI,NBS1,NBS2,FCJ)
            CALL CALFCI(JD,I1,J0+1,AI,BI,NBS1,NBS2,FCJ)
          ENDIF
        ENDIF
        IF (NDIF.EQ.-NBS1) THEN
          JD = 1
          I1 = MAX(1,I0-3)
          IF (I1.NE.1) THEN
            CALL CALFCI(JD,I1,J0,AI,BI,NBS1,NBS2,FCJ)
            CALL CALFCI(JD,I1,J0-1,AI,BI,NBS1,NBS2,FCJ)
          ENDIF
        ENDIF
  100 CONTINUE
C***********************************************************************
C                 LA ZONE PRECEDENTE AYANT ETE TRAITEE
C             RECHERCHE SI NECESSAIRE D UNE NOUVELLE ZONE
C***********************************************************************
      IF (NBRE.LT.NBELTR) THEN
        DO 140 II=NELDEB+1,NBELTR
          IF (ITRAV(II).NE.0) THEN
            NELDEB = II
            GOTO 70
          ENDIF
  140   CONTINUE
      ENDIF
  150 CONTINUE
      RETURN
      END
