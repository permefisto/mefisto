      SUBROUTINE JOTHE4( NUSD , NBJOIN , MXNOJO ,
     S                   LINUOB , LINUNO , NBNOEU , COOR ,
     S                   JOIM , NBJOIM , JOIE , NBJOIE , MNMAJO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER ET FACTORISER LA MATRICE Q ASSOCIEE A
C ---    CHAQUE COTE DE CHACUN DES JOINTS
C
C ENTREES :
C ---------
C NUSD   : NUMERO D'OBJET DU SOUS-DOMAINE TRAITE
C NBJOIN : NOMBRE DE JOINTS
C LINUOB : LES NUMEROS DES OBJETS ASSOCIES AUX LIGNES
C LINUNO : POUR CHAQUE NOEUD DE CHAQUE JOINT,
C          LE NUMERO DU NOEUD DANS L'OBJET ASSOCIE
C NBNOEU : NOMBRE DE NOEUDS DE CHAQUE JOINT
C COOR   : LES COORDONNEES DE CES NOEUDS
C JOIM   : LISTE DES NBJOIM JOINTS MAITRES
C JOIE   : LISTE DES NBJOIE JOINTS ESCLAVES
C MNMAJO : ADRESSE MCN DES TABLEAUX Q ET LU ASSOCIES AUX JOINTS
C
C SORTIES :
C ---------
C QM     : LA MATRICE ASSOCIEE AU JOINT COTE MAITRE
C LUM    : SA FACTORISEE              (COTE MAITRE)
C SM     : LE SECOND MEMBRE           (COTE MAITRE)
C QE     : LA MATRICE ASSOCIEE AU JOINT (COTE ESCLAVE)
C LUE    : SA FACTORISEE                (COTE ESCLAVE)
C SE     : LE SECOND MEMBRE             (COTE ESCLAVE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / MSIMTA / NOIMPR
      INTEGER JOIM(NBJOIM),JOIE(NBJOIE)
      INTEGER  NBNOEU(NBJOIN,2),LINUOB(NBJOIN,4)
      INTEGER  LINUNO(NBJOIN,2,MXNOJO)
      REAL  COOR(NBJOIN,2,MXNOJO,3),XYZ(3),XYZE(3)
      EQUIVALENCE       (MCN(1),RMCN(1),DMCN(1))

C     LES JOINTS MAITRES DU DOMAINE
C     -----------------------------
      NBK = 0
      NBJ = 0
      IJK = 0
      IKJ = 0

      DO 1 J = 1 , NBJOIM
C        LE NUMERO DU JOINT MAITRE
         NJ = JOIM(J)
C====================================================
         IF (NOIMPR.GT.10) WRITE(IMPRIM,9001) NJ
 9001    FORMAT(/,' JOINT MAITRE',I6/,1X,12(1H-))
C====================================================
         NJ10 = (NJ-1)*10
C        LE NOMBRE DE NOEUDS DU JOINT
         NBJ1=NBNOEU(NJ,1)
         NBJ2=NBNOEU(NJ,2)
C        LES SOUS-DOMAINES LIES AU JOINT
         NUSD1 = LINUOB(NJ,2)
         NUSD2 = LINUOB(NJ,4)
C        LE NUMERO DU SOUS-DOMAINE
         IF ( NUSD .EQ. NUSD1 ) THEN
           IJK = 1
           IKJ = 2
           NBJ = NBJ1
           NBK = NBJ2
         ELSE IF ( NUSD .EQ. NUSD2 ) THEN
           IJK = 2
           IKJ = 1
           NBJ = NBJ2
           NBK = NBJ1
         ENDIF
         MNQM = MCN(MNMAJO+NJ10)
         IAQM = ( MNQM - 1 ) / 2 + 1
         IAQ  = IAQM
         MNSM = MCN(MNMAJO+NJ10+2)
         IASM = ( MNSM - 1 ) / 2 + 1
         IAS  = IASM
         IAQ1 = IAQ + 3
         IAS1 = IAS
C        LE PREMIER NOEUD DU JOINT
         X1 = COOR(NJ,IJK,1,1)
         Y1 = COOR(NJ,IJK,1,2)
         X2 = COOR(NJ,IJK,2,1)
         Y2 = COOR(NJ,IJK,2,2)
         XY = X2 - X1 + Y2 - Y1
         XY1 = XY
C====================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,9002) LINUNO(NJ,IJK,1)
         WRITE(IMPRIM,9003) X1,Y1
 9002    FORMAT(' NOEUD',I6)
 9003    FORMAT(' COORDONNEES',2G15.6)
         ENDIF
C====================================================
         DMCN(IAQ)   = 0.D0
         DMCN(IAQ+1) = 1.D0
         DMCN(IAQ+2) = 0.D0
C*       DMCN(IAS)   =  XY / 2.D0
         DMCN(IAS)   =  - XY / 2.D0
         DO 2 NBJA = 2 , NBJ - 1
            X1 = COOR(NJ,IJK,NBJA,1)
            Y1 = COOR(NJ,IJK,NBJA,2)
            X2 = COOR(NJ,IJK,NBJA+1,1)
            Y2 = COOR(NJ,IJK,NBJA+1,2)
            XYA = XY
            XY = X2 - X1 + Y2 - Y1
C====================================================
            IF (NOIMPR.GT.10)  THEN
            WRITE(IMPRIM,9002) LINUNO(NJ,IJK,NBJA)
            WRITE(IMPRIM,9003) X1,Y1
            ENDIF
C====================================================
            IAQ = IAQ + 3
            IAS = IAS + 1
            DMCN(IAQ) = DMCN(IAQ-1)
            DMCN(IAQ+1) = ( XYA + XY ) / 3.D0
C*          DMCN(IAQ+1) = 2.D0 * XY / 3.D0
            DMCN(IAQ+2) = XY / 6.D0
            DMCN(IAS) = 0.D0
 2       CONTINUE
C        LE DERNIER NOEUD DU JOINT
         X1 = COOR(NJ,IJK,NBJ-1,1)
         Y1 = COOR(NJ,IJK,NBJ-1,2)
         X2 = COOR(NJ,IJK,NBJ,1)
         Y2 = COOR(NJ,IJK,NBJ,2)
         XY = X2 - X1 + Y2 - Y1
         XYN = XY
C====================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,9002) LINUNO(NJ,IJK,NBJ)
         WRITE(IMPRIM,9003) X2,Y2
         ENDIF
C====================================================
         IAQ2 = IAQ
         IAS2 = IAS + 1
         IAQ = IAQ + 3
         IAS = IAS + 1
         DMCN(IAQ)   = 0.D0
         DMCN(IAQ+1) = 1.D0
         DMCN(IAQ+2) = 0.D0
C*       DMCN(IAS)   = XY / 2.D0
         DMCN(IAS)   = - XY / 2.D0
C        CORRECTION DE LA DEUXIEME LIGNE DE Q
C*       DMCN(IAQ1)   = DMCN(IAS1)
         DMCN(IAQ1)   = 0.D0
C*       DMCN(IAQ1+1) = DMCN(IAQ1+1) + DMCN(IAS1)
         DMCN(IAQ1+1) = DMCN(IAQ1+1) + XY1 / 6.D0
C        CORRECTION DE L'AVANT DERNIERE LIGNE DE Q
C*       DMCN(IAQ2+1) = DMCN(IAQ2+1) + DMCN(IAS2)
         DMCN(IAQ2+1) = DMCN(IAQ2+1) + XYN / 6.D0
C*       DMCN(IAQ2+2) = DMCN(IAS2)
         DMCN(IAQ2+2) = 0.D0
C        FACTORISATION LU DE LA MATRICE QM
         NBLU = NBJ-2
         IAQM = ( MNQM - 1 ) / 2 + 1
         IALU = ( MCN(MNMAJO+NJ10+1) - 1 ) / 2 + 1
         IF (NBLU.NE.0) CALL SDRES10(DMCN(IAQM+3),DMCN(IALU),NBLU)
C==================================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,4425) ((DMCN(IAQM+JJ*3+K),K=0,2),
     S                        DMCN(IASM+JJ),JJ=0,NBJ-1)
         WRITE(IMPRIM,4426) ((DMCN(IALU+JJ*3+K),K=0,2),JJ=0,NBJ-3)
 4425    FORMAT(' MATRICE QM  ',3D15.6,' SM',D15.6)
 4426    FORMAT(' MATRICE LUM ',3D15.6)
         ENDIF
C==================================================================
C        LA MATRICE QQ D'INTERPOLATION
         MNQQM = MCN(MNMAJO+NJ10+3)
         IAQQ = ( MNQQM - 1 ) / 2 + 1
         IAQ = IAQQ
         MNNUM = MCN(MNMAJO+NJ10+4)
         MNN = MNNUM
         DO 11 NBJM = 1 , NBJ
            XM = COOR(NJ,IJK,NBJM,1)
            YM = COOR(NJ,IJK,NBJM,2)
            XYZE(1) = XM
            XYZE(2) = YM
            XYZE(3) = 0.D0
            DO 12 NBJA = 1  , NBK
               X1 = COOR(NJ,IKJ,NBJA,1)
               Y1 = COOR(NJ,IKJ,NBJA,2)
               XYZ(1) = X1
               XYZ(2) = Y1
               XYZ(3) = 0.D0
               CALL XYZIDE(XYZE,XYZ,IDENT)
               IF (IDENT.EQ.1) THEN
C              LES POINTS SONT CONFONDUS !
                  DMCN(IAQ)   = 1.D0
                  DMCN(IAQ+1) = 0.D0
                  IAQ = IAQ + 2
                  MCN(MNN)   = NBJA
                  MCN(MNN+1) = NBJA
                  MNN = MNN + 2
C==================================================================
                 IF (NOIMPR.GT.10)  THEN
                 WRITE(IMPRIM,9004) LINUNO(NJ,IJK,NBJM),XM,YM,
     S                              LINUNO(NJ,IKJ,NBJA),X1,Y1
 9004            FORMAT(' NOEUD',I6/,' COORDONNEES',2G15.6/,
     S           ' CONFONDU AVEC LE NOEUD',I6/,
     S           ' COORDONNEES',2G15.6)
                 ENDIF
C==================================================================
                  GOTO 11
               ELSE
                  X2 = COOR(NJ,IKJ,NBJA+1,1)
                  Y2 = COOR(NJ,IKJ,NBJA+1,2)
                  XY  = X2 - X1 + Y2 - Y1
                  XY1 = XM - X1 + YM - Y1
                  XY2 = X2 - XM + Y2 - YM
                  RAP = XY1 / XY
                  IF (RAP.GE.0.D0 .AND. RAP.LT.1.D0) THEN
                     DMCN(IAQ)   = XY2 / XY
                     DMCN(IAQ+1) = XY1 / XY
                     IAQ = IAQ + 2
                     MCN(MNN)   = NBJA
                     MCN(MNN+1) = NBJA+1
                     MNN = MNN + 2
C==================================================================
                   IF (NOIMPR.GT.10)  THEN
                   WRITE(IMPRIM,9005) LINUNO(NJ,IJK,NBJM),XM,YM,
     S                              LINUNO(NJ,IKJ,NBJA),X1,Y1,
     S                              LINUNO(NJ,IKJ,NBJA+1),X2,Y2
 9005              FORMAT(' NOEUD',I6/,' COORDONNEES',2G15.6/,
     S             ' COMPRIS ENTRE LE NOEUD',I6/,
     S             ' COORDONNEES',2G15.6/,
     S             ' ET LE NOEUD',I6/,
     S             ' COORDONNEES',2G15.6)
                   ENDIF
C==================================================================
                     GOTO 11
                  ENDIF
               ENDIF
 12         CONTINUE
 11      CONTINUE
C=====================================================================
      IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,4525) (LINUNO(NJ,IJK,JJ+1),(DMCN(IAQQ+JJ*2+K),
     S                       LINUNO(NJ,IKJ,MCN(MNNUM+JJ*2+K)),
     S                       K=0,1),JJ=0,NBJ-1)
 4525    FORMAT(' NOEUD',I5,' : QI  ',D15.6,I6,D15.6,I6)
      ENDIF
C=====================================================================
C
 1    CONTINUE
C
C     LES JOINTS ESCLAVES DU DOMAINE
C     ------------------------------
C
      DO 3 J = 1 , NBJOIE
C        LE NUMERO DU JOINT ESCLAVE
         NJ = JOIE(J)
C====================================================
         IF (NOIMPR.GT.10) WRITE(IMPRIM,9006) NJ
 9006    FORMAT(/,' JOINT ESCLAVE',I6/,1X,13(1H-))
C====================================================
         NJ10 = (NJ-1)*10
C        LE NOMBRE DE NOEUDS DU JOINT
         NBJ1=NBNOEU(NJ,1)
         NBJ2=NBNOEU(NJ,2)
C        LES SOUS-DOMAINES LIES AU JOINT
         NUSD1 = LINUOB(NJ,2)
         NUSD2 = LINUOB(NJ,4)
C        LE NUMERO DU SOUS-DOMAINE
         IF ( NUSD .EQ. NUSD1 ) THEN
           IJK = 1
           IKJ = 2
           NBJ = NBJ1
           NBK = NBJ2
         ELSE IF ( NUSD .EQ. NUSD2 ) THEN
           IJK = 2
           IKJ = 1
           NBJ = NBJ2
           NBK = NBJ1
         ENDIF
         MNQE = MCN(MNMAJO+NJ10+5)
         IAQE = ( MNQE - 1 ) / 2 + 1
         IAQ  = IAQE
         MNSE = MCN(MNMAJO+NJ10+7)
         IASE = ( MNSE - 1 ) / 2 + 1
         IAS  = IASE
         IAQ1 = IAQ + 3
         IAS1 = IAS
C        LE PREMIER NOEUD DU JOINT
         X1 = COOR(NJ,IJK,1,1)
         Y1 = COOR(NJ,IJK,1,2)
         X2 = COOR(NJ,IJK,2,1)
         Y2 = COOR(NJ,IJK,2,2)
         XY = X2 - X1 + Y2 - Y1
         XY1 = XY
C====================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,9002) LINUNO(NJ,IJK,1)
         WRITE(IMPRIM,9003) X1,Y1
         ENDIF
C====================================================
         DMCN(IAQ)   = 0.D0
         DMCN(IAQ+1) = 1.D0
         DMCN(IAQ+2) = 0.D0
         DMCN(IAS)   = XY / 2.D0
         DO 4 NBJE = 2 , NBJ - 1
            X1 = COOR(NJ,IJK,NBJE,1)
            Y1 = COOR(NJ,IJK,NBJE,2)
            X2 = COOR(NJ,IJK,NBJE+1,1)
            Y2 = COOR(NJ,IJK,NBJE+1,2)
            XYA = XY
            XY = X2 - X1 + Y2 - Y1
C====================================================
            IF (NOIMPR.GT.10)  THEN
            WRITE(IMPRIM,9002) LINUNO(NJ,IJK,NBJE)
            WRITE(IMPRIM,9003) X1,Y1
            ENDIF
C====================================================
            IAQ = IAQ + 3
            IAS = IAS +1
            DMCN(IAQ)   = DMCN(IAQ-1)
C8          DMCN(IAQ+1) = 2.D0 * XY / 3.D0
            DMCN(IAQ+1) = ( XY + XYA ) / 3.D0
            DMCN(IAQ+2) = XY / 6.D0
            DMCN(IAS) = 0.D0
 4       CONTINUE
C        LE DERNIER NOEUD DU JOINT
         X1 = COOR(NJ,IJK,NBJ-1,1)
         Y1 = COOR(NJ,IJK,NBJ-1,2)
         X2 = COOR(NJ,IJK,NBJ,1)
         Y2 = COOR(NJ,IJK,NBJ,2)
         XY = X2 - X1 + Y2 - Y1
         XYN = XY
C====================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,9002) LINUNO(NJ,IJK,NBJ)
         WRITE(IMPRIM,9003) X2,Y2
         ENDIF
C====================================================
         IAQ2 = IAQ
         IAS2 = IAS +1
         IAQ = IAQ + 3
         IAS = IAS +1
         DMCN(IAQ)   = 0.D0
         DMCN(IAQ+1) = 1.D0
         DMCN(IAQ+2) = 0.D0
C*       DMCN(IAS)   = XY / 2.D0
         DMCN(IAS)   = - XY / 2.D0
C        CORRECTION DE LA DEUXIEME LIGNE DE Q
C*       DMCN(IAQ1)   = DMCN(IAS1)
         DMCN(IAQ1)   = 0.D0
C*       DMCN(IAQ1+1) = DMCN(IAQ1+1) + DMCN(IAS1)
         DMCN(IAQ1+1) = DMCN(IAQ1+1) + XY1 / 6.D0
C        CORRECTION DE L'AVANT DERNIERE LIGNE DE Q
C*       DMCN(IAQ2+1) = DMCN(IAQ2+1) + DMCN(IAS2)
         DMCN(IAQ2+1) = DMCN(IAQ2+1) + XYN / 6.D0
C*       DMCN(IAQ2+2) = DMCN(IAS2)
         DMCN(IAQ2+2) = 0.D0
C        FACTORISATION LU DE LA MATRICE QE
         NBLU = NBJ-2
         IAQE = ( MNQE - 1 ) / 2 + 1
         IALU = ( MCN(MNMAJO+NJ10+6) - 1 ) / 2 + 1
         IF (NBLU.NE.0) CALL SDRES10(DMCN(IAQE+3),DMCN(IALU),NBLU)
C==================================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,4427) ((DMCN(IAQE+JJ*3+K),K=0,2),
     S                        DMCN(IASE+JJ),JJ=0,NBJ-1)
         WRITE(IMPRIM,4428) ((DMCN(IALU+JJ*3+K),K=0,2),JJ=0,NBJ-3)
 4427    FORMAT(' MATRICE QE  ',3D15.6,' SM',D15.6)
 4428    FORMAT(' MATRICE LUE ',3D15.6)
         ENDIF
C==================================================================
C        LA MATRICE QQ D'INTERPOLATION
         MNQQE = MCN(MNMAJO+NJ10+8)
         IAQQ = ( MNQQE - 1 ) / 2 + 1
         IAQ = IAQQ
         MNNUE = MCN(MNMAJO+NJ10+9)
         MNN = MNNUE
         DO 13 NBJE = 1 , NBJ
            XE = COOR(NJ,IJK,NBJE,1)
            YE = COOR(NJ,IJK,NBJE,2)
            XYZE(1) = XE
            XYZE(2) = YE
            XYZE(3) = 0.D0
            DO 14 NBJA =  1  , NBK
               X1 = COOR(NJ,IKJ,NBJA,1)
               Y1 = COOR(NJ,IKJ,NBJA,2)
               XYZ(1) = X1
               XYZ(2) = Y1
               XYZ(3) = 0.D0
               CALL XYZIDE(XYZE,XYZ,IDENT)
               IF (IDENT.EQ.1) THEN
C              LES POINTS SONT CONFONDUS !
                  DMCN(IAQ)   = 1.D0
                  DMCN(IAQ+1) = 0.D0
                  IAQ = IAQ + 2
                  MCN(MNN)   = NBJA
                  MCN(MNN+1) = NBJA
                  MNN = MNN + 2
C==================================================================
                  IF (NOIMPR.GT.10)  THEN
                  WRITE(IMPRIM,9004) LINUNO(NJ,IJK,NBJE),XE,YE,
     S                               LINUNO(NJ,IKJ,NBJA),X1,Y1
                  ENDIF
C==================================================================
                  GOTO 13
               ELSE
                  X2 = COOR(NJ,IKJ,NBJA+1,1)
                  Y2 = COOR(NJ,IKJ,NBJA+1,2)
                  XY  = X2 - X1 + Y2 - Y1
                  XY1 = XE - X1 + YE - Y1
                  XY2 = X2 - XE + Y2 - YE
                  RAP = XY1 / XY
                  IF (RAP.GE.0.D0 .AND. RAP.LT.1.D0) THEN
                     DMCN(IAQ)   = XY2 / XY
                     DMCN(IAQ+1) = XY1 / XY
                     IAQ = IAQ + 2
                     MCN(MNN)   = NBJA
                     MCN(MNN+1) = NBJA+1
                     MNN = MNN + 2
C==================================================================
                 IF (NOIMPR.GT.10)  THEN
                 WRITE(IMPRIM,9005) LINUNO(NJ,IJK,NBJE),XE,YE,
     S                              LINUNO(NJ,IKJ,NBJA),X1,Y1,
     S                              LINUNO(NJ,IKJ,NBJA+1),X2,Y2
                 ENDIF
C==================================================================
                     GOTO 13
                  ENDIF
               ENDIF
 14         CONTINUE
 13      CONTINUE
C====================================================================
         IF (NOIMPR.GT.10)  THEN
         WRITE(IMPRIM,4525) (LINUNO(NJ,IJK,JJ+1),(DMCN(IAQQ+JJ*2+K),
     S                       LINUNO(NJ,IKJ,MCN(MNNUE+JJ*2+K)),
     S                       K=0,1),JJ=0,NBJ-1)
         ENDIF
C====================================================================
C
 3    CONTINUE
C
      RETURN
      END
