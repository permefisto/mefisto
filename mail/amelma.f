      SUBROUTINE AMELMA(NBS1,NBS2,NBS3,COSO,NITC,
     S                  FCI,FCJ,FCK,FADIST,DISTID,LGMAX,PROF)
C***********************************************************************
C BUT: CALCUL DES FONCTIONS DE CONTROLE POUR L AMELIORATIONS DES
C      MAILLES AU BORD
C***********************************************************************
C ENTREES:
C           NBS1    : NOMBRE DE POINTS SUR LIGNES X
C           NBS2    : NOMBRE DE POINTS SUR LIGNES Y
C           NBS3    : NOMBRE DE POINTS SUR LIGNES Z
C           NITC    : NUMERO DE L ITERATION DE CORRECTION
C
C ENTREES ET RESULTAT:
C           COSO  : COORDONNEES DES POINTS DU MAILLAGE
C
C TRAVAIL:
C           FCI   : VECTEUR FONCTIONS DE CONTROLE EN I
C           FCJ   : VECTEUR FONCTIONS DE CONTROLE EN J
C           FCK   : VECTEUR FONCTIONS DE CONTROLE EN K
C***********************************************************************
C AUTEUR : DOURSAT CHRISTOPHE  ANALYSE NUMERIQUE PARIS  OCTOBRE 1989
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C***********************************************************************
      DIMENSION COSO(3,LGMAX)
      DIMENSION DISTID(LGMAX),FCI(LGMAX),FCJ(LGMAX),FCK(LGMAX)
      DIMENSION V1(3),V2(3),V(3),VN(3),V3(3),V4(3),FADIST(6)
      NBASE = NBS1*NBS2
C***********************************************************************
C         CALCUL DES FONCTIONS DE CONTROLE SUR CHAQUE FACE
C***********************************************************************
C                          SURFACE DU BAS
C***********************************************************************
      DO 80 J=2,NBS2-1
        DO 70 I=2,NBS1-1
          NEU = (J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 10 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBS1)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBS1)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU+NBASE)-COSO(II,NEU)
            DD     = DD+V(II)**2
   10     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(1).GE.0 ) THEN
              DISTID(NEU) = FADIST(1)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(1).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(1).EQ.-2) THEN
              NI1B =          (J-1)*NBS1+1
              NI1H = NBASE   +(J-1)*NBS1+1
              NI2B =          (J-1)*NBS1+NBS1
              NI2H = NBASE   +(J-1)*NBS1+NBS1
              NJ1B =                     I
              NJ1H = NBASE              +I
              NJ2B =       (NBS2-1)*NBS1+I
              NJ2H = NBASE+(NBS2-1)*NBS1+I
              DI1 = 0.
              DI2 = 0.
              DJ1 = 0.
              DJ2 = 0.
              DO 15 II=1,3
                DI1 = DI1 + (COSO(II,NI1H)-COSO(II,NI1B))**2
                DI2 = DI2 + (COSO(II,NI2H)-COSO(II,NI2B))**2
                DJ1 = DJ1 + (COSO(II,NJ1H)-COSO(II,NJ1B))**2
                DJ2 = DJ2 + (COSO(II,NJ2H)-COSO(II,NJ2B))**2
   15         CONTINUE
              DI1 = SQRT (DI1)
              DI2 = SQRT (DI2)
              DJ1 = SQRT (DJ1)
              DJ2 = SQRT (DJ2)
              DISTID(NEU) = 0.5 * ( (NBS1-I)/(NBS1-1.) * DI1
     S                             +(I-1.)  /(NBS1-1.) * DI2
     S                             +(NBS2-J)/(NBS2-1.) * DJ1
     S                             +(J-1.)  /(NBS2-1.) * DJ2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 20 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
   20     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 30 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
   30     CONTINUE
          VN(1) = V1(2)*V2(3)-V1(3)*V2(2)
          VN(2) = V1(3)*V2(1)-V1(1)*V2(3)
          VN(3) = V1(1)*V2(2)-V1(2)*V2(1)
          D = 0
          DO 35 II=1,3
            D = D+VN(II)**2
   35     CONTINUE
          D = SQRT(D)
          DO 37 II=1,3
            VN(II) = VN(II)/D
   37     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 40 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
   40     CONTINUE
          D3 = 0
          D4 = 0
          DO 50 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
   50     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 60 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
   60     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCI(NEU) = ATAN(ALP1)
          FCJ(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(1).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCK(NEU) = (DD-DISTI)/DISTI
            ELSE
              FCK(NEU) = ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
   70   CONTINUE
   80 CONTINUE
C***********************************************************************
C                          SURFACE DU HAUT
C***********************************************************************
      K = NBS3
      DO 180 J=2,NBS2-1
        DO 170 I=2,NBS1-1
          NEU = (K-1)*NBASE+(J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 110 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBS1)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBS1)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU-NBASE)-COSO(II,NEU)
            DD     = DD+V(II)**2
  110     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(4).GE.0 ) THEN
              DISTID(NEU) = FADIST(4)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(4).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(4).EQ.-2) THEN
              NI1B = (K-1)*NBASE +    (J-1)*NBS1 +1
              NI1H = (K-2)*NBASE +    (J-1)*NBS1 +1
              NI2B = (K-1)*NBASE +    (J-1)*NBS1 +NBS1
              NI2H = (K-2)*NBASE +    (J-1)*NBS1 +NBS1
              NJ1B = (K-1)*NBASE                 +I
              NJ1H = (K-2)*NBASE                 +I
              NJ2B = (K-1)*NBASE + (NBS2-1)*NBS1 +I
              NJ2H = (K-2)*NBASE + (NBS2-1)*NBS1 +I
              DI1 = 0.
              DI2 = 0.
              DJ1 = 0.
              DJ2 = 0.
              DO 115 II=1,3
                DI1 = DI1 + (COSO(II,NI1H)-COSO(II,NI1B))**2
                DI2 = DI2 + (COSO(II,NI2H)-COSO(II,NI2B))**2
                DJ1 = DJ1 + (COSO(II,NJ1H)-COSO(II,NJ1B))**2
                DJ2 = DJ2 + (COSO(II,NJ2H)-COSO(II,NJ2B))**2
  115         CONTINUE
              DI1 = SQRT (DI1)
              DI2 = SQRT (DI2)
              DJ1 = SQRT (DJ1)
              DJ2 = SQRT (DJ2)
              DISTID(NEU) = 0.5 * ( (NBS1-I)/(NBS1-1.) * DI1
     S                             +(I-1.)  /(NBS1-1.) * DI2
     S                             +(NBS2-J)/(NBS2-1.) * DJ1
     S                             +(J-1.)  /(NBS2-1.) * DJ2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 120 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
  120     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 130 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
  130     CONTINUE
          VN(1) = -(V1(2)*V2(3)-V1(3)*V2(2))
          VN(2) = -(V1(3)*V2(1)-V1(1)*V2(3))
          VN(3) = -(V1(1)*V2(2)-V1(2)*V2(1))
          D = 0
          DO 135 II=1,3
            D = D+VN(II)**2
  135     CONTINUE
          D = SQRT(D)
          DO 137 II=1,3
            VN(II) = VN(II)/D
  137     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 140 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
  140     CONTINUE
          D3 = 0
          D4 = 0
          DO 150 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
  150     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 160 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
  160     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCI(NEU) = ATAN(ALP1)
          FCJ(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(4).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCK(NEU) = -(DD-DISTI)/DISTI
            ELSE
              FCK(NEU) = -ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
  170   CONTINUE
  180 CONTINUE
C***********************************************************************
C                          SURFACE DE DEVANT
C***********************************************************************
      J = 1
      DO 280 K=2,NBS3-1
        DO 270 I=2,NBS1-1
          NEU = (K-1)*NBASE+(J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 210 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBASE)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBASE)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU+NBS1)-COSO(II,NEU)
            DD     = DD+V(II)**2
  210     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(2).GE.0 ) THEN
              DISTID(NEU) = FADIST(2)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(2).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(2).EQ.-2) THEN
              NI1B =    (K-1)*NBASE        + 1
              NI1H =    (K-1)*NBASE + NBS1 + 1
              NI2B =    (K-1)*NBASE        + NBS1
              NI2H =    (K-1)*NBASE + NBS1 + NBS1
              NK1B =                         I
              NK1H =                  NBS1 + I
              NK2B = (NBS3-1)*NBASE        + I
              NK2H = (NBS3-1)*NBASE + NBS1 + I
              DI1 = 0.
              DI2 = 0.
              DK1 = 0.
              DK2 = 0.
              DO 215 II=1,3
                DI1 = DI1 + (COSO(II,NI1H)-COSO(II,NI1B))**2
                DI2 = DI2 + (COSO(II,NI2H)-COSO(II,NI2B))**2
                DK1 = DK1 + (COSO(II,NK1H)-COSO(II,NK1B))**2
                DK2 = DK2 + (COSO(II,NK2H)-COSO(II,NK2B))**2
  215         CONTINUE
              DI1 = SQRT (DI1)
              DI2 = SQRT (DI2)
              DK1 = SQRT (DK1)
              DK2 = SQRT (DK2)
              DISTID(NEU) = 0.5 * ( (NBS1-I)/(NBS1-1.) * DI1
     S                             +(I-1.)  /(NBS1-1.) * DI2
     S                             +(NBS3-K)/(NBS3-1.) * DK1
     S                             +(K-1.)  /(NBS3-1.) * DK2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 220 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
  220     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 230 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
  230     CONTINUE
          VN(1) = -(V1(2)*V2(3)-V1(3)*V2(2))
          VN(2) = -(V1(3)*V2(1)-V1(1)*V2(3))
          VN(3) = -(V1(1)*V2(2)-V1(2)*V2(1))
          D = 0
          DO 235 II=1,3
            D = D+VN(II)**2
  235     CONTINUE
          D = SQRT(D)
          DO 237 II=1,3
            VN(II) = VN(II)/D
  237     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 240 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
  240     CONTINUE
          D3 = 0
          D4 = 0
          DO 250 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
  250     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 260 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
  260     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCI(NEU) = ATAN(ALP1)
          FCK(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(2).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCJ(NEU) = (DD-DISTI)/DISTI
            ELSE
              FCJ(NEU) = ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
  270   CONTINUE
  280 CONTINUE
C***********************************************************************
C                          SURFACE DE DERRIERE
C***********************************************************************
      J = NBS2
      DO 380 K=2,NBS3-1
        DO 370 I=2,NBS1-1
          NEU = (K-1)*NBASE+(J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 310 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBASE)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBASE)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU-NBS1)-COSO(II,NEU)
            DD     = DD+V(II)**2
  310     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(5).GE.0 ) THEN
              DISTID(NEU) = FADIST(5)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(5).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(5).EQ.-2) THEN
              NI1B =    (K-1)*NBASE + (NBS2-1)*NBS1 + 1
              NI1H =    (K-1)*NBASE + (NBS2-2)*NBS1 + 1
              NI2B =    (K-1)*NBASE + (NBS2-1)*NBS1 + NBS1
              NI2H =    (K-1)*NBASE + (NBS2-2)*NBS1 + NBS1
              NK1B =                  (NBS2-1)*NBS1 + I
              NK1H =                  (NBS2-2)*NBS1 + I
              NK2B = (NBS3-1)*NBASE + (NBS2-1)*NBS1 + I
              NK2H = (NBS3-1)*NBASE + (NBS2-2)*NBS1 + I
              DI1 = 0.
              DI2 = 0.
              DK1 = 0.
              DK2 = 0.
              DO 315 II=1,3
                DI1 = DI1 + (COSO(II,NI1H)-COSO(II,NI1B))**2
                DI2 = DI2 + (COSO(II,NI2H)-COSO(II,NI2B))**2
                DK1 = DK1 + (COSO(II,NK1H)-COSO(II,NK1B))**2
                DK2 = DK2 + (COSO(II,NK2H)-COSO(II,NK2B))**2
  315         CONTINUE
              DI1 = SQRT (DI1)
              DI2 = SQRT (DI2)
              DK1 = SQRT (DK1)
              DK2 = SQRT (DK2)
              DISTID(NEU) = 0.5 * ( (NBS1-I)/(NBS1-1.) * DI1
     S                             +(I-1.)  /(NBS1-1.) * DI2
     S                             +(NBS3-K)/(NBS3-1.) * DK1
     S                             +(K-1.)  /(NBS3-1.) * DK2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 320 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
  320     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 330 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
  330     CONTINUE
          VN(1) = V1(2)*V2(3)-V1(3)*V2(2)
          VN(2) = V1(3)*V2(1)-V1(1)*V2(3)
          VN(3) = V1(1)*V2(2)-V1(2)*V2(1)
          D = 0
          DO 335 II=1,3
            D = D+VN(II)**2
  335     CONTINUE
          D = SQRT(D)
          DO 337 II=1,3
            VN(II) = VN(II)/D
  337     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 340 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
  340     CONTINUE
          D3 = 0
          D4 = 0
          DO 350 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
  350     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 360 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
  360     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCI(NEU) = ATAN(ALP1)
          FCK(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(5).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCJ(NEU) = -(DD-DISTI)/DISTI
            ELSE
              FCJ(NEU) = -ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
  370   CONTINUE
  380 CONTINUE
C***********************************************************************
C                          SURFACE DE GAUCHE
C***********************************************************************
      I = 1
      DO 480 K=2,NBS3-1
        DO 470 J=2,NBS2-1
          NEU = (K-1)*NBASE+(J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 410 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-NBS1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+NBS1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBASE)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBASE)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU+1)-COSO(II,NEU)
            DD     = DD+V(II)**2
  410     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(3).GE.0 ) THEN
              DISTID(NEU) = FADIST(3)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(3).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(3).EQ.-2) THEN
              NJ1B =    (K-1)*NBASE                 + 1
              NJ1H =    (K-1)*NBASE                 + 2
              NJ2B =    (K-1)*NBASE + (NBS2-1)*NBS1 + 1
              NJ2H =    (K-1)*NBASE + (NBS2-1)*NBS1 + 2
              NK1B =                     (J-1)*NBS1 + 1
              NK1H =                     (J-1)*NBS1 + 2
              NK2B = (NBS3-1)*NBASE +    (J-1)*NBS1 + 1
              NK2H = (NBS3-1)*NBASE +    (J-1)*NBS1 + 2
              DJ1 = 0.
              DJ2 = 0.
              DK1 = 0.
              DK2 = 0.
              DO 415 II=1,3
                DJ1 = DJ1 + (COSO(II,NJ1H)-COSO(II,NJ1B))**2
                DJ2 = DJ2 + (COSO(II,NJ2H)-COSO(II,NJ2B))**2
                DK1 = DK1 + (COSO(II,NK1H)-COSO(II,NK1B))**2
                DK2 = DK2 + (COSO(II,NK2H)-COSO(II,NK2B))**2
  415         CONTINUE
              DJ1 = SQRT (DJ1)
              DJ2 = SQRT (DJ2)
              DK1 = SQRT (DK1)
              DK2 = SQRT (DK2)
              DISTID(NEU) = 0.5 * ( (NBS2-J)/(NBS2-1.) * DJ1
     S                             +(J-1.)  /(NBS2-1.) * DJ2
     S                             +(NBS3-K)/(NBS3-1.) * DK1
     S                             +(K-1.)  /(NBS3-1.) * DK2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 420 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
  420     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 430 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
  430     CONTINUE
          VN(1) = V1(2)*V2(3)-V1(3)*V2(2)
          VN(2) = V1(3)*V2(1)-V1(1)*V2(3)
          VN(3) = V1(1)*V2(2)-V1(2)*V2(1)
          D = 0
          DO 435 II=1,3
            D = D+VN(II)**2
  435     CONTINUE
          D = SQRT(D)
          DO 437 II=1,3
            VN(II) = VN(II)/D
  437     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 440 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
  440     CONTINUE
          D3 = 0
          D4 = 0
          DO 450 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
  450     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 460 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
  460     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCJ(NEU) = ATAN(ALP1)
          FCK(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(3).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCI(NEU) = (DD-DISTI)/DISTI
            ELSE
              FCI(NEU) = ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
  470   CONTINUE
  480 CONTINUE
C***********************************************************************
C                          SURFACE DE DROITE
C***********************************************************************
      I = NBS1
      DO 580 K=2,NBS3-1
        DO 570 J=2,NBS2-1
          NEU = (K-1)*NBASE+(J-1)*NBS1+I
C-----------------------------------------------------------------------
C                CALCUL DU VECTEUR NORMAL A LA SURFACE
C-----------------------------------------------------------------------
          D1 = 0
          D2 = 0
          D3 = 0
          D4 = 0
          DD = 0
          DO 510 II=1,3
            V1(II) = COSO(II,NEU)-COSO(II,NEU-NBS1)
            D1     = D1+V1(II)**2
            V2(II) = COSO(II,NEU+NBS1)-COSO(II,NEU)
            D2     = D2+V2(II)**2
            V3(II) = COSO(II,NEU)-COSO(II,NEU-NBASE)
            D3     = D3+V3(II)**2
            V4(II) = COSO(II,NEU+NBASE)-COSO(II,NEU)
            D4     = D4+V4(II)**2
            V(II)  = COSO(II,NEU-1)-COSO(II,NEU)
            DD     = DD+V(II)**2
  510     CONTINUE
          D1 = SQRT(D1)
          D2 = SQRT(D2)
          D3 = SQRT(D3)
          D4 = SQRT(D4)
          DD = SQRT(DD)
          IF (NITC.EQ.1) THEN
            IF ( FADIST(6).GE.0 ) THEN
              DISTID(NEU) = FADIST(6)*(D1+D2+D3+D4)/4
            ELSE IF (FADIST(6).EQ.-1) THEN
              DISTID(NEU) = DD
            ELSE IF (FADIST(6).EQ.-2) THEN
              NJ1B =    (K-1)*NBASE                 + NBS1
              NJ1H =    (K-1)*NBASE                 + NBS1-1
              NJ2B =    (K-1)*NBASE + (NBS2-1)*NBS1 + NBS1
              NJ2H =    (K-1)*NBASE + (NBS2-1)*NBS1 + NBS1-1
              NK1B =                     (J-1)*NBS1 + NBS1
              NK1H =                     (J-1)*NBS1 + NBS1-1
              NK2B = (NBS3-1)*NBASE +    (J-1)*NBS1 + NBS1
              NK2H = (NBS3-1)*NBASE +    (J-1)*NBS1 + NBS1-1
              DJ1 = 0.
              DJ2 = 0.
              DK1 = 0.
              DK2 = 0.
              DO 515 II=1,3
                DJ1 = DJ1 + (COSO(II,NJ1H)-COSO(II,NJ1B))**2
                DJ2 = DJ2 + (COSO(II,NJ2H)-COSO(II,NJ2B))**2
                DK1 = DK1 + (COSO(II,NK1H)-COSO(II,NK1B))**2
                DK2 = DK2 + (COSO(II,NK2H)-COSO(II,NK2B))**2
  515         CONTINUE
              DJ1 = SQRT (DJ1)
              DJ2 = SQRT (DJ2)
              DK1 = SQRT (DK1)
              DK2 = SQRT (DK2)
              DISTID(NEU) = 0.5 * ( (NBS2-J)/(NBS2-1.) * DJ1
     S                             +(J-1.)  /(NBS2-1.) * DJ2
     S                             +(NBS3-K)/(NBS3-1.) * DK1
     S                             +(K-1.)  /(NBS3-1.) * DK2 )
            ENDIF
          ENDIF
          DD1 = 0
          DD2 = 0
          DO 520 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
            V3(II) = V3(II)/D3
            V4(II) = V4(II)/D4
            V(II)  = V(II)/DD
            V1(II) = V1(II)+V2(II)
            V2(II) = V3(II)+V4(II)
            DD1    = DD1+V1(II)**2
            DD2    = DD2+V2(II)**2
  520     CONTINUE
          D1 = SQRT(DD1)
          D2 = SQRT(DD2)
          DO 530 II=1,3
            V1(II) = V1(II)/D1
            V2(II) = V2(II)/D2
  530     CONTINUE
          VN(1) = -(V1(2)*V2(3)-V1(3)*V2(2))
          VN(2) = -(V1(3)*V2(1)-V1(1)*V2(3))
          VN(3) = -(V1(1)*V2(2)-V1(2)*V2(1))
          D = 0
          DO 535 II=1,3
            D = D+VN(II)**2
  535     CONTINUE
          D = SQRT(D)
          DO 537 II=1,3
            VN(II) = VN(II)/D
  537     CONTINUE
C-----------------------------------------------------------------------
C        CALCUL DES ERREURS SUR LES ANGLES ET DISTANCES
C-----------------------------------------------------------------------
          V3(1) = V1(2)*VN(3)-V1(3)*VN(2)
          V3(2) = V1(3)*VN(1)-V1(1)*VN(3)
          V3(3) = V1(1)*VN(2)-V1(2)*VN(1)
          V4(1) = V2(2)*VN(3)-V2(3)*VN(2)
          V4(2) = V2(3)*VN(1)-V2(1)*VN(3)
          V4(3) = V2(1)*VN(2)-V2(2)*VN(1)
          PS1 = 0
          PS2 = 0
          DO 540 II=1,3
            PS1 = PS1+V(II)*V3(II)
            PS2 = PS2+V(II)*V4(II)
  540     CONTINUE
          D3 = 0
          D4 = 0
          DO 550 II=1,3
            V3(II) = V(II)-PS1*V3(II)
            V4(II) = V(II)-PS2*V4(II)
            D3    = D3+V3(II)**2
            D4    = D4+V4(II)**2
  550     CONTINUE
          D3  = SQRT(D3)
          D4  = SQRT(D4)
          PS1 = 0
          PS2 = 0
          PSENS1 = 0
          PSENS2 = 0
          DO 560 II=1,3
            V3(II)  = V3(II)/D3
            V4(II)  = V4(II)/D4
            PS1 = PS1+V3(II)*VN(II)
            PS2 = PS2+V4(II)*VN(II)
            PSENS1 = PSENS1+V3(II)*V1(II)
            PSENS2 = PSENS2+V4(II)*V2(II)
  560     CONTINUE
          PS1 = MAX(PS1,-1.0)
          PS1 = MIN(PS1,1.0)
          PS2 = MAX(PS2,-1.0)
          PS2 = MIN(PS2,1.0)
          ALP1 = ACOS(PS1)
          ALP2 = ACOS(PS2)
          IF (PSENS1.LT.0) THEN
            ALP1 = -ALP1
          ENDIF
          IF (PSENS2.LT.0) THEN
            ALP2 = -ALP2
          ENDIF
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES ANGLES
C-----------------------------------------------------------------------
          FCJ(NEU) = ATAN(ALP1)
          FCK(NEU) = ATAN(ALP2)
C-----------------------------------------------------------------------
C                   CORRECTION SUR LES LONGUEURS
C-----------------------------------------------------------------------
          IF ( FADIST(6).NE.0 ) THEN
            DISTI = DISTID(NEU)
            IF (DD.LT.DISTI) THEN
              FCI(NEU) = -(DD-DISTI)/DISTI
            ELSE
              FCI(NEU) = -ATAN((DD-DISTI)/DD)
            ENDIF
          ENDIF
  570   CONTINUE
  580 CONTINUE
C***********************************************************************
C                 INTERPOLATION PAR EXPONENTIELLES
C***********************************************************************
      COD = 1.5/2.
      COA = 1.2/2.
      DO 1030 K=2,NBS3-1
        DO 1020 J=2,NBS2-1
          DO 1010 I=2,NBS1-1
            NEU  =    (K-1)*NBASE   +(J-1)*NBS1   +I
            NEUB =                   (J-1)*NBS1   +I
            NEUH = (NBS3-1)*NBASE   +(J-1)*NBS1   +I
            NEUG =    (K-1)*NBASE   +(J-1)*NBS1   +1
            NEUD =    (K-1)*NBASE   +(J-1)*NBS1+NBS1
            NEUF =    (K-1)*NBASE                 +I
            NEUT =    (K-1)*NBASE+(NBS2-1)*NBS1   +I
C-----------------------------------------------------------------------
C       CALCUL DES COEFFICIENTS POIDS DEPENDANT DE LA PROXIMITE
C            DU NOEUD PAR RAPPORT AUX  FACES DU CUBE
C-----------------------------------------------------------------------
            A1 = ABS((I-1.)/(NBS1-1)-0.5)
            A2 = ABS((J-1.)/(NBS2-1)-0.5)
            A3 = ABS((K-1.)/(NBS3-1)-0.5)
            DIDIST = 1600*(MIN(A2,A3))**4+PROF
            DIANGJ = 1600*(MIN(A1,A2))**4+PROF
            DIANGK = 1600*(MIN(A1,A3))**4+PROF
            DJDIST = 1600*(MIN(A1,A3))**4+PROF
            DJANGI = 1600*(MIN(A1,A2))**4+PROF
            DJANGK = 1600*(MIN(A2,A3))**4+PROF
            DKDIST = 1600*(MIN(A1,A2))**4+PROF
            DKANGI = 1600*(MIN(A1,A3))**4+PROF
            DKANGJ = 1600*(MIN(A2,A3))**4+PROF
            CID  = COD*EXP((DIDIST-PROF)/(NBS1-1.))
            CIAJ = COA*EXP((DIANGJ-PROF)/(NBS2-1.))
            CIAK = COA*EXP((DIANGK-PROF)/(NBS3-1.))
            CJD  = COD*EXP((DJDIST-PROF)/(NBS2-1.))
            CJAI = COA*EXP((DJANGI-PROF)/(NBS1-1.))
            CJAK = COA*EXP((DJANGK-PROF)/(NBS3-1.))
            CKD  = COD*EXP((DKDIST-PROF)/(NBS3-1.))
            CKAI = COA*EXP((DKANGI-PROF)/(NBS1-1.))
            CKAJ = COA*EXP((DKANGJ-PROF)/(NBS2-1.))
C            DIDIST = 1600*(MIN(A2,A3))**4+5.
C            DIANGJ = 1600*(MIN(A1,A2))**4+5.
C            DIANGK = 1600*(MIN(A1,A3))**4+5.
C            DJDIST = 1600*(MIN(A1,A3))**4+5.
C            DJANGI = 1600*(MIN(A1,A2))**4+5.
C            DJANGK = 1600*(MIN(A2,A3))**4+5.
C            DKDIST = 1600*(MIN(A1,A2))**4+5.
C            DKANGI = 1600*(MIN(A1,A3))**4+5.
C            DKANGJ = 1600*(MIN(A2,A3))**4+5.
C            CID  = COD*EXP((DIDIST-5.)/(NBS1-1.))
C            CIAJ = COA*EXP((DIANGJ-5.)/(NBS2-1.))
C            CIAK = COA*EXP((DIANGK-5.)/(NBS3-1.))
C            CJD  = COD*EXP((DJDIST-5.)/(NBS2-1.))
C            CJAI = COA*EXP((DJANGI-5.)/(NBS1-1.))
C            CJAK = COA*EXP((DJANGK-5.)/(NBS3-1.))
C            CKD  = COD*EXP((DKDIST-5.)/(NBS3-1.))
C            CKAI = COA*EXP((DKANGI-5.)/(NBS1-1.))
C            CKAJ = COA*EXP((DKANGJ-5.)/(NBS2-1.))
C-----------------------------------------------------------------------
C                        INTERPOLATION
C-----------------------------------------------------------------------
            FCI(NEU) = FCI(NEU)
     S                +CID *EXP(-DIDIST*(I-1.)/(NBS1-1.))  *FCI(NEUG)
     S                +CID *EXP(-DIDIST*(NBS1-I)/(NBS1-1.))*FCI(NEUD)
     S                +CIAJ*EXP(-DIANGJ*(J-1.)/(NBS2-1.))  *FCI(NEUF)
     S                +CIAJ*EXP(-DIANGJ*(NBS2-J)/(NBS2-1.))*FCI(NEUT)
     S                +CIAK*EXP(-DIANGK*(K-1.)/(NBS3-1.))  *FCI(NEUB)
     S                +CIAK*EXP(-DIANGK*(NBS3-K)/(NBS3-1.))*FCI(NEUH)
            FCJ(NEU) = FCJ(NEU)
     S                +CJAI*EXP(-DJANGI*(I-1.)/(NBS1-1.))  *FCJ(NEUG)
     S                +CJAI*EXP(-DJANGI*(NBS1-I)/(NBS1-1.))*FCJ(NEUD)
     S                +CJD *EXP(-DJDIST*(J-1.)/(NBS2-1.))  *FCJ(NEUF)
     S                +CJD *EXP(-DJDIST*(NBS2-J)/(NBS2-1.))*FCJ(NEUT)
     S                +CJAK*EXP(-DJANGK*(K-1.)/(NBS3-1.))  *FCJ(NEUB)
     S                +CJAK*EXP(-DJANGK*(NBS3-K)/(NBS3-1.))*FCJ(NEUH)
            FCK(NEU) = FCK(NEU)
     S                +CKAI*EXP(-DKANGI*(I-1.)/(NBS1-1.))  *FCK(NEUG)
     S                +CKAI*EXP(-DKANGI*(NBS1-I)/(NBS1-1.))*FCK(NEUD)
     S                +CKAJ*EXP(-DKANGJ*(J-1.)/(NBS2-1.))  *FCK(NEUF)
     S                +CKAJ*EXP(-DKANGJ*(NBS2-J)/(NBS2-1.))*FCK(NEUT)
     S                +CKD *EXP(-DKDIST*(K-1.)/(NBS3-1.))  *FCK(NEUB)
     S                +CKD *EXP(-DKDIST*(NBS3-K)/(NBS3-1.))*FCK(NEUH)
 1010     CONTINUE
 1020   CONTINUE
 1030 CONTINUE
      RETURN
      END
