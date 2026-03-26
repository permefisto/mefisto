      SUBROUTINE AMELMA2D( NBS1,NBS2,NITC,COSO,ALPHA,
     S                     SENSFR,FCI,FCJ,DISTID,FADIST,  TRAV )
C***********************************************************************
C BUT :    PROCEDURE DE CALCUL DES FONCTIONS DE CONTROLE
C          POUR AMELIORER LA GEOMETRIE DES MAILLES FRONTIERES
C***********************************************************************
C
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           NITC   : NOMBRE D'ITERATIONS DE CORRECTION
C           FCI    : VECTEUR DES FONCTIONS DE CONTROLE EN I
C           FCJ    : VECTEUR DES FONCTIONS DE CONTROLE EN J
C           FADIST : COEFFICIENT DE REDUCTION POUR LES MAILLES VOISINES
C                    DU BORD, POUR CHAQUE COTE
C TRAVAIL:
C           DISTID : TABLEAU OU SONT STOCKEES LES DISTANCES IDEALES
C                    AU BORD DES MAILLES
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
      DIMENSION COSO(3,NBS1*NBS2),ALPHA(2*(NBS1+NBS2)-3)
      DIMENSION FCI(NBS1*NBS2),FCJ(NBS1*NBS2),TRAV(NBS1*NBS2)
      DIMENSION V1(2),V2(2),FADIST(4)
      DIMENSION DISTID(2*(NBS1+NBS2-2))
      INTEGER   SENSFR
C
      PI = 4.E0*ATAN(1.E0)
CCC      WRITE (IMPRIM,*) ' '
CCC      WRITE (IMPRIM,1000)
CCC      WRITE (IMPRIM,1010)
CCC      WRITE (IMPRIM,1020)
CCC      WRITE (IMPRIM,1060) NITC
CCC      WRITE (IMPRIM,1000)
CCC 1000 FORMAT('**************************************')
CCC 1010 FORMAT('*    AMELIORATION DE LA GEOMETRIE    *')
CCC 1020 FORMAT('*       DES MAILLES FRONTIERES       *')
CCC 1060 FORMAT('*          ITERATION',I4,'             *')
      IF (NITC.EQ.1) THEN
C***********************************************************************
C             CALCUL DE XI ET ET ABSCISSES CURVILIGNES
C***********************************************************************
C-----------------------------------------------------------------------
C                            BORD DU BAS
C-----------------------------------------------------------------------
      SOMTOT = 0.E0
      DO 200 I=2,NBS1
        DIS    = (COSO(1,I)-COSO(1,I-1))**2 + (COSO(2,I)-COSO(2,I-1))**2
        SOMTOT = SOMTOT + SQRT(DIS)
  200 CONTINUE
      SOM      = 0.E0
      DO 210 I=2,NBS1-1
        DIS    = (COSO(1,I)-COSO(1,I-1))**2 + (COSO(2,I)-COSO(2,I-1))**2
        SOM    = SOM + SQRT(DIS)
        TRAV(I)  = SOM / SOMTOT
  210 CONTINUE
C-----------------------------------------------------------------------
C                            BORD DROIT
C-----------------------------------------------------------------------
      SOMTOT = 0.E0
      DO 220 J=2,NBS2
        DIS    =  (COSO(1,J*NBS1)-COSO(1,(J-1)*NBS1))**2
     S          + (COSO(2,J*NBS1)-COSO(2,(J-1)*NBS1))**2
        SOMTOT = SOMTOT + SQRT(DIS)
  220 CONTINUE
      SOM      = 0.E0
      DO 230 J=2,NBS2-1
        DIS    =  (COSO(1,J*NBS1)-COSO(1,(J-1)*NBS1))**2
     S          + (COSO(2,J*NBS1)-COSO(2,(J-1)*NBS1))**2
        SOM    = SOM + SQRT(DIS)
        TRAV(J*NBS1)  = SOM / SOMTOT
  230 CONTINUE
C-----------------------------------------------------------------------
C                            BORD DU HAUT
C-----------------------------------------------------------------------
      SOMTOT = 0.E0
      DO 240 I=2,NBS1
        NEU    = (NBS2-1)*NBS1+I
        DIS    =  (COSO(1,NEU)-COSO(1,NEU-1))**2
     S          + (COSO(2,NEU)-COSO(2,NEU-1))**2
        SOMTOT = SOMTOT + SQRT(DIS)
  240 CONTINUE
      SOM      = 0.E0
      DO 250 I=2,NBS1-1
        NEU    = (NBS2-1)*NBS1+I
        DIS    =  (COSO(1,NEU)-COSO(1,NEU-1))**2
     S          + (COSO(2,NEU)-COSO(2,NEU-1))**2
        SOM    =  SOM + SQRT(DIS)
        TRAV((NBS2-1)*NBS1+I)  =  SOM / SOMTOT
  250 CONTINUE
C-----------------------------------------------------------------------
C                            BORD GAUCHE
C-----------------------------------------------------------------------
      SOMTOT = 0.E0
      DO 260 J=2,NBS2
        DIS    =  (COSO(1,(J-1)*NBS1+1)-COSO(1,(J-2)*NBS1+1))**2
     S          + (COSO(2,(J-1)*NBS1+1)-COSO(2,(J-2)*NBS1+1))**2
        SOMTOT = SOMTOT+SQRT(DIS)
  260 CONTINUE
      SOM      = 0.E0
      DO 270 J=2,NBS2-1
        DIS    =  (COSO(1,(J-1)*NBS1+1)-COSO(1,(J-2)*NBS1+1))**2
     S          + (COSO(2,(J-1)*NBS1+1)-COSO(2,(J-2)*NBS1+1))**2
        SOM    =  SOM + SQRT(DIS)
        TRAV((J-1)*NBS1+1)  =  SOM / SOMTOT
  270 CONTINUE
      ENDIF
C***********************************************************************
C                            BORD DU BAS
C***********************************************************************
      DO 20 I=2,NBS1-1
        NEUB = I
        NEU  = I
        V1(1) = COSO(1,NEU-1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU-1)-COSO(2,NEU)
        D1    = SQRT(V1(1)**2+V1(2)**2)
        V1(1) = COSO(1,NEU+1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU+1)-COSO(2,NEU)
        D2    = SQRT(V1(1)**2+V1(2)**2)
        DISTD = (D1+D2)/2
        V1(1) = V1(1)/D2
        V1(2) = V1(2)/D2
        BET   = ACOS(V1(1))
        IF (V1(2).LE.0.) THEN
          BET = -BET
        ENDIF
        ALPD  = ALPHA(NEU)/2
        V2(1) = COSO(1,NEU+NBS1)-COSO(1,NEU)
        V2(2) = COSO(2,NEU+NBS1)-COSO(2,NEU)
        DIST  = SQRT(V2(1)**2+V2(2)**2)
        V2(1) = V2(1)/DIST
        V2(2) = V2(2)/DIST
C-----------------------------------------------------------------------
C    SI LE NOEUD CONSIDERE EST VOISIN D'UN COIN ON CALCULE L'ANGLE
C         MAXIMUM DUQUEL PEUT ETRE DEPLACE LE NOEUD INTERIEUR
C-----------------------------------------------------------------------
        IF (I.EQ.2) THEN
          ALPH = ALPHA(NEUB-1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D2**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D1**2+DIST**2-DISTI**2)/(2*D1*DIST))
              ALPM = ALPHA(NEUB)-ALPM
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
        IF (I.EQ.NBS1-1) THEN
          ALPH = ALPHA(NEUB+1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D2**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D2**2+DIST**2-DISTI**2)/(2*D2*DIST))
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C               CALCUL DE L'ANGLE REEL A CORRIGER
C-----------------------------------------------------------------------
        ALP   = BET+SENSFR*ALPD
        V1(1) = COS(ALP)
        V1(2) = SIN(ALP)
        PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
        PSCAL = AMIN1(PSCAL,1.)
        PSCAL = AMAX1(PSCAL,-1.)
        ALP   = ACOS(PSCAL)
C-----------------------------------------------------------------------
C     A LA PREMIERE ITERATION DE CORRECTION ON GARDE EN MEMOIRE LES
C                         DISTANCES IDEALES
C-----------------------------------------------------------------------
        IF (NITC.EQ.1) THEN
          IF (FADIST(1).GT.0) THEN
            DISTID(NEUB) = FADIST(1)*DISTD
          ELSE IF (FADIST(1).EQ.-1) THEN
            DISTID(NEUB) = DIST
            DMOY = 0.
            DO 11 JJ=1,NBS2-1
              NEU1 = NEU+(JJ-1)*NBS1
              NEU2 = NEU+JJ*NBS1
              DMOY = DMOY+SQRT( (COSO(1,NEU2)-COSO(1,NEU1))**2
     S                         +(COSO(2,NEU2)-COSO(2,NEU1))**2
     S                         +(COSO(3,NEU2)-COSO(3,NEU1))**2 )
   11       CONTINUE
            DMOY = DMOY/(NBS2-1)
            DISTID(NEUB) = MIN(DISTD,DMOY)
          ELSE IF (FADIST(1).EQ.-2) THEN
            D1 = 0.
            D2 = 0.
            DO 10 II=1,3
              D1 = D1+(COSO(II,NBS1+1)-COSO(II,1))**2
              D2 = D2+(COSO(II,2*NBS1)-COSO(II,NBS1))**2
   10       CONTINUE
            D1 = SQRT(D1)
            D2 = SQRT(D2)
C            DISTID(NEUB) = (NBS1-I)/(NBS1-1.)*D1+(I-1.)/(NBS1-1.)*D2
            DISTID(NEUB) = (1.E0-TRAV(NEU))*D1+TRAV(NEU)*D2
          ENDIF
        ENDIF
        PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
        IF (PSCAL.LT.0.) THEN
          ALP = -ALP
        ENDIF
C-----------------------------------------------------------------------
C                 CORRECTION SUR LES ANGLES AU BORD
C-----------------------------------------------------------------------
        IF (ABS(ALP).LT.ALPD) THEN
          FCI(NEU)= -SENSFR*0.5*ATAN(ALP)
        ENDIF
C-----------------------------------------------------------------------
C                CORRECTION SUR LES DISTANCES AU BORD
C-----------------------------------------------------------------------
        IF (FADIST(1).NE.0) THEN
          DISTD = DISTID(NEUB)
          IF (DIST.GT.DISTD) THEN
            FCJ(NEU)= -ATAN((DISTD-DIST)/DISTD)
          ELSE
            FCJ(NEU)= -(DISTD-DIST)/DISTD
          ENDIF
        ENDIF
   20 CONTINUE
C***********************************************************************
C                            BORD DROIT
C***********************************************************************
      DO 40 J=2,NBS2-1
        NEUB = NBS1+J-1
        NEU  = J*NBS1
        V1(1) = COSO(1,NEU-NBS1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU-NBS1)-COSO(2,NEU)
        D1    = SQRT(V1(1)**2+V1(2)**2)
        V1(1) = COSO(1,NEU+NBS1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU+NBS1)-COSO(2,NEU)
        D2    = SQRT(V1(1)**2+V1(2)**2)
        DISTD = (D1+D2)/2
        V1(1) = V1(1)/D2
        V1(2) = V1(2)/D2
        BET   = ACOS(V1(1))
        IF (V1(2).LE.0.) THEN
          BET = -BET
        ENDIF
        ALPD  = ALPHA(NEUB)/2
        V2(1) = COSO(1,NEU-1)-COSO(1,NEU)
        V2(2) = COSO(2,NEU-1)-COSO(2,NEU)
        DIST  = SQRT(V2(1)**2+V2(2)**2)
        V2(1) = V2(1)/DIST
        V2(2) = V2(2)/DIST
C-----------------------------------------------------------------------
C    SI LE NOEUD CONSIDERE EST VOISIN D'UN COIN ON CALCULE L'ANGLE
C         MAXIMUM DUQUEL PEUT ETRE DEPLACE LE NOEUD INTERIEUR
C-----------------------------------------------------------------------
        IF (J.EQ.2) THEN
          ALPH = ALPHA(NEUB-1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D1**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D1**2+DIST**2-DISTI**2)/(2*D1*DIST))
              ALPM = ALPHA(NEUB)-ALPM
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
        IF (J.EQ.NBS2-1) THEN
          ALPH = ALPHA(NEUB+1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D2**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D2**2+DIST**2-DISTI**2)/(2*D2*DIST))
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C               CALCUL DE L'ANGLE REEL A CORRIGER
C-----------------------------------------------------------------------
        ALP   = BET+SENSFR*ALPD
        V1(1) = COS(ALP)
        V1(2) = SIN(ALP)
        PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
        PSCAL = AMIN1(PSCAL,1.)
        PSCAL = AMAX1(PSCAL,-1.)
        ALP   = ACOS(PSCAL)
C-----------------------------------------------------------------------
C     A LA PREMIERE ITERATION DE CORRECTION ON GARDE EN MEMOIRE LES
C                         DISTANCES IDEALES
C-----------------------------------------------------------------------
        IF (NITC.EQ.1) THEN
          IF (FADIST(2).GT.0) THEN
            DISTID(NEUB) = FADIST(2)*DISTD
          ELSE IF (FADIST(2).EQ.-1) THEN
            DISTID(NEUB) = DIST
            DMOY = 0.
            DO 31 JJ=1,NBS1-1
              NEU1 = NEU-JJ+1
              NEU2 = NEU-JJ
              DMOY = DMOY+SQRT( (COSO(1,NEU2)-COSO(1,NEU1))**2
     S                         +(COSO(2,NEU2)-COSO(2,NEU1))**2
     S                         +(COSO(3,NEU2)-COSO(3,NEU1))**2 )
   31       CONTINUE
            DMOY = DMOY/(NBS1-1)
            DISTID(NEUB) = MIN(DISTD,DMOY)
          ELSE IF (FADIST(2).EQ.-2) THEN
            D1 = 0.
            D2 = 0.
            DO 30 II=1,3
              D1 = D1+(COSO(II,NBS1-1)-COSO(II,NBS1))**2
              D2 = D2+(COSO(II,NBS2*NBS1-1)-COSO(II,NBS2*NBS1))**2
   30       CONTINUE
            D1 = SQRT(D1)
            D2 = SQRT(D2)
C            DISTID(NEUB) = (NBS2-J)/(NBS2-1.)*D1+(J-1.)/(NBS2-1.)*D2
            DISTID(NEUB) = (1.E0-TRAV(NEU))*D1+TRAV(NEU)*D2
          ENDIF
        ENDIF
        PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
        IF (PSCAL.LT.0.) THEN
          ALP = -ALP
        ENDIF
C-----------------------------------------------------------------------
C                 CORRECTION SUR LES ANGLES AU BORD
C-----------------------------------------------------------------------
        IF (ABS(ALP).LT.ALPD) THEN
          FCJ(NEU)= -SENSFR*0.5*ATAN(ALP)
        ENDIF
C-----------------------------------------------------------------------
C                CORRECTION SUR LES DISTANCES AU BORD
C-----------------------------------------------------------------------
        IF (FADIST(2).NE.0) THEN
          DISTD = DISTID(NEUB)
          IF (DIST.GT.DISTD) THEN
            FCI(NEU)= ATAN((DISTD-DIST)/DISTD)
          ELSE
            FCI(NEU)= (DISTD-DIST)/DISTD
          ENDIF
        ENDIF
   40 CONTINUE
C***********************************************************************
C                            BORD DU HAUT
C***********************************************************************
      DO 60 I=NBS1-1,2,-1
        NEUB = 2*NBS1+NBS2-I-1
        NEU  = (NBS2-1)*NBS1+I
        V1(1) = COSO(1,NEU+1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU+1)-COSO(2,NEU)
        D1 = SQRT(V1(1)**2+V1(2)**2)
        V1(1) = COSO(1,NEU-1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU-1)-COSO(2,NEU)
        D2  = SQRT(V1(1)**2+V1(2)**2)
        DISTD = (D1+D2)/2
        V1(1) = V1(1)/D2
        V1(2) = V1(2)/D2
        BET   = ACOS(V1(1))
        IF (V1(2).LE.0.) THEN
          BET = -BET
        ENDIF
        ALPD  = ALPHA(NEUB)/2
        V2(1) = COSO(1,NEU-NBS1)-COSO(1,NEU)
        V2(2) = COSO(2,NEU-NBS1)-COSO(2,NEU)
        DIST  = SQRT(V2(1)**2+V2(2)**2)
        V2(1) = V2(1)/DIST
        V2(2) = V2(2)/DIST
C-----------------------------------------------------------------------
C    SI LE NOEUD CONSIDERE EST VOISIN D'UN COIN ON CALCULE L'ANGLE
C         MAXIMUM DUQUEL PEUT ETRE DEPLACE LE NOEUD INTERIEUR
C-----------------------------------------------------------------------
        IF (I.EQ.2) THEN
          ALPH = ALPHA(NEUB+1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D2**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D2**2+DIST**2-DISTI**2)/(2*D2*DIST))
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
        IF (I.EQ.NBS1-1) THEN
          ALPH = ALPHA(NEUB-1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D1**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D1**2+DIST**2-DISTI**2)/(2*D1*DIST))
              ALPM = ALPHA(NEUB)-ALPM
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C               CALCUL DE L'ANGLE REEL A CORRIGER
C-----------------------------------------------------------------------
        ALP   = BET+SENSFR*ALPD
        V1(1) = COS(ALP)
        V1(2) = SIN(ALP)
        PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
        PSCAL = AMIN1(PSCAL,1.)
        PSCAL = AMAX1(PSCAL,-1.)
        ALP   = ACOS(PSCAL)
C-----------------------------------------------------------------------
C     A LA PREMIERE ITERATION DE CORRECTION ON GARDE EN MEMOIRE LES
C                         DISTANCES IDEALES
C-----------------------------------------------------------------------
        IF (NITC.EQ.1) THEN
          IF (FADIST(3).GT.0) THEN
            DISTID(NEUB) = FADIST(3)*DISTD
          ELSE IF (FADIST(3).EQ.-1) THEN
            DISTID(NEUB) = DIST
            DMOY = 0.
            DO 51 JJ=1,NBS2-1
              NEU1 = NEU-(JJ-1)*NBS1
              NEU2 = NEU-JJ*NBS1
              DMOY = DMOY+SQRT( (COSO(1,NEU2)-COSO(1,NEU1))**2
     S                         +(COSO(2,NEU2)-COSO(2,NEU1))**2
     S                         +(COSO(3,NEU2)-COSO(3,NEU1))**2 )
   51       CONTINUE
            DMOY = DMOY/(NBS2-1)
            DISTID(NEUB) = MIN(DISTD,DMOY)
          ELSE IF (FADIST(3).EQ.-2) THEN
            D1 = 0.
            D2 = 0.
            DO 50 II=1,3
              D1 = D1+(COSO(II,(NBS2-2)*NBS1+1)
     S                -COSO(II,(NBS2-1)*NBS1+1))**2
              D2 = D2+(COSO(II,(NBS2-1)*NBS1)-COSO(II,NBS2*NBS1))**2
   50       CONTINUE
            D1 = SQRT(D1)
            D2 = SQRT(D2)
C            DISTID(NEUB) = (NBS1-I)/(NBS1-1.)*D1+(I-1.)/(NBS1-1.)*D2
            DISTID(NEUB) = (1.E0-TRAV(NEU))*D1+TRAV(NEU)*D2
          ENDIF
        ENDIF
        PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
        IF (PSCAL.LT.0.) THEN
          ALP = -ALP
        ENDIF
C-----------------------------------------------------------------------
C                 CORRECTION SUR LES ANGLES AU BORD
C-----------------------------------------------------------------------
        IF (ABS(ALP).LT.ALPD) THEN
          FCI(NEU)= SENSFR*0.5*ATAN(ALP)
        ENDIF
C-----------------------------------------------------------------------
C                CORRECTION SUR LES DISTANCES AU BORD
C-----------------------------------------------------------------------
        IF (FADIST(3).NE.0) THEN
          DISTD = DISTID(NEUB)
          IF (DIST.GT.DISTD) THEN
            FCJ(NEU)= ATAN((DISTD-DIST)/DISTD)
          ELSE
            FCJ(NEU)= (DISTD-DIST)/DISTD
          ENDIF
        ENDIF
   60 CONTINUE
C***********************************************************************
C                            BORD GAUCHE
C***********************************************************************
      DO 80 J=NBS2-1,2,-1
        NEUB = 2*NBS1+2*NBS2-J-2
        NEU  = (J-1)*NBS1+1
        V1(1) = COSO(1,NEU+NBS1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU+NBS1)-COSO(2,NEU)
        D1 = SQRT(V1(1)**2+V1(2)**2)
        V1(1) = COSO(1,NEU-NBS1)-COSO(1,NEU)
        V1(2) = COSO(2,NEU-NBS1)-COSO(2,NEU)
        D2  = SQRT(V1(1)**2+V1(2)**2)
        DISTD = (D1+D2)/2
        V1(1) = V1(1)/D2
        V1(2) = V1(2)/D2
        BET   = ACOS(V1(1))
        IF (V1(2).LE.0.) THEN
          BET = -BET
        ENDIF
        ALPD  = ALPHA(NEUB)/2
        V2(1) = COSO(1,NEU+1)-COSO(1,NEU)
        V2(2) = COSO(2,NEU+1)-COSO(2,NEU)
        DIST  = SQRT(V2(1)**2+V2(2)**2)
        V2(1) = V2(1)/DIST
        V2(2) = V2(2)/DIST
C-----------------------------------------------------------------------
C    SI LE NOEUD CONSIDERE EST VOISIN D'UN COIN ON CALCULE L'ANGLE
C         MAXIMUM DUQUEL PEUT ETRE DEPLACE LE NOEUD INTERIEUR
C-----------------------------------------------------------------------
        IF (J.EQ.2) THEN
          ALPH = ALPHA(1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D2**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D2**2+DIST**2-DISTI**2)/(2*D2*DIST))
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
        IF (J.EQ.NBS2-1) THEN
          ALPH = ALPHA(NEUB-1)
          IF (ALPH.LT.PI/2) THEN
            DELTA = DIST**2-D1**2*(SIN(ALPH))**2
            IF (DELTA.GT.0) THEN
              DISTI= COS(ALPH)*D2+SQRT(DELTA)
              ALPM = ACOS((D1**2+DIST**2-DISTI**2)/(2*D1*DIST))
              ALPM = ALPHA(NEUB)-ALPM
              ALPD = AMIN1(ALPD,ALPM)
            ENDIF
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C               CALCUL DE L'ANGLE REEL A CORRIGER
C-----------------------------------------------------------------------
        ALP   = BET+SENSFR*ALPD
        V1(1) = COS(ALP)
        V1(2) = SIN(ALP)
        PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
        PSCAL = AMIN1(PSCAL,1.)
        PSCAL = AMAX1(PSCAL,-1.)
        ALP   = ACOS(PSCAL)
C-----------------------------------------------------------------------
C     A LA PREMIERE ITERATION DE CORRECTION ON GARDE EN MEMOIRE LES
C                         DISTANCES IDEALES
C-----------------------------------------------------------------------
        IF (NITC.EQ.1) THEN
          IF (FADIST(4).GT.0) THEN
            DISTID(NEUB) = FADIST(4)*DISTD
          ELSE IF (FADIST(4).EQ.-1) THEN
            DISTID(NEUB) = DIST
            DMOY = 0.
            DO 71 JJ=1,NBS1-1
              NEU1 = NEU+JJ-1
              NEU2 = NEU+JJ
              DMOY = DMOY+SQRT( (COSO(1,NEU2)-COSO(1,NEU1))**2
     S                         +(COSO(2,NEU2)-COSO(2,NEU1))**2
     S                         +(COSO(3,NEU2)-COSO(3,NEU1))**2 )
   71       CONTINUE
            DMOY = DMOY/(NBS1-1)
            DISTID(NEUB) = MIN(DISTD,DMOY)
          ELSE IF (FADIST(4).EQ.-2) THEN
            D1 = 0.
            D2 = 0.
            DO 70 II=1,3
              D1 = D1+(COSO(II,2)-COSO(II,1))**2
              D2 = D2+(COSO(II,(NBS2-1)*NBS1+2)
     S                -COSO(II,(NBS2-1)*NBS1+1))**2
   70       CONTINUE
            D1 = SQRT(D1)
            D2 = SQRT(D2)
C            DISTID(NEUB) = (NBS2-J)/(NBS2-1.)*D1+(J-1.)/(NBS2-1.)*D2
            DISTID(NEUB) = (1.E0-TRAV(NEU))*D1+TRAV(NEU)*D2
          ENDIF
        ENDIF
        PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
        IF (PSCAL.LT.0.) THEN
          ALP = -ALP
        ENDIF
C-----------------------------------------------------------------------
C                 CORRECTION SUR LES ANGLES AU BORD
C-----------------------------------------------------------------------
        IF (ABS(ALP).LT.ALPD) THEN
          FCJ(NEU)= SENSFR*0.5*ATAN(ALP)
        ENDIF
C-----------------------------------------------------------------------
C                CORRECTION SUR LES DISTANCES AU BORD
C-----------------------------------------------------------------------
        IF (FADIST(4).NE.0) THEN
          DISTD = DISTID(NEUB)
          IF (DIST.GT.DISTD) THEN
            FCI(NEU)= -ATAN((DISTD-DIST)/DISTD)
          ELSE
            FCI(NEU)= -(DISTD-DIST)/DISTD
          ENDIF
        ENDIF
   80 CONTINUE
C***********************************************************************
C      CALCUL DES FONCTIONS DE CONTROLE DANS TOUT LE DOMAINE
C***********************************************************************
       RAPP = 0.1
       COD = 0.1
       COA = 0.4
C-----------------------------------------------------------------------
C                INTERPOLATION PAR EXPONENTIELLES
C-----------------------------------------------------------------------
      DO 100 J=2,NBS2-1
        DO 90 I=2,NBS1-1
          NEU = (J-1)*NBS1+I
          NEUH = (NBS2-1)*NBS1+I
          NEUB = I
          NEUG = (J-1)*NBS1+1
          NEUD = J*NBS1
          SI     = AMIN1(0.25,((J-0.5)/(NBS1-1.)))
          SI     = AMIN1(SI,((NBS2-J+0.5)/(NBS1-1.)))
          SJ     = AMIN1(0.25,((I-0.5)/(NBS2-1.)))
          SJ     = AMIN1(SJ,((NBS1-I+0.5)/(NBS2-1.)))
          DLONDI = -ALOG(RAPP)/SI
          DLONAI = -ALOG(RAPP)/SJ
          DLONDJ = -ALOG(RAPP)/SJ
          DLONAJ = -ALOG(RAPP)/SI
          FCI(NEU) = FCI(NEU)
     S              +COD*EXP(-DLONDI*(I-1.)/(NBS1-1.))  *FCI(NEUG)
     S              +COD*EXP(-DLONDI*(NBS1-I)/(NBS1-1.))*FCI(NEUD)
     S              +COA*EXP(-DLONAI*(J-1.)/(NBS2-1.))  *FCI(NEUB)
     S              +COA*EXP(-DLONAI*(NBS2-J)/(NBS2-1.))*FCI(NEUH)
          FCJ(NEU) = FCJ(NEU)
     S              +COA*EXP(-DLONAJ*(I-1.)/(NBS1-1.))  *FCJ(NEUG)
     S              +COA*EXP(-DLONAJ*(NBS1-I)/(NBS1-1.))*FCJ(NEUD)
     S              +COD*EXP(-DLONDJ*(J-1.)/(NBS2-1.))  *FCJ(NEUB)
     S              +COD*EXP(-DLONDJ*(NBS2-J)/(NBS2-1.))*FCJ(NEUH)
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
