      SUBROUTINE VERGEO(NBS1,NBS2,COSO,SENSFR,NBELTR,NBELTD,NELT)
C***********************************************************************
C BUT :  VERIFICATION DE LA GEOMETRIE DE L ENSEMBLE DES MAILLES
C***********************************************************************
C
C ENTREES:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DE TOUS LES POINTS DU MAILLAGE
C           SENSFR : SENS DE ROTATION DE LA FRONTIERE
C                    -1 : SENS INDIRECT
C                     0 : FRONTIERE DEGENEREE
C                     1 : SENS DIRECT
C           NBELTR : -1 N'ETABLIT PAS LA LISTE DES ELTS RET (NELT)
C
C SORTIES:
C           NBELTR : NOMBRE D'ELEMENTS RETOURNES DU MAILLAGE
C           NBELTD : NOMBRE D'ELEMENTS DEGENERES DU MAILLAGE
C           NELT   : LISTE DES ELEMENTS RETOURNES DU MAILLAGE
C                    SEULEMENT SI NBELTR <> -1 EN ENTREE
C                    (1 A NBELTR) SUIVIE DE LA LISTE DES ELEMENTS
C                    DEGENERES (NBPEXT+1 A NBPEXT+NBPFRT)
C
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DIMENSION COSO(3,NBS1*NBS2),NELT(NBS1*NBS2)
      DIMENSION COIN(6,2),V1(2),V2(2)
      INTEGER SENSRO,SENSFR
C
      PI = 4.0E0*ATAN(1.E0)
C***********************************************************************
C                RECHERCHE DES ELEMENTS RETOURNES
C      <QUI NE TOURNENT PAS DANS LE MEME SENS QUE LA FRONTIERE>
C          ET DES ELEMENTS DEGENERES QUI SONT DU TYPE
C                            *-----*
C                             \   /
C                              \ /
C                               X
C                              / \
C                             *---*
C***********************************************************************
      IF (NBELTR.EQ.-1) THEN
        IREMP = 0
      ELSE
        IREMP = 1
      ENDIF
      NBELTR = 0
      NBELTD = 0
      SENSRO = SENSFR
      DO 40 J=1,NBS2-1
        DO 30 I=1,NBS1-1
          NUMELT = (J-1)*(NBS1-1)+I
          NEU    = (J-1)*NBS1+I
C-----------------------------------------------------------------------
C     REMPLISSAGE DU TABLEAU COIN POUR APPELER LA FONCTION DE CALCUL
C               DU SENS DE ROTATION DE L ELEMENT
C-----------------------------------------------------------------------
          COIN(1,1) = COSO(1,NEU+NBS1)
          COIN(1,2) = COSO(2,NEU+NBS1)
          COIN(2,1) = COSO(1,NEU)
          COIN(2,2) = COSO(2,NEU)
          COIN(3,1) = COSO(1,NEU+1)
          COIN(3,2) = COSO(2,NEU+1)
          COIN(4,1) = COSO(1,NEU+NBS1+1)
          COIN(4,2) = COSO(2,NEU+NBS1+1)
          COIN(5,1) = COSO(1,NEU+NBS1)
          COIN(5,2) = COSO(2,NEU+NBS1)
          COIN(6,1) = COSO(1,NEU)
          COIN(6,2) = COSO(2,NEU)
C-----------------------------------------------------------------------
C CALCUL DE LA SOMME DES ANGLES AUX QUATRES COINS. SI ELLE EST EGALE A
C    *  4*PI/2   LE BORD DE L ELEMENT TOURNE DANS LE SENS DIRECT
C    *  6*PI/2   LE BORD DE L ELEMENT EST DEGENERE
C    *  8*PI/2   LE BORD DE L ELEMENT TOURNE DANS LE SENS INDIRECT.
C-----------------------------------------------------------------------
          SOMALP = 0.0E0
          DO 10 II=2,5
            V1(1) = COIN(II+1,1)-COIN(II,1)
            V1(2) = COIN(II+1,2)-COIN(II,2)
            V2(1) = COIN(II-1,1)-COIN(II,1)
            V2(2) = COIN(II-1,2)-COIN(II,2)
            DIST  = SQRT(V1(1)**2+V1(2)**2)
            IF (DIST.NE.0.0E0) THEN
              V1(1) = V1(1)/DIST
              V1(2) = V1(2)/DIST
            ELSE
              NBELTD = NBELTD+1
              IF (IREMP.EQ.1) THEN
                NELT(NBELTR+NBELTD) = NUMELT
              ENDIF
              GOTO 30
            ENDIF
            DIST  = SQRT(V2(1)**2+V2(2)**2)
            IF (DIST.NE.0.0E0) THEN
              V2(1) = V2(1)/DIST
              V2(2) = V2(2)/DIST
            ELSE
              NBELTD = NBELTD+1
              IF (IREMP.EQ.1) THEN
                NELT(NBELTR+NBELTD) = NUMELT
              ENDIF
              GOTO 30
            ENDIF
            PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
            PSCAL = AMIN1(PSCAL,1.0E0)
            PSCAL = AMAX1(PSCAL,-1.0E0)
            ALP   = ACOS(PSCAL)
            PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
            IF (PSCAL.LT.0.0E0) THEN
              ALP = 2.0E0*PI-ALP
            ENDIF
            SOMALP = SOMALP+ALP
   10     CONTINUE
          IF (ABS(SOMALP-2.0E0*PI).LT.PI) THEN
            SENSRO = 1
          ELSE
            IF (ABS(SOMALP-6.0E0*PI).LT.PI) THEN
              SENSRO = -1
            ELSE
              SENSRO = 0
            ENDIF
          ENDIF
          IF (SENSRO.NE.0) THEN
            IF (SENSRO.NE.SENSFR) THEN
              IF (IREMP.EQ.1) THEN
                DO 20 K=NBELTD,1,-1
                  NELT(NBELTR+K+1) = NELT(NBELTR+K)
   20           CONTINUE
              ENDIF
              NBELTR = NBELTR+1
              IF (IREMP.EQ.1) THEN
                NELT(NBELTR) = NUMELT
              ENDIF
            ENDIF
          ELSE
            NBELTD = NBELTD+1
            IF (IREMP.EQ.1) THEN
              NELT(NBELTR+NBELTD) = NUMELT
            ENDIF
          ENDIF
   30   CONTINUE
   40 CONTINUE
C***********************************************************************
C           IMPRESSION DU NOMBRE D ELEMENTS 'ANORMAUX'
C***********************************************************************
      IF ((NBELTR.NE.0).OR.(NBELTD.NE.0)) THEN
        WRITE (IMPRIM,*) ' '
        WRITE (IMPRIM,1000)
        WRITE (IMPRIM,1001)
        WRITE (IMPRIM,1002) NBELTR
        WRITE (IMPRIM,1003) NBELTD
        WRITE (IMPRIM,1000)
        WRITE (IMPRIM,*) ' '
      ENDIF
1000  FORMAT('**************************************')
1001  FORMAT('*       TOPOLOGIE DES MAILLES        *')
1002  FORMAT('*   ',I4,' ELEMENT(S) RETOURNE(S)      *')
1003  FORMAT('*   ',I4,' ELEMENT(S) DEGENERE(S)      *')
      RETURN
      END
