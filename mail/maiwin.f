      SUBROUTINE MAIWIN(NBS1,NBS2,NBS3,COSO,DISTID,FCI,FCJ,FCK,NITCOR,
     S                  FADIST,LGMAX,TCPU,IERR)
C***********************************************************************
C BUT: MAILLER UN VOLUME L ENSEMBLE DES COORDONNEES ETANT INITIALISEES
C***********************************************************************
C
C ENTREES:
C           NBS1    : NOMBRE DE POINTS SUR LIGNES X
C           NBS2    : NOMBRE DE POINTS SUR LIGNES Y
C           NBS3    : NOMBRE DE POINTS SUR LIGNES Z
C           LGMAX   : LONGUEUR MAXIMALE DES VECTEURS
C
C ENTREES ET RESULTAT:
C           COSO  : COORDONNEES DES POINTS DU MAILLAGE
C
C TRAVAIL:
C           DISTID: VECTEUR DE TRAVAIL DES DISTANCES IDEALES
C           FCI   : VECTEUR DES FONCTIONS DE CONTROLE EN I
C           FCJ   : VECTEUR DES FONCTIONS DE CONTROLE EN J
C           FCK   : VECTEUR DES FONCTIONS DE CONTROLE EN K
C***********************************************************************
C AUTEUR : DOURSAT CHRISTOPHE  ANALYSE NUMERIQUE PARIS  OCTOBRE 1989
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C***********************************************************************
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DIMENSION  COSO(3,LGMAX)
      DIMENSION  DISTID(LGMAX),FCI(LGMAX),FCJ(LGMAX),FCK(LGMAX)
      DIMENSION  DEL(3,3),FADIST(6)
      REAL       MINEUR(3,3)
      DOUBLE PRECISION DINFO,TCPU(2,100),TCPUT
      DIMENSION  NBIT(100)
C***********************************************************************
C                   INITIALISATION DES VARIABLES
C***********************************************************************
      DATA TESTL1/0.01/
      TESTL2 = 0.05
      NBASE  = NBS1*NBS2
      NTCPU = 0
      TCPUT = 0
C***********************************************************************
C                ENTREE DES DONNES DU CALCUL
C***********************************************************************
      IF (NITCOR.NE.0) THEN
C        WRITE (*,*) 'PROFONDEUR DE CORRECTION ? ENTIER > 5'
C        READ  (*,*) PROF
         PROF = 25
      ENDIF
C-----------------------------------------------------------------------
C         DEFINITION DES FORMATS NECESSAIRES AUX IMPRESSIONS
C-----------------------------------------------------------------------
 1010 FORMAT('**************************************')
C***********************************************************************
C                        INITIALISATION
C***********************************************************************
      NITC   = 1
      DO 30 K=1,NBS3
        DO 20 J=1,NBS2
          DO 10 I=1,NBS1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            FCI(NEU) = 0.
            FCJ(NEU) = 0.
            FCK(NEU) = 0.
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C***********************************************************************
C         CALCUL DE LA LONGUEUR MOYENNE MINIMALE DES MAILLES
C                  SUR LE BORD POUR UN PREMIER CALCUL
C***********************************************************************
      DMOY = 10000.
C-----------------------------------------------------------------------
C                       FACES A K CONSTANT
C-----------------------------------------------------------------------
      DO 60 K=1,NBS3,NBS3-1
        DO 50 J=1,NBS2
          DO 40 I=1,NBS1-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+1)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+1)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+1)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
   40     CONTINUE
   50   CONTINUE
        DO 55 I=1,NBS1
          DO 45 J=1,NBS2-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+NBS1)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+NBS1)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+NBS1)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
   45     CONTINUE
   55   CONTINUE
   60 CONTINUE
C-----------------------------------------------------------------------
C                       FACES A J CONSTANT
C-----------------------------------------------------------------------
      DO 90 J=1,NBS2,NBS2-1
        DO 80 K=1,NBS3
          DO 70 I=1,NBS1-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+1)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+1)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+1)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
   70     CONTINUE
   80   CONTINUE
        DO 85 I=1,NBS1
          DO 75 K=1,NBS3-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+NBASE)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+NBASE)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+NBASE)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
   75     CONTINUE
   85   CONTINUE
   90 CONTINUE
C-----------------------------------------------------------------------
C                       FACES A I CONSTANT
C-----------------------------------------------------------------------
      DO 120 I=1,NBS1,NBS1-1
        DO 110 K=1,NBS3
          DO 100 J=1,NBS2-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+NBS1)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+NBS1)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+NBS1)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
  100     CONTINUE
  110   CONTINUE
        DO 115 J=1,NBS2
          DO 105 K=1,NBS3-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            DIST   = SQRT( (COSO(1,NEU+NBASE)-COSO(1,NEU))**2
     S                    +(COSO(2,NEU+NBASE)-COSO(2,NEU))**2
     S                    +(COSO(3,NEU+NBASE)-COSO(3,NEU))**2 )
            DMOY = AMIN1(DMOY,DIST)
  105     CONTINUE
  115   CONTINUE
  120 CONTINUE
      NDMOY = 0
C     WRITE (*,*) 'DMOY CALCULEE AU DEPART ',DMOY
      DO 533 K=2,NBS3-1
        DO 532 J=2,NBS2-1
          DO 531 I=2,NBS1-1
            NEU    = (K-1)*NBASE+(J-1)*NBS1+I
            NEUB   =                (J-1)*NBS1+I
            NEUH   = (NBS3-1)*NBASE +(J-1)*NBS1   +I
            NEUG   =    (K-1)*NBASE               +I
            NEUD   =    (K-1)*NBASE +(NBS2-1)*NBS1+I
            NEUF   =    (K-1)*NBASE +(J-1)*NBS1
            NEUC   =    (K-1)*NBASE +(J-1)*NBS1   +NBS1
            DO 534 KK=1,3
              COSO(KK,NEU) =    (K-1.)/(NBS3-1.)*COSO(KK,NEUH)
     S                       +(NBS3-K)/(NBS3-1.)*COSO(KK,NEUB)
     S                       +  (J-1.)/(NBS2-1.)*COSO(KK,NEUD)
     S                       +(NBS2-J)/(NBS2-1.)*COSO(KK,NEUG)
     S                       +  (I-1.)/(NBS1-1.)*COSO(KK,NEUC)
     S                       +(NBS1-I)/(NBS1-1.)*COSO(KK,NEUF)
              COSO(KK,NEU) = COSO(KK,NEU) / 3.0
  534       CONTINUE
  531     CONTINUE
  532   CONTINUE
  533 CONTINUE
      NTCPU         = NTCPU+1
      TCPU(1,NTCPU) = DINFO('DELTA CPU')
      TCPUT         = TCPU(1,NTCPU)+TCPUT
      TCPU(2,NTCPU) = TCPUT
C***********************************************************************
C
C                              RESOLUTION
C
C***********************************************************************
      LM = 0
      DO 200 L=1,400
        TEST = 0.
        DO 160 K=2,NBS3-1
          DO 150 J=2,NBS2-1
            DO 140 I=2,NBS1-1
              NE = (K-1)*NBASE+(J-1)*NBS1+I
C***********************************************************************
C           CALCUL INTERMEDIAIRE DU JACOBIEN APPELE DEL(3,3)
C***********************************************************************
              DEL(1,1)=COSO(1,NE+1)-COSO(1,NE-1)
              DEL(2,1)=COSO(1,NE+NBS1)-COSO(1,NE-NBS1)
              DEL(3,1)=COSO(1,NE+NBASE)-COSO(1,NE-NBASE)
              DEL(1,2)=COSO(2,NE+1)-COSO(2,NE-1)
              DEL(2,2)=COSO(2,NE+NBS1)-COSO(2,NE-NBS1)
              DEL(3,2)=COSO(2,NE+NBASE)-COSO(2,NE-NBASE)
              DEL(1,3)=COSO(3,NE+1)-COSO(3,NE-1)
              DEL(2,3)=COSO(3,NE+NBS1)-COSO(3,NE-NBS1)
              DEL(3,3)=COSO(3,NE+NBASE)-COSO(3,NE-NBASE)
C***********************************************************************
C                 CALCUL INTERMEDIAIRE DES MINEURS
C***********************************************************************
              MINEUR(1,1)=DEL(2,2)*DEL(3,3)-DEL(3,2)*DEL(2,3)
              MINEUR(2,1)=DEL(1,2)*DEL(3,3)-DEL(3,2)*DEL(1,3)
              MINEUR(3,1)=DEL(1,2)*DEL(2,3)-DEL(2,2)*DEL(1,3)
              MINEUR(1,2)=DEL(2,1)*DEL(3,3)-DEL(3,1)*DEL(2,3)
              MINEUR(2,2)=DEL(1,1)*DEL(3,3)-DEL(3,1)*DEL(1,3)
              MINEUR(3,2)=DEL(1,1)*DEL(2,3)-DEL(2,1)*DEL(1,3)
              MINEUR(1,3)=DEL(2,1)*DEL(3,2)-DEL(3,1)*DEL(2,2)
              MINEUR(2,3)=DEL(1,1)*DEL(3,2)-DEL(3,1)*DEL(1,2)
              MINEUR(3,3)=DEL(1,1)*DEL(2,2)-DEL(2,1)*DEL(1,2)
C***********************************************************************
C                 CALCUL DES COEFFICIENTS DE LA MATRICE
C***********************************************************************
              CX2=(MINEUR(1,1)**2+MINEUR(1,2)**2+MINEUR(1,3)**2)/16
              CY2=(MINEUR(2,1)**2+MINEUR(2,2)**2+MINEUR(2,3)**2)/16
              CZ2=(MINEUR(3,1)**2+MINEUR(3,2)**2+MINEUR(3,3)**2)/16
              CXY=MINEUR(1,1)*MINEUR(2,1)+MINEUR(1,2)*MINEUR(2,2)+
     %            MINEUR(1,3)*MINEUR(2,3)
              CXY=CXY/32
              CXZ=MINEUR(1,1)*MINEUR(3,1)+MINEUR(1,2)*MINEUR(3,2)+
     %            MINEUR(1,3)*MINEUR(3,3)
              CXZ=CXZ/32
              CYZ=MINEUR(2,1)*MINEUR(3,1)+MINEUR(2,2)*MINEUR(3,2)+
     %            MINEUR(2,3)*MINEUR(3,3)
              CYZ=CYZ/32
              DIAG = -2*(CX2+CY2+CZ2)
C***********************************************************************
C                 CALCUL DES VALEURS INTERMEDIAIRES
C***********************************************************************
              XINT= CYZ*COSO(1,NE-NBASE-NBS1)-CXZ*COSO(1,NE-NBASE-1)
     S             -CZ2*COSO(1,NE-NBASE)   +CXZ*COSO(1,NE-NBASE+1)
     S             -CYZ*COSO(1,NE-NBASE+NBS1)
     S             +CXY*COSO(1,NE-NBS1-1)    -CY2*COSO(1,NE-NBS1)
     S             -CXY*COSO(1,NE-NBS1+1)    -CX2*COSO(1,NE-1)
     S             -CX2*COSO(1,NE+1)         -CXY*COSO(1,NE+NBS1-1)
     S             -CY2*COSO(1,NE+NBS1)      +CXY*COSO(1,NE+NBS1+1)
     S             -CYZ*COSO(1,NE+NBASE-NBS1)
     S             +CXZ*COSO(1,NE+NBASE-1)   -CZ2*COSO(1,NE+NBASE)
     S             -CXZ*COSO(1,NE+NBASE+1)   +CYZ*COSO(1,NE+NBASE+NBS1)
     S             +CX2*(COSO(1,NE+1)-COSO(1,NE-1))        *FCI(NE)/2
     S             +CY2*(COSO(1,NE+NBS1)-COSO(1,NE-NBS1))  *FCJ(NE)/2
     S             +CZ2*(COSO(1,NE+NBASE)-COSO(1,NE-NBASE))*FCK(NE)/2
              YINT= CYZ*COSO(2,NE-NBASE-NBS1)-CXZ*COSO(2,NE-NBASE-1)
     S             -CZ2*COSO(2,NE-NBASE)     +CXZ*COSO(2,NE-NBASE+1)
     S             -CYZ*COSO(2,NE-NBASE+NBS1)
     S             +CXY*COSO(2,NE-NBS1-1)    -CY2*COSO(2,NE-NBS1)
     S             -CXY*COSO(2,NE-NBS1+1)    -CX2*COSO(2,NE-1)
     S             -CX2*COSO(2,NE+1)         -CXY*COSO(2,NE+NBS1-1)
     S             -CY2*COSO(2,NE+NBS1)      +CXY*COSO(2,NE+NBS1+1)
     S             -CYZ*COSO(2,NE+NBASE-NBS1)
     S             +CXZ*COSO(2,NE+NBASE-1)   -CZ2*COSO(2,NE+NBASE)
     S             -CXZ*COSO(2,NE+NBASE+1)   +CYZ*COSO(2,NE+NBASE+NBS1)
     S             +CX2*(COSO(2,NE+1)-COSO(2,NE-1))        *FCI(NE)/2
     S             +CY2*(COSO(2,NE+NBS1)-COSO(2,NE-NBS1))  *FCJ(NE)/2
     S             +CZ2*(COSO(2,NE+NBASE)-COSO(2,NE-NBASE))*FCK(NE)/2
              ZINT= CYZ*COSO(3,NE-NBASE-NBS1)-CXZ*COSO(3,NE-NBASE-1)
     S             -CZ2*COSO(3,NE-NBASE)     +CXZ*COSO(3,NE-NBASE+1)
     S             -CYZ*COSO(3,NE-NBASE+NBS1)
     S             +CXY*COSO(3,NE-NBS1-1)    -CY2*COSO(3,NE-NBS1)
     S             -CXY*COSO(3,NE-NBS1+1)    -CX2*COSO(3,NE-1)
     S             -CX2*COSO(3,NE+1)         -CXY*COSO(3,NE+NBS1-1)
     S             -CY2*COSO(3,NE+NBS1)      +CXY*COSO(3,NE+NBS1+1)
     S             -CYZ*COSO(3,NE+NBASE-NBS1)
     S             +CXZ*COSO(3,NE+NBASE-1)   -CZ2*COSO(3,NE+NBASE)
     S             -CXZ*COSO(3,NE+NBASE+1)   +CYZ*COSO(3,NE+NBASE+NBS1)
     S             +CX2*(COSO(3,NE+1)-COSO(3,NE-1))        *FCI(NE)/2
     S             +CY2*(COSO(3,NE+NBS1)-COSO(3,NE-NBS1))  *FCJ(NE)/2
     S             +CZ2*(COSO(3,NE+NBASE)-COSO(3,NE-NBASE))*FCK(NE)/2
              XINT = XINT/DIAG
              YINT = YINT/DIAG
              ZINT = ZINT/DIAG
              TEST=MAX(TEST,ABS(XINT-COSO(1,NE)))
              COSO(1,NE) = XINT
              TEST=MAX(TEST,ABS(YINT-COSO(2,NE)))
              COSO(2,NE) = YINT
              TEST=MAX(TEST,ABS(ZINT-COSO(3,NE)))
              COSO(3,NE) = ZINT
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
C***********************************************************************
C            AIGUILLAGE SELON LE TEST DE SORTIE DU CALCUL
C***********************************************************************
        TEST = TEST/DMOY
CCC        WRITE (IMPRIM,1030) L,TEST
        IF (TEST.LE.TESTL2) THEN
C-----------------------------------------------------------------------
C                   ON EST A LA PRECISION SOUHAITEE
C-----------------------------------------------------------------------
          IF (NDMOY.EQ.0) THEN
            NDMOY=1
C***********************************************************************
C         CALCUL DE LA LONGUEUR MOYENNE MINIMALE DES MAILLES
C                    SUR LE MAILLAGE COOK
C***********************************************************************
            DMOY = 10000.
C-----------------------------------------------------------------------
C                           LIGNES DES I
C-----------------------------------------------------------------------
            DO 760 K=2,NBS3-1
              DO 750 J=2,NBS2-1
                DO 740 I=1,NBS1-1
                  NEU    = (K-1)*NBASE+(J-1)*NBS1+I
                  DIST   = SQRT( (COSO(1,NEU+1)-COSO(1,NEU))**2
     S                          +(COSO(2,NEU+1)-COSO(2,NEU))**2
     S                          +(COSO(3,NEU+1)-COSO(3,NEU))**2 )
                  DMOY = AMIN1(DMOY,DIST)
  740           CONTINUE
  750         CONTINUE
  760       CONTINUE
C-----------------------------------------------------------------------
C                           LIGNES DES J
C-----------------------------------------------------------------------
            DO 790 K=2,NBS3-1
              DO 780 I=2,NBS1-1
                DO 770 J=1,NBS2-1
                  NEU    = (K-1)*NBASE+(J-1)*NBS1+I
                  DIST   = SQRT( (COSO(1,NEU+NBS1)-COSO(1,NEU))**2
     S                          +(COSO(2,NEU+NBS1)-COSO(2,NEU))**2
     S                          +(COSO(3,NEU+NBS1)-COSO(3,NEU))**2 )
                  DMOY = AMIN1(DMOY,DIST)
  770           CONTINUE
  780         CONTINUE
  790       CONTINUE
C-----------------------------------------------------------------------
C                           LIGNES DES K
C-----------------------------------------------------------------------
            DO 720 J=2,NBS2-1
              DO 710 I=2,NBS1-1
                DO 700 K=1,NBS3-1
                  NEU    = (K-1)*NBASE+(J-1)*NBS1+I
                  DIST   = SQRT( (COSO(1,NEU+NBASE)-COSO(1,NEU))**2
     S                          +(COSO(2,NEU+NBASE)-COSO(2,NEU))**2
     S                          +(COSO(3,NEU+NBASE)-COSO(3,NEU))**2 )
                  DMOY = AMIN1(DMOY,DIST)
  700           CONTINUE
  710         CONTINUE
  720       CONTINUE
C           WRITE (*,*) ' '
C           WRITE (*,*) 'CHANGEMENT DE TEST '
C           WRITE (*,*) '    DMOY CALCULEE ',DMOY
C           WRITE (*,*) ' '
            NTCPU         = NTCPU+1
            TCPU(1,NTCPU) = DINFO('DELTA CPU')
            TCPUT         = TCPU(1,NTCPU)+TCPUT
            TCPU(2,NTCPU) = TCPUT
            NBIT(NTCPU) = L
            NBITT       = L
          ELSE
            WRITE (IMPRIM,1031) L-LM
            LM = L
            NTCPU         = NTCPU+1
            TCPU(1,NTCPU) = DINFO('DELTA CPU')
            TCPUT         = TCPU(1,NTCPU)+TCPUT
            TCPU(2,NTCPU) = TCPUT
            NBIT(NTCPU)   = L-NBITT
            NBITT         = L
C             WRITE (*,*) 'VOULEZ VOUS SORTIR ? 0 NON 1 OUI'
C             READ (*,*) LL
C             IF (LL.EQ.1) THEN
C               RETURN
C             ENDIF
            IF (NITC.LE.NITCOR) THEN
C-----------------------------------------------------------------------
C            ON AMELIORE LA GEOMETRIE DES MAILLES VOISINES DU BORD
C-----------------------------------------------------------------------
              WRITE (IMPRIM,1032) NITC
              CALL AMELMA(NBS1,NBS2,NBS3,COSO,NITC,
     S                    FCI,FCJ,FCK,FADIST,DISTID,LGMAX,PROF)
              NITC = NITC+1
            ELSE
              IF (TESTL2.NE.TESTL1) THEN
C-----------------------------------------------------------------------
C       L AMELIORATION DES MAILLES VOISINES DU BORD EST EFFECTUEE
C       ON RENVOIE LE CALCUL AVEC UNE PRECISION FINALE PLUS FAIBLE
C-----------------------------------------------------------------------
                TESTL2 = TESTL1
              ELSE
C-----------------------------------------------------------------------
C           AFFICHAGE A L'ECRAN DES RESULATS DE CONVERGENCE
C-----------------------------------------------------------------------
CCC                WRITE (IMPRIM,1010)
                WRITE (IMPRIM,1040) L
CCC                WRITE (IMPRIM,1050)
CCC                WRITE (IMPRIM,1041)
                WRITE (IMPRIM,1042) TCPUT
CCC                WRITE (IMPRIM,1010)
C                DO 175 NC=0,2
C                 IF (NITCOR.EQ.0) THEN
C                   METHOD = 1
C                 ELSE
C                   METHOD = 2
C                 ENDIF
C                 CALL QUAL3D(PROF,METHOD,NC,NBS1,NBS2,NBS3,
C    S                        COSO,FCI,FCJ,FCK)
C 175           CONTINUE
                TCPU(2,1) = TCPUT
                RETURN
              ENDIF
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
C***********************************************************************
C         DEFINITION DES FORMATS NECESSAIRES AUX IMPRESSIONS
C***********************************************************************
CCC 1030 FORMAT('   BOUCLE ',I5,' TEST ',1PE10.3)
 1031 FORMAT(' NOMBRE DE PROBLEMES LINEARISES RESOLUS',T45,':',I5)
 1032 FORMAT(' ITERATION DE CORRECTION DE L''ORTHOGONALITE',T45,
     S             ':',I5)
 1040 FORMAT('CONVERGENCE EN ',I4,' ITERATIONS')
 1042 FORMAT('TEMPS CPU : ',1PE11.3)
CCC 1040 FORMAT('*   CONVERGENCE EN',I4,' ITERATIONS    *')
CCC 1041 FORMAT('*         TEMPS CPU UTILISE          *')
CCC 1042 FORMAT('*           ',1PE11.3,'              *')
CCC 1043 FORMAT('* PHASE  * ITERA *   TCPU  * T.CUMUL.*')
CCC 1044 FORMAT('*  COOK  *       *',1PE8.1,' *',1PE8.1,' *')
CCC 1045 FORMAT('* INITIA *       *',1PE8.1,' *',1PE8.1,' *')
CCC 1046 FORMAT('* CALC.1 * ',I4,'  *',1PE8.1,' *',1PE8.1,' *')
CCC 1047 FORMAT('* CALC.2 * ',I4,'  *',1PE8.1,' *',1PE8.1,' *')
CCC 1048 FORMAT('* AMELIO * ',I4,'  *',1PE8.1,' *',1PE8.1,' *')
CCC 1049 FORMAT('* CAL FI * ',I4,'  *',1PE8.1,' *',1PE8.1,' *')
CCC 1050 FORMAT('*                                    *')
C***********************************************************************
C TOUTES LES BOUCLES ONT ETE EFFECTUEES DONC LA METHODE N A PAS CONVERGE
C                     SORTIE AVEC MESSAGE D ERREUR
C***********************************************************************
      IERR = 1
      WRITE (IMPRIM,1010)
      WRITE (IMPRIM,1060)
      WRITE (IMPRIM,1061)
      WRITE (IMPRIM,1062)
      WRITE (IMPRIM,1063)
      WRITE (IMPRIM,1010)
 1060 FORMAT('*             ATTENTION              *')
 1061 FORMAT('*  PROBLEME DANS LA ROUTINE MAIWIN   *')
 1062 FORMAT('*   TOUTES LES ITERATIONS ONT ETE    *')
 1063 FORMAT('*    EFFECTUEES SANS CONVERGENCE     *')
      RETURN
      END
