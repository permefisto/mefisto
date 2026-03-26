        SUBROUTINE SUEX07( NBS1,NBS2,LADEFI,RADEFI,NBRPRO,
     S                     LIMITE,CEN,VECN,RAYON,COSO, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INITIALISATION PROJECTION SUR CERCLE
C -----
C ENTREE:
C -------
C NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C TRAVAIL:
C --------
C COSO   : COORDONNEES DES NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C IERR  : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C....+7..............................................................012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      DIMENSION COSO(3,NBS1*NBS2)
      DIMENSION LIMITE(4,NBRPRO)
      DIMENSION CEN(3,NBRPRO),VECN(3,NBRPRO),RAYON(NBRPRO)
C***********************************************************************
C     INITIALISER WBRPRO POUR ASSURER LA BONNE ADRESSE DANS LADEFI
C      NAD  = WBRPRO   VALEUR AU DESSOUS 0 EST FAUSSE ...
      NAD = 0
      IERR = 0
C***********************************************************************
C         RECHERCHE DES LIMITES TOPOLOGIQUES DES TRANSFORMATIONS
C***********************************************************************
      DO 30 N=1,NBRPRO
        DO 20 I=1,4
          LIMITE(I,N) = LADEFI(NAD+I)
   20   CONTINUE
        NAD = NAD+4
   30 CONTINUE
C***********************************************************************
C                RECHERCHE DES CENTRES DES PROJECTIONS
C           ET DES POINTS DEFINISSANT LES VECTEURS NORMAUX
C***********************************************************************
      DO 70 N=1,NBRPRO
        NOPO1 = LADEFI(NAD+2*N-1)
        NOPO2 = LADEFI(NAD+2*N)
        CALL XYZPOI(NOPO1,MNXYZ1,IERR)
        IF (IERR.NE.0) THEN
           NBLGRC(NRERR) = 1
           KERR(1) = 'ERREUR SUEX07: POINT NON RETROUVE'
           CALL LEREUR
           IERR = 1
           RETURN
        ENDIF
        CALL XYZPOI(NOPO2,MNXYZ2,IERR)
        IF (IERR.NE.0) THEN
          NBLGRC(NRERR) = 1
          KERR(1) =  'ERREUR SUEX07: POINT NON RETROUVE'
          CALL LEREUR
          IERR = 1
          RETURN
        ENDIF
        DO 60 I=1,3
          CEN(I,N)  = RMCN(MNXYZ1+I-1)
          VECN(I,N) = RMCN(MNXYZ2+I-1)-RMCN(MNXYZ1+I-1)
   60   CONTINUE
   70 CONTINUE
C***********************************************************************
C                RECHERCHE DES RAYONS DES CYLINDRES ET CERCLES
C***********************************************************************
      NAD = NAD+2*NBRPRO
      DO 80 N=1,NBRPRO
        RAYON(N) = RADEFI(NAD+N)
   80 CONTINUE
      CALL PROJEC(NBS1,NBS2,COSO,NBRPRO,LIMITE,CEN,VECN,RAYON,IERR)
      RETURN
      END
      SUBROUTINE PROJEC(NBS1,NBS2,COSO,NBRPRO,LIMITE,
     S                  CENTRE,VEN,RAYON,IERR)
C***********************************************************************
C BUT :    PROJECTION SUR CYLINDRE OU CERCLE
C***********************************************************************
C
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C TRAVAIL:
C           COSO   : COORDONNEES DES NOEUDS DU MAILLAGE
C
C SORTIES:
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      DIMENSION COSO(3,NBS1*NBS2)
      DIMENSION NDESCR(9),NDES(9)
      DIMENSION CEN(3),V1(3),V2(3),VECN(3),VECN1(3),VECN2(3),VN(3)
      DIMENSION LIMITE(4,NBRPRO)
      DIMENSION CENTRE(3,NBRPRO),VEN(3,NBRPRO),RAYON(NBRPRO)
      DO 2000 N=1,NBRPRO
C***********************************************************************
C                         ENTREE DES DONNEES
C***********************************************************************
C-----------------------------------------------------------------------
C                        ENTREE DE LA RUSTINE
C-----------------------------------------------------------------------
      IINF = LIMITE(1,N)
      ISUP = LIMITE(2,N)
      JINF = LIMITE(3,N)
      JSUP = LIMITE(4,N)
      IF (IINF.EQ.0) THEN
        IINF = 1
      ENDIF
      IF (ISUP.EQ.0) THEN
        ISUP = NBS1
      ENDIF
      IF (JINF.EQ.0) THEN
        JINF = 1
      ENDIF
      IF (JSUP.EQ.0) THEN
        JSUP = NBS2
      ENDIF
      IF (ISUP-IINF.LT.2) THEN
        WRITE (IMPRIM,*) 'MAUVAISE DEFINITION DES I'
        NBLGRC(NRERR) = 1
        KERR(1) =  'ERREUR SUEX07: MAUVAISE DEFINITION DES I'
        CALL LEREUR
        IERR = 1
        RETURN
      ENDIF
      IF (JSUP-JINF.LT.2) THEN
        WRITE (IMPRIM,*) 'MAUVAISE DEFINITION DES J'
        NBLGRC(NRERR) = 1
        KERR(1) =  'ERREUR SUEX07: MAUVAISE DEFINITION DES J'
        CALL LEREUR
        IERR = 1
        RETURN
      ENDIF
C-----------------------------------------------------------------------
C             ENTREE DES COORDONNEES GEOMETRIQUES DU CENTRE
C              ET DU VECTEUR DIRECTEUR DE L AXE DU CYLINDRE
C-----------------------------------------------------------------------
      DIST = 0.E0
      DO 10 I=1,3
        CEN(I) = CENTRE(I,N)
        VECN(I) = VEN(I,N)
        DIST = DIST+VECN(I)**2
   10 CONTINUE
      DIST = SQRT(DIST)
      VECN(1) = VECN(1)/DIST
      VECN(2) = VECN(2)/DIST
      VECN(3) = VECN(3)/DIST
C-----------------------------------------------------------------------
C                 ENTREE DU RAYON DU CYLINDRE OU CERCLE
C-----------------------------------------------------------------------
      RAY = RAYON(N)
C***********************************************************************
C          RECHERCHE AUTOMATIQUE DE CORRESPONDANCE DU CENTRE
C***********************************************************************
      DISREF = 1000000.E0
      I0 = 0
      J0 = 0
      DO 40 J=1,NBS2
        DO 30 I=1,NBS1
          NEU = (J-1)*NBS1+I
          DIST = 0
          DD   = 0
          PSCAL = 0
          DO 20 K=1,3
            DIST = DIST +(CEN(K)-COSO(K,NEU))**2
            PSCAL = PSCAL +((CEN(K)-COSO(K,NEU))*VECN(K))
   20     CONTINUE
          DO 21 K=1,3
            DD = DD + (PSCAL*VECN(K))**2
   21     CONTINUE
          DIST = DIST-DD
          DIST = SQRT(DIST)
          IF (DIST.LT.DISREF) THEN
            DISREF = DIST
            I0 = I
            J0 = J
          ENDIF
   30   CONTINUE
   40 CONTINUE
C***********************************************************************
C            CALCUL DE LA DISTANCE TOPOLOGIQUE DE PROJECTION
C***********************************************************************
      NN = 1
      SURFRE = 3.14156*RAY**2
   50 CONTINUE
      SURFT = 0
      II = MAX(I0-NN,1)
      IS = MIN(I0+NN-1,NBS1-1)
      JI = MAX(J0-NN,1)
      JS = MIN(J0+NN-1,NBS2-1)
      DO 110 J=JI,JS
        DO 100 I=II,IS
          NEU = (J-1)*NBS1+I
          SURF1 = 0
          DO 60 K=1,3
            V1(K) = COSO(K,NEU+1)-COSO(K,NEU)
            V2(K) = COSO(K,NEU+NBS1)-COSO(K,NEU)
   60     CONTINUE
          VECN1(1) = V1(2)*V2(3)-V1(3)*V2(2)
          VECN1(2) = V1(3)*V2(1)-V1(1)*V2(3)
          VECN1(3) = V1(1)*V2(2)-V1(2)*V2(1)
          DO 70 K=1,3
            SURF1 = SURF1+VECN1(K)**2
   70     CONTINUE
          SURF1 = SQRT(SURF1)/2
          NEU = J*NBS1+I+1
          SURF2 = 0.
          DO 80 K=1,3
            V1(K) = COSO(K,NEU-1)-COSO(K,NEU)
            V2(K) = COSO(K,NEU-NBS1)-COSO(K,NEU)
   80     CONTINUE
          VECN2(1) = V1(2)*V2(3)-V1(3)*V2(2)
          VECN2(2) = V1(3)*V2(1)-V1(1)*V2(3)
          VECN2(3) = V1(1)*V2(2)-V1(2)*V2(1)
          SIGN = 0.
          DO 90 K=1,3
            SURF2 = SURF2+VECN2(K)**2
            SIGN = SIGN+VECN1(K)*VECN2(K)
   90     CONTINUE
          SURF2 = SQRT(SURF2)/2
          IF (SIGN.GT.0) THEN
            SURFT = SURFT+SURF1+SURF2
          ELSE
            WRITE (IMPRIM,*) 'PROBLEME DE CALCUL DE SURFACE'
            NBLGRC(NRERR) = 1
            KERR(1) =  'ERREUR SUEX07: MAUVAIS CALCUL DES SURFACES'
            CALL LEREUR
            IERR = 1
            RETURN
          ENDIF
  100   CONTINUE
  110 CONTINUE
      XNBRET = 4.*NN**2
      XNBRE  = (JS-JI+1)*(IS-II+1)
      COR   = XNBRET/XNBRE
      SURFT = SURFT*COR
      IF (SURFT.LT.SURFRE) THEN
        NN = NN+1
        GOTO 50
      ENDIF
      N0 = NN
      NIINF = MAX(I0-N0,1)
      NISUP = MIN(I0+N0,NBS1)
      NJINF = MAX(J0-N0,1)
      NJSUP = MIN(J0+N0,NBS2)
C***********************************************************************
C                    CALCUL DU TYPE DE PROJECTION
C              INTERIEURE OU TOUCHANT UN OU DES BORDS
C***********************************************************************
C      N0I  = MAX(NISUP-NIINF,1)
C      N0J  = MAX(NJSUP-NJINF,1)
      NTYPE = 0
      NBRC  = 4
      IF (NJINF.EQ.1) THEN
        NTYPE = 1
        NBRC = NBRC-1
      ENDIF
      IF (NISUP.EQ.NBS1) THEN
        NTYPE = 2
        NBRC = NBRC-1
      ENDIF
      IF (NJSUP.EQ.NBS2) THEN
        NTYPE = 3
        NBRC = NBRC-1
      ENDIF
      IF (NIINF.EQ.1) THEN
        IF (NTYPE.NE.1) THEN
          NTYPE = 4
        ENDIF
        NBRC = NBRC-1
      ENDIF
C***********************************************************************
C                        PREPARATION DU CARRE
C***********************************************************************
      DO 500 II=1,9
        NDESCR(II) = 0
  500 CONTINUE
      NEUL            = (NJINF-1)*NBS1+I0
      NAD             = MOD(4-NTYPE,4)+1
      NDESCR(2*NAD-1) = NEUL
      NEUL            = (J0-1)*NBS1+NISUP
      NAD             = MOD(5-NTYPE,4)+1
      NDESCR(2*NAD-1) = NEUL
      NEUL            = (NJSUP-1)*NBS1+I0
      NAD             = MOD(6-NTYPE,4)+1
      NDESCR(2*NAD-1) = NEUL
      NEUL            = (J0-1)*NBS1+NIINF
      NAD             = MOD(7-NTYPE,4)+1
      NDESCR(2*NAD-1) = NEUL
      NDESCR(2*NBRC+1) = NDESCR(1)
      IF ((I0.NE.1).AND.(J0.NE.1).AND.
     S    (I0.NE.NBS1).AND.(J0.NE.NBS2)) THEN
        NEU0  = (J0-1)*NBS1+I0
        DO 510 KK=1,3
          V1(KK) = CEN(KK)-COSO(KK,NEU0)
  510   CONTINUE
        DO 540 I=NIINF,NISUP
          DO 530 J=NJINF,NJSUP
            NEU = (J-1)*NBS1+I
            DO 520 KK=1,3
              COSO(KK,NEU) = COSO(KK,NEU)+V1(KK)
  520       CONTINUE
  530     CONTINUE
  540   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C                    PROJECTION SUR LE CERCLE
C-----------------------------------------------------------------------
      DO 580 I=1,2*NBRC-1,2
        NEU   = NDESCR(I)
        J1   = (NEU-1)/NBS1+1
        I1   = NEU-(J1-1)*NBS1
        IF ((I1.NE.1).AND.(I1.NE.NBS1).AND.
     S      (J1.NE.1).AND.(J1.NE.NBS2)) THEN
          PSCAL = 0
          DO 550 K=1,3
            V1(K) = COSO(K,NEU)-CEN(K)
            PSCAL = PSCAL+V1(K)*VECN(K)
  550     CONTINUE
          PNORM = 0
          DO 560 K=1,3
            V2(K) = -PSCAL*VECN(K)+V1(K)
            PNORM = PNORM+V2(K)**2
  560     CONTINUE
          PNORM = SQRT(PNORM)
          DO 570 K=1,3
            COSO(K,NEU) = CEN(K)+RAY/PNORM*V2(K)
  570     CONTINUE
        ENDIF
  580 CONTINUE
C***********************************************************************
C      PROJECTION DES AUTRES POINTS SUR LE CARRE PAR HOMOGENEISATION
C                  DES LONGUEURS D ARRETES SUR LE CERCLE
C***********************************************************************
C              TRANSFORMATION DU TABLEAU DESCRIPTEUR
C-----------------------------------------------------------------------
      IF (NTYPE.NE.0) THEN
        NDES(1) = (NJINF-1)*NBS1 + NISUP
        NDES(2) = (NJSUP-1)*NBS1 + NISUP
        NDES(3) = (NJSUP-1)*NBS1 + NIINF
        NDES(4) = (NJINF-1)*NBS1 + NIINF
        NAD1 = MOD(NTYPE-1,4)+1
        NEU1 = NDES(NAD1)
        IF (NEU1.EQ.NDESCR(1)) THEN
          DO 800 I=2,NBRC
            NDESCR(I) = NDESCR(2*I-1)
  800     CONTINUE
        ELSE
          DO 810 I=1,NBRC
            NDESCR(I+1) = NDESCR(2*I-1)
  810     CONTINUE
          NDESCR(1) = NEU1
          NBRC = NBRC+1
        ENDIF
        NAD2 = MOD(NTYPE+NBRC-2,4)+1
        NEU2 = NDES(NAD2)
        IF (NEU2.NE.NDESCR(NBRC)) THEN
          NBRC = NBRC+1
          NDESCR(NBRC) = NEU2
        ENDIF
      ELSE
        NBRC = NBRC+1
        DO 820 I=1,NBRC
          NDESCR(I) = NDESCR(2*I-1)
  820   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C                    PROJECTION DES AUTRES POINTS
C-----------------------------------------------------------------------
      DO 900 I=1,NBRC-1
C-----------------------------------------------------------------------
C         CALCUL DE L ANGLE ENTRE DEUX POINTS DU TABLEAU DESCRIPTEUR
C-----------------------------------------------------------------------
        NEU1 = NDESCR(I)
        NEU2 = NDESCR(I+1)
        PSCAL = 0.
        PS1   = 0.
        PS2   = 0.
        DO 830 K=1,3
          V1(K) = COSO(K,NEU1)-CEN(K)
          V2(K) = COSO(K,NEU2)-CEN(K)
          PS1   = PS1+V1(K)**2
          PS2   = PS2+V2(K)**2
          PSCAL = PSCAL+V1(K)*V2(K)
  830   CONTINUE
        PS1   = SQRT(PS1)
        PS2   = SQRT(PS2)
        DO 835 K=1,3
          V1(K) = V1(K)/PS1
  835   CONTINUE
        PSCAL = PSCAL/(PS1*PS2)
        TETA  = ACOS(PSCAL)
        PS = 0
        DO 840 K=1,3
          VN(K) = V2(K)-PSCAL*PS2*V1(K)
          PS    = PS+VN(K)**2
  840   CONTINUE
        PS = SQRT(PS)
        DO 850 K=1,3
          VN(K) = VN(K)/PS
  850   CONTINUE
        J1   = (NEU1-1)/NBS1+1
        I1   = NEU1-(J1-1)*NBS1
        J2   = (NEU2-1)/NBS1+1
        I2   = NEU2-(J2-1)*NBS1
        NI = I2-I1
        NJ = J2-J1
        IF (NI.NE.0) THEN
          ID = ABS(NI)/NI
        ELSE
          ID = 0
        ENDIF
        IF (NJ.NE.0) THEN
          JD = ABS(NJ)/NJ
        ELSE
          JD = 0
        ENDIF
        NBRE = ABS(NI)+ABS(NJ)
        TETA = TETA/NBRE
        NEU = NEU1
C-----------------------------------------------------------------------
C    CALCUL DES NOUVELLES COORDONNEES DES POINTS ENTRE LES DEUX POINTS
C-----------------------------------------------------------------------
        DO 870 NN=1,NBRE-1
          IF (NI*NJ.GE.0) THEN
            IF (ABS(NI).GT.0) THEN
              NEU = NEU+ID
              NI  = NI-ID
            ELSE
              NEU = NEU+JD*NBS1
            ENDIF
          ELSE
            IF (ABS(NJ).GT.0) THEN
              NEU = NEU+JD*NBS1
              NJ  = NJ-JD
            ELSE
              NEU = NEU+ID
            ENDIF
          ENDIF
          DO 860 K=1,3
            V2(K) = RAY*(COS(NN*TETA)*V1(K)+SIN(NN*TETA)*VN(K))
            COSO(K,NEU) = CEN(K)+V2(K)
  860     CONTINUE
  870   CONTINUE
  900 CONTINUE
C***********************************************************************
C                    REMAILLAGE INTERNE DU TROU
C***********************************************************************
      CALL REMAIL(NBS1,NBS2,NIINF,NISUP,NJINF,NJSUP,COSO)
C***********************************************************************
C                    REMAILLAGE EXTERNE DU TROU
C***********************************************************************
C           REMAILLAGE DE LA ZONE SUPERIEURE DU TROU
C-----------------------------------------------------------------------
      CALL REMAIL(NBS1,NBS2,IINF,ISUP,NJSUP,JSUP,COSO)
C-----------------------------------------------------------------------
C           REMAILLAGE DE LA ZONE INFERIEURE DU TROU
C-----------------------------------------------------------------------
      CALL REMAIL(NBS1,NBS2,IINF,ISUP,JINF,NJINF,COSO)
C-----------------------------------------------------------------------
C           REMAILLAGE DE LA ZONE DROITE DU TROU
C-----------------------------------------------------------------------
      CALL REMAIL(NBS1,NBS2,NISUP,ISUP,JINF,JSUP,COSO)
C-----------------------------------------------------------------------
C           REMAILLAGE DE LA ZONE GAUCHE DU TROU
C-----------------------------------------------------------------------
      CALL REMAIL(NBS1,NBS2,IINF,NIINF,JINF,JSUP,COSO)
 2000 CONTINUE
      RETURN
      END
      SUBROUTINE REMAIL(NBS1,NBS2,NII,NIS,NJI,NJS,COSO)
C***********************************************************************
C BUT :    PROJECTION SUR CYLINDRE OU CERCLE
C***********************************************************************
C
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C SORTIES:
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C****+7**************************************************************012
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      DIMENSION COSO(3,NBS1*NBS2)
      N1 = NIS-NII+1
      N2 = NJS-NJI+1
      IF ((N1.LE.1).OR.(N2.LE.1)) THEN
        RETURN
      ENDIF
      LCOS   = 3*N1*N2
      NADCOS = 0
      NADPO1 = 0
      NADPO2 = 0
      NADPO3 = 0
      NADPO4 = 0
      CALL TNMCDC( 'REEL' , LCOS , NADCOS )
      CALL TNMCDC( 'REEL' , N1   , NADPO1 )
      CALL TNMCDC( 'REEL' , N2   , NADPO2 )
      CALL TNMCDC( 'REEL' , N1   , NADPO3 )
      CALL TNMCDC( 'REEL' , N2   , NADPO4 )
      DO 450 J=NJI,NJS
        DO 440 I=NII,NIS
          NEU  =  (J-1)*NBS1+I
          NEUL =  (J-NJI)*N1+I-NII+1
          DO 430 K=1,3
            NAD = NADCOS+K-1+3*(NEUL-1)
            RMCN(NAD) = COSO(K,NEU)
  430     CONTINUE
  440   CONTINUE
  450 CONTINUE
C     ATTENTION IL Y A DES TGS MAINTENANT!!!
      CALL SUEXQ222(N1,N2,RMCN(NADPO1),RMCN(NADPO2),RMCN(NADPO3),
     S            RMCN(NADPO4),RMCN(NADCOS))
      DO 480 J=NJI,NJS
        DO 470 I=NII,NIS
          NEU  =  (J-1)*NBS1+I
          NEUL =  (J-NJI)*N1+I-NII+1
          DO 460 K=1,3
            NAD = NADCOS+K-1+3*(NEUL-1)
            COSO(K,NEU) = RMCN(NAD)
  460     CONTINUE
  470   CONTINUE
  480 CONTINUE
      CALL TNMCDS( 'REEL' , LCOS , NADCOS )
      CALL TNMCDS( 'REEL' , N1   , NADPO1 )
      CALL TNMCDS( 'REEL' , N2   , NADPO2 )
      CALL TNMCDS( 'REEL' , N1   , NADPO3 )
      CALL TNMCDS( 'REEL' , N2   , NADPO4 )
      RETURN
      END
