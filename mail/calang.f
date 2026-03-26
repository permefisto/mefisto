      SUBROUTINE CALANG(NBS,COSOBO,IMPRES,SENSFR,ALPHA)
C***********************************************************************
C BUT : CALCUL DE TOUS LES ANGLES SUR LE BORD EN PASSANT PAR L INTERIEUR
C       ET DU SENS DE ROTATION DE LA FRONTIERE
C***********************************************************************
C
C ENTREES:
C    NBS    : NOMBRE DE POINTS SUR FORMANT LE CONTOUR FERME
C                                                     *****
C    COSOBO : MATRICE DES COORDONNEES DE LA FRONTIERE
C    IMPRES : 0 : LA PROCEDURE DEVIENT MUETTE (PAS DE MESSAGE D'ERREUR)
C             1 : AFFICHAGE DES MESSAGES D'ERREUR A L'ECRAN
C
C SORTIES:
C    SENSFR : SENS DE ROTATION DE LA FRONTIERE
C                  -1 : SENS INDIRECT
C                   0 : FRONTIERE DEGENEREE
C                   1 : SENS DIRECT
C             -2 OU 2 : UN DES COINS NON CONVEXE
C    ALPHA  : ANGLE DES DIFFERENTS POINTS DU BORD EN PASSANT PAR
C                    L'INTERIEUR
C
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      DIMENSION         COSOBO(3,NBS)
      DIMENSION         ALPHA(NBS)
      DIMENSION         V1(2), V2(2)
      REAL              XYZ1(3), XYZ2(3)
      INTEGER           SENSFR
C
      PI = 4.E0*ATAN(1.E0)
C***********************************************************************
C  VERIFICATION QUE LE CONTOUR EST BIEN FERME
C***********************************************************************
      DO 1 I=1,2
         XYZ1(I) = COSOBO(I,  1)
         XYZ2(I) = COSOBO(I,NBS)
 1    CONTINUE
      XYZ1(3) = 0.0
      XYZ2(3) = 0.0
      CALL XYZIDE( XYZ1, XYZ2, I )
      IF ( I .EQ. 0 ) THEN
         IF (IMPRES.EQ.1) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'ERREUR CALANG :'
            KERR(2) = 'LE CONTOUR N''EST PAS FERME'
            CALL LEREUR
         ENDIF
         RETURN
      ENDIF
C***********************************************************************
C  CALCUL DES ANGLES ENTRE LES SEGMENTS(M(I),M(I+1)) ET (M(I),M(I-1)
C           COMPTES DANS LES SENS TRIGONOMETRIQUE
C***********************************************************************
      NBAR = NBS-1
      NBRALP = 0
      SOMALP = 0.
      NDECAL = 0
      IDEB   = 0
      DO 10 I=1,NBAR
        V1(1) = COSOBO(1,I+1)-COSOBO(1,I)
        V1(2) = COSOBO(2,I+1)-COSOBO(2,I)
        DIST  = SQRT(V1(1)**2+V1(2)**2)
        IF (DIST.NE.0) THEN
          IDEB = I
          GOTO 20
        ENDIF
 10   CONTINUE
C
 20   IFIN = 0
      DO 30 I=NBAR,1,-1
        V1(1) = COSOBO(1,I+1)-COSOBO(1,I)
        V1(2) = COSOBO(2,I+1)-COSOBO(2,I)
        DIST  = SQRT(V1(1)**2+V1(2)**2)
        IF (DIST.NE.0) THEN
          IFIN = I
          GOTO 40
        ENDIF
 30   CONTINUE
C
 40   IF ((COSOBO(1,IDEB).NE.COSOBO(1,IFIN+1))
     S    .OR.(COSOBO(2,IDEB).NE.COSOBO(2,IFIN+1))) THEN
        WRITE (IMPRIM,*) 'ERREUR CALANG : JONCTION '
        SENSFR = 0
        RETURN
      ENDIF
      DO 50 I=IDEB,IFIN
        V1(1) = COSOBO(1,I+1)-COSOBO(1,I)
        V1(2) = COSOBO(2,I+1)-COSOBO(2,I)
        IF (I.EQ.IDEB) THEN
          V2(1) = COSOBO(1,IFIN)-COSOBO(1,IFIN+1)
          V2(2) = COSOBO(2,IFIN)-COSOBO(2,IFIN+1)
        ELSE
          V2(1) = COSOBO(1,I-1+NDECAL)-COSOBO(1,I+NDECAL)
          V2(2) = COSOBO(2,I-1+NDECAL)-COSOBO(2,I+NDECAL)
        ENDIF
        DIST  = SQRT(V1(1)**2+V1(2)**2)
        IF (DIST.NE.0) THEN
          V1(1) = V1(1)/DIST
          V1(2) = V1(2)/DIST
        ELSE
          NDECAL = NDECAL-1
          GOTO 50
        ENDIF
        DIST  = SQRT(V2(1)**2+V2(2)**2)
        IF (DIST.NE.0) THEN
          V2(1) = V2(1)/DIST
          V2(2) = V2(2)/DIST
        ELSE
          WRITE (IMPRIM,*) 'ERREUR CALANG : UNE BARRE NULLE'
          SENSFR = 0
          RETURN
        ENDIF
        NDECAL= 0
        PSCAL = V1(1)*V2(1)+V1(2)*V2(2)
        PSCAL = AMIN1(PSCAL,1.)
        PSCAL = AMAX1(PSCAL,-1.)
        ALP   = ACOS(PSCAL)
        PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
        IF (PSCAL.LT.0.) THEN
          ALP = 2*PI-ALP
        ENDIF
        ALPHA(I) = ALP
        NBRALP   = NBRALP+1
        SOMALP   = SOMALP+ALP
   50 CONTINUE
C***********************************************************************
C     RECHERCHE DU SENS DE ROTATION SUIVANT LA SOMME DES ANGLES
C***********************************************************************
      SOMDIR = 2*PI+(NBRALP-4)*PI
      SOMIND = SOMDIR+4*PI
      IF (ABS(SOMALP-SOMIND).LT.PI) THEN
        DO 60 I=1,NBAR
          ALPHA(I) = 2*PI-ALPHA(I)
   60   CONTINUE
        SENSFR = -1
      ELSE
        IF (ABS(SOMALP-SOMDIR).LT.PI) THEN
          SENSFR = 1
        ELSE
          SENSFR = 0
          IF (IMPRES.EQ.1) THEN
            NBLGRC(NRERR) = 2
            KERR(1) =  'ERREUR CALANG :'
            KERR(2) =  'LE CONTOUR EST DEGENERE'
            CALL LEREUR
          ENDIF
        ENDIF
      ENDIF
      ALPHA(NBS) = ALPHA(1)
      RETURN
      END
