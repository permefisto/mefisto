      SUBROUTINE QUAQUA( NOSOEL, XYZSOM, QUALIT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA QUALITE D'UN QUADRANGLE DANS R**3
C -----   QUALITE: 1 - ( MAX( | Pi/2 - ANGLE AU SOMMET I | ) / Pi/2
C                      I=1,...,4
C
C ENTREES:
C --------
C NOSOEL : NUMERO DES 4 SOMMETS DANS XYZSOM DES SOMMETS DU QUADRANGLE
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE
C
C SORTIE :
C --------
C QUALIT : QUALITE DU QUADRANGLE (VALEUR COMPRISE ENTRE 0 ET 1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : CHRISTOPHE DOURSAT  ANALYSE NUMERIQUE UPMC PARIS    MARS 1991
C2345X7..............................................................012
      DOUBLE PRECISION  PI, PI2
      PARAMETER        (PI  = 3.1415926535897932D0,
     %                  PI2 = 1.5707963267948966D0 )
C
      INTEGER           NOSOEL(1:4)
      REAL              XYZSOM(3,*), QUALIT
      DOUBLE PRECISION  ALP, SOMALP, TETA12
C
C     TABLEAUX AUXILIAIRES
      INTEGER           NVOIS(4,2)
      REAL              V1(3), V2(3)
C
      DO 10 I=1,4
         NVOIS(I,1) =  1
         NVOIS(I,2) = -1
  10  CONTINUE
      NVOIS(1,2) =  3
      NVOIS(4,1) = -3
C
      QUALIT = 1.
      SENSFR = 0
C
      SOMALP = 0
      SOMPSC = 0
C
      DO 60 NC=1,4
C
         NEU  = NOSOEL( NC )
         NEU1 = NOSOEL( NC+NVOIS(NC,1) )
         NEU2 = NOSOEL( NC+NVOIS(NC,2) )
         XLG1  = 0.
         XLG2  = 0.
         XPS12 = 0.
         DO 50 II=1,3
            V1(II) = XYZSOM(II,NEU1) - XYZSOM(II,NEU)
            V2(II) = XYZSOM(II,NEU2) - XYZSOM(II,NEU)
            XLG1   = XLG1 +  V1(II) ** 2
            XLG2   = XLG2 +  V2(II) ** 2
            XPS12  = XPS12 + V1(II) * V2(II)
  50     CONTINUE
         XLG1 = SQRT(XLG1)
         XLG2 = SQRT(XLG2)
C
         IF ((XLG1*XLG2).NE.0) THEN
            XPS12 = XPS12/(XLG1*XLG2)
            XPS12 = MAX ( XPS12 , -1.E0 )
            XPS12 = MIN ( XPS12 ,  1.E0 )
            TETA12 = ACOS(XPS12)
            IF(SENSFR.NE.0) THEN
              PSCAL = -V1(2)*V2(1)+V1(1)*V2(2)
              IF (PSCAL.LT.0.) THEN
                ALP = 2*PI-TETA12
              ELSE
                ALP = TETA12
              ENDIF
              SOMALP = SOMALP+ALP
              IF (PSCAL.GT.0) THEN
                SOMPSC = SOMPSC+1
              ELSE
                SOMPSC = SOMPSC-1
              ENDIF
            ENDIF
            TETA12 = 1.-(ABS(TETA12-PI2)/PI2)
C
C           CRITERE D ORTHOGONALITE PAR LE MINIMUM
            IF( QUALIT .GT. TETA12 ) THEN
               QUALIT = REAL( TETA12 )
            ENDIF
C
         ELSE
C
            QUALIT = 0
C
         ENDIF
  60  CONTINUE
C
      IF (SENSFR.NE.0) THEN
        IF (ABS(SOMALP-2*PI).LT.PI) THEN
           SENSRO = 1
        ELSE IF (ABS(SOMALP-6*PI).LT.PI) THEN
           SENSRO = -1
        ELSE
            SENSRO = 0
        ENDIF
        IF (SENSRO.NE.SENSFR) THEN
           QUALIT = 0
        ENDIF
        IF ((SOMPSC.NE.4).AND.(SOMPSC.NE.-4)) THEN
           QUALIT = 0
        ENDIF
      ENDIF
C
      RETURN
      END
