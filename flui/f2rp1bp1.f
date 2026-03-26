      SUBROUTINE F2RP1BP1( X,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE PROFIL DU TRIANGLE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C
C ENTREES:	
C --------
C X      : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES FLUIDES
C
C SORTIES:
C --------
C AE     : MATRICE PROFIL ELEMENTAIRE 11x11  AVANT GAUSS FRONTALE
C          SUR LES DL VITESSE DU BARYCENTRE  4 8
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Novembre 2008
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
C
      DOUBLE PRECISION   EPSILON
      PARAMETER         (EPSILON=1D-7)
C ATTENTION: VALEUR CRUCIALE 1D-10 pour la cavity2d DONNE DE MAUVAIS RESULTATS!
C ATTENTION: VALEUR CRUCIALE 1D-14 pour le tube2d   DONNE DE MAUVAIS RESULTATS!
C ATTENTION: VALEUR CRUCIALE 1D-16 pour le tube2d   DONNE DE MAUVAIS RESULTATS!
C ATTENTION: VALEUR PIRE     1D-20 pour le tube2d   DONNE DE MAUVAIS RESULTATS!
C
      REAL              X(3,2)
      DOUBLE PRECISION  AE(66)
      DOUBLE PRECISION  VISCOS, VISDEL, COEFPRES, COEFPRED, COPRES
      INTEGER           NOOBSF, NUMISU, NUMASU
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      DOUBLE PRECISION  DELTA, X21, Y21, X31, Y31
      DOUBLE PRECISION  TDFDF(2,2)
      INTRINSIC         ABS
C
      DOUBLE PRECISION  TDPDP(4,2,4,2)
C     TDPDP(i,k,j,l) = integrale dPi/dxk DPj/dxl dX
C
      DATA TDPDP/
     %  0.95D+00,   -0.5D-01,
     %  0.45D+00,   -0.135D+01,
     %  0.725D+00,   0.225D+00,
     % -0.275D+00,  -0.675D+00,
     % -0.5D-01,     0.95D+00,
     %  0.45D+00,   -0.135D+01,
     % -0.275D+00,   0.225D+00,
     %  0.725D+00,  -0.675D+00,
     %  0.45D+00,    0.45D+00,
     %  0.45D+00,   -0.135D+01,
     %  0.225D+00,   0.225D+00,
     %  0.225D+00,  -0.675D+00,
     % -0.135D+01,  -0.135D+01,
     % -0.135D+01,   0.405D+01,
     % -0.675D+00,  -0.675D+00,
     % -0.675D+00,   0.2025D+01,
     %  0.725D+00,  -0.275D+00,
     %  0.225D+00,  -0.675D+00,
     %  0.95D+00,    0.45D+00,
     % -0.5D-01,    -0.135D+01,
     %  0.225D+00,   0.225D+00,
     %  0.225D+00,  -0.675D+00,
     %  0.45D+00,    0.45D+00,
     %  0.45D+00,   -0.135D+01,
     % -0.275D+00,   0.725D+00,
     %  0.225D+00,  -0.675D+00,
     % -0.5D-01,     0.45D+00,
     %  0.95D+00,   -0.135D+01,
     % -0.675D+00,  -0.675D+00,
     % -0.675D+00,   0.2025D+01,
     % -0.135D+01,  -0.135D+01,
     % -0.135D+01,   0.405D+01 /
C
      include"./incl/dp1bla2d.inc"
C     INTEGRALES P1 DP1B SUR LE TRIANGLE de REFERENCE
C     DOUBLE PRECISION  DP1Bla2d(4,2,3)
C     DP1Bla2d(i,k,j) = integrale DP1Bi/dxk Lambdaj dX
C
C     CALCUL DE LA VISCOSITE AU BARYCENTRE DU TRIANGLE
C     ------------------------------------------------
      X21 = (X(1,1) + X(2,1) + X(3,1)) / 3.D0
      Y21 = (X(1,2) + X(2,2) + X(3,2)) / 3.D0
      CALL REVISC( 3, NOOBSF, X21, Y21, 0D0,
     %             LTDESU(LPVISC,NOOBSF), VISCOS )
C
C     RECHERCHE DU COEFFICIENT SUR LE GRADIENT DE LA PRESSION
C     -------------------------------------------------------
      IF( LTDESU(LPCPRE,NOOBSF) .GT. 0 ) THEN
C        IL EXISTE UN COEFFICIENT DEVANT LA PRESSION
         CALL RECPRE( 3, NOOBSF, X21, Y21, 0D0,
     %                LTDESU(LPCPRE,NOOBSF), COPRES )
      ELSE
C        IL N'EXISTE PAS DE COEFFICIENT DEVANT LA PRESSION
         COPRES = 1D0
      ENDIF
C
C     CALCULS INITIAUX SUR Fe: e ref -> e
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
C
C     CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
      DELTA = ABS( X21*Y31 - X31*Y21 )
C
C     CALCUL DES COEFFICIENTS DE LA SOUS MATRICE t DF-1 * DF-1 * VISCOS / DELTA
      VISDEL = VISCOS / DELTA
      TDFDF(1,1) = VISDEL * (  X31 * X31 + Y31 * Y31 )
      TDFDF(1,2) = VISDEL * ( -X21 * X31 - Y21 * Y31 )
      TDFDF(2,1) = TDFDF(1,2)
      TDFDF(2,2) = VISDEL * (  X21 * X21 + Y21 * Y21 )
C
C     LES 2 BLOCS DIAGONAUX 4x4 Vitesse Vitesse DE LA MATRICE ELEMENTAIRE
C     STOCKAGE DE LA PARTIE TRIANGULAIRE SUPERIEURE COLONNE PAR COLONNE
      M11 = 0
      M22 = 10
      DO J=1,4
         M22 = M22 + 4
         DO I=1,J
            M11 = M11 + 1
            AE(M11) = 0D0
            DO K=1,2
               DO L=1,2
                  AE(M11) = AE(M11) + TDPDP(I,K,J,L) * TDFDF(K,L)
               ENDDO
            ENDDO
            M22 = M22 + 1
            AE(M22) = AE(M11)
         ENDDO
      ENDDO
C
C     LE BLOC A12 VITESSExVITESSE est NUL
      M12 = 10
      DO J=1,4
         DO I=1,4
            M12 = M12 + 1
            AE(M12) = 0D0
         ENDDO
         M12 = M12 + J
      ENDDO
C
C     LES BLOCS  A31 4x3  et  A32 4x3 DE LA MATRICE ELEMENTAIRE
      M31 = 36
      M32 = 40
      DO J=1,3
         DO I=1,4
            M31 = M31 + 1
            AE(M31) = (-Y31*Dp1bla2d(I,1,J) +Y21*Dp1bla2d(I,2,J))*COPRES
            M32 = M32 + 1
            AE(M32) = ( X31*Dp1bla2d(I,1,J) -X21*Dp1bla2d(I,2,J))*COPRES
         ENDDO
         M31 = M32 + J
         M32 = M31 + 4
      ENDDO
C
C
C     PENALISATION RELATIVE DU BLOC DIAGONAL EN PRESSION
C     --------------------------------------------------
C     COEFFICIENT DIAGONAL DU BLOC
      COEFPRED = EPSILON * DELTA * VISCOS / 12D0
C     COEFFICIENT NON  DIAGONAL
      COEFPRES = COEFPRED / 2D0
C
C     LE BLOC INTEGRALE EPSILON P Q dx dy
      AE(45) = COEFPRED
C
      AE(54) = COEFPRES
      AE(55) = COEFPRED
C
      AE(64) = COEFPRES
      AE(65) = COEFPRES
      AE(66) = COEFPRED
C
      RETURN
      END
