      SUBROUTINE F2MP1BP1( X, NOOBSF, NUMISU, NUMASU, LTDESU,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE LA MATRICE AE DU TRIANGLE BREZZI-FORTIN
C -----    LAGRANGE DE DEGRE 1 + BULLE P3 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1            POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA DENSITE DE MASSE EST SUPPOSEE CONSTANTE
C
C ENTREES:	
C --------
C X      : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE DE L'OBJET
C          (ICI LA MASSE)
C
C SORTIES:
C --------
C AE : MATRICE ELEMENTAIRE 11x11 STOCKEE SYMETRIQUE PLEINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Juin 2007
C MODIF : Alain PERRONNET LJLL UPMC et St Pierre du Perray      Mai 2010
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
C
      REAL               X(3,2)
      DOUBLE PRECISION   AE(66)
      DOUBLE PRECISION   VMASSE
      INTEGER            LTDESU( 1:MXDOFL, NUMISU:NUMASU )
C
      DOUBLE PRECISION   X21, Y21, X31, Y31, DELTA
      DOUBLE PRECISION   A, B, C, D, COEF
C
C     RECHERCHE DE LA DENSITE DE MASSE AU BARYCENTRE DU TRIANGLE
C     ==========================================================
      X21 = (X(1,1) + X(2,1) + X(3,1)) / 3.D0
      Y21 = (X(1,2) + X(2,2) + X(3,2)) / 3.D0
      CALL REMASS2( 3, NOOBSF, X21, Y21, 0D0,
     %              LTDESU(LPMASS,NOOBSF), VMASSE )
C
C     CALCUL DU JACOBIEN de Fe: e ref -> e  P1 par composante
C     =======================================================
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
C
      DELTA = ABS( X21*Y31 - X31*Y21 )
C
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
C     ====================================================
      COEF = VMASSE * DELTA
      A    = COEF * 83D0 / 1680D0
      B    = COEF * 13D0 / 1680D0
      C    = COEF *  3D0 /  112D0
      D    = COEF * 81D0 /  560D0
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
C     ==================================================================
C     BLOC A[NTDL + 2*NUEF1,1]  integrale Pi Pj * rho moyen * delta
      AE(1) = A
C
      AE(2) = B
      AE(3) = A
C
      AE(4) = B
      AE(5) = B
      AE(6) = A
C
      AE(7)  = C
      AE(8)  = C
      AE(9)  = C
      AE(10) = D
C
C     BLOC A[2,1]
      AE(11) = 0D0
      AE(12) = 0D0
      AE(13) = 0D0
      AE(14) = 0D0
C     BLOC A[2,2]
      AE(15) = A
C
C     BLOC A[2,1]
      AE(16) = 0D0
      AE(17) = 0D0
      AE(18) = 0D0
      AE(19) = 0D0
C     BLOC A[2,2]
      AE(20) = B
      AE(21) = A
C
C     BLOC A[2,1]
      AE(22) = 0D0
      AE(23) = 0D0
      AE(24) = 0D0
      AE(25) = 0D0
C     BLOC A[2,2]
      AE(26) = B
      AE(27) = B
      AE(28) = A
C
C     BLOC A[2,1]
      AE(29) = 0D0
      AE(30) = 0D0
      AE(31) = 0D0
      AE(32) = 0D0
C     BLOC A[2,2]
      AE(33) = C
      AE(34) = C
      AE(35) = C
      AE(36) = D
C
C     BLOC A[3,1] A[3,2] A[3,3]
      DO I = 37, 66
        AE(I) = 0D0
      ENDDO
C
      RETURN
      END
