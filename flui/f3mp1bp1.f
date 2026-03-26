      SUBROUTINE F3MP1BP1( X, NOVOLU, NUMIVO, NUMAVO, LTDEVO,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE MASSE AE DU TETRAEDRE BREZZI-FORTIN
C -----    LAGRANGE DE DEGRE 1 + BULLE P4 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1            POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA DENSITE DE MASSE EST SUPPOSEE CONSTANTE
C
C ENTREES:	
C --------
C X      : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C NOVOLU : NUMERO DE L'OBJET VOLUME DE CET ELEMENT FINI
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE DE L'OBJET
C          (ICI SEULE LA DENSITE VOLUMIQUE DE MASSE EST UTILISEE)
C
C SORTIES:
C --------
C AE : MATRICE ELEMENTAIRE 19x19 STOCKEE SOUS FORME SYMETRIQUE PLEINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Novembre 2008
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
C
      REAL               X(4,3)
      DOUBLE PRECISION   AE(190)
      DOUBLE PRECISION   VMASSE
      INTEGER            NOVOLU,NUMIVO,NUMAVO
      INTEGER            LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER            I
      DOUBLE PRECISION   X1, Y1, Z1, ABS, DETM33, DELTA
      DOUBLE PRECISION   A, B, C, D, COEF
C
C     RECHERCHE DE LA DENSITE DE MASSE AU BARYCENTRE DU TETRAEDRE
C     ===========================================================
      X1 = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      Y1 = ( X(1,2) + X(2,2) + X(3,2) + X(4,2) ) * 0.25D0
      Z1 = ( X(1,3) + X(2,3) + X(3,3) + X(4,3) ) * 0.25D0
      CALL REMASS2( 4, NOVOLU, X1, Y1, Z1,
     %              LTDEVO(LPMASS,NOVOLU), VMASSE )
C
C     CALCUL DU JACOBIEN
C     ==================
      X1 = X(1,1)
      Y1 = X(1,2)
      Z1 = X(1,3)
      DELTA = ABS( DETM33( X(2,1)-X1, X(3,1)-X1, X(4,1)-X1,
     %                     X(2,2)-Y1, X(3,2)-Y1, X(4,2)-Y1,
     %                     X(2,3)-Z1, X(3,3)-Z1, X(4,3)-Z1 ) )
C
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
C     ====================================================
      COEF = VMASSE * DELTA
      A = 29836D0 / 2494800D0 * COEF
      B =  9046D0 / 2494800D0 * COEF
      C =   956D0 /  155925D0 * COEF
      D =  4096D0 /  155925D0 * COEF
C
CCC     A=  0.11959275292608627D-01
CCC     B=  0.36259419592752926D-02
CCC     C=  0.61311527978194641D-02
CCC     D=  0.26269039602372937D-01
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
C     ==================================================================
C     BLOC A[1,1]  integrale Pi Pj * rho moyen * delta
      AE(1) = A
C
      AE(2) = B
      AE(3) = A
C
      AE(4) = B
      AE(5) = B
      AE(6) = A
C
      AE(7) = B
      AE(8) = B
      AE(9) = B
      AE(10)= A
C
      AE(11) = C
      AE(12) = C
      AE(13) = C
      AE(14) = C
      AE(15) = D
C
C     RESTE DE LA MATRICE ELEMENTAIRE MISE A ZERO
      DO I = 16, 190
        AE(I) = 0D0
      ENDDO
C
C     MODIFICATION DU BLOC A[2,2]
      AE(21) = A
C
      AE(27) = B
      AE(28) = A
C
      AE(34) = B
      AE(35) = B
      AE(36) = A
C
      AE(42) = B
      AE(43) = B
      AE(44) = B
      AE(45) = A
C
      AE(51) = C
      AE(52) = C
      AE(53) = C
      AE(54) = C
      AE(55) = D
C
C     MODIFICATION DU BLOC A[3,3]
      AE(66) = A
C
      AE(77) = B
      AE(78) = A
C
      AE(89) = B
      AE(90) = B
      AE(91) = A
C
      AE(102) = B
      AE(103) = B
      AE(104) = B
      AE(105) = A
C
      AE(116) = C
      AE(117) = C
      AE(118) = C
      AE(119) = C
      AE(120) = D
C
      RETURN
      END
