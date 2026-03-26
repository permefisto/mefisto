      SUBROUTINE F2MP2P1( X, NOOBSF, NUMISU, NUMASU, LTDESU,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE LA MATRICE AE DU TRIANGLE  TAYLOR HOOD
C -----    LAGRANGE DE DEGRE 2 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1 POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA DENSITE DE MASSE EST SUPPOSEE CONSTANTE
C
C ENTREES:	
C --------
C X      : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          (ICI LA MASSE)
C
C SORTIES:
C --------
C AE : MATRICE ELEMENTAIRE 15x15 STOCKEE SYMETRIQUE PLEINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Mai 2007
C MODIF  : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Mars 2010
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
C
      REAL               X(6,2)
      DOUBLE PRECISION   AE(120)
      DOUBLE PRECISION   VMASSE
      INTEGER            LTDESU(1:MXDOFL,NUMISU:NUMASU)
C
      DOUBLE PRECISION   DELTA, XB, YB, ZB
      DOUBLE PRECISION   X21, Y21, X31, Y31
      DOUBLE PRECISION   A, B, C, D, E, COEF
C
C     CALCUL DU JACOBIEN
C     ==================
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
C
      DELTA = ABS( X21*Y31 - X31*Y21 )
C
C     RECHERCHE DE LA DENSITE DE MASSE AU BARYCENTRE DU TRIANGLE
C     ==========================================================
      XB = (X(1,1) + X(2,1) + X(3,1)) / 3.D0
      YB = (X(1,2) + X(2,2) + X(3,2)) / 3.D0
      ZB = 0D0
      CALL REMASS2( 3, NOOBSF, XB, YB, ZB,
     %              LTDESU(LPMASS,NOOBSF), VMASSE )
C
C     LES COEFFICIENTS DE LA MATRICE DE MASSE
C     =======================================
      COEF = VMASSE  * DELTA
      A    =  COEF /  60.D0
      B    = -COEF / 360.D0
      C    =  COEF * (2D0/45D0)
      D    =  COEF * (4D0/45D0)
      E    = -COEF /  90.D0
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
C     ==================================================================
C     BLOC A[1,1]  integrale P2i P2j * rho moyen * delta
      AE(1)= A
C
      AE(2)= B
      AE(3)= A
C
      AE(4)= B
      AE(5)= B
      AE(6)= A
C
      AE(7)= 0.D0
      AE(8)= 0.D0
      AE(9)= E
      AE(10)=D
C
      AE(11)=E
      AE(12)=0.D0
      AE(13)=0.D0
      AE(14)=C
      AE(15)=D
C
      AE(16)=0.D0
      AE(17)=E
      AE(18)=0.D0
      AE(19)=C
      AE(20)=C
      AE(21)=D
C
C     BLOC A[2,1]
      AE(22)=0D0
      AE(23)=0D0
      AE(24)=0D0
      AE(25)=0D0
      AE(26)=0D0
      AE(27)=0D0
C     BLOC A[2,2]
      AE(28)=A
C
C     BLOC A[2,1]
      AE(29)=0D0
      AE(30)=0D0
      AE(31)=0D0
      AE(32)=0D0
      AE(33)=0D0
      AE(34)=0D0
C
C     BLOC A[2,2]
      AE(35)= B
      AE(36)= A
C
C     BLOC A[2,1]
      AE(37)=0D0
      AE(38)=0D0
      AE(39)=0D0
      AE(40)=0D0
      AE(41)=0D0
      AE(42)=0D0
C     BLOC A[2,2]
      AE(43)= B
      AE(44)= B
      AE(45)= A
C
C     BLOC A[2,1]
      AE(46)=0D0
      AE(47)=0D0
      AE(48)=0D0
      AE(49)=0D0
      AE(50)=0D0
      AE(51)=0D0
C     BLOC A[2,2]
      AE(52)= 0.D0
      AE(53)= 0.D0
      AE(54)= E
      AE(55)= D
C
C     BLOC A[2,1]
      AE(56)=0D0
      AE(57)=0D0
      AE(58)=0D0
      AE(59)=0D0
      AE(60)=0D0
      AE(61)=0D0
C     BLOC A[2,2]
      AE(62)=E
      AE(63)=0.D0
      AE(64)=0.D0
      AE(65)=C
      AE(66)=D
C
C     BLOC A[2,1]
      AE(67)=0D0
      AE(68)=0D0
      AE(69)=0D0
      AE(70)=0D0
      AE(71)=0D0
      AE(72)=0D0
C     BLOC A[2,2]
      AE(73)=0.D0
      AE(74)=E
      AE(75)=0.D0
      AE(76)=C
      AE(77)=C
      AE(78)=D
C
C     BLOC A[3,1] A[3,2] A[3,3]
      DO I = 79, 120
        AE(I)=0.D0
      ENDDO
C
      RETURN
      END
