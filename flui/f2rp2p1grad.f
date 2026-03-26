      SUBROUTINE F2RP2P1GRAD( X, NOOBSF, NUMISU, NUMASU, LTDESU,   AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE VISCOSITE AE DU TRIANGLE  TAYLOR HOOD
C -----    LAGRANGE DE DEGRE 2 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1 POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA VISCOSITE EST SUPPOSEE CONSTANTE
C          SUR CHAQUE EF
C VARIANTE: Integrale Epsilon Grad Q grad P dx avec Q,P Polynomes P1
C           pour le BLOC A(3,3)
C
C ENTREES:	
C --------
C X      : LES 2 COORDONNEES DES 3 SOMMETS 3 MILIEUX DES ARETES DU TRIANGLE
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE DE L'OBJET
C
C SORTIES:
C --------
C AE : MATRICE ELEMENTAIRE 15x15 STOCKEE SYMETRIQUE PLEINE (6+6+3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Auteur : ALAIN PERRONNET LJLL UPMC et Saint Pierre du Perray  Mai 2011
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
C
      DOUBLE PRECISION   EPSILON
ccc      PARAMETER         (EPSILON=1D-7)
      PARAMETER         (EPSILON=1D-9)
C ATTENTION: VALEUR CRUCIALE 1D-14 pour le tube DONNE DE MAUVAIS RESULTATS!
C ATTENTION: VALEUR CRUCIALE 1D-16 pour le tube DONNE DE MAUVAIS RESULTATS!
C ATTENTION: VALEUR PIRE     1D-20 pour le tube DONNE DE MAUVAIS RESULTATS!
C
      REAL               X(6,2)
      DOUBLE PRECISION   AE(120)
      DOUBLE PRECISION   VISCOS, COPRES
      INTEGER            NOOBSF, NUMISU, NUMASU
      INTEGER            LTDESU(1:MXDOFL,NUMISU:NUMASU)
C
      DOUBLE PRECISION   DELTA, XD, YD
      DOUBLE PRECISION   X21, Y12, X13, Y31, X32, Y32
      DOUBLE PRECISION   A, B, D, COEFF1, COEFF2
C
C     RECHERCHE DE LA VISCOSITE AU BARYCENTRE DU TRIANGLE
C     ===================================================
      XD = (X(1,1) + X(2,1) + X(3,1)) / 3.D0
      YD = (X(1,2) + X(2,2) + X(3,2)) / 3.D0
      CALL REVISC( 3, NOOBSF, XD, YD, 0D0,
     %             LTDESU(LPVISC,NOOBSF), VISCOS )
C
C     RECHERCHE DU COEFFICIENT SUR LE GRADIENT DE LA PRESSION
C     =======================================================
      IF( LTDESU(LPCPRE,NOOBSF) .GT. 0 ) THEN
C        IL EXISTE UN COEFFICIENT DEVANT LA PRESSION
         CALL RECPRE( 3, NOOBSF, XD, YD, 0D0,
     %                LTDESU(LPCPRE,NOOBSF), COPRES )
      ELSE
C        IL N'EXISTE PAS DE COEFFICIENT DEVANT LA PRESSION
         COPRES = 1D0
      ENDIF
C
C     CALCUL DE QUELQUES CONSTANTES UTILES
C     ====================================
      X21 = X(2,1) - X(1,1)
      X13 = X(1,1) - X(3,1)
C
      Y12 = X(1,2) - X(2,2)
      Y31 = X(3,2) - X(1,2)
C
      X32 = X(3,1) - X(2,1)
      Y32 = X(3,2) - X(2,2)
C
      A = Y31 * Y31 + X13 * X13
      B = Y31 * Y12 + X13 * X21
      D = Y12 * Y12 + X21 * X21
C
      DELTA = ABS( X21*Y31 - X13*Y12 )
C
C     COEFF1 = VISCOS / (ABS( X21*Y31 - X13*Y12 ))
      COEFF1 = VISCOS  / DELTA
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
C     ==================================================================
C     BLOC A[1,1]
      AE(1)=COEFF1*(0.5D0*(A+D)+B)
C
      AE(2)=COEFF1/6D0*(A+B)
      AE(3)=COEFF1*0.5D0*A
C
      AE(4)= COEFF1*(B+D)/6D0
      AE(5)=-COEFF1/6D0*B
      AE(6)= COEFF1*0.5D0*D
C
      AE(7)= COEFF1*(-2D0/3D0*(A+B))
      AE(8)= COEFF1*(-2D0/3D0*(A+B))
      AE(9)= 0D0
      AE(10)=COEFF1*(4D0/3D0*(A+B+D))
C
      AE(11)=0D0
      AE(12)=COEFF1*( 2D0/3D0*B)
      AE(13)=COEFF1*( 2D0/3D0*B)
      AE(14)=COEFF1*(-4D0/3D0*(B+D))
      AE(15)=COEFF1*( 4D0/3D0*(A+B+D))
C
      AE(16)=COEFF1*(-2D0/3D0*(B+D))
      AE(17)=0D0
      AE(18)=COEFF1*(-2D0/3D0*(B+D))
      AE(19)=COEFF1*( 4D0/3D0*B)
      AE(20)=COEFF1*(-4D0/3D0*(A+B))
      AE(21)=COEFF1*( 4D0/3D0*(A+B+D))
C
C     BLOC A[2,1]
      AE(22)=0D0
      AE(23)=0D0
      AE(24)=0D0
      AE(25)=0D0
      AE(26)=0D0
      AE(27)=0D0
C     BLOC A[2,2]
      AE(28)=COEFF1*(0.5D0*(A+D)+B)
C
C     BLOC A[2,1]
      AE(29)=0D0
      AE(30)=0D0
      AE(31)=0D0
      AE(32)=0D0
      AE(33)=0D0
      AE(34)=0D0
C     BLOC A[2,2]
      AE(35)=COEFF1/6D0*(A+B)
      AE(36)=COEFF1*0.5D0*A
C
C     BLOC A[2,1]
      AE(37)=0D0
      AE(38)=0D0
      AE(39)=0D0
      AE(40)=0D0
      AE(41)=0D0
      AE(42)=0D0
C     BLOC A[2,2]
      AE(43)= COEFF1/6D0*(B+D)
      AE(44)=-COEFF1/6D0*B
      AE(45)= COEFF1*0.5D0*D
C
C     BLOC A[2,1]
      AE(46)=0D0
      AE(47)=0D0
      AE(48)=0D0
      AE(49)=0D0
      AE(50)=0D0
      AE(51)=0D0
C     BLOC A[2,2]
      AE(52)=COEFF1*(-2D0/3D0*(A+B))
      AE(53)=COEFF1*(-2D0/3D0*(A+B))
      AE(54)=0D0
      AE(55)=COEFF1*( 4D0/3D0*(A+B+D))
C
C     BLOC A[2,1]
      AE(56)=0D0
      AE(57)=0D0
      AE(58)=0D0
      AE(59)=0D0
      AE(60)=0D0
      AE(61)=0D0
C     BLOC A[2,2]
      AE(62)=0D0
      AE(63)=COEFF1*( 2D0/3D0*B)
      AE(64)=COEFF1*( 2D0/3D0*B)
      AE(65)=COEFF1*(-4D0/3D0*(B+D))
      AE(66)=COEFF1*( 4D0/3D0*(A+B+D))
C
C     BLOC A[2,1]
      AE(67)=0D0
      AE(68)=0D0
      AE(69)=0D0
      AE(70)=0D0
      AE(71)=0D0
      AE(72)=0D0
C     BLOC A[2,2]
      AE(73)=COEFF1*(-2D0/3D0*(B+D))
      AE(74)=0D0
      AE(75)=COEFF1*(-2D0/3D0*(B+D))
      AE(76)=COEFF1*( 4D0/3D0*B)
      AE(77)=COEFF1*(-4D0/3D0*(A+B))
      AE(78)=COEFF1*( 4D0/3D0*(A+B+D))
C
C     BLOC A[3,1]
      AE(79)= (Y31+Y12)/6D0 * COPRES
      AE(80)=  0D0
      AE(81)=  0D0
      AE(82)= (Y12-Y31)/6D0 * COPRES
      AE(83)=-(Y31+Y12)/6D0 * COPRES
      AE(84)= (Y31-Y12)/6D0 * COPRES
C     BLOC A[3,2]
      AE(85)= (X13+X21)/6D0 * COPRES
      AE(86)=  0D0
      AE(87)=  0D0
      AE(88)= (X21-X13)/6D0 * COPRES
      AE(89)=-(X13+X21)/6D0 * COPRES
      AE(90)= (X13-X21)/6D0 * COPRES
C
C     BLOC A[3,3]
C     Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) dX sur l'EF NEF
      COEFF2 = EPSILON * VISCOS / ( 2D0 * DELTA )
      AE(91) = ( Y32 * Y32 + X32 * X32 ) * COEFF2
C
C     BLOC A[3,1]
      AE(92)= 0D0
      AE(93)= -Y31/6D0 * COPRES
      AE(94)= 0D0
      AE(95)=( Y31/6D0 + Y12/3D0) * COPRES
      AE(96)=(-Y31/6D0 - Y12/3D0) * COPRES
      AE(97)=( Y31/6D0 ) * COPRES
C     BLOC A[3,2]
      AE(98) = 0D0
      AE(99) = -X13/6D0 * COPRES
      AE(100)= 0D0
      AE(101)=( X13/6D0 + X21/3D0) * COPRES
      AE(102)=(-X13/6D0 - X21/3D0) * COPRES
      AE(103)=  X13/6D0 * COPRES
C
C     BLOC A[3,3]
      AE(104) = (-Y32 * Y31 + X32 * X13 ) * COEFF2
      AE(105) = ( Y31 * Y31 + X13 * X13 ) * COEFF2
C
C     BLOC A[3,1]
      AE(106)= 0D0
      AE(107)= 0D0
      AE(108)= -Y12/6D0 * COPRES
      AE(109)=  Y12/6D0 * COPRES
      AE(110)=(-Y12/6D0 - Y31/3D0) * COPRES
      AE(111)=( Y12/6D0 + Y31/3D0) * COPRES
C     BLOC A[3,2]
      AE(112)= 0D0
      AE(113)= 0D0
      AE(114)= -X21/6D0 * COPRES
      AE(115)=  X21/6D0 * COPRES
      AE(116)=(-X21/6D0 - X13/3D0) * COPRES
      AE(117)=( X21/6D0 + X13/3D0) * COPRES
C
C     BLOC A[3,3]
      AE(118) = (-Y32 * Y12 + X32 * X21 ) * COEFF2
      AE(119) = ( Y31 * Y12 + X13 * X21 ) * COEFF2
      AE(120) = ( Y12 * Y12 + X21 * X21 ) * COEFF2
C
      RETURN
      END
