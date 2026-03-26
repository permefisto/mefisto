      SUBROUTINE F3MP2P1( X, NOVOLU, NUMIVO, NUMAVO, LTDEVO,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE MASSE AE DU TRIANGLE  TAYLOR HOOD
C -----    LAGRANGE DE DEGRE 2 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1 POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA DENSITE DE MASSE EST SUPPOSEE CONSTANTE
C
C ENTREES:	
C --------
C X      : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C NOVOLU : NUMERO DE VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          (ICI LA MASSE)
C
C SORTIES:
C --------
C AE : MATRICE ELEMENTAIRE 34x34 STOCKEE SYMETRIQUE PLEINE  (10+10+10+4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    Juin 2007
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/donflu.inc"
C
      REAL               X(10,3)
      DOUBLE PRECISION   AE(595)
      DOUBLE PRECISION   VMASSE
      INTEGER            NOVOLU, NUMIVO, NUMAVO
      INTEGER            LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER            I
      DOUBLE PRECISION   X1, Y1, Z1
      DOUBLE PRECISION   ABS, DETM33, DELTA
      DOUBLE PRECISION   XYD(3)
      DOUBLE PRECISION   A, B, C, D, E, F, G, COEF
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
C     RECHERCHE DE LA DENSITE DE MASSE AU BARYCENTRE DU TETRAEDRE
C     ===========================================================
      XYD(1) = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      XYD(2) = ( X(1,2) + X(2,2) + X(3,2) + X(4,2) ) * 0.25D0
      XYD(3) = ( X(1,3) + X(2,3) + X(3,3) + X(4,3) ) * 0.25D0
      CALL REMASS2( 4, NOVOLU, XYD(1), XYD(2), XYD(3),
     %              LTDEVO(LPMASS,NOVOLU), VMASSE )
C
C     LES COEFFICIENTS DE LA MATRICE DE MASSE
C     =======================================
      COEF = VMASSE * DELTA
C
      A =  COEF / 420.D0
      G =  COEF / 315.D0
      B =  4.D0 * G
      C =  COEF / 2520.D0
      D = -COEF /  630.D0
      E = -A
      F =  2.D0 * G
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
C     ==================================================================
C     BLOC A[1,1]  integrale P2i P2j * rho moyen * delta
      AE(1)  = A
C
      AE(2)  = C
      AE(3)  = A
C
      AE(4)  = C
      AE(5)  = C
      AE(6)  = A
C
      AE(7)  = C
      AE(8)  = C
      AE(9)  = C
      AE(10) = A
C
      AE(11) = D
      AE(12) = D
      AE(13) = E
      AE(14) = E
      AE(15) = B
C
      AE(16) = E
      AE(17) = D
      AE(18) = D
      AE(19) = E
      AE(20) = F
      AE(21) = B
C
      AE(22) = D
      AE(23) = E
      AE(24) = D
      AE(25) = E
      AE(26) = F
      AE(27) = F
      AE(28) = B
C
      AE(29) = D
      AE(30) = E
      AE(31) = E
      AE(32) = D
      AE(33) = F
      AE(34) = G
      AE(35) = F
      AE(36) = B
C
      AE(37) = E
      AE(38) = D
      AE(39) = E
      AE(40) = D
      AE(41) = F
      AE(42) = F
      AE(43) = G
      AE(44) = F
      AE(45) = B
C
      AE(46) = E
      AE(47) = E
      AE(48) = D
      AE(49) = D
      AE(50) = G
      AE(51) = F
      AE(52) = F
      AE(53) = F
      AE(54) = F
      AE(55) = B
C
C     MISE A ZERO GENERALE RESTANTE
      DO I = 56, 595
         AE(I) = 0.D0
      ENDDO
C
C     BLOC A[2,2]  integrale P2i P2j * rho moyen * delta
      AE( 66) = A
C
      AE( 77) = C
      AE( 78) = A
C
      AE( 89) = C
      AE( 90) = C
      AE( 91) = A
C
      AE(102) = C
      AE(103) = C
      AE(104) = C
      AE(105) = A
C
      AE(116) = D
      AE(117) = D
      AE(118) = E
      AE(119) = E
      AE(120) = B
C
      AE(131) = E
      AE(132) = D
      AE(133) = D
      AE(134) = E
      AE(135) = F
      AE(136) = B
C
      AE(147) = D
      AE(148) = E
      AE(149) = D
      AE(150) = E
      AE(151) = F
      AE(152) = F
      AE(153) = B
C
      AE(164) = D
      AE(165) = E
      AE(166) = E
      AE(167) = D
      AE(168) = F
      AE(169) = G
      AE(170) = F
      AE(171) = B
C
      AE(182) = E
      AE(183) = D
      AE(184) = E
      AE(185) = D
      AE(186) = F
      AE(187) = F
      AE(188) = G
      AE(189) = F
      AE(190) = B
C
      AE(201) = E
      AE(202) = E
      AE(203) = D
      AE(204) = D
      AE(205) = G
      AE(206) = F
      AE(207) = F
      AE(208) = F
      AE(209) = F
      AE(210) = B
C
C     BLOC A[3,3]  integrale P2i P2j * rho moyen * delta
      AE(231) = A
C
      AE(252) = C
      AE(253) = A
C
      AE(274) = C
      AE(275) = C
      AE(276) = A
C
      AE(297) = C
      AE(298) = C
      AE(299) = C
      AE(300) = A
C
      AE(321) = D
      AE(322) = D
      AE(323) = E
      AE(324) = E
      AE(325) = B
C
      AE(346) = E
      AE(347) = D
      AE(348) = D
      AE(349) = E
      AE(350) = F
      AE(351) = B
C
      AE(372) = D
      AE(373) = E
      AE(374) = D
      AE(375) = E
      AE(376) = F
      AE(377) = F
      AE(378) = B
C
      AE(399) = D
      AE(400) = E
      AE(401) = E
      AE(402) = D
      AE(403) = F
      AE(404) = G
      AE(405) = F
      AE(406) = B
C
      AE(427) = E
      AE(428) = D
      AE(429) = E
      AE(430) = D
      AE(431) = F
      AE(432) = F
      AE(433) = G
      AE(434) = F
      AE(435) = B
C
      AE(456) = E
      AE(457) = E
      AE(458) = D
      AE(459) = D
      AE(460) = G
      AE(461) = F
      AE(462) = F
      AE(463) = F
      AE(464) = F
      AE(465) = B
C
      RETURN
      END
