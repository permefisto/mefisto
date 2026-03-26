      SUBROUTINE F3RP2P1( X, DP2DP2, NOOBVO, NUMIVO, NUMAVO, LTDEVO, AE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE VISCOSITE AE DU TETRAEDRE TAYLOR HOOD
C -----    LAGRANGE DE DEGRE 2 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1 POUR LA PRESSION
C          INTEGRATION EXACTE CAR LA VISCOSITE EST SUPPOSEE CONSTANTE
C
C ENTREES:	
C --------
C X      : LES 3 COORDONNEES DES 10 POINTS DU TETRAEDRE
C DP2DP2 : INTEGRALE SUR L'EF REFERENCE des dpi/dxk dpj/dxl avec p=P2
C NOOBVO : NUMERO DE L'OBJET VOLUME DE CE FLUIDE
C NUMIVO : NUMERO MINIMAL DES VOLUMES
C NUMAVO : NUMERO MAXIMAL DES VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C
C SORTIES:
C --------
C AE  : MATRICE ELEMENTAIRE 34x34 STOCKEE SYMETRIQUE PLEINE (10+10+10+4)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    Juin 2007
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray Avril 2010
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
C
      DOUBLE PRECISION   EPSILON
      PARAMETER         (EPSILON=1D-7)
ccc      PARAMETER         (EPSILON=1D-10)
ccc      PARAMETER         (EPSILON=0D0)  boucle infinie dans crout!
C
      REAL               X(10,3)
      DOUBLE PRECISION   DP2DP2(3,10,3,10)
      DOUBLE PRECISION   AE(595)
      DOUBLE PRECISION   VISCOS, VISDEL, COPRES
      INTEGER            NOOBVO, NUMIVO, NUMAVO
      INTEGER            LTDEVO( 1:MXDOFL, NUMIVO:NUMAVO )
C
      INTEGER            I, J, M1, M2, M3
      DOUBLE PRECISION   DELTA, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION   XD, YD, ZD, C1, C2, S, F1, F2, F3
      DOUBLE PRECISION   D15, D30, D40, D120
C
C     MISE A ZERO GENERALE EXCEPTE LE PREMIER BLOC DIAGONAL
C     =====================================================
      DO I = 56, 591
         AE(I)=0D0
      ENDDO
C
C     RECHERCHE DE LA VISCOSITE AU BARYCENTRE DU TETRAEDRE
C     ====================================================
      XD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      YD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      ZD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      CALL REVISC( 4, NOOBVO, XD, YD, ZD,
     %             LTDEVO(LPVISC,NOOBVO), VISCOS )
C
C     RECHERCHE DU COEFFICIENT SUR LE GRADIENT DE LA PRESSION
C     =======================================================
      IF( LTDEVO(LPCPRE,NOOBVO) .GT. 0 ) THEN
C        IL EXISTE UN COEFFICIENT DEVANT LA PRESSION
         CALL RECPRE( 4, NOOBVO, XD, YD, ZD,
     %                LTDEVO(LPCPRE,NOOBVO), COPRES )
      ELSE
C        IL N'EXISTE PAS DE COEFFICIENT DEVANT LA PRESSION
         COPRES = 1D0
      ENDIF
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
      XD = X(1,1)
      YD = X(1,2)
      ZD = X(1,3)
      DF(1,1) = X(2,1) - XD
      DF(1,2) = X(2,2) - YD
      DF(1,3) = X(2,3) - ZD
C
      DF(2,1) = X(3,1) - XD
      DF(2,2) = X(3,2) - YD
      DF(2,3) = X(3,3) - ZD
C
      DF(3,1) = X(4,1) - XD
      DF(3,2) = X(4,2) - YD
      DF(3,3) = X(4,3) - ZD
C
C     [DF]-1
      CALL M33INV( DF, DELTA, DFM1 )
C
C     t[DF]-1 [DF]-1  EST UNE MATRICE SYMETRIQUE
      TDFDF(1,1) = DFM1(1,1) * DFM1(1,1)
     %           + DFM1(2,1) * DFM1(2,1)
     %           + DFM1(3,1) * DFM1(3,1)
C
      TDFDF(2,1) = DFM1(1,2) * DFM1(1,1)
     %           + DFM1(2,2) * DFM1(2,1)
     %           + DFM1(3,2) * DFM1(3,1)
      TDFDF(1,2) = TDFDF(2,1)
C
      TDFDF(3,1) = DFM1(1,3) * DFM1(1,1)
     %           + DFM1(2,3) * DFM1(2,1)
     %           + DFM1(3,3) * DFM1(3,1)
      TDFDF(1,3) = TDFDF(3,1)
C
      TDFDF(2,2) = DFM1(1,2) * DFM1(1,2)
     %           + DFM1(2,2) * DFM1(2,2)
     %           + DFM1(3,2) * DFM1(3,2)
C
      TDFDF(3,2) = DFM1(1,3) * DFM1(1,2)
     %           + DFM1(2,3) * DFM1(2,2)
     %           + DFM1(3,3) * DFM1(3,2)
      TDFDF(2,3) = TDFDF(3,2)
C
      TDFDF(3,3) = DFM1(1,3) * DFM1(1,3)
     %           + DFM1(2,3) * DFM1(2,3)
     %           + DFM1(3,3) * DFM1(3,3)
C
C     LES 3 BLOCS DIAGONAUX [Vi,Vi] SONT MULTIPLIES PAR VISCOS JACOBIEN
C     =================================================================
C     VISCOSITE * DELTA
      VISDEL = VISCOS * DELTA
      M1 =   0
      M2 =  55
      M3 = 210
      DO 40 I=1,10
         M2 = M2 + 10
         M3 = M3 + 20
         DO 30 J=1,I
C
C           LE BLOC DIAGONAL SYMETRIQUE [V1,V1]
            M1 = M1 + 1
            S  = 0D0
            DO 20 K=1,3
               DO 10 L=1,3
                  S = S + DP2DP2(K,I,L,J) * TDFDF(K,L)
 10            CONTINUE
 20         CONTINUE
            AE(M1) = S * VISDEL
C
C           LE BLOC DIAGONAL SYMETRIQUE [V2,V2]
            M2 = M2 + 1
            AE(M2) = AE(M1)
C
C           LE BLOC DIAGONAL SYMETRIQUE [V3,V3]
            M3 = M3 + 1
            AE(M3) = AE(M1)
C
 30      CONTINUE
 40   CONTINUE
C
C     PRISE EN COMPTE DU COEFFICIENT DU GRADIENT DE LA PRESSION
C     =========================================================
      C1   = DELTA * COPRES
      D15  = 1D0 /  15D0
      D30  = 1D0 /  30D0
      D40  = 1D0 /  40D0
      D120 = 1D0 / 120D0
C
C     BLOC [V1,P] LIGNE P1
C     --------------------
      F1 = - DFM1(1,1) * C1
      F2 = - DFM1(1,2) * C1
      F3 = - DFM1(1,3) * C1
C
      AE(466) = - ( F1 + F2 + F3 ) * D40
      AE(467) = - F1 * D120
      AE(468) = - F2 * D120
      AE(469) = - F3 * D120
C
      AE(470) =  ( F1 - F2 - F3 ) * D30
      AE(471) =  ( F1 + F2      ) * D30
      AE(472) =  (-F1 + F2 - F3 ) * D30
      AE(473) =  (-F1 - F2 + F3 ) * D30
      AE(474) =  ( F1      + F3 ) * D30
      AE(475) =  (      F2 + F3 ) * D30
C
C     BLOC [V1,P] LIGNE P2
      AE(497) = ( F1 + F2 + F3 ) * D120
      AE(498) =   F1             * D40
      AE(499) =      - F2        * D120
      AE(500) =            - F3  * D120
C
      AE(501) = - F1 * D30 - (F2+F3) * D15
      AE(502) =   F1 * D30 + F2* D15
      AE(503) = -(F1 + F3) * D30
      AE(504) = -(F1 + F2) * D30
      AE(505) =   F1 * D30 + F3 * D15
      AE(506) =  (F2 + F3) * D30
C
C     BLOC [V1,P] LIGNE P3
      AE(529) = ( F1 + F2 + F3 ) * D120
      AE(530) = - F1 * D120
      AE(531) =   F2 * D40
      AE(532) = - F3 * D120
C
      AE(533) = -( F2 + F3 ) * D30
      AE(534) =    F1 * D15 + F2 * D30
      AE(535) = -( F1+ F3 ) * D15 - F2 * D30
      AE(536) = -( F1+ F2 ) * D30
      AE(537) =  ( F1+ F3 ) * D30
      AE(538) =    F2 * D30 + F3 * D15
C
C     BLOC [V1,P] LIGNE P4
      AE(562) = ( F1 + F2 + F3 ) * D120
      AE(563) = - F1 * D120
      AE(564) = - F2 * D120
      AE(565) =   F3 * D40
C
      AE(566) = -( F2 + F3 ) * D30
      AE(567) =  ( F1 + F2 ) * D30
      AE(568) = -( F1 + F3 ) * D30
      AE(569) = -( F1 + F2 ) * D15 - F3 * D30
      AE(570) =    F1 * D15 + F3 * D30
      AE(571) =    F2 * D15 + F3 * D30
C
C
C     BLOC [V2,P] LIGNE P1
C     --------------------
      F1 = - DFM1(2,1) * C1
      F2 = - DFM1(2,2) * C1
      F3 = - DFM1(2,3) * C1
C
      AE(476) = - ( F1 + F2 + F3 ) * D40
      AE(477) = - F1 * D120
      AE(478) = - F2 * D120
      AE(479) = - F3 * D120
C
      AE(480) =  ( F1 - F2 - F3 ) * D30
      AE(481) =  ( F1 + F2      ) * D30
      AE(482) =  (-F1 + F2 - F3 ) * D30
      AE(483) =  (-F1 - F2 + F3 ) * D30
      AE(484) =  ( F1      + F3 ) * D30
      AE(485) =  (      F2 + F3 ) * D30
C
C     BLOC [V2,P] LIGNE P2
      AE(507) = ( F1 + F2 + F3 ) * D120
      AE(508) =   F1             * D40
      AE(509) =      - F2        * D120
      AE(510) =            - F3  * D120
C
      AE(511) = - F1 * D30 - (F2+F3) * D15
      AE(512) =   F1 * D30 + F2* D15
      AE(513) = -(F1 + F3) * D30
      AE(514) = -(F1 + F2) * D30
      AE(515) =   F1 * D30 + F3 * D15
      AE(516) =  (F2 + F3) * D30
C
C     BLOC [V2,P] LIGNE P3
      AE(539) = ( F1 + F2 + F3 ) * D120
      AE(540) = - F1 * D120
      AE(541) =   F2 * D40
      AE(542) = - F3 * D120
C
      AE(543) = -( F2 + F3 ) * D30
      AE(544) =    F1 * D15 + F2 * D30
      AE(545) = -( F1+ F3 ) * D15 - F2 * D30
      AE(546) = -( F1+ F2 ) * D30
      AE(547) =  ( F1+ F3 ) * D30
      AE(548) =    F2 * D30 + F3 * D15
C
C     BLOC [V2,P] LIGNE P4
      AE(572) = ( F1 + F2 + F3 ) * D120
      AE(573) = - F1 * D120
      AE(574) = - F2 * D120
      AE(575) =   F3 * D40
C
      AE(576) = -( F2 + F3 ) * D30
      AE(577) =  ( F1 + F2 ) * D30
      AE(578) = -( F1 + F3 ) * D30
      AE(579) = -( F1 + F2 ) * D15 - F3 * D30
      AE(580) =    F1 * D15 + F3 * D30
      AE(581) =    F2 * D15 + F3 * D30
C
C
C     BLOC [V3,P] LIGNE P1
C     --------------------
      F1 = - DFM1(3,1) * C1
      F2 = - DFM1(3,2) * C1
      F3 = - DFM1(3,3) * C1
C
      AE(486) = - ( F1 + F2 + F3 ) * D40
      AE(487) = - F1 * D120
      AE(488) = - F2 * D120
      AE(489) = - F3 * D120
C
      AE(490) =  ( F1 - F2 - F3 ) * D30
      AE(491) =  ( F1 + F2      ) * D30
      AE(492) =  (-F1 + F2 - F3 ) * D30
      AE(493) =  (-F1 - F2 + F3 ) * D30
      AE(494) =  ( F1      + F3 ) * D30
      AE(495) =  (      F2 + F3 ) * D30
C
C     BLOC [V3,P] LIGNE P2
      AE(517) = ( F1 + F2 + F3 ) * D120
      AE(518) =   F1             * D40
      AE(519) =      - F2        * D120
      AE(520) =            - F3  * D120
C
      AE(521) = - F1 * D30 - (F2+F3) * D15
      AE(522) =   F1 * D30 + F2* D15
      AE(523) = -(F1 + F3) * D30
      AE(524) = -(F1 + F2) * D30
      AE(525) =   F1 * D30 + F3 * D15
      AE(526) =  (F2 + F3) * D30
C
C     BLOC [V3,P] LIGNE P3
      AE(549) = ( F1 + F2 + F3 ) * D120
      AE(550) = - F1 * D120
      AE(551) =   F2 * D40
      AE(552) = - F3 * D120
C
      AE(553) = -( F2 + F3 ) * D30
      AE(554) =    F1 * D15 + F2 * D30
      AE(555) = -( F1+ F3 ) * D15 - F2 * D30
      AE(556) = -( F1+ F2 ) * D30
      AE(557) =  ( F1+ F3 ) * D30
      AE(558) =    F2 * D30 + F3 * D15
C
C     BLOC [V3,P] LIGNE P4
      AE(582) = ( F1 + F2 + F3 ) * D120
      AE(583) = - F1 * D120
      AE(584) = - F2 * D120
      AE(585) =   F3 * D40
C
      AE(586) = -( F2 + F3 ) * D30
      AE(587) =  ( F1 + F2 ) * D30
      AE(588) = -( F1 + F3 ) * D30
      AE(589) = -( F1 + F2 ) * D15 - F3 * D30
      AE(590) =    F1 * D15 + F3 * D30
      AE(591) =    F2 * D15 + F3 * D30
C
C     PENALISATION RELATIVE DU BLOC DIAGONAL PRESSION
C     -----------------------------------------------
      C2 = EPSILON * VISDEL / 60D0
      C1 = C2 / 2D0
C
C     LE BLOC DIAGONAL SYMETRIQUE FINAL [P,P]
      AE(496) = C2
C
      AE(527) = C1
      AE(528) = C2
C
      AE(559) = C1
      AE(560) = C1
      AE(561) = C2
C
      AE(592) = C1
      AE(593) = C1
      AE(594) = C1
      AE(595) = C2
C
      RETURN
      END
