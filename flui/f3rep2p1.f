      SUBROUTINE F3REP2P1( DT,  X, TP2P2,  TDP2DP2,
     %                     NOVOLU, NUMIVO, NUMAVO, LTDEVO,
     %                     AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE VISCOSITE AE DU TETRAEDRE TAYLOR HOOD
C -----    Ae = ( [Me(1)] + dt ([Ke(Mu/Rho)]+[Re(omega)]))
C          LAGRANGE DE DEGRE 2 POUR LA VITESSE.
C          LAGRANGE DE DEGRE 1 POUR LA PRESSION
C          INTEGRATION EN XYZ EXACTE CAR LA VISCOSITE EST SUPPOSEE CONSTANTE
C          INTEGRATION EN TEMPS PAR SCHEMA D'EULER LE PLUS STABLE (THETA=1)
C
C ENTREES:	
C --------
C DT     : PAS DE TEMPS DU SCHEMA D'EULER D'INTEGRATION EN TEMPS
C X      : LES 3 COORDONNEES DES 10 POINTS DU TETRAEDRE
C TP2P2  : INTEGRALE SUR L'EF REFERENCE des pi pj dX avec p=P2
C TDP2DP2: INTEGRALE SUR L'EF REFERENCE des dpi/dxk dpj/dxl dX avec p=P2
C NOVOLU : NUMERO DE L'OBJET VOLUME DE CE FLUIDE
C NUMIVO : NUMERO MINIMAL DES VOLUMES
C NUMAVO : NUMERO MAXIMAL DES VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C
C SORTIES:
C --------
C AE  : MATRICE ELEMENTAIRE 34x34 STOCKEE SYMETRIQUE PLEINE (10+10+10+4)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Mai 2010
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
C
      DOUBLE PRECISION   EPSILON
      PARAMETER         (EPSILON=1D-7)
C
      REAL               DT, X(10,3)
      DOUBLE PRECISION   TP2P2(10,10), TDP2DP2(3,10,3,10)
      DOUBLE PRECISION   AE(34,34)
      INTEGER            NOVOLU, NUMIVO, NUMAVO
      INTEGER            LTDEVO( 1:MXDOFL, NUMIVO:NUMAVO )
C
      INTEGER            I, J, K, L
      DOUBLE PRECISION   VISCOS, VISDEL, VMASSE, VITANG(3), V
      DOUBLE PRECISION   DELTA, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION   XD, YD, ZD, C1, C2, S, F1, F2, F3
      DOUBLE PRECISION   D15, D30, D40, D120, DTDELTA
C
C     LES 3 COORDONNEES DU BARYCENTRE DU TETRAEDRE
      XD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      YD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
      ZD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
C
C     COEFFICIENT DE LA DENSITE DE MASSE AU BARYCENTRE DU TETRAEDRE
      CALL REMASS2( 4, NOVOLU, XD, YD, ZD,
     %              LTDEVO(LPMASS,NOVOLU), VMASSE )
C
C     COEFFICIENT DE LA VISCOSITE AU BARYCENTRE DU TETRAEDRE
      CALL REVISC( 4, NOVOLU, XD, YD, ZD,
     %             LTDEVO(LPVISC,NOVOLU), VISCOS )
C
C     VECTEUR(3) DE VITESSE ANGULAIRE AU BARYCENTRE DU TETRAEDRE
      CALL REVIAN( 4, NOVOLU, XD, YD, ZD,
     %             LTDEVO(LPVIAN,NOVOLU), VITANG )
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
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
C     CALCUL DES COEFFICIENTS DE LA SOUS MATRICE
C     tDF-1 * DF-1 * DELTA * VISCOS / VMASSE
      DTDELTA = DT * DELTA
      VISDEL  = VISCOS * DTDELTA / VMASSE
      TDFDF(1,1) = ( DFM1(1,1) * DFM1(1,1)
     %             + DFM1(2,1) * DFM1(2,1)
     %             + DFM1(3,1) * DFM1(3,1) ) * VISDEL
C
      TDFDF(2,1) = ( DFM1(1,2) * DFM1(1,1)
     %             + DFM1(2,2) * DFM1(2,1)
     %             + DFM1(3,2) * DFM1(3,1) ) * VISDEL
      TDFDF(1,2) = TDFDF(2,1)
C
      TDFDF(3,1) = ( DFM1(1,3) * DFM1(1,1)
     %             + DFM1(2,3) * DFM1(2,1)
     %             + DFM1(3,3) * DFM1(3,1) ) * VISDEL
      TDFDF(1,3) = TDFDF(3,1)
C
      TDFDF(2,2) = ( DFM1(1,2) * DFM1(1,2)
     %             + DFM1(2,2) * DFM1(2,2)
     %             + DFM1(3,2) * DFM1(3,2) ) * VISDEL
C
      TDFDF(3,2) = ( DFM1(1,3) * DFM1(1,2)
     %             + DFM1(2,3) * DFM1(2,2)
     %             + DFM1(3,3) * DFM1(3,2) ) * VISDEL
      TDFDF(2,3) = TDFDF(3,2)
C
      TDFDF(3,3) = ( DFM1(1,3) * DFM1(1,3)
     %             + DFM1(2,3) * DFM1(2,3)
     %             + DFM1(3,3) * DFM1(3,3) ) * VISDEL
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE GAUCHE A DROITE
C     ==================================================================
      DO I = 1, 10
         DO J = 1, I
C
C           LE  BLOC DIAGONAL 10x10 V1 x V1 de la MATRICE de MASSE
C           BLOC A[1,1]  Integrale Pi Pj * delta
            S = TP2P2(I,J) * DELTA
C
C           LE  BLOC DIAGONAL 10x10 V1 x V1 de la MATRICE de VISCOSITE
            DO K=1,3
               DO L=1,3
                  S = S + TDP2DP2(K,I,L,J) * TDFDF(K,L)
               ENDDO
            ENDDO
C
C           LES BLOCS DIAGONAUX DES VITESSES V2 V2 et V3 V3
            AE(I   ,J   ) = S
            AE(I+10,J+10) = S
            AE(I+20,J+20) = S
C
         ENDDO
      ENDDO
C
C     LE BLOC V2 V1 EST NUL
C     =====================
      DO I = 11, 20
         DO J = 1, 10
            AE(I,J) = 0D0
         ENDDO
      ENDDO
C
C     LES BLOCS V3 V1 et V3 V2 SONT NULS
C     ==================================
      DO I = 21, 30
         DO J = 1, 20
            AE(I,J) = 0D0
         ENDDO
      ENDDO
C
C     LES BLOCS  A41 A42 A43  MATRICES 4x10 DE LA MATRICE de VISCOSITE
C     ================================================================
      C1   = DTDELTA / VMASSE
      D15  = 1D0 /  15D0
      D30  = 1D0 /  30D0
      D40  = 1D0 /  40D0
      D120 = 1D0 / 120D0
C
C     BLOC [P,V1] LIGNE PRESSION ST1
C     ------------------------------
      F1 = - DFM1(1,1) * C1
      F2 = - DFM1(1,2) * C1
      F3 = - DFM1(1,3) * C1
C
      AE(31, 1) = - ( F1 + F2 + F3 ) * D40
      AE(31, 2) = - F1 * D120
      AE(31, 3) = - F2 * D120
      AE(31, 4) = - F3 * D120
C
      AE(31, 5) =  ( F1 - F2 - F3 ) * D30
      AE(31, 6) =  ( F1 + F2      ) * D30
      AE(31, 7) =  (-F1 + F2 - F3 ) * D30
      AE(31, 8) =  (-F1 - F2 + F3 ) * D30
      AE(31, 9) =  ( F1      + F3 ) * D30
      AE(31,10) =  (      F2 + F3 ) * D30
C
C     BLOC [P,V1] LIGNE PRESSION ST2
      AE(32, 1) = ( F1 + F2 + F3 ) * D120
      AE(32, 2) =   F1             * D40
      AE(32, 3) =      - F2        * D120
      AE(32, 4) =            - F3  * D120
C
      AE(32, 5) = - F1 * D30 - (F2+F3) * D15
      AE(32, 6) =   F1 * D30 + F2* D15
      AE(32, 7) = -(F1 + F3) * D30
      AE(32, 8) = -(F1 + F2) * D30
      AE(32, 9) =   F1 * D30 + F3 * D15
      AE(32,10) =  (F2 + F3) * D30
C
C     BLOC [P,V1] LIGNE PRESSION ST3
      AE(33, 1) = ( F1 + F2 + F3 ) * D120
      AE(33, 2) = - F1 * D120
      AE(33, 3) =   F2 * D40
      AE(33, 4) = - F3 * D120
C
      AE(33, 5) = -( F2 + F3 ) * D30
      AE(33, 6) =    F1 * D15 + F2 * D30
      AE(33, 7) = -( F1+ F3 ) * D15 - F2 * D30
      AE(33, 8) = -( F1+ F2 ) * D30
      AE(33, 9) =  ( F1+ F3 ) * D30
      AE(33,10) =    F2 * D30 + F3 * D15
C
C     BLOC [P,V1] LIGNE PRESSION ST4
      AE(34, 1) = ( F1 + F2 + F3 ) * D120
      AE(34, 2) = - F1 * D120
      AE(34, 3) = - F2 * D120
      AE(34, 4) =   F3 * D40
C
      AE(34, 5) = -( F2 + F3 ) * D30
      AE(34, 6) =  ( F1 + F2 ) * D30
      AE(34, 7) = -( F1 + F3 ) * D30
      AE(34, 8) = -( F1 + F2 ) * D15 - F3 * D30
      AE(34, 9) =    F1 * D15 + F3 * D30
      AE(34,10) =    F2 * D15 + F3 * D30
C
C
C     BLOC [P,V2] LIGNE PRESSION ST1
C     ------------------------------
      F1 = - DFM1(2,1) * C1
      F2 = - DFM1(2,2) * C1
      F3 = - DFM1(2,3) * C1
C
      AE(31,11) = - ( F1 + F2 + F3 ) * D40
      AE(31,12) = - F1 * D120
      AE(31,13) = - F2 * D120
      AE(31,14) = - F3 * D120
C
      AE(31,15) =  ( F1 - F2 - F3 ) * D30
      AE(31,16) =  ( F1 + F2      ) * D30
      AE(31,17) =  (-F1 + F2 - F3 ) * D30
      AE(31,18) =  (-F1 - F2 + F3 ) * D30
      AE(31,19) =  ( F1      + F3 ) * D30
      AE(31,20) =  (      F2 + F3 ) * D30
C
C     BLOC [P,V2] LIGNE  PRESSION ST2
      AE(32,11) = ( F1 + F2 + F3 ) * D120
      AE(32,12) =   F1             * D40
      AE(32,13) =      - F2        * D120
      AE(32,14) =            - F3  * D120
C
      AE(32,15) = - F1 * D30 - (F2+F3) * D15
      AE(32,16) =   F1 * D30 + F2* D15
      AE(32,17) = -(F1 + F3) * D30
      AE(32,18) = -(F1 + F2) * D30
      AE(32,19) =   F1 * D30 + F3 * D15
      AE(32,20) =  (F2 + F3) * D30
C
C     BLOC [P,V2] LIGNE PRESSION ST3
      AE(33,11) = ( F1 + F2 + F3 ) * D120
      AE(33,12) = - F1 * D120
      AE(33,13) =   F2 * D40
      AE(33,14) = - F3 * D120
C
      AE(33,15) = -( F2 + F3 ) * D30
      AE(33,16) =    F1 * D15 + F2 * D30
      AE(33,17) = -( F1+ F3 ) * D15 - F2 * D30
      AE(33,18) = -( F1+ F2 ) * D30
      AE(33,19) =  ( F1+ F3 ) * D30
      AE(33,20) =    F2 * D30 + F3 * D15
C
C     BLOC [P,V2] LIGNE PRESSION ST4
      AE(34,11) = ( F1 + F2 + F3 ) * D120
      AE(34,12) = - F1 * D120
      AE(34,13) = - F2 * D120
      AE(34,14) =   F3 * D40
C
      AE(34,15) = -( F2 + F3 ) * D30
      AE(34,16) =  ( F1 + F2 ) * D30
      AE(34,17) = -( F1 + F3 ) * D30
      AE(34,18) = -( F1 + F2 ) * D15 - F3 * D30
      AE(34,19) =    F1 * D15 + F3 * D30
      AE(34,20) =    F2 * D15 + F3 * D30
C
C
C     BLOC [P,V3] LIGNE PRESSION ST1
C     ------------------------------
      F1 = - DFM1(3,1) * C1
      F2 = - DFM1(3,2) * C1
      F3 = - DFM1(3,3) * C1
C
      AE(31,21) = - ( F1 + F2 + F3 ) * D40
      AE(31,22) = - F1 * D120
      AE(31,23) = - F2 * D120
      AE(31,24) = - F3 * D120
C
      AE(31,25) =  ( F1 - F2 - F3 ) * D30
      AE(31,26) =  ( F1 + F2      ) * D30
      AE(31,27) =  (-F1 + F2 - F3 ) * D30
      AE(31,28) =  (-F1 - F2 + F3 ) * D30
      AE(31,29) =  ( F1      + F3 ) * D30
      AE(31,30) =  (      F2 + F3 ) * D30
C
C     BLOC [P,V3] LIGNE PRESSION ST2
      AE(32,21) = ( F1 + F2 + F3 ) * D120
      AE(32,22) =   F1             * D40
      AE(32,23) =      - F2        * D120
      AE(32,24) =            - F3  * D120
C
      AE(32,25) = - F1 * D30 - (F2+F3) * D15
      AE(32,26) =   F1 * D30 + F2* D15
      AE(32,27) = -(F1 + F3) * D30
      AE(32,28) = -(F1 + F2) * D30
      AE(32,29) =   F1 * D30 + F3 * D15
      AE(32,30) =  (F2 + F3) * D30
C
C     BLOC [P,V3] LIGNE PRESSION ST3
      AE(33,21) = ( F1 + F2 + F3 ) * D120
      AE(33,22) = - F1 * D120
      AE(33,23) =   F2 * D40
      AE(33,24) = - F3 * D120
C
      AE(33,25) = -( F2 + F3 ) * D30
      AE(33,26) =    F1 * D15 + F2 * D30
      AE(33,27) = -( F1+ F3 ) * D15 - F2 * D30
      AE(33,28) = -( F1+ F2 ) * D30
      AE(33,29) =  ( F1+ F3 ) * D30
      AE(33,30) =    F2 * D30 + F3 * D15
C
C     BLOC [P,V3] LIGNE PRESSION ST4
      AE(34,21) = ( F1 + F2 + F3 ) * D120
      AE(34,22) = - F1 * D120
      AE(34,23) = - F2 * D120
      AE(34,24) =   F3 * D40
C
      AE(34,25) = -( F2 + F3 ) * D30
      AE(34,26) =  ( F1 + F2 ) * D30
      AE(34,27) = -( F1 + F3 ) * D30
      AE(34,28) = -( F1 + F2 ) * D15 - F3 * D30
      AE(34,29) =    F1 * D15 + F3 * D30
      AE(34,30) =    F2 * D15 + F3 * D30
C
C     PENALISATION RELATIVE DU BLOC INTEGRALE EPSILON P Q dx dy
C     DIAGONAL 4x4 DE LA PRESSION
C     =========================================================
C     COEFFICIENT DIAGONAL DU BLOC
      C2 = EPSILON * VISCOS * C1 / 60D0
C     COEFFICIENT NON  DIAGONAL
      C1 = C2 / 2D0
C
C     LE BLOC DIAGONAL SYMETRIQUE FINAL [P,P]
      AE(31,31) = C2
C
      AE(32,31) = C1
      AE(32,32) = C2
C
      AE(33,31) = C1
      AE(33,32) = C1
      AE(33,33) = C2
C
      AE(34,31) = C1
      AE(34,32) = C1
      AE(34,33) = C1
      AE(34,34) = C2
C
C     LA PARTIE SYMETRIQUE EST TRANSPOSEE
C     ===================================
      DO I = 2, 34
         DO J = 1, I-1
            AE(J,I) = AE(I,J)
         ENDDO
      ENDDO
C
C     LES 6 BLOCS NON SYMETRIQUES DUS A LA ROTATION DU DOMAINE
C     LES 3 BLOCS DIAGONAUX SONT NULS
C     ========================================================
      C1 = 2D0 * DTDELTA
      DO J = 1, 10
         DO I = 1, 10
            V = TP2P2(I,J) * C1
C
C           BLOC V1 x V2
            AE(I,J+10)    = AE(I,J+10)    - VITANG(3) * V
C           BLOC V1 x V3
            AE(I,J+20)    = AE(I,J+20)    + VITANG(2) * V
C
C           BLOC V2 x V1
            AE(I+10,J)    = AE(I+10,J)    + VITANG(3) * V
C           BLOC V2 x V3
            AE(I+10,J+20) = AE(I+10,J+20) - VITANG(1) * V
C
C           BLOC V3 x V1
            AE(I+20,J)    = AE(I+20,J)    - VITANG(2) * V
C           BLOC V3 x V2
            AE(I+20,J+10) = AE(I+20,J+10) + VITANG(1) * V
         ENDDO
      ENDDO
C
      RETURN
      END
