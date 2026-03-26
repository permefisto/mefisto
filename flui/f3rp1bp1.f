      SUBROUTINE F3RP1BP1( X,
     %                     NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                     AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE VISCOSITE + PRESSION
C -----    DU TETRAEDRE BREZZI FORTIN AVEC INTERPOLATIONS
C          P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          VISCOSITE et COEFFICIENT de la PRESSION CALCULES
C          AU BARYCENTRE DES 4 SOMMETS DU TETRAEDRE
C          LES INTEGRALES DES PRODUITS DE POLYNOMES SONT INTEGRES
C          EXACTEMENT SUR LE TETRAEDRE DE REFERENCE
C
C ENTREES:	
C --------
C X      : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C NOOBVO : NUMERO DE VOLUME DE CET ELEMENT FINI
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES FLUIDES
C
C SORTIES:
C --------
C AE     : MATRICE PROFIL ELEMENTAIRE 19x19  AVANT GAUSS FRONTALE
C          SUR LES DL VITESSE DU BARYCENTRE  5 10 15
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Decembre 2008
C MODIFS: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray    Mai 2010
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/donflu.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      DOUBLE PRECISION  EPSILON
      PARAMETER        (EPSILON=1D-7)
C     ATTENTION: VALEUR CRUCIALE 1D-10 pour la cavity3d => MAUVAIS RESULTATS!
C
      INTEGER           NOOBVO, NUMIVO, NUMAVO
      INTEGER           LTDEVO( 1:MXDOFL, NUMIVO:NUMAVO )
      REAL              X(4,3)
      DOUBLE PRECISION  AE(190)
C
      INTEGER           I, J, K, L, M1, M2, M3
      DOUBLE PRECISION  VISCOS, VISDEL, COEFPRES, COEFPRED, COPRES
      DOUBLE PRECISION  DELTA, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION  XD, YD, ZD, S
C
C     DOUBLE PRECISION  TDPDP(5,3,5,3)
C     TDPDP(i,k,j,l) = integrale dP1Bi/dxk dP1Bj/dxl dX
      include"./incl/tdp53dp53.inc"
C
C     DOUBLE PRECISION  DP1Bla3d(5,3,4)
C     DP1Bla3d(i,k,j)   = integrale dP1Bi/dxk Lambdaj dX
      include"./incl/dp1bla3d.inc"
C
C     MISE A ZERO GENERALE EXCEPTE LE PREMIER BLOC DIAGONAL
C     =====================================================
      DO I = 16, 186
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
C     LES 3 BLOCS DIAGONAUX 5x5 Vitesse Vitesse DE LA MATRICE ELEMENTAIRE
C     STOCKAGE DE LA PARTIE TRIANGULAIRE INFERIEURE LIGNE PAR LIGNE
C     ===================================================================
C     VISCOSITE * DELTA
      VISDEL = VISCOS * DELTA
      M1 =  0
      M2 = 15
      M3 = 55
      DO I=1,5
         M2 = M2 +  5
         M3 = M3 + 10
         DO J=1,I
C
C           LE BLOC DIAGONAL SYMETRIQUE [V1,V1]
            M1 = M1 + 1
C           LE BLOC DIAGONAL SYMETRIQUE [V2,V2]
            M2 = M2 + 1
C           LE BLOC DIAGONAL SYMETRIQUE [V3,V3]
            M3 = M3 + 1
C
C           CALCUL DU COEFFICIENT Ae(I,J)
            S = 0D0
            DO K=1,3
               DO L=1,3
                  S = S + TDPDP(I,K,J,L) * TDFDF(K,L)
               ENDDO
            ENDDO
            S = S * VISDEL
C
            AE(M1) = S
            AE(M2) = S
            AE(M3) = S
C
         ENDDO
      ENDDO
C
C     LES BLOCS A41 A42 A43  MATRICES 4x5 DE LA MATRICE ELEMENTAIRE
C     =============================================================
      COEFPRES = DELTA * COPRES
      M1 = 120
      M2 = 125
      M3 = 130
      DO I=1,4
         DO J=1,5
C
C           LE COEFFICIENT DU BLOC NON DIAGONAL [P, V1]
            M1 = M1 + 1
C           LE COEFFICIENT DU BLOC NON DIAGONAL [P, V2]
            M2 = M2 + 1
C           LE COEFFICIENT DU BLOC NON DIAGONAL [P, V3]
            M3 = M3 + 1
C
            DO K=1,3
C              LE MOINS POUR LE TERME - INTEGRALE div v  coef p dx
               S = DP1Bla3d(J,K,I) * COEFPRES
               AE(M1) = AE(M1) - DFM1(1,K) * S
               AE(M2) = AE(M2) - DFM1(2,K) * S
               AE(M3) = AE(M3) - DFM1(3,K) * S
            ENDDO
C
         ENDDO
         M1 = M3 + I
         M2 = M1 + 5
         M3 = M2 + 5
      ENDDO
C
C     PENALISATION RELATIVE DU BLOC DIAGONAL 4x4 EN PRESSION
C     ------------------------------------------------------
C     COEFFICIENT DIAGONAL DU BLOC
      COEFPRED = EPSILON * VISDEL / 60D0
C     COEFFICIENT NON  DIAGONAL
      COEFPRES = COEFPRED / 2D0
C
C     LE BLOC INTEGRALE EPSILON P Q dx dy
      AE(136) = COEFPRED
C
      AE(152) = COEFPRES
      AE(153) = COEFPRED
C
      AE(169) = COEFPRES
      AE(170) = COEFPRES
      AE(171) = COEFPRED
C
      AE(187) = COEFPRES
      AE(188) = COEFPRES
      AE(189) = COEFPRES
      AE(190) = COEFPRED
C
      RETURN
      END
