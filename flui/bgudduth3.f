      SUBROUTINE BGUDDUTH3( NBSOM,  NBNOVI, XYZNOE,
     %                      NBNOEF, NBEF,   NONOEF, NONOSO,
     %                      NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                      VXVYVZ, CoForc, CoTeNL,   BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU SYSTEME
C ----     -CoGrPr LAPLACIEN P = Div( -CoForc Fomega + CoTeNL u. Grad u )
C          DIRECTEMENT SANS UTILISER LA FORMULE DE GREEN POUR REPORTER
C          LA DIVERGENCE SUR L'INTERPOLATION P1 DE LA FONCTION TEST tQ
C          POUR LE TETRAEDRE de TAYLOR-HOOD
C          PISO Version 2: Bg avec Div calculee sur la Vitesse  => DDP2
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C NOOBSF : NUMERO DE SURFACE DU FLUIDE
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES SURFACES
C CoForc : Coefficient du terme des FORCES
C CoTeNL : DENSITE DE MASSE also COEFFICIENT OF THE NON LINEAR TERM
C NOOBSF : NUMERO DE SURFACE DU FLUIDE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C VXVYVZ   : COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en X puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Y
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
C
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOVI, NBNOEF, NBEF,
     %                  NONOEF(NBEF,NBNOEF), NONOSO(NBNOVI),
     %                  NOOBVC, NUMIVO, NUMAVO
      DOUBLE PRECISION  CoForc, CoTeNL, VXVYVZ(NBNOVI,3), BG(NBSOM)
C
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      DOUBLE PRECISION  DELTAe, DF(3,3), DFM1(3,3), BeI, SJ, SL, SK, FL,
     %                  VKJ, VLJJ, X, Y, Z
      INTEGER           NOSOTE(4), NS1, NS2, NS3, NS4
      EQUIVALENCE      (NOSOTE(1),NS1),(NOSOTE(2),NS2),(NOSOTE(3),NS3),
     %                 (NOSOTE(4),NS4)
      DOUBLE PRECISION  FORCE(3,10)
      INTEGER           I, J, JJ, K, L, N, NN, NONOJ, NEF, IEFORC
C
C     INTEGRALES P1 DP2 SUR LE TETRAEDRE de REFERENCE
C     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy sur TETRAEDRE UNITE
ccc   DOUBLE PRECISION  P1DP23D(4,3,10)
      include"./incl/p1dp23d.inc"
C
C     INTEGRALES P1 DP2 DP2 SUR LE TETRAEDRE de REFERENCE
C     P1DP2DP2(i,k,j,l,m) = integrale P1i DP2j/dxk DP2m/dxl dx dy
ccc   DOUBLE PRECISION  P1DP2DP2(4,3,10,3,10)
      include"./incl/p1dp2dp23d.inc"
C
C     INTEGRALES P1 P2 DDP2 SUR LE TETRAEDRE de REFERENCE
C     P1DP2DP2(i,j,k,l,m) = integrale P1i P2j DDP2m/dxkdxl dx dy
ccc   DOUBLE PRECISION  P1P2DDP2(4,10,3,3,10)
      include"./incl/p1p2ddp23d.inc"
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 4 SOMMETS DU TETRAEDRE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
         NS4 = NONOEF(NEF,4)
C
C        CONSTRUCTION DE LA MATRICE DF
         X = XYZNOE(1,NS1)
         DF(1,1) = XYZNOE(1,NS2) - X
         DF(2,1) = XYZNOE(1,NS3) - X
         DF(3,1) = XYZNOE(1,NS4) - X
C
         Y = XYZNOE(2,NS1)
         DF(1,2) = XYZNOE(2,NS2) - Y
         DF(2,2) = XYZNOE(2,NS3) - Y
         DF(3,2) = XYZNOE(2,NS4) - Y
C
         Z = XYZNOE(3,NS1)
         DF(1,3) = XYZNOE(3,NS2) - Z
         DF(2,3) = XYZNOE(3,NS3) - Z
         DF(3,3) = XYZNOE(3,NS4) - Z
C
C        LE DETERMINANT DE DF
         DELTAe = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %          + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %          + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1
         DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) / DELTAe
         DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) ) / DELTAe
         DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) ) / DELTAe
C
         DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) ) / DELTAe
         DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) ) / DELTAe
         DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) ) / DELTAe
C
         DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) / DELTAe
         DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) ) / DELTAe
         DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) ) / DELTAe
C
C        CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES TAYLOR-HOOD
C        -----------------------------------------------------------
         IEFORC = LTDEVO(LPFORC,NOOBVC)
         IF( IEFORC .GT. 0 ) THEN
C           VALEUR DES EFFORTS VOLUMIQUES AUX 10 NOEUDS DU TETRAEDRE
            DO J=1,10
               NONOJ = NONOEF(NEF,J)
               X = XYZNOE(1,NONOJ)
               Y = XYZNOE(2,NONOJ)
               Z = XYZNOE(3,NONOJ)
               CALL REFORC( 4, NOOBVC, 3, X,Y,Z,  0D0,0D0,0D0,
     %                      LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
            ENDDO
         ENDIF
C
         DO I=1,4
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            BeI = 0D0
C
            DO L=1,3
C
C              COEFFICIENTS DES EFFORTS SURFACIQUES L
C              Integrale  P1 DivP2 Fel  dx
               FL = 0D0
               IF( IEFORC .GT. 0 ) THEN
                  DO J=1,10
                     DO N=1,3
                        FL = FL + DFM1(L,N) * P1DP23D(I,N,J)
     %                          * FORCE(L,J) * CoForc
                     ENDDO
                  ENDDO
               ENDIF
C
C              Integrale  P1  d/dxl Somme  uk. d/dxk  ul dx
C                                   k=1,3
               SL = 0D0
               DO K=1,3
C
                  SK = 0D0
                  DO J=1,10
C
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NONOEF(NEF,J)
                     VKJ = VXVYVZ( NONOEF(NEF,J), K )
C
C                    LE 1-ER TERME DERIVE
                     SJ = 0D0
                     DO N=1,3
                        DO JJ=1,10
C
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NONOEF(NEF,JJ)
                           VLJJ = VXVYVZ( NONOEF(NEF,JJ), L )
C
                           DO NN=1,3
                              SJ = SJ + DFM1(L,N) * DFM1(K,NN)
     %                                * P1DP2DP2(I,N,J,NN,JJ)
     %                                * VLJJ
                           ENDDO
C
                        ENDDO
                     ENDDO
C
C                    LE SECOND TERME DERIVE
                     DO JJ=1,10
C
C                       COMPOSANTE L DE LA VITESSE AU NOEUD NONOEF(NEF,JJ)
                        VLJJ = VXVYVZ( NONOEF(NEF,JJ), L )
C
                        DO N=1,3
                           DO NN=1,3
                              SJ = SJ + DFM1(L,N)
     %                                * P1P2DDP2(I,J,N,NN,JJ)
     %                                * DFM1(K,NN) * VLJJ
                           ENDDO
                        ENDDO
                     ENDDO
C
                     SK = SK + SJ * VKJ
C
                  ENDDO
C
                  SL = SL + SK
C
               ENDDO
C
               BeI = BeI - FL + CoTeNL * SL
C
            ENDDO
C
C           ASSEMBLAGE DE Be(I) DANS BG( NONOSO( NOSOTE(I) ) )
            NONOJ = NONOSO( NOSOTE(I) )
            BG(NONOJ) = BG(NONOJ) + BeI * DELTAe
C
         ENDDO
C
 100  CONTINUE
C
      RETURN
      END
