      SUBROUTINE BGUDDUTH2( NBSOM,  NBNOVI, XYZNOE,
     %                      NBNOEF, NBEF,   NONOEF, NONOSO,
     %                      NOOBSF, NUMISU, NUMASU, LTDESU,
     %                      VXVY,   CoForc, CoTeNL,   BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU SYSTEME
C ----     -CoGrPr LAPLACIEN P = Div( -CoForc Fomega + CoTeNL u. Grad u )
C          DIRECTEMENT SANS UTILISER LA FORMULE DE GREEN POUR REPORTER
C          LA DIVERGENCE SUR L'INTERPOLATION P1 DE LA FONCTION TEST tQ
C          POUR LE TRIANGLE de TAYLOR-HOOD
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
C CoTeNL : DENSITE DE MASSE also COEFFICIENT OF THE NON LINEAR TERM
C NOOBSF : NUMERO DE SURFACE DU FLUIDE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C VXVY   : COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en X puis
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
     %                  NOOBSF, NUMISU, NUMASU
      DOUBLE PRECISION  CoForc, CoTeNL, VXVY(NBNOVI,2), BG(NBSOM)
C
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      DOUBLE PRECISION  DELTAe, DFM1(2,2), BeI, SJ, SL, SK, FL,
     %                  VKJ, VLJJ, X, Y, X21, Y21, X31, Y31
      INTEGER           NOSOTR(3), NS1, NS2, NS3
      EQUIVALENCE      (NOSOTR(1),NS1),(NOSOTR(2),NS2),(NOSOTR(3),NS3)
      DOUBLE PRECISION  FORCE(2,6)
      INTEGER           I, J, JJ, K, L, N, NN, NONOJ, NEF, IEFORC
C
C     INTEGRALES P1 DP2 SUR LE TRIANGLE de REFERENCE
C     P1DP22D(i,k,j) = integrale P1i DP2j/dxk dx dy sur TRIANGLE UNITE
ccc   DOUBLE PRECISION  P1DP22D(3,2,6)
      include"./incl/p1dp22d.inc"
C
C     INTEGRALES P1 DP2 DP2 SUR LE TRIANGLE de REFERENCE
C     P1DP2DP2(i,k,j,l,m) = integrale P1i DP2j/dxk DP2m/dxl dx dy
ccc   DOUBLE PRECISION  P1DP2DP2(3,2,6,2,6)
      include"./incl/p1dp2dp22d.inc"
C
C     INTEGRALES P1 P2 DDP2 SUR LE TRIANGLE de REFERENCE
C     P1DP2DP2(i,j,k,l,m) = integrale P1i P2j DDP2m/dxkdxl dx dy
ccc   DOUBLE PRECISION  P1P2DDP2(3,6,2,2,6)
      include"./incl/p1p2ddp22d.inc"
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
C
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X21 = XYZNOE(1,NS2) - XYZNOE(1,NS1)
         X31 = XYZNOE(1,NS3) - XYZNOE(1,NS1)
C
         Y21 = XYZNOE(2,NS2) - XYZNOE(2,NS1)
         Y31 = XYZNOE(2,NS3) - XYZNOE(2,NS1)
C
C        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DELTAe = ABS( X21*Y31 - X31*Y21 )
C
C        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
C        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1
         DFM1(1,1) =  Y31 / DELTAe
         DFM1(2,1) = -X31 / DELTAe
C
         DFM1(1,2) = -Y21 / DELTAe
         DFM1(2,2) =  X21 / DELTAe
C
C        CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES TAYLOR-HOOD
C        -----------------------------------------------------------
         IEFORC = LTDESU(LPFORC,NOOBSF)
         IF( IEFORC .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS SURFACIQUES AUX 3 MILIEUX
C           DES ARETES DU TRIANGLE (INTEGRALE=0 AUX SOMMETS)
            DO J=1,6
               NONOJ = NONOEF(NEF,J)
               X = XYZNOE(1,NONOJ)
               Y = XYZNOE(2,NONOJ)
ccc            Z = 0D0
               CALL REFORC( 3,NOOBSF, 2, X,Y,0D0,  0D0,0D0,0D0,
     %                      LTDESU(LPFORC,NOOBSF), FORCE(1,J) )
            ENDDO
C
         ENDIF
C
         DO I=1,3
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            BeI = 0D0
C
            DO L=1,2
C
C              COEFFICIENTS DES EFFORTS SURFACIQUES L
C              Integrale  P1 DivP2 Fel  dx
               FL = 0D0
               IF( IEFORC .GT. 0 ) THEN
                  DO J=1,6
                     DO N=1,2
                        FL = FL + DFM1(L,N) * P1DP22D(I,N,J)
     %                          * FORCE(L,J) * CoForc
                     ENDDO
                  ENDDO
               ENDIF
C
C              Integrale  P1  d/dxl Somme  uk. d/dxk  ul dx
C                                   k=1,2
               SL = 0D0
               DO K=1,2
C
                  SK = 0D0
                  DO J=1,6
C
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NONOEF(NEF,J)
                     VKJ = VXVY( NONOEF(NEF,J), K )
C
C                    LE 1-ER TERME DERIVE
                     SJ = 0D0
                     DO N=1,2
                        DO JJ=1,6
C
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NONOEF(NEF,JJ)
                           VLJJ = VXVY( NONOEF(NEF,JJ), L )
C
                           DO NN=1,2
                              SJ = SJ + DFM1(L,N) * DFM1(K,NN)
     %                                * P1DP2DP2(I,N,J,NN,JJ)
     %                                * VLJJ
                           ENDDO
C
                        ENDDO
                     ENDDO
C
C                    LE SECOND TERME DERIVE
                     DO JJ=1,6
C
C                       COMPOSANTE L DE LA VITESSE AU NOEUD NONOEF(NEF,JJ)
                        VLJJ = VXVY( NONOEF(NEF,JJ), L )
C
                        DO N=1,2
                           DO NN=1,2
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
C           ASSEMBLAGE DE Be(I) DANS BG( NONOSO( NOSOTR(I) ) )
            NONOJ = NONOSO( NOSOTR(I) )
            BG(NONOJ) = BG(NONOJ) + BeI * DELTAe
C
         ENDDO
C
 100  CONTINUE
C
      RETURN
      END
