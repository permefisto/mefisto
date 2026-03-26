      SUBROUTINE CORPRETH3( NBSOM,  NBNOVI, MNXYZN, MNPGEL, NONOSO,
     %                      MNELE,  NBEF,   NOOBSF,
     %                      Omega, Wtn1,
     %                      BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
C ----     POUR CORRIGER LA PRESSION APRES LE CALCUL DE W(tn+1)
C          -1/Rho LAPLACIEN P = Div( 2(Omega(tn+1) x W(tn+1))) -2|Omega|**2
C          -1/Rho dP/dn = n .(2(Omega(tn+1) x W(tn+1) + Omega x (Omega x r))
C          POUR L'ELEMENT FINI TETRAEDRE TAYLOR-HOOD
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C MNXYZN : ADRESSE MCN DU TABLEAU DES COORDONNEES DES NOEUDS
C MNPGEL : ADRESSE MCN DU TABLEAU DES COORDONNEES DES POINTS
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C
C MNELE  : ADRESSE MCN DU TABLEAU NPEF DES TETRAEDRES
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NOOBSF : NUMERO DES OBJETS SURFACES DES FACES DE L'ELEMENT FINI
C
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C Mhu    : Coefficient du Laplacien  Integrale Mhu DP2 DP2 dX
C Omega  : 3 COMPOSANTES DE LA VITESSE ANGULAIRE
C
C Wtn1   : 3 COMPOSANTES DE LA VITESSE W(tn+1) EN LES NBNOVI NOEUDS
C          VITESSE DE LA TETRAEDRISATION
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      INTEGER           MNXYZN, MNPGEL, MNELE,
     %                  NBSOM, NBEF, NBNOVI,
     %                  NBNOEF, NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOEF(10), NONOSO(NBNOVI)
      INTEGER           I, J, K, K1, K2, L, NEF, NS, NSJ
      INTEGER           NS1, NS2, NS3, NSF(3)
      EQUIVALENCE      (NS1,NSF(1)), (NS2,NSF(2)), (NS3,NSF(3))
C
      REAL              XYZEF(10,3)
      DOUBLE PRECISION  Omega(3), OmegaS, BG(NBSOM),
     %                  Wtn1(NBNOVI,3)
      DOUBLE PRECISION  DELTA, DF(3,3), DFM1(3,3), DFM1DLA(3,4), VN(3),
     %                  DGL(2,3), X, Y, Z, D, DD,
     %                  X1, X2, X3,  Y1, Y2, Y3,  Z1, Z2, Z3,
     %                  PI1, PI2, PI3
C
      INTEGER           NONOFAP2(6,4)
      DATA              NONOFAP2 / 1, 3, 2, 7,  6, 5,
     %                             1, 4, 3, 8, 10, 7,
     %                             1, 2, 4, 5,  9, 8,
     %                             2, 3, 4, 6, 10, 9 /
C
C     INTEGRALE P1i P1j dx SUR LE TRIANGLE UNITE
      DOUBLE PRECISION  P1P1(3,3)
      DATA              P1P1    /
     %  0.83333333333333333D-01, 0.41666666666666667D-01,
     %  0.41666666666666667D-01,
     %  0.41666666666666667D-01, 0.83333333333333333D-01,
     %  0.41666666666666667D-01,
     %  0.41666666666666667D-01, 0.41666666666666667D-01,
     %  0.83333333333333333D-01 /
C
C     INTEGRALE P1i P2j dx SUR LE TRIANGLE UNITE
      DOUBLE PRECISION  P1P2(3,6)
      DATA              P1P2    /
     %  0.16666666666666667D-01, -0.83333333333333333D-02,
     % -0.83333333333333333D-02,
     % -0.83333333333333333D-02,  0.16666666666666667D-01,
     % -0.83333333333333333D-02,
     % -0.83333333333333333D-02, -0.83333333333333333D-02,
     %  0.16666666666666667D-01,
     %  0.66666666666666667D-01,  0.66666666666666667D-01,
     %  0.33333333333333333D-01,
     %  0.33333333333333333D-01,  0.66666666666666667D-01,
     %  0.66666666666666667D-01,
     %  0.66666666666666667D-01,  0.33333333333333333D-01,
     %  0.66666666666666667D-01 /
C
C     DOUBLE PRECISION  P1DP23D(4,3,10)
C     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
      include"./incl/p1dp23d.inc"
C
C     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
C     CONTRIBUTION DE LA VOLUME AU SECOND MEMBRE
C     ===========================================
      DO NEF = 1, NBEF
C
C        NO DES NOEUDS DE L'ELEMENT FINI NEF
         CALL EFNOEU( MNELE, NEF, NBNOEF, NONOEF )
C
C        NO DE  POINTS  LIGNES SURFACES VOLUME
C           DES SOMMETS FACES  FACES    VOLUME DU TETRAEDRE NEF
         CALL EFPLSV( MNELE , NEF,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C        COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
         CALL EFXYZP( 3, MNXYZN, NBEF, NEF, MNPGEL, NBNOEF,  XYZEF )
C
C        CONSTRUCTION DE LA MATRICE DF
         X = XYZEF(1,1)
         DF(1,1) = XYZEF(2,1) - X
         DF(2,1) = XYZEF(3,1) - X
         DF(3,1) = XYZEF(4,1) - X
C
         Y = XYZEF(1,2)
         DF(1,2) = XYZEF(2,2) - Y
         DF(2,2) = XYZEF(3,2) - Y
         DF(3,2) = XYZEF(4,2) - Y
C
         Z = XYZEF(1,3)
         DF(1,3) = XYZEF(2,3) - Z
         DF(2,3) = XYZEF(3,3) - Z
         DF(3,3) = XYZEF(4,3) - Z
C
C        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3))
     %              +DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3))
     %              +DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1
         DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) / DELTA
         DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) ) / DELTA
         DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) ) / DELTA
C
         DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) ) / DELTA
         DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) ) / DELTA
         DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) ) / DELTA
C
         DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) / DELTA
         DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) ) / DELTA
         DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) ) / DELTA
C
C        [DFM1] [DLa]
         DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
         DFM1DLa(1,2) =  DFM1(1,1)
         DFM1DLa(1,3) =  DFM1(1,2)
         DFM1DLa(1,4) =  DFM1(1,3)
C
         DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
         DFM1DLa(2,2) =  DFM1(2,1)
         DFM1DLa(2,3) =  DFM1(2,2)
         DFM1DLa(2,4) =  DFM1(2,3)
C
         DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
         DFM1DLa(3,2) =  DFM1(3,1)
         DFM1DLa(3,3) =  DFM1(3,2)
         DFM1DLa(3,4) =  DFM1(3,3)
C
C        BE(I) COMPOSANTE I DU SECOND MEMBRE
C        -----------------------------------
         OmegaS = Omega(1)**2 + Omega(2)**2 + Omega(3)**2
C
         DELTA = DELTA * 2D0
         DO I=1,4
C
C           Integrale tP1 Div{ 2 Omega(tn+1) x W(tn+1) } dX
            D = 0D0
            DO K=1,3
               IF( K .NE. 3 ) THEN
                  K1 = K + 1
               ELSE
                  K1 = 1
               ENDIF
C
               IF( K1 .NE. 3 ) THEN
                  K2 = K1 + 1
               ELSE
                  K2 = 1
               ENDIF
C
               DO J=1,10
C
C                 NO GLOBAL DU NOEUD J DU TETRAEDRE TH
                  NS = NONOEF(J)
                  DD = Omega(K1) * Wtn1(NS,K2) - Omega(K2) * Wtn1(NS,K1)
C
                  DO L=1,3
C                    Integrale tP1 Div{ 2 Omega(tn+1) x (Wn+1,m-W(tn)) } dX
C                    P1DP23D(i,l,j) = integrale Lambdai dPj/dxl dx dy dz
                     D = D + DD * DFM1(K,L) * P1DP23D(I,L,J)
                  ENDDO
C
               ENDDO
C
            ENDDO
C
C           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NoSOMMET i ) )
            NS = NONOSO( NONOEF(I) )
            BG(NS) = BG(NS) + DELTA * ( D - OmegaS / 24D0 )
C
         ENDDO
C
C        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
C        -1/Rho dP/dn = n .(2(Omega(tn+1) x W(tn+1) + Omega x (Omega x r))
C        =================================================================
         DO K=1,4
C
C           LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE K
            NS = NOOBSF(K)
            IF( NS .GT. 0 ) THEN
C
C              TOUTE FACE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
C              LE NUMERO DANS LE TETRAEDRE DES 3 SOMMETS DE LA FACE
C              K DU TETRAEDRE
               NS1 = NONOFAP2(1,K)
               NS2 = NONOFAP2(2,K)
               NS3 = NONOFAP2(3,K)
C
C              CALCUL DU JACOBIEN DE G
               X = XYZEF( NS1, 1 )
               DGL(1,1) = XYZEF( NS2, 1 ) - X
               DGL(2,1) = XYZEF( NS3, 1 ) - X
C
               Y = XYZEF( NS1, 2 )
               DGL(1,2) = XYZEF( NS2, 2 ) - Y
               DGL(2,2) = XYZEF( NS3, 2 ) - Y
C
               Z = XYZEF( NS1, 3 )
               DGL(1,3) = XYZEF( NS2, 3 ) - Z
               DGL(2,3) = XYZEF( NS3, 3 ) - Z
C
C              CALCUL DE LA NORMALE: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
               VN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
               VN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
               VN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)
C
cccC           LE JACOBIEN DE G
ccc            DELTA = SQRT( VN(1) ** 2 + VN(2) ** 2 + VN(3) ** 2 )
C
C              LE VECTEUR NORMAL UNITAIRE
ccc            ( NON /DELTA COMPENSEE PAR LE JACOBIEN )
ccc            VN(1) = VN(1) / DELTA
ccc            VN(2) = VN(2) / DELTA
ccc            VN(3) = VN(3) / DELTA
C
C              LES 3 COORDONNEES DES 3 SOMMETS DE LA FACE K
               X1 = XYZEF(NS1,1)
               Y1 = XYZEF(NS1,2)
               Z1 = XYZEF(NS1,3)
C
               X2 = XYZEF(NS2,1)
               Y2 = XYZEF(NS2,2)
               Z2 = XYZEF(NS2,3)
C
               X3 = XYZEF(NS3,1)
               Y3 = XYZEF(NS3,2)
               Z3 = XYZEF(NS3,3)
C
C              COMPOSANTE NSi DU SECOND MEMBRE BE
               DO I = 1, 3
C
                  PI1 = P1P1(I,1)
                  PI2 = P1P1(I,2)
                  PI3 = P1P1(I,3)
C
                  D = 0D0
                  DO J=1,6
C                   NO GLOBAL DU NOEUD J DE LA FACE K DU TETRAEDRE
                    NSJ = NONOFAP2( J, K )
                    NS  = NONOEF( NSJ )
                    D   = D + P1P2(I,J) * 2D0 *
     %                 ( VN(1)*(Omega(2)*Wtn1(NS,3)-Omega(3)*Wtn1(NS,2))
     %                 + VN(2)*(Omega(3)*Wtn1(NS,1)-Omega(1)*Wtn1(NS,3))
     %                 + VN(3)*(Omega(1)*Wtn1(NS,2)-Omega(2)*Wtn1(NS,1))
     %                 )
                  ENDDO
C
                  D = D
     %              + VN(1) * ( Omega(2) * ( Omega(1) * ( PI1 * Y1
     %                                                  + PI2 * Y2
     %                                                  + PI3 * Y3 )
     %                                     - Omega(2) * ( PI1 * X1
     %                                                  + PI2 * X2
     %                                                  + PI3 * X3 ) )
     %                        - Omega(3) * ( Omega(3) * ( PI1 * X1
     %                                                  + PI2 * X2
     %                                                  + PI3 * X3 )
     %                                     - Omega(1) * ( PI1 * Z1
     %                                                  + PI2 * Z2
     %                                                  + PI3 * Z3 ) ) )
     %              + VN(2) * ( Omega(3) * ( Omega(2) * ( PI1 * Z1
     %                                                  + PI2 * Z2
     %                                                  + PI3 * Z3 )
     %                                     - Omega(3) * ( PI1 * Y1
     %                                                  + PI2 * Y2
     %                                                  + PI3 * Y3 ) )
     %                        - Omega(1) * ( Omega(1) * ( PI1 * Y1
     %                                                  + PI2 * Y2
     %                                                  + PI3 * Y3 )
     %                                     - Omega(2) * ( PI1 * X1
     %                                                  + PI2 * X2
     %                                                  + PI3 * X3 ) ) )
     %              + VN(3) * ( Omega(1) * ( Omega(3) * ( PI1 * X1
     %                                                  + PI2 * X2
     %                                                  + PI3 * X3 )
     %                                     - Omega(1) * ( PI1 * Z1
     %                                                  + PI2 * Z2
     %                                                  + PI3 * Z3 ) )
     %                        - Omega(2) * ( Omega(2) * ( PI1 * Z1
     %                                                  + PI2 * Z2
     %                                                  + PI3 * Z3 )
     %                                     - Omega(3) * ( PI1 * Y1
     %                                                  + PI2 * Y2
     %                                                  + PI3 * Y3 ) ) )
C
C              BG = BG + Be(NSi)
C              ASSEMBLAGE DE Be(NS3) DANS BG( NONOSO( NoSOMMETi ) )
               NS = NONOSO( NONOEF( NSF(I) ) )
               BG(NS) = BG(NS) + D
C
               ENDDO
C
            ENDIF
C
C           FIN TRAITEMENT FACE K
         ENDDO
C
C        FIN TRAITEMENT DE L'EF NEF
      ENDDO
C
      RETURN
      END
