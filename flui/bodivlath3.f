      SUBROUTINE BODIVLATH3( NBSOM,  NBNOVI, MNXYZN, MNPGEL, NONOSO,
     %                       MNELE,  NBEF,   NOOBSF,
     %                       DELTAT, Rho, Mhu, Omega, WSTAR, WM, WTN,
     %                       BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
C ----     -1/Rho LAPLACIEN P = -(Id/dt - Mhu/Rho Laplacien) Div  WSTAR*
C                       + 2 Div ( Omega(tn+1) x ( Wn+1,m - W(tn) ) )
C                -1/Rho dP/dn = Mhu/Rho n . Laplacien WSTAR*
C          POUR L'ELEMENT FINI TETRAEDRE TAYLOR-HOOD
C          Div WSTAR* EST INTERPOLEE P1 EXACTEMENT SUR LE TETRAEDRE
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
C DELTAT : PAS DE TEMPS
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C Mhu    : Coefficient du Laplacien  Integrale Mhu DP2 DP2 dX
C Omega  : 3 COMPOSANTES DE LA VITESSE ANGULAIRE
C
C WSTAR  : 3 COMPOSANTES DE LA VITESSE W*    EN LES NBNOVI NOEUDS VITESSE
C WM     : 3 COMPOSANTES DE LA VITESSE WM    EN LES NBNOVI NOEUDS VITESSE
C WTN    : 3 COMPOSANTES DE LA VITESSE W(tn) EN LES NBNOVI NOEUDS VITESSE
C          DE LA TETRAEDRISATION
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      INTEGER           MNXYZN, MNPGEL, MNELE, NBSOM, NBEF, NBNOVI,
     %                  NBNOEF, NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOEF(10), NONOSO(NBNOVI)
      INTEGER           I, J, K, K1, K2, L, M, NS1, NS2, NS3, NEF, NS
C
      REAL              XYZEF(10,3)
      DOUBLE PRECISION  DELTAT, Rho, Mhu, Omega(3), BG(NBSOM),
     %                  WSTAR(NBNOVI,3), WM(NBNOVI,3), WTN(NBNOVI,3)
      DOUBLE PRECISION  DIVWP2(4), TDLDL(4,4)
      DOUBLE PRECISION  DELTA, DF(3,3), DFM1(3,3), DFM1DLA(3,4), VN(3),
     %                  DGL(2,3), X, Y, Z, V(4), S, D, DD,
     %                  CRho, CMhu
C
      INTEGER           NOSOFAP1(3,4)
      DATA              NOSOFAP1 / 1,3,2,  1,4,3,  1,2,4,  2,3,4 /
C
C     DOUBLE PRECISION  P1DP23D(4,3,10)
C     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
      include"./incl/p1dp23d.inc"
C
C     [DP2(4 Sommets)]  SUR LE TETRAEDRE de REFERENCE
      DOUBLE PRECISION  DP2S(3,10,4)
      DATA              DP2S /
     % -0.3D+01, -0.3D+01,
     % -0.3D+01, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     %  0.4D+01,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.4D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.4D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.1D+01,  0.1D+01,
     %  0.1D+01,  0.3D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     % -0.4D+01, -0.4D+01,
     % -0.4D+01,  0.0D+00,
     %  0.4D+01,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.4D+01,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.1D+01,  0.1D+01,
     %  0.1D+01, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.3D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.4D+01,
     %  0.0D+00,  0.0D+00,
     % -0.4D+01, -0.4D+01,
     % -0.4D+01,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.4D+01,
     %  0.1D+01,  0.1D+01,
     %  0.1D+01, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.1D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.3D+01,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.0D+00, -0.4D+01,
     % -0.4D+01, -0.4D+01,
     %  0.4D+01,  0.0D+00,
     %  0.0D+00,  0.0D+00,
     %  0.4D+01,  0.0D+00 /
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
C        CALCUL de -Div [P2] {W*(4 Sommets)}
C        -----------------------------------
         DO M = 1, 4
C           AU SOMMET M
            S = 0D0
            DO I = 1, 3
C              dUi/dxi
               DO J = 1, 3
                  DO L = 1, 10
C                    NO GLOBAL DU NOEUD J DU TETRAEDRE NEF
                     NS = NONOEF(L)
                     S = S - DFM1(I,J) * DP2S(J,L,M) * WSTAR(NS,I)
                  ENDDO
               ENDDO
            ENDDO
            DIVWP2( M ) = S
         ENDDO
C
C        tDLambda tDFM1  DFM1 DLambda
C        ----------------------------
         DO I=1,4
            DO J=1,4
               S = 0D0
               DO K=1,3
                  S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
               ENDDO
               TDLDL(I,J) = S
            ENDDO
         ENDDO
C
C        BE(I) COMPOSANTE I DU SECOND MEMBRE
C        -----------------------------------
C        Coefficient de Integrale P1 P1 dX
         CRho = 1D0 / ( DELTAT * 120D0 )
         CMhu = Mhu / Rho
C
         DO I=1,4
C
C           Integrale tP1 Div{ - (Id/dt - Mhu/Rho Laplacien) WSTAR
            V(1) = CRho + CMhu * TDLDL(I,1)
            V(2) = CRho + CMhu * TDLDL(I,2)
            V(3) = CRho + CMhu * TDLDL(I,3)
            V(4) = CRho + CMhu * TDLDL(I,4)
            V(I) = V(I) + CRho
C
C           Integrale tP1 Div{ 2 Omega(tn+1) x (Wn+1,m-W(tn)) } dX
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
                  DD = DELTA* 2D0*( Omega(K1) * (WM(NS,K2)- WTN(NS,K2))
     %                             -Omega(K2) * (WM(NS,K1)- WTN(NS,K1)))
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
            BG(NS) = BG(NS) + D
     %             + ( V(1) * DIVWP2(1) + V(2) * DIVWP2(2)
     %               + V(3) * DIVWP2(3) + V(4) * DIVWP2(4) ) * DELTA
C
         ENDDO
C
C        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
C        =====================================================
         DO K=1,4
C
C           LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE
            NS = NOOBSF(K)
            IF( NS .GT. 0 ) THEN
C
C              TOUTE FACE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
C              LE NUMERO DANS LE TETRAEDRE DES 3 SOMMETS DE LA FACE
C              K DU TETRAEDRE
               NS1 = NOSOFAP1(1,K)
               NS2 = NOSOFAP1(2,K)
               NS3 = NOSOFAP1(3,K)
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
C              COMPOSANTE NS1=NS2=NS3 DU SECOND MEMBRE BE
               D = CMhu / 6D0 *
     %            ( VN(1) * ( DFM1DLA(1,1) * DIVWP2(1)
     %                      + DFM1DLA(1,2) * DIVWP2(2)
     %                      + DFM1DLA(1,3) * DIVWP2(3)
     %                      + DFM1DLA(1,4) * DIVWP2(4) )
     %            + VN(2) * ( DFM1DLA(2,1) * DIVWP2(1)
     %                      + DFM1DLA(2,2) * DIVWP2(2)
     %                      + DFM1DLA(2,3) * DIVWP2(3)
     %                      + DFM1DLA(2,4) * DIVWP2(4) )
     %            + VN(3) * ( DFM1DLA(3,1) * DIVWP2(1)
     %                      + DFM1DLA(3,2) * DIVWP2(2)
     %                      + DFM1DLA(3,3) * DIVWP2(3)
     %                      + DFM1DLA(3,4) * DIVWP2(4) ) )
C
C              ASSEMBLAGE DE Be(NS1) DANS BG( NONOSO( NoSOMMET1 ) )
               NS = NONOSO( NS1 )
               BG(NS) = BG(NS) - D
C
C              BG = BG + Be(NS2)
C              ASSEMBLAGE DE Be(NS2) DANS BG( NONOSO( NoSOMMET2 ) )
               NS = NONOSO( NS2 )
               BG(NS) = BG(NS) - D
C
C              BG = BG + Be(NS3)
C              ASSEMBLAGE DE Be(NS3) DANS BG( NONOSO( NoSOMMET3 ) )
               NS = NONOSO( NS3 )
               BG(NS) = BG(NS) - D
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
