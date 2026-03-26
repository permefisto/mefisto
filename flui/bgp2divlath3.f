      SUBROUTINE BGP2DIVLATH3( NBNOVI, MNXYZN, MNPGEL, MNELE, NBEF,
     %                         VITXYZ, Rho,    dtMhu,   BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
C ----   - dt CoGrPr LAPLACIEN P = -(Rho - dtMhu Laplacien) Div VITXYZ*
C        - dt CoGrPr dP/dn       =     n . dtMhu Laplacien VITXYZ*
C          POUR L'ELEMENT FINI TETRAEDRE TAYLOR-HOOD
C          Div VITXYZ* EST INTERPOLEE P1 EXACTEMENT
C          LES FONCTIONS TESTS SONT P2 et la PRESSION P2
C
C ENTREES:
C --------
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C MNXYZN : ADRESSE MCN DU TABLEAU DES COORDONNEES DES NOEUDS
C MNPGEL : ADRESSE MCN DU TABLEAU DES COORDONNEES DES POINTS
C MNELE  : ADRESSE MCN DU TABLEAU NPEF DES TETRAEDRES
C NBEF   : NOMBRE D'EF DU MAILLAGE
C
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C dtMhu  : Coefficient du Laplacien  Integrale dtMhu DP2 DP2 dX
C VITXYZ : 3 COMPOSANTES X Y Z DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
C          DE LA TETRAEDRISATION
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC Paris & St Pierre du Perray juin 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      INTEGER           MNXYZN, MNPGEL, MNELE, NBEF, NBNOVI,
     %                  NBNOEF, NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOEF(10), I, J, K, L, M, NS1, NS2, NS3, NEF,NS
C
      REAL              XYZEF(10,3)
      DOUBLE PRECISION  Rho, dtMhu, VITXYZ(NBNOVI,3), BG(NBNOVI)
      DOUBLE PRECISION  DIVUP2(4), DELTA, DF(3,3), DFM1(3,3),
     %                  DFM1DLA(3,4), VN(3), DGL(2,3), S, D, dtMhu6
C
      DOUBLE PRECISION  IntP2La(10,4)
      DATA              IntP2La /
     %  0.00000000000000000D+00, -0.27777777777777778D-02,
     % -0.27777777777777778D-02, -0.27777777777777778D-02,
     %  0.11111111111111111D-01,  0.55555555555555556D-02,
     %  0.11111111111111111D-01,  0.11111111111111111D-01,
     %  0.55555555555555556D-02,  0.55555555555555556D-02,
     % -0.27777777777777778D-02,  0.00000000000000000D+00,
     % -0.27777777777777778D-02, -0.27777777777777778D-02,
     %  0.11111111111111111D-01,  0.11111111111111111D-01,
     %  0.55555555555555556D-02,  0.55555555555555556D-02,
     %  0.11111111111111111D-01,  0.55555555555555556D-02,
     % -0.27777777777777778D-02, -0.27777777777777778D-02,
     %  0.00000000000000000D+00, -0.27777777777777778D-02,
     %  0.55555555555555556D-02,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.55555555555555556D-02,
     %  0.55555555555555556D-02,  0.11111111111111111D-01,
     % -0.27777777777777778D-02, -0.27777777777777778D-02,
     % -0.27777777777777778D-02,  0.00000000000000000D+00,
     %  0.55555555555555556D-02,  0.55555555555555556D-02,
     %  0.55555555555555556D-02,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.11111111111111111D-01 /
C
      DOUBLE PRECISION IntDP2(3,10)
      DATA             IntDP2 /
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.00000000000000000D+00, -0.16666666666666667D+00,
     % -0.16666666666666667D+00,  0.16666666666666667D+00,
     %  0.16666666666666667D+00,  0.00000000000000000D+00,
     % -0.16666666666666667D+00,  0.00000000000000000D+00,
     % -0.16666666666666667D+00, -0.16666666666666667D+00,
     % -0.16666666666666667D+00,  0.00000000000000000D+00,
     %  0.16666666666666667D+00,  0.00000000000000000D+00,
     %  0.16666666666666667D+00,  0.00000000000000000D+00,
     %  0.16666666666666667D+00,  0.16666666666666667D+00 /
C
      INTEGER           NONOFAP2(6,4)
      DATA              NONOFAP2 / 1, 3, 2,  7,  6, 5,
     %                             1, 4, 3,  8, 10, 7,
     %                             1, 2, 4,  5,  9, 8,
     %                             2, 3, 4,  6, 10, 9  /
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
      DO I = 1, NBNOVI
         BG( I ) = 0D0
      ENDDO
C
C     CONTRIBUTION DU VOLUME AU SECOND MEMBRE
C     =======================================
      dtMhu6 = dtMhu / 6D0
C
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
         S = XYZEF(1,1)
         DF(1,1) = XYZEF(2,1) - S
         DF(2,1) = XYZEF(3,1) - S
         DF(3,1) = XYZEF(4,1) - S
C
         S = XYZEF(1,2)
         DF(1,2) = XYZEF(2,2) - S
         DF(2,2) = XYZEF(3,2) - S
         DF(3,2) = XYZEF(4,2) - S
C
         S = XYZEF(1,3)
         DF(1,3) = XYZEF(2,3) - S
         DF(2,3) = XYZEF(3,3) - S
         DF(3,3) = XYZEF(4,3) - S
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
C        CALCUL de -Div [P2] {U(4 Sommets)}
C        ----------------------------------
         DO M = 1, 4
C           AU SOMMET M
            S = 0D0
            DO I = 1, 3
C              dUi/dxi
               DO J = 1, 3
                  DO L = 1, 10
C                    NO GLOBAL DU NOEUD L DU TETRAEDRE NEF
                     NS = NONOEF(L)
                     S = S - DFM1(I,J) * DP2S(J,L,M) * VITXYZ(NS,I)
                  ENDDO
               ENDDO
            ENDDO
            DIVUP2( M ) = S
         ENDDO
C
C        (int P2 Lambda dx + int tDP2 dx tDFM1  DFM1 DLambda) (-div U(se) )
C        ------------------------------------------------------------------
         DO I=1,10
C
C           int P2 Lambda dx (-div U(se) )
            D = 0D0
            DO M=1,4
               D = D + IntP2La(I,M) * DIVUP2(M)
            ENDDO
C
C           int tDP2 dx tDFM1  DFM1 DLambda (-div U(se) )
            S = 0D0
            DO K=1,3
               DO L=1,3
                  DO M=1,4
                     S = S + IntDP2(K,I) * DFM1(K,L) * DFM1DLa(L,M)
     %                                   * DIVUP2(M)
                  ENDDO
               ENDDO
            ENDDO
C           ASSEMBLAGE DE BE(I) DANS BG( NONOEF(I) )
            NS = NONOEF(I)
            BG(NS) = BG(NS) + ( Rho * D + dtMhu * S ) * DELTA
C
         ENDDO
C
C        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
C        - dt CoGrPr dP/dn = n . dtMhu Laplacien VITXYZ*
C        int dt CoGrPr tP2 dP2/dn dG = -int tP2 dLambda/dn dG {-div VITXYZ*}
C        ===================================================================
         DO K=1,4
C
C           LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE
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
               S = XYZEF( NS1, 1 )
               DGL(1,1) = XYZEF( NS2, 1 ) - S
               DGL(2,1) = XYZEF( NS3, 1 ) - S
C
               S = XYZEF( NS1, 2 )
               DGL(1,2) = XYZEF( NS2, 2 ) - S
               DGL(2,2) = XYZEF( NS3, 2 ) - S
C
               S = XYZEF( NS1, 3 )
               DGL(1,3) = XYZEF( NS2, 3 ) - S
               DGL(2,3) = XYZEF( NS3, 3 ) - S
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
ccc            ( NON /DELTA car COMPENSEE PAR LE JACOBIEN )
ccc            VN(1) = VN(1) / DELTA
ccc            VN(2) = VN(2) / DELTA
ccc            VN(3) = VN(3) / DELTA
C
C              COMPOSANTE NS1=NS2=NS3 DU SECOND MEMBRE BE
               D = dtMhu6 *
     %            ( VN(1) * ( DFM1DLA(1,1) * DIVUP2(1)
     %                      + DFM1DLA(1,2) * DIVUP2(2)
     %                      + DFM1DLA(1,3) * DIVUP2(3)
     %                      + DFM1DLA(1,4) * DIVUP2(4) )
     %            + VN(2) * ( DFM1DLA(2,1) * DIVUP2(1)
     %                      + DFM1DLA(2,2) * DIVUP2(2)
     %                      + DFM1DLA(2,3) * DIVUP2(3)
     %                      + DFM1DLA(2,4) * DIVUP2(4) )
     %            + VN(3) * ( DFM1DLA(3,1) * DIVUP2(1)
     %                      + DFM1DLA(3,2) * DIVUP2(2)
     %                      + DFM1DLA(3,3) * DIVUP2(3)
     %                      + DFM1DLA(3,4) * DIVUP2(4) ) )
C
C              L'integrale de P2 sur le triangle = 0   pour i=1,2,3
C                                                = 1/6 pour i=4,5,6
C              SEULS LES MILIEUX DES 3 ARETES DE LA FACE K ONT UNE CONTRIBUTION
               DO M=4,6
                  NS = NONOEF( NONOFAP2(M,K) )
                  BG(NS) = BG(NS) - D
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
