      SUBROUTINE BGDIVLATH3( Rho,    dtMhu,  VITXYZ, NBNOVI, XYZNOE,
     %                       NONOSO, NBELEM, NUNOEF, NUSFEF,
     %                       NBSOM,  BG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
! ----     -dt CoGrPr LAPLACIEN P = -(Rho - dtMhu Laplacien) Div  VITXYZ
!          -dt CoGrPr       dP/dn = dtMhu n . Laplacien VITXYZ
!          POUR L'ELEMENT FINI TETRAEDRE TAYLOR-HOOD
!          Div VITXYZ EST INTERPOLEE P1 EXACTEMENT

! ENTREES:
! --------
! Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
! dtMhu  : Coefficient du Laplacien  Integrale dtMhu DP2 DP2 dX
! VITXYZ : 3 COMPOSANTES XYZ DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
!          DE LA TETRAEDRISATION
! NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
! XYZNOE : (3,NBNOVI) TABLEAU DES 3 COORDONNEES DES NBNOVI NOEUDS
! NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,10) NO DES 10 NOEUDS DES NBELEM EF
! NUSFEF : TABLEAU DU NUMERO DE SURFACE DES 4 FACES DE CHAQUE TETRAEDRE
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE ou NOMBRE DE COMPOSANTES DE BG

! SORTIE :
! --------
! BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE

! REMARQUE: SI UN EF EST DE SURFACE NULLE, SA CONTRIBUTION A BG EST NULLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      INTEGER           NBSOM, NBELEM, NBNOVI
      INTEGER           NUNOEF(NBELEM,10), NUSFEF(4,NBELEM),
     %                  NONOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  Rho, dtMhu, VITXYZ(NBNOVI,3), BG(NBSOM)

      INTEGER           I, J, K, L, M, NS1, NS2, NS3, NUELEM, NS
      DOUBLE PRECISION  DIVUP2(4), TDLDL(4,4), DELTA,
     %                  DF(3,3), DFM1(3,3), DFM1DLA(3,4), VN(3),
     %                  DGL(2,3), X, Y, Z, V1, V2, V3, V4, S,
     %                  Rho120, dtMhu6

      INTEGER           NOSOFAP1(3,4)
      DATA              NOSOFAP1 / 1,3,2,  1,4,3,  1,2,4,  2,3,4 /

!     [DP2(4 Sommets)]  SUR LE TETRAEDRE de REFERENCE
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

!     MISE A ZERO DU VECTEUR GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO

!     CONTRIBUTION DU VOLUME AU SECOND MEMBRE
!     =======================================
      Rho120 = Rho / 120D0
      dtMhu6 = dtMhu / 6D0

      DO 100 NUELEM = 1, NBELEM

!        CONSTRUCTION DE LA MATRICE DF
         NS = NUNOEF(NUELEM,1)
         I  = NUNOEF(NUELEM,2)
         J  = NUNOEF(NUELEM,3)
         K  = NUNOEF(NUELEM,4)

         X = XYZNOE(1,NS)
         DF(1,1) = XYZNOE(1,I) - X
         DF(2,1) = XYZNOE(1,J) - X
         DF(3,1) = XYZNOE(1,K) - X

         Y = XYZNOE(2,NS)
         DF(1,2) = XYZNOE(2,I) - Y
         DF(2,2) = XYZNOE(2,J) - Y
         DF(3,2) = XYZNOE(2,K) - Y

         Z = XYZNOE(3,NS)
         DF(1,3) = XYZNOE(3,I) - Z
         DF(2,3) = XYZNOE(3,J) - Z
         DF(3,3) = XYZNOE(3,K) - Z

!        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3))
     %              +DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3))
     %              +DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )
         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT *, 'BGDIVLATH3: ATTENTION EF',NUELEM,
     %                  ' de VOLUME*6=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT *, 'BGDIVLATH3: ATTENTION FE',NUELEM,
     %                  ' of VOLUME*6=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

!        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
!        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1
         DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) / DELTA
         DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) ) / DELTA
         DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) ) / DELTA

         DFM1(1,2) = (  DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) ) / DELTA
         DFM1(2,2) = (  DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) ) / DELTA
         DFM1(3,2) = (  DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) ) / DELTA

         DFM1(1,3) = (  DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) / DELTA
         DFM1(2,3) = (  DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) ) / DELTA
         DFM1(3,3) = (  DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) ) / DELTA

!        [DFM1] [DLa]
         DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
         DFM1DLa(1,2) =  DFM1(1,1)
         DFM1DLa(1,3) =  DFM1(1,2)
         DFM1DLa(1,4) =  DFM1(1,3)

         DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
         DFM1DLa(2,2) =  DFM1(2,1)
         DFM1DLa(2,3) =  DFM1(2,2)
         DFM1DLa(2,4) =  DFM1(2,3)

         DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
         DFM1DLa(3,2) =  DFM1(3,1)
         DFM1DLa(3,3) =  DFM1(3,2)
         DFM1DLa(3,4) =  DFM1(3,3)

!        CALCUL de -Div [P2] {U(4 Sommets)}
!        ----------------------------------
         DO M = 1, 4
!           AU SOMMET M
            S = 0D0
            DO I = 1, 3
!              dUi/dxi
               DO J = 1, 3
                  DO L = 1, 10
!                    NO GLOBAL DU NOEUD L DU TETRAEDRE NUELEM
                     NS = NUNOEF(NUELEM,L)
                     S = S - DFM1(I,J) * DP2S(J,L,M) * VITXYZ(NS,I)
                  ENDDO
               ENDDO
            ENDDO
            DIVUP2( M ) = S
         ENDDO

!        tDLambda tDFM1  DFM1 DLambda
!        ----------------------------
         DO I=1,4
            DO J=1,4
               S = 0D0
               DO K=1,3
                  S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
               ENDDO
               TDLDL(I,J) = S
            ENDDO
         ENDDO

!        PREMIERE COMPOSANTE DU SECOND MEMBRE BE(1)
         V1 = Rho120 * 2D0 + dtMhu6 * TDLDL(1,1)
         V2 = Rho120       + dtMhu6 * TDLDL(1,2)
         V3 = Rho120       + dtMhu6 * TDLDL(1,3)
         V4 = Rho120       + dtMhu6 * TDLDL(1,4)
!        ASSEMBLAGE DE BE(1) DANS BG( NONOSO( NoSOMMET1 ) )
         NS = NONOSO( NUNOEF(NUELEM,1) )
         BG(NS) = BG(NS)
     %          + ( V1*DIVUP2(1) + V2*DIVUP2(2)
     %            + V3*DIVUP2(3) + V4*DIVUP2(4) ) * DELTA

!        SECONDE COMPOSANTE DU SECOND MEMBRE BE(2)
         V1 = Rho120       + dtMhu6 * TDLDL(2,1)
         V2 = Rho120 * 2D0 + dtMhu6 * TDLDL(2,2)
         V3 = Rho120       + dtMhu6 * TDLDL(2,3)
         V4 = Rho120       + dtMhu6 * TDLDL(2,4)
!        ASSEMBLAGE DE BE(2) DANS BG( NONOSO( NoSOMMET2 ) )
         NS = NONOSO( NUNOEF(NUELEM,2) )
         BG(NS) = BG(NS)
     %          + ( V1*DIVUP2(1) + V2*DIVUP2(2)
     %            + V3*DIVUP2(3) + V4*DIVUP2(4) ) * DELTA

!        TROISIEME COMPOSANTE DU SECOND MEMBRE BE(3)
         V1 = Rho120       + dtMhu6 * TDLDL(3,1)
         V2 = Rho120       + dtMhu6 * TDLDL(3,2)
         V3 = Rho120 * 2D0 + dtMhu6 * TDLDL(3,3)
         V4 = Rho120       + dtMhu6 * TDLDL(3,4)
!        ASSEMBLAGE DE BE(3) DANS BG( NONOSO( NoSOMMET3 ) )
         NS = NONOSO( NUNOEF(NUELEM,3) )
         BG(NS) = BG(NS)
     %          + ( V1*DIVUP2(1) + V2*DIVUP2(2)
     %            + V3*DIVUP2(3) + V4*DIVUP2(4) ) * DELTA

!        QUATRIEME COMPOSANTE DU SECOND MEMBRE BE(4)
         V1 = Rho120       + dtMhu6 * TDLDL(4,1)
         V2 = Rho120       + dtMhu6 * TDLDL(4,2)
         V3 = Rho120       + dtMhu6 * TDLDL(4,3)
         V4 = Rho120 * 2D0 + dtMhu6 * TDLDL(4,4)
!        ASSEMBLAGE DE BE(4) DANS BG( NONOSO( NoSOMMET4 ) )
         NS = NONOSO( NUNOEF(NUELEM,4) )
         BG(NS) = BG(NS)
     %          + ( V1*DIVUP2(1) + V2*DIVUP2(2)
     %            + V3*DIVUP2(3) + V4*DIVUP2(4) ) * DELTA

!        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
!        - dt CoGrPr dP/dn = n . dtMhu Laplacien VITXYZ
!        int dt CoGrPr tP1 dP1/dn dG = -int tP1 dLambda/dn dG {-div VITXYZ}
!        ==================================================================
         DO K=1,4

!           LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE K DE NUELEM
            NS = NUSFEF( K, NUELEM )
            IF( NS .GT. 0 ) THEN

!              TOUTE FACE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
!              LE NUMERO DANS LE TETRAEDRE DES 3 SOMMETS DE LA FACE
!              K DU TETRAEDRE
               NS1 = NUNOEF( NUELEM, NOSOFAP1(1,K) )
               NS2 = NUNOEF( NUELEM, NOSOFAP1(2,K) )
               NS3 = NUNOEF( NUELEM, NOSOFAP1(3,K) )

!              CALCUL DU JACOBIEN DE G
               X = XYZNOE( 1, NS1 )
               DGL(1,1) = XYZNOE( 1, NS2 ) - X
               DGL(2,1) = XYZNOE( 1, NS3 ) - X
!
               Y = XYZNOE( 2, NS1 )
               DGL(1,2) = XYZNOE( 2, NS2 ) - Y
               DGL(2,2) = XYZNOE( 2, NS3 ) - Y
!
               Z = XYZNOE( 3, NS1 )
               DGL(1,3) = XYZNOE( 3, NS2 ) - Z
               DGL(2,3) = XYZNOE( 3, NS3 ) - Z

!              CALCUL DE LA NORMALE: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
               VN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
               VN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
               VN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)

!!!!           LE JACOBIEN DE G
!!!            DELTA = SQRT( VN(1) ** 2 + VN(2) ** 2 + VN(3) ** 2 )

!!!            LE VECTEUR NORMAL UNITAIRE
!!!            ( NON /DELTA COMPENSEE PAR LE JACOBIEN )
!!!            VN(1) = VN(1) / DELTA
!!!            VN(2) = VN(2) / DELTA
!!!            VN(3) = VN(3) / DELTA
!
!              COMPOSANTE NS1=NS2=NS3 DU SECOND MEMBRE BE
               S = dtMhu6 *
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

!              ASSEMBLAGE DE Be(NS1) DANS BG( NONOSO( NoSOMMET1 ) )
               NS = NONOSO( NS1 )
               BG(NS) = BG(NS) - S

!              BG = BG + Be(NS2)
!              ASSEMBLAGE DE Be(NS2) DANS BG( NONOSO( NoSOMMET2 ) )
               NS = NONOSO( NS2 )
               BG(NS) = BG(NS) - S

!              BG = BG + Be(NS3)
!              ASSEMBLAGE DE Be(NS3) DANS BG( NONOSO( NoSOMMET3 ) )
               NS = NONOSO( NS3 )
               BG(NS) = BG(NS) - S

            ENDIF

!           FIN TRAITEMENT FACE K
         ENDDO

!        FIN TRAITEMENT DE L'EF NUELEM
 100     CONTINUE

      RETURN
      END
