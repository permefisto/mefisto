      SUBROUTINE F3BPGTH( UnSRho,
     %                    NBSOM,  NBNOVI, MNXYZN, MNPGEL, NONOSO,
     %                    MNELE,  NBEF,
     %                    NOOBSF,
     %                    NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                    VITm,   P2DP2,  BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:  CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU SYSTEME (2)
C ----
C (2) - 1/Rho Laplacien Pm+1(X) = Div { -Fg(tn+1)/Rho
C                                     + Vm . Grad Vm } dans Omega
C     Condition aux limites
C            -1/Rho    dPm+1/dn = n . ( -Fg(tn+1)/Rho
C                                     +  Vm . Grad Vm } sur Gamma
C     POUR LE TETRAEDRE TAYLOR-HOOD
C     PISO Version 1: Bg avec Div reportee sur la Pression en VOLUME
C
C ENTREES:
C --------
C UnSRho : 1/DENSITE VOLUMIQUE DE MASSE DU FLUIDE
C
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C MNXYZN : ADRESSE MCN DU TABLEAU DES COORDONNEES DES NOEUDS
C MNPGEL : ADRESSE MCN DU TABLEAU DES COORDONNEES DES POINTS
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C MNELE  : ADRESSE MCN DU TABLEAU NPEF DES TETRAEDRES
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NOOBSF : NUMERO DES OBJETS SURFACES DES FACES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C NOOBVC : NUMERO DE VOLUME DU FLUIDE
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES VOLUMES
C
C VITm   : COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en X puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Y puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Z
C P2DP2  : P2DP2(i,k,j) = Integrale P2i DP2j/Dxk dX SUR LE TETRAEDRE UNITE
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C          AUX NBSOM SOMMETS DE LA TETRAEDRISATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      DOUBLE PRECISION  UnSRho,
     %                  VITm(NBNOVI,3), BG(NBSOM)
      INTEGER           NBSOM, NBNOVI, NBNOEF, NBEF, MNXYZN, MNPGEL,
     %                  NONOSO(NBNOVI), MNELE,
     %                  NUMIVO, NUMAVO
C
      INTEGER           NONOTE(10),
     %                  LTDEVO(1:MXDOFL,NUMIVO:NUMAVO),
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8),
     %                  NOVCEL, NOSFEL, NOLAEL, NOPSEL
C
      REAL              XYZEF(10,3)
      DOUBLE PRECISION  FORCE(3,10), P2DP2(10,3,10), VJ, VJJ, Coef,
     %                  SF, SP,
     %                  DELTAe, DF(3,3), DFM1(3,3), DFM1DLa(3,4),
     %                  D, S, X, Y, Z, VN(3), DGL(2,3)
C
      INTEGER           NF, I, J, JJ, K, L, M, N, NONOJ, NONOJJ, NS,
     %                  NEF, IEFORC, NS1, NS2, NS3, NS4, NS5, NS6,NSF(6)
      EQUIVALENCE       (NSF(1),NS1),(NSF(2),NS2),(NSF(3),NS3),
     %                  (NSF(4),NS4),(NSF(5),NS5),(NSF(6),NS6)
C
      INTEGER           NONOFAP2(6,4)
      DATA              NONOFAP2 / 1, 3, 2, 7,  6, 5,
     %                             1, 4, 3, 8, 10, 7,
     %                             1, 2, 4, 5,  9, 8,
     %                             2, 3, 4, 6, 10, 9 /
C
C     [DP2(4 Sommets)]  SUR LE TETRAEDRE de REFERENCE
     include"./incl/dp2s3d.inc"
ccc   DOUBLE PRECISION  DP2S3D(3,10,4)
C
C     INTEGRALE P2i dx SUR LE TETRAEDRE UNITE
      DOUBLE PRECISION  INTP23D(10)
      DATA              INTP23D/
     %  -0.8333333333333333333D-2,  -0.8333333333333333333D-2,
     %  -0.8333333333333333333D-2,  -0.8333333333333333333D-2,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1 /
C
C     INTEGRALE P1i P2j dx SUR LE TRIANGLE UNITE
      DOUBLE PRECISION  P1P22D(3,6)
      DATA              P1P22D    /
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
C     INTEGRALE P1i P2j P1k dx dy SUR LE TRIANGLE UNITE
      DOUBLE PRECISION  P1P2P12D(3,6,3)
      DATA              P1P2P12D  /
     %  0.16666666666666667D-01,  0.00000000000000000D+00,
     %  0.00000000000000000D+00, -0.55555555555555556D-02,
     %  0.00000000000000000D+00, -0.27777777777777778D-02,
     % -0.55555555555555556D-02, -0.27777777777777778D-02,
     %  0.00000000000000000D+00,  0.33333333333333333D-01,
     %  0.22222222222222222D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.33333333333333333D-01,
     %  0.11111111111111111D-01,  0.22222222222222222D-01,
     %  0.00000000000000000D+00, -0.55555555555555556D-02,
     % -0.27777777777777778D-02,  0.00000000000000000D+00,
     %  0.16666666666666667D-01,  0.00000000000000000D+00,
     % -0.27777777777777778D-02, -0.55555555555555556D-02,
     %  0.00000000000000000D+00,  0.22222222222222222D-01,
     %  0.33333333333333333D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.33333333333333333D-01,
     %  0.22222222222222222D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.11111111111111111D-01,
     %  0.00000000000000000D+00, -0.27777777777777778D-02,
     % -0.55555555555555556D-02, -0.27777777777777778D-02,
     %  0.00000000000000000D+00, -0.55555555555555556D-02,
     %  0.00000000000000000D+00,  0.00000000000000000D+00,
     %  0.16666666666666667D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.11111111111111111D-01,
     %  0.11111111111111111D-01,  0.22222222222222222D-01,
     %  0.33333333333333333D-01,  0.22222222222222222D-01,
     %  0.11111111111111111D-01,  0.33333333333333333D-01 /
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NO DES NOEUDS DE L'ELEMENT FINI NEF
         CALL EFNOEU( MNELE, NEF, NBNOEF, NONOTE )
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
C        LE DETERMINANT DE DF
         DELTAe = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %          + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %          + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
         DFM1(1,1) = DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)
         DFM1(2,1) = DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1)
         DFM1(3,1) = DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2)
C
         DFM1(1,2) = DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3)
         DFM1(2,2) = DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1)
         DFM1(3,2) = DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2)
C
         DFM1(1,3) = DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)
         DFM1(2,3) = DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1)
         DFM1(3,3) = DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2)
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
C        CONTRIBUTION DES EFFORTS VOLUMIQUES INTERPOLES TAYLOR-HOOD
C        ----------------------------------------------------------
         IEFORC = LTDEVO(LPFORC,NOOBVC)
         IF( IEFORC .GT. 0 ) THEN
            DO J=1,10
               X = XYZEF(J,1)
               Y = XYZEF(J,2)
               Z = XYZEF(J,3)
C              VALEUR DES EFFORTS VOLUMIQUES AU NOEUD J DU TETRAEDRE
               CALL REFORC( 4, NOOBVC, 3, X,Y,Z,  0D0,0D0,0D0,
     %                      LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
            ENDDO
         ENDIF
C
C        REDUCTION DU NOMBRE DE DIVISIONS
         Coef = 1D0 / DELTAe
         DO I=1,4
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            S = 0D0
C
            DO L=1,3
C
C              COEFFICIENTS DES EFFORTS VOLUMIQUES COMPOSANTE L
C              Integrale  tDP1 ( 1/Rho P2 Fg(tn+1,Se) dX
               SF = 0D0
               IF( IEFORC .GT. 0 ) THEN
                  DO J=1,10
                     SF = SF + INTP23D(J) * FORCE( L,J) * UnSRho
                  ENDDO
               ENDIF
C
C              - Integrale  Grad P1 ( Vm. Grad ) Vm dx
               DO K=1,3
                  DO J=1,10
C
C                    NO GLOBAL DU NOEUD J DU TETRAEDRE TH
                     NONOJ = NONOTE(J)
C
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NONOJ
                     VJ = VITm(NONOJ,K) * Coef
C
                     DO N=1,3
C
                        D = DFM1(K,N) * VJ
C
                        DO JJ=1,10
C
C                          NO GLOBAL DU NOEUD JJ DU TETRAEDRE TH
                           NONOJJ = NONOTE(JJ)
C
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NONOJJ
                           VJJ = VITm(NONOJJ,L)
C
C                          P2DP2(J,N,JJ) = integrale P2j dP2jj/dxn dx dy dz
                           SF = SF - P2DP2(J,N,JJ) * D * VJJ
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               S = S + DFM1DLa(L,I) * SF
C
            ENDDO
C
C           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTE(I) ) )
C           AU SOMMET DU DL NONOJ DE PRESSION
            NONOJ = NONOSO( NONOTE(I) )
            BG(NONOJ) = BG(NONOJ) + S
C
         ENDDO
C
C        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
C        TOUTE LA SURFACE FRONTIERE DOIT ETRE PRESENTE DANS L'OBJET
C        ==========================================================
         DO NF=1,4
C
C           LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE
            NS = NOOBSF(NF)
            IF( NS .GT. 0 ) THEN
C
C              TOUTE FACE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
C              LE NUMERO DANS LE TETRAEDRE DES 3 SOMMETS DE LA FACE
C              NF DU TETRAEDRE
               NS1 = NONOFAP2(1,NF)
               NS2 = NONOFAP2(2,NF)
               NS3 = NONOFAP2(3,NF)
               NS4 = NONOFAP2(4,NF)
               NS5 = NONOFAP2(5,NF)
               NS6 = NONOFAP2(6,NF)
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
               DO I=1,3
C
C                 COEFFICIENT NONOFAP2(I,NF) DE BE
                  S = 0D0
                  DO L=1,3
C
C                    LA FORCE EXTERIEURE SUR LA SURFACE DE CETTE FACE NF
                     SF = 0D0
                     IF( IEFORC .GT. 0 ) THEN
                        DO J=1,6
                           SF = SF - P1P22D(I,J) * FORCE(L,NSF(J))
                        ENDDO
                        SF = SF * UnSRho
                     ENDIF
C
C                    LE TERME DE TRANSPORT NON LINEAIRE   Vm . Grad Vm
                     SP = 0D0
                     DO K=1,3
                        DO J=1,6
C                          COMPOSANTE K DE LA VITESSEm AU NOEUD NONOTE(NSF(J))
                           VJ = VITm( NONOTE(NSF(J)), K )
                           DO N=1,3
                              D = DFM1(K,N) * VJ
                              DO JJ=1,10
C                                COMPOSANTE L DE LA VITESSE AU NOEUD NONOTE(JJ)
                                 VJJ = VITm(NONOTE(JJ),L)
                                 DO M=1,3
                                    SP = SP + P1P2P12D(I,J,M)
     %                                 * D * DP2S3D(N,JJ,NSF(M)) *VJJ
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
C
                     S = S + VN(L) * ( SF + SP / DELTAe )
C
                  ENDDO
C
C                 ASSEMBLAGE DE Be(NSi) DANS BG( NONOSO( NoSOMMETi ) )
                  NS = NONOSO( NONOTE( NSF(I) ) )
                  BG(NS) = BG(NS) + S
C
               ENDDO
C
            ENDIF
C
C           FIN TRAITEMENT FACE NF
         ENDDO
C
C        FIN TRAITEMENT DE L'EF NEF
 100  CONTINUE
C
      RETURN
      END
