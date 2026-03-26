      SUBROUTINE F2BPGTH( UnSRho,
     %                    NBSOM,  NBNOVI, MNXYZN, MNPGEL, NONOSO,
     %                    MNELE,  NBEF,
     %                    NOOBLA,
     %                    NOOBSF, NUMISU, NUMASU, LTDESU,
     %                    VITm,   P2DP2,  BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:  CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU SYSTEME (2)
C ----
C (2) - 1/Rho Laplacien Pm+1(X) = Div { -Fg(tn+1)/Rho
C                                     + Vm . Grad Vm } dans Omega
C     Condition aux limites
C            -1/Rho    dPm+1/dn = n . ( -Fg(tn+1)/Rho
C                                     +  Vm . Grad Vm } sur Gamma
C     POUR LE TRIANGLE TAYLOR-HOOD
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
C
C MNELE  : ADRESSE MCN DU TABLEAU NPEF DES TRIANGLES
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NOOBSF : NUMERO DES OBJETS SURFACES DES FACES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C
C VITm   : COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en X puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Y puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Z
C P2DP2  : P2DP2(i,k,j) = Integrale P2i DP2j/Dxk dX SUR LE TRIANGLE UNITE
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
     %                  VITm(NBNOVI,2), BG(NBSOM)
      INTEGER           NBSOM, NBNOVI, NBNOEF, NBEF, MNXYZN, MNPGEL,
     %                  NONOSO(NBNOVI), MNELE,
     %                  NUMISU, NUMASU
C
      INTEGER           NONOTR(6)
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU),
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8),
     %                  NOVCEL, NOSFEL, NOLAEL, NOPSEL, NOOBS
C
      REAL              XYEF(6,2)
      DOUBLE PRECISION  FORCE(2,6), P2DP2(6,2,6), VJ, Coef, SF, SP,
     %                  X21, Y21, X31, Y31, DELTAe, DFM1(2,2),
     %                  DFM1DLa(2,3), S, X, Y, VN(2)
C
      INTEGER           NA, I, J, JJ, K, L, M, N, NONOJ, NONOJJ, NS,
     %                  NEF, IEFORC, NS1, NS2, NS3, NSF(3)
      EQUIVALENCE      (NSF(1),NS1),(NSF(2),NS2),(NSF(3),NS3)
C
C     INTEGRALE P2i dx SUR LE TRIANGLE UNITE  INTP2(1:3)=0D0
      DOUBLE PRECISION  INTP2(4:6)
      DATA              INTP2/ 0.1666666666666666667D0,
     %                         0.1666666666666666667D0,
     %                         0.1666666666666666667D0 /
C
C     [DP2(3Sommets)]   SUR LE TRIANGLE de REFERENCE
      DOUBLE PRECISION  DP2S(2,6,3)
      DATA              DP2S /
     %             -0.3D+01, -0.3D+01,
     %             -0.1D+01,  0.0D+00,
     %              0.0D+00, -0.1D+01,
     %              0.4D+01,  0.0D+00,
     %              0.0D+00,  0.0D+00,
     %              0.0D+00,  0.4D+01,
     %              0.1D+01,  0.1D+01,
     %              0.3D+01,  0.0D+00,
     %              0.0D+00, -0.1D+01,
     %             -0.4D+01, -0.4D+01,
     %              0.0D+00,  0.4D+01,
     %              0.0D+00,  0.0D+00,
     %              0.1D+01,  0.1D+01,
     %             -0.1D+01,  0.0D+00,
     %              0.0D+00,  0.3D+01,
     %              0.0D+00,  0.0D+00,
     %              0.4D+01,  0.0D+00,
     %             -0.4D+01, -0.4D+01  /
C
C     INTEGRALE P1i P2j dx SUR L'ARETE UNITE
      DOUBLE PRECISION  P1P2(2,3)
      DATA              P1P2    /
     %     0.16666666666666667D0,   0D0,
     %     0.33333333333333333D0,   0.33333333333333333D0,
     %     0D0,                     0.16666666666666667D0  /
C
C     INTEGRALE P1i P2j P1k dx dy SUR L'ARETE UNITE
      DOUBLE PRECISION  P1P2P1(2,3,2)
      DATA              P1P2P1  /
     %     0.15D0,                   0.16666666666666667D-1,
     %    -0.1666666666666667D-1,    0.16666666666666667D-1,
     %     0.2D0,                    0.13333333333333333D0,
     %     0.16666666666666667D-1,  -0.16666666666666667D-1,
     %     0.16666666666666667D-1,   0.15D0,
     %     0.13333333333333333D0,    0.2D0     /
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NO DES NOEUDS DE L'ELEMENT FINI NEF
         CALL EFNOEU( MNELE, NEF, NBNOEF, NONOTR )
C
C        NO DE  POINTS  LIGNES SURFACES VOLUME
C           DES SOMMETS FACES  FACES    VOLUME DU TRIANGLE NEF
         CALL EFPLSV( MNELE , NEF,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C        COORDONNEES DES NBNOEF SOMMETS=POINTS=NOEUDS DE L'EF
         CALL EFXYZP( 2, MNXYZN, NBEF, NEF, MNPGEL, NBNOEF,  XYEF )
C
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X21 = XYEF(2,1) - XYEF(1,1)
         X31 = XYEF(3,1) - XYEF(1,1)
C
         Y21 = XYEF(2,2) - XYEF(1,2)
         Y31 = XYEF(3,2) - XYEF(1,2)
C
C        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DELTAe = ABS( X21*Y31 - X31*Y21 )
C
C        LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
C        LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
         DFM1(1,1) =  Y31
         DFM1(2,1) = -X31
C
         DFM1(1,2) = -Y21
         DFM1(2,2) =  X21
C
C        PRODUIT des MATRICES DFM1(2,2) DLambda(2,3)
         DFM1DLa(1,1) = Y21 - Y31
         DFM1DLa(2,1) = X31 - X21
C
         DFM1DLa(1,2) = Y31
         DFM1DLa(2,2) =-X31
C
         DFM1DLa(1,3) =-Y21
         DFM1DLa(2,3) = X21
C
C        CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES TAYLOR-HOOD
C        -----------------------------------------------------------
         NOOBS  = NOOBSF(1)
         IEFORC = LTDESU(LPFORC,NOOBS)
         IF( IEFORC .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS SURFACIQUES AUX 3 SOMMETS ET MILIEUX
C           DES ARETES DU TRIANGLE  (INTEGRALE=0 AUX SOMMETS)
C           MAIS VALEUR UTILE POUR LES FORCES SUR LES ARETES
            DO J=1,6
               X = XYEF(J,1)
               Y = XYEF(J,2)
               CALL REFORC( 3,NOOBS, 2, X,Y,0D0,  0D0,0D0,0D0,
     %                      LTDESU(LPFORC,NOOBS), FORCE(1,J) )
            ENDDO
C
         ENDIF
C
C        REDUCTION DU NOMBRE DE DIVISIONS
         Coef = 1D0 / DELTAe
         DO I=1,3
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            S = 0D0
C
            DO L=1,2
C
C              COEFFICIENTS DES EFFORTS VOLUMIQUES COMPOSANTE L
C              Integrale  tDP1 ( 1/Rho P2 Fg(tn+1,Se) dX
               SF = 0D0
               IF( IEFORC .GT. 0 ) THEN
                  DO J=4,6
                     SF = SF + INTP2(J) * FORCE(L,J) * UnSRho
                  ENDDO
               ENDIF
C
C              - Integrale  Grad P1 ( Vm. Grad ) Vm dx
               DO K=1,2
                  DO J=1,6
C
C                    NO GLOBAL DU NOEUD J DU TRIANGLE TH
                     NONOJ = NONOTR(J)
C
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NONOJ
                     VJ = VITm(NONOJ,K)
C
                     DO N=1,2
                        DO JJ=1,6
C
C                          NO GLOBAL DU NOEUD JJ DU TRIANGLE TH
                           NONOJJ = NONOTR(JJ)
C
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NONOJJ
ccc                        VJJ = VITm(NONOJJ,L)
C
C                          P2DP2(J,N,JJ) = integrale P2j dP2jj/dxn dx dy dz
                           SF = SF - Coef * P2DP2(J,N,JJ) * VJ
     %                                    * DFM1(K,N) * VITm(NONOJJ,L)
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
            NONOJ = NONOSO( NONOTR(I) )
            BG(NONOJ) = BG(NONOJ) + S
C
         ENDDO
C
C        CONTRIBUTION DES FACES FRONTALIERES AU SECOND MEMBRE
C        TOUTE LA SURFACE FRONTIERE DOIT ETRE PRESENTE DANS L'OBJET
C        ==========================================================
         DO NA=1,3
C
C           LE NUMERO EVENTUEL DE LA LIGNE DE CETTE ARETE
            NOOBS = NOOBLA(NA)
            IF( NOOBS .GT. 0 ) THEN
C
C              TOUTE ARETE FRONTIERE EST TRAITEE SANS PLUS DE DONNEES
C              LE NUMERO DANS LE TRIANGLE DES 2 SOMMETS DE L'ARETE NA
               NS1 = NA
               IF( NS1 .NE. 3 ) THEN
                  NS2 = NS1 + 1
               ELSE
                  NS2 = 1
               ENDIF
               NS3 = NS1 + 3
C
C              LE VECTEUR NORMAL UNITAIRE A L'ARETE
               VN(1) = XYEF(NS2,2) - XYEF(NS1,2)
               VN(2) = XYEF(NS1,1) - XYEF(NS2,1)
C
ccc            DELTAK= SQRT( VN(1)**2 + VN(2)**2 )
ccc            VN(1) = VN(1) / DELTAK
ccc            VN(2) = VN(2) / DELTAK
ccc          ( NON /DELTA CAR COMPENSEE PAR LE JACOBIEN )
C
               DO I=1,2
C
                  S = 0D0
                  DO L=1,3
C
C                    LA FORCE EXTERIEURE SUR LA LIGNE DE CETTE ARETE NA
                     SF = 0D0
                     IF( IEFORC .GT. 0 ) THEN
                        DO J=1,3
                           SF = SF -P1P2(I,J) * FORCE(L,NSF(J)) * UnSRho
                        ENDDO
                     ENDIF
C
C                    LE TERME DE TRANSPORT NON LINEAIRE   Vm . Grad Vm
                     SP = 0D0
                     DO K=1,2
                        DO J=1,3
                           DO N=1,2
                              DO JJ=1,6
                                 DO M=1,2
                               SP = SP + VITm(NONOTR(NSF(J)),K)
     %                                 * P1P2P1(NSF(I),NSF(J),NSF(M))
     %                                 * DFM1(K,N) * DP2S(N,JJ,NSF(M))
     %                                 * VITm(NONOTR(JJ),L)
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
                  NS = NONOSO( NONOTR( NSF(I) ) )
                  BG(NS) = BG(NS) + S
C
               ENDDO
C
            ENDIF
C
C           FIN TRAITEMENT ARETE NA
         ENDDO
C
C        FIN TRAITEMENT DE L'EF NEF
 100  CONTINUE
C
      RETURN
      END
