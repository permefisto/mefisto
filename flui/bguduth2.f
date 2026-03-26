      SUBROUTINE BGUDUTH2( CoForc, CoTeNL,
     %                     NBSOM,  NBNOVI, MNXYZN, MNPGEL, NONOSO,
     %                     MNELE,  NBEF,
     %                     NOOBLA,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     VXVY,   BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU SYSTEME
C ----     -CoGrPr LAPLACIEN P = Div( -CoForc Fomega + CoTeNL u. Grad u )
C          -CoGrPr       dP/dn = n .( -CoForc Fomega + CoTeNL u. Grad u )
C          POUR LE TRIANGLE TAYLOR-HOOD
C          PISO Version 1: Bg avec Div reportee sur la Pression en VOLUME
C
C ENTREES:
C --------
C CoForc : Coefficient du terme des FORCES
C CoTeNL : Coefficient du TERME NON LINEAIRE DE NAVIER-STOKES (TRANSPORT)
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C MNXYZN : ADRESSE MCN DU TABLEAU DES COORDONNEES DES NOEUDS
C MNPGEL : ADRESSE MCN DU TABLEAU DES COORDONNEES DES POINTS
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C
C MNELE  : ADRESSE MCN DU TABLEAU NPEF DES TRIANGLES
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NOOBLA : NUMERO DES OBJETS LIGNES  DES ARETES DE L'ELEMENT FINI
C NOOBSF : NUMERO DES OBJETS SURFACES DES FACES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES VOLUMES
C
C VXVY   : COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en X puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Y puis
C          COMPOSANTES AUX SOMMETS ET MILIEUX DE LA VITESSE en Z
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
      INTEGER           NBSOM,  NBNOVI, NBNOEF, NBEF, MNELE,
     %                  MNXYZN, MNPGEL, NONOTR(6), NONOSO(NBNOVI),
     %                  NUMISU, NUMASU
      DOUBLE PRECISION  CoForc, CoTeNL, VXVY(NBNOVI,2), BG(NBSOM)
C
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU),
     %                  NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8),
     %                  NOVCEL, NOSFEL, NOLAEL, NOPSEL, NOOBS
      REAL              XYEF(6,2)
      DOUBLE PRECISION  DELTAe, DFM1(2,2), DFM1DLa(2,3), SP, Coef,
     %                  D, S, SF, X, Y, X21, Y21, X31, Y31, VJ, VJJ,
     %                  FORCE(2,6), VN(2)
      INTEGER           I, J, JJ, K, L, M, N, NEF, IEFORC, NS, NA,
     %                  NONOJ, NONOJJ,
     %                  NS1, NS2, NS3, NSF(3)
      EQUIVALENCE      (NS1,NSF(1)), (NS2,NSF(2)), (NS3,NSF(3))
C
C     INTEGRALES P2 DP2 SUR LE TRIANGLE de REFERENCE
C     P2DP22D(i,k,j) = integrale P2i DP2j/dxk dx dy sur TRIANGLE UNITE
C     DOUBLE PRECISION  P2DP22D(6,2,6)
      include"./incl/p2dp22d.inc"
C
C     [DP2(4 Sommets)]  SUR LE TRIANGLE UNITE
C     DOUBLE PRECISION  DP2S2D(2,6,3)
      include"./incl/dp2s2d.inc"
C
C     INTEGRALE P2i dx SUR LE TRIANGLE UNITE  INTP22D(1:3)=0D0
      DOUBLE PRECISION  INTP22D(4:6)
      DATA              INTP22D/ 0.1666666666666666667D0,
     %                           0.1666666666666666667D0,
     %                           0.1666666666666666667D0 /
C
C     INTEGRALE P1i P2j dx SUR L'ARETE UNITE
      DOUBLE PRECISION  P1P21D(2,3)
      DATA              P1P21D    /
     %     0.16666666666666667D0,   0D0,
     %     0.33333333333333333D0,   0.33333333333333333D0,
     %     0D0,                     0.16666666666666667D0  /
C
C     INTEGRALE P1i P2j P1k dx dy SUR L'ARETE UNITE
      DOUBLE PRECISION  P1P2P11D(2,3,2)
      DATA              P1P2P11D  /
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
C           DES ARETES DU TRIANGLE (INTEGRALE=0 AUX SOMMETS)
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
         Coef = CoTeNL / DELTAe
         DO I=1,3
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            S = 0D0
C
            DO L=1,2
C
C              COEFFICIENTS DES EFFORTS SURFACIQUES L
               SF = 0D0
               IF( IEFORC .GT. 0 ) THEN
                  DO J=4,6
                     SF = SF + INTP22D(J) * FORCE(L,J) * CoForc
                  ENDDO
               ENDIF
C
C              - Integrale  Grad P1  CoTeNL ( u. Grad u) dx
               DO K=1,2
                  DO J=1,6
C
C                    NO GLOBAL DU NOEUD J DU TRIANGLE
                     NONOJ = NONOTR(J)
C
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NONOJ
                     VJ = VXVY(NONOJ,K) * Coef
C
                     DO N=1,2
C
                        D = DFM1(K,N) * VJ
C
                        DO JJ=1,6
C
C                          NO GLOBAL DU NOEUD JJ DU TRIANGLE
                           NONOJJ = NONOTR(JJ)
C
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NONOJJ
                           VJJ = VXVY(NONOJJ,L)
C
C                          P2DP22D(i,k,j) = integrale P2j dP2jj/dxn dx dy
                           SF = SF - P2DP22D(J,N,JJ) * D * VJJ
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
C           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NONOTR(I) ) )
            N = NONOSO( NONOTR(I) )
            BG(N) = BG(N) + S
C
         ENDDO
C
C        CONTRIBUTION DES ARETES FRONTALIERES AU SECOND MEMBRE
C        =====================================================
         DO NA=1,3
C
C           LE NUMERO EVENTUEL DE LA LIGNE DE CETTE ARETE
            NS = NOOBLA(NA)
            IF( NS .GT. 0 ) THEN
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
ccc          ( NON /DELTA COMPENSEE PAR LE JACOBIEN )
C
               DO I=1,2
C
                  S = 0D0
                  DO L=1,2
C
C                    LA FORCE
                     SF = 0D0
                     IF( IEFORC .GT. 0 ) THEN
                        DO J=1,3
                           SF= SF - P1P21D(I,J) *FORCE(L,NSF(J))
                        ENDDO
                        SF = SF * CoForc
                     ENDIF
C
C                    LE TERME DE TRANSPORT NON LINEAIRE   Wm . Grad Wm
                     SP = 0D0
                     DO K=1,2
                        DO J=1,3
C                          COMPOSANTE K DE LA WITESSEm AU NOEUD NONOTR(NSF(J))
                           VJ = VXVY( NONOTR(NSF(J)), K )
                           DO N=1,2
                              D = DFM1(K,N) * VJ
                              DO JJ=1,6
C                                COMPOSANTE L DE LA VITESSE AU NOEUD NONOTR(JJ)
                                 VJJ = VXVY( NONOTR(JJ), L )
                                 DO M=1,2
                                    SP = SP + P1P2P11D(I,J,M)
     %                                      *D *DP2S2D(N,JJ,NSF(M)) *VJJ
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
C
                     S = S + VN(L) * ( SF + SP * CoTeNL / DELTAe )
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
