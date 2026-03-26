      SUBROUTINE F2BVETH( CoefVit,CoefFo, CoefNL, CoGrPr,
     %                    XYEF,   NONOTR, NBNOVI, NBSOM,  NONOSO,
     %                    NOOBLA, NUMILI, NUMALI, LTDELI, IEFOCL,
     %                    NOOBSF, NUMISU, NUMASU, LTDESU, IEFOIN,
     %                    VITtn,  VITm,   PRESSm,
     %                    VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE ELEMENTAIRE EOL PISO + (V.Grad)V
C -----    P2 CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION D'UN TRIANGLE TAYLOR-HOOD
C          CoefVit  Integrale  tP2 P2 dX V(tn,Se)
C        - CoefNL Integrale  tP2 P2 Vem(Se) . Grad P2 Vem(Se) dX
C        + CoefFo Integrale  tP2 P2 Force(tn+1,Se) dX
C        - CoGrPr CoefFo Integrale  tP2 Grad P1 Pm(Se) )  dX

C ENTREES:
C --------
C CoefVit: COEFFICIENT DES DERIVEES EN TEMPS
C CoefFo : COEFFICIENT DU GRADIENT DE PRESSION et DES FORCES
C CoefNL : COEFFICIENT DU TERME NON LINEAIRE DE CONVECTION
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C
C XYEF   : 2 COORDONNEES DES 6 NOEUDS DU TRIANGLE
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME NOEUD DU TRIANGLE I=1,...,6
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES TRIANGLES
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION = NOMBRE DE DL DE LA PRESSION
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL
C
C NOOBSF : NUMERO DES OBJETS SURFACES DES ARETES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C IEFOCL :=0 SI PAS DE FORCES DONNEES SUR LES FRONTIERES
C         >0 NOMBRE DE TMS FORCE "AUX LIMITES" DES PLS DE L'OBJET RETROUVES
C
C NOOBVC : NUMERO DE VOLUME DE CET ELEMENT FINI
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS VOLUMES
C IEFOIN :=0 SI PAS DE FORCES INTERNES AU VOLUME DONNEES
C         >0 NOMBRE DE TMS FORCE INTERNE     DES SV DE L'OBJET RETROUVES
C
C VITtn  : DL DES 2 COMPOSANTES DE LA VITESSE A L'INSTANT tn
C VITm   : DL DES 2 COMPOSANTES DE LA VITESSE A L'INSTANT tn+1,m DU POINT FIXE
C PRESSm : DL DE LA PRESSION(tn+1,m) AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : (6,2) SECOND MEMBRE ELEMENTAIRE DE LA VITESSE A L'INSTANT tn+1,m+1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR      Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
C
      REAL              XYEF(6,2)
      INTEGER           NBSOM, NBNOVI, NONOTR(6), NONOSO(NBNOVI)
      DOUBLE PRECISION  CoefVit, CoefFo, CoefNL, CoGrPr,
     %                  VITtn(NBNOVI,2), VITm(NBNOVI,2),
     %                  PRESSm(NBSOM), VE(6,2)
      INTEGER           NOOBSF(6),  NUMISU, NUMASU,
     %                  NOOBLA(12), NUMILI, NUMALI
      INTEGER           LTDELI(1:MXDOFL,NUMILI:NUMALI), IEFOCL
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU), IEFOIN
C
      INTEGER           I, J, K, L, M, N, NOOB, NOE, NOEN, NA
      DOUBLE PRECISION  FORCE(2,6)
      DOUBLE PRECISION  DFM1(2,2), DELTAe, DFM1DLa(2,3),
     %                  DGLN, VN(2)
      DOUBLE PRECISION  X, Y, X21, Y21, X31, Y31, A, B, C, D, S, COEF
      INTEGER           NS1, NS2, NS3, NSA(3)
      EQUIVALENCE      (NS1,NSA(1)), (NS2,NSA(2)), (NS3,NSA(3))
      INTRINSIC         ABS
C
      include "./incl/p2p22d.inc"
C     INTEGRALES P2 P2 SUR LE TRIANGLE P2 de REFERENCE
C     DOUBLE PRECISION  P2P22D(6,6)
C                       P2P22D(i,j) = integrale P2i P2j dX
C
      include "./incl/p2p2dp22d.inc"
C     INTEGRALES P2 P2 DP2 SUR LE TRIANGLE de REFERENCE
C     P2P2DP2(i,j,k,l) = integrale P2i P2j DP2l/dxk dz dx dy
C     DOUBLE PRECISION  P2P2DP2(6,6,2,6)
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      X = XYEF(1,1)
      Y = XYEF(1,2)
C
      X21 = XYEF(2,1) - X
      Y21 = XYEF(2,2) - Y
      X31 = XYEF(3,1) - X
      Y31 = XYEF(3,2) - Y
C
C     CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
      DELTAe = ABS( X21*Y31 - X31*Y21 )
C
C     L'ARETE EST SUPPOSE DE SURFACE NON NULLE
C     LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1
      DFM1(1,1) =  Y31 / DELTAe
      DFM1(2,1) = -X31 / DELTAe
C
      DFM1(1,2) = -Y21 / DELTAe
      DFM1(2,2) =  X21 / DELTAe
C
C     DFM1 DLambda (2,3)
      DFM1DLa(1,1) = Y21 - Y31
      DFM1DLa(2,1) = X31 - X21
C
      DFM1DLa(1,2) = Y31
      DFM1DLa(2,2) =-X31
C
      DFM1DLa(1,3) =-Y21
      DFM1DLa(2,3) = X21
C
C     CoefVit Integrale  tP2 1/dt P2 dX  V(tn,Se)
C     -------------------------------------------
      IF( CoefVit .EQ. 0D0 ) THEN
C
C        CAS STATIONNAIRE
         CALL AZEROD( 12, VE )
C
      ELSE
C
C        CAS INSTATIONNAIRE
         COEF = CoefVit * DELTAe
         DO I=1,6
            DO K=1,2
               S = 0D0
               DO J=1,6
                  S = S + P2P22D(I,J) * VITtn(NONOTR(J),K)
               ENDDO
               VE(I,K) = S * COEF
            ENDDO
         ENDDO
C
      ENDIF
C
C     - CoefNL Integrale  tP2 (P2 Vem(Se) . Grad) P2 Vem(Se) dX  COMPOSANTE k
C     -----------------------------------------------------------------------
      COEF = DELTAe * CoefNL
      DO K=1,2
         DO I=1,6
            S = 0D0
            DO L=1,2
               DO M=1,6
                  NOE = NONOTR( M )
                  A   = Vitm(NOE,L)
                  DO N=1,6
                     NOEN = NONOTR( N )
                     B    = Vitm(NOEN,K) * A
                     DO J=1,2
                        S = S + DFM1(L,J) * P2P2DP2(I,M,J,N) * B
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            VE(I,K) = VE(I,K) - COEF * S
         ENDDO
      ENDDO
C
C     + Integrale  CoefFo tP2 P2 Force(tn+1,Se) dX
C     INTERPOLATION DE FOMEGA AUX 6 NOEUDS DU TRIANGLE
C     ------------------------------------------------
      IF( IEFOIN .GT. 0 ) THEN
C
C        IL EXISTE UNE SURFACE AVEC FORCE INTERNE
         NOOB = NOOBSF(1)
         IF( LTDESU(LPFORC,NOOB) .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS VOLUMIQUES AUX 6 NOEUDS DU TRIANGLE
            DO J=1,6
               X = XYEF(J,1)
               Y = XYEF(J,2)
               CALL REFORC( 3,NOOB, 2, X,Y,0D0,  0D0,0D0,0D0,
     %                      LTDESU(LPFORC,NOOB), FORCE(1,J) )
            ENDDO
C
C           COEFFICIENTS DU SECOND MEMBRE LIES AUX FORCES SURFACIQUES
            COEF = CoefFo * DELTAe
            DO I=1,6
               DO J=1,6
                  S = P2P22D(I,J) * COEF
                  VE(I,1) = VE(I,1) + S * FORCE(1,J)
                  VE(I,2) = VE(I,2) + S * FORCE(2,J)
               ENDDO
            ENDDO
C
         ENDIF
      ENDIF
C
C     CONTRIBUTION DES FORCES SUR LES ARETES DU TRIANGLE
C     INTERPOLATION DE FGAMMA AUX 3 NOEUDS DE L'ARETE
C     --------------------------------------------------
      IF( IEFOCL .LE. 0 ) GOTO 20
C
C     Integrales P2i P2j dX sur l'ARETE UNITE
      A =  2D0 / 15D0
      B =  8D0 / 15D0
      C = -1D0 / 30D0
      D =  1D0 / 15D0
C
      DO NA=1,3
C
C        NO DE LIGNE DE L'ARETE NA
         NOOB = NOOBLA(NA)
         IF( NOOB .GT. 0 ) THEN
            IF( LTDELI( LPFORC, NOOB ) .GT. 0 ) THEN
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
C              LE VECTEUR NORMAL A L'ARETE NA
               VN(1) = XYEF(NS2,2) - XYEF(NS1,2)
               VN(2) = XYEF(NS1,1) - XYEF(NS2,1)
C
C              LE JACOBIEN DE L'ARETE NA
               DGLN = SQRT( VN(1)**2 + VN(2)**2 )
C
C              LES 2 COMPOSANTES DE LA FORCE ELEMENTAIRE AUX 3 NOEUDS
C              DE L'ARETE NA
               DO J=1,3
                  K = NSA(J)
                  X = XYEF(K,1)
                  Y = XYEF(K,2)
                  CALL REFORC( 2, NOOB, 2,
     %                         X, Y, 0D0,  VN(1), VN(2), 0D0,
     %                         LTDELI(LPFORC,NOOB), FORCE(1,J) )
               ENDDO
C
               DO K=1,2
C                 COMPOSANTE K AU NOEUD NS1
                  VE(NS1,K) = VE(NS1,K) + DGLN *
     %          ( A*FORCE(K,1) + C*FORCE(K,2) + D*FORCE(K,3) )
               ENDDO
C
               DO K=1,2
C                 COMPOSANTE K AU NOEUD NS2
                  VE(NS2,K) = VE(NS2,K) + DGLN *
     %          ( C*FORCE(K,1) + A*FORCE(K,2) + D*FORCE(K,3) )
               ENDDO
C
               DO K=1,2
C                 COMPOSANTE K AU NOEUD NS3
                  VE(NS3,K) = VE(NS3,K) + DGLN *
     %          ( D*(FORCE(K,1) + FORCE(K,2)) + B*FORCE(K,3) )
               ENDDO
C
            ENDIF
         ENDIF
C
C     FIN TRAITEMENT ARETE NA
      ENDDO
C
C     - CoGrPr Integrale  tP2 Grad P1 Pm(Se) )  dX / Rho
C     -------------------------------------------
C     Integrale P2i dX = 0D0 i=1,3 et 1/6 pour i=4,6
 20   A = CoGrPr * CoefFo * DELTAe / 6D0
C
      DO K=1,2
         DO I=4,6
            S = 0D0
            DO J=1,3
               S = S + DFM1DLa(K,J) * PRESSm(NONOSO(NONOTR(J)))
            ENDDO
            VE(I,K) = VE(I,K) - A * S
         ENDDO
      ENDDO
C
      RETURN
      END
