      SUBROUTINE F3BVETH( CoefVit,CoefFo, CoefNL, CoGrPr,
     %                    XYZEF,  NONOTE, NBNOVI, NBSOM,  NONOSO,
     %                    NOOBSF, NUMISU, NUMASU, LTDESU, IEFOCL,
     %                    NOOBVC, NUMIVO, NUMAVO, LTDEVO, IEFOIN,
     %                    VITtn,  VITm,   PRESSm,
     %                    VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE ELEMENTAIRE PISO + (V.Grad)V
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION D'UN TETRAEDRE TAYLOR-HOOD
C          CoefVit Integrale  tP2 P2 dX V(tn,Se)
C        - CoefNL  Integrale  tP2 P2 Vem(Se) . Grad P2 Vem(Se) dX
C        + CoefFo  Integrale  tP2 P2 Force(tn+1,Se) dX
C        - CoefFo CoGrPr Integrale  tP2 Grad P1 Pm(Se) )  dX
C
C ENTREES:
C --------
C CoefVit: COEFFICIENT DES DERIVEES EN TEMPS
C CoefFo : COEFFICIENT DU GRADIENT DE PRESSION et DES FORCES
C CoefNL : COEFFICIENT DU TERME NON LINEAIRE DE CONVECTION
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES TETRAEDRES
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
C VITtn  : DL DES 3 COMPOSANTES DE LA VITESSE A L'INSTANT tn
C VITm   : DL DES 3 COMPOSANTES DE LA VITESSE A L'INSTANT tn+1,m DU POINT FIXE
C PRESSm : DL DE LA PRESSION(tn+1,m) AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : (10,3) SECOND MEMBRE ELEMENTAIRE DE LA VITESSE A L'INSTANT tn+1,m+1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR      Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
C
      REAL              XYZEF(10,3)
      INTEGER           NBSOM, NBNOVI, NONOTE(10), NONOSO(NBNOVI)
      DOUBLE PRECISION  CoefVit, CoefFo, CoefNL, CoGrPr,
     %                  VITtn(NBNOVI,3), VITm(NBNOVI,3),
     %                  PRESSm(NBSOM), VE(10,3)
      INTEGER           NOOBSF(6), NUMISU, NUMASU,
     %                  NOOBVC,    NUMIVO, NUMAVO
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU), IEFOCL
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO), IEFOIN
C
      INTEGER           I, J, K, L, M, N, NF, NOOB, NOE, L1, N1
      DOUBLE PRECISION  FORCE(3,10)
      DOUBLE PRECISION  DFM1(3,3), DF(3,3), DELTAe, DFM1DLa(3,4),
     %                  GL(3), DGL(2,3), DGLN, VN(3)
      DOUBLE PRECISION  X, Y, Z, S, COEF, A, B
      INTRINSIC         ABS
C
      include "./incl/p2p22d.inc"
C     INTEGRALES P2 P2 SUR LE TRIANGLE P2 de REFERENCE
C     DOUBLE PRECISION  P2P22D(6,6)
C                       P2P22D(i,j) = integrale P2i P2j dX
C
      include "./incl/p2p23d.inc"
C     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
C     DOUBLE PRECISION  P2P23D(10,10)
C                       P2P23D(i,j) = integrale P2i P2j dX
C
      include "./incl/p2p2dp23d.inc"
C     INTEGRALES P2 P2 DP2 SUR LE TETRAEDRE de REFERENCE
C     P2P2DP2(i,j,k,l) = integrale P2i P2j DP2l/dxk dz dx dy
C     DOUBLE PRECISION  P2P2DP2(10,10,3,10)
C
C     NO DES 3 SOMMETS ET 3 MILIEUX DES ARETES DES 4 FACES DU TETRAEDRE
      INTEGER  NOSFTE(6,4)
      DATA     NOSFTE / 1,3,2, 7, 6,5,
     %                  1,4,3, 8,10,7,
     %                  1,2,4, 5, 9,8,
     %                  2,3,4, 6,10,9 /
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      X = XYZEF(1,1)
      Y = XYZEF(1,2)
      Z = XYZEF(1,3)
C
      DF(1,1) = XYZEF(2,1) - X
      DF(1,2) = XYZEF(2,2) - Y
      DF(1,3) = XYZEF(2,3) - Z
C
      DF(2,1) = XYZEF(3,1) - X
      DF(2,2) = XYZEF(3,2) - Y
      DF(2,3) = XYZEF(3,3) - Z
C
      DF(3,1) = XYZEF(4,1) - X
      DF(3,2) = XYZEF(4,2) - Y
      DF(3,3) = XYZEF(4,3) - Z
C
C     LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1
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
C     [DFM1] [DLa] (3,4)
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
C     CoefVit Integrale  tP2 1/dt P2 dX  V(tn,Se)
C     -------------------------------------------
      IF( CoefVit .EQ. 0D0 ) THEN
C
C        CAS STATIONNAIRE
         CALL AZEROD( 30, VE )
C
      ELSE
C
C        CAS INSTATIONNAIRE
         COEF = CoefVit * DELTAe
         DO I=1,10
            DO K=1,3
               S = 0D0
               DO J=1,10
                  S = S + P2P23D(I,J) * VITtn(NONOTE(J),K)
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
      DO K=1,3
         DO I=1,10
            S = 0D0
            DO L=1,3
               DO M=1,10
                  NOE = NONOTE( M )
                  A   = VITm(NOE,L)
                  DO N=1,10
                     B    = VITm(NONOTE(N),K) * A
                     DO J=1,3
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
C     INTERPOLATION DE FOMEGA AUX 10 NOEUDS DU TETRAEDRE
C     --------------------------------------------------
      IF( IEFOIN .GT. 0 ) THEN
C
C        IL EXISTE UN VOLUME AVEC FORCE INTERNE
         IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS VOLUMIQUES AUX 10 NOEUDS DU TETRAEDRE
            DO J=1,10
               X = XYZEF(J,1)
               Y = XYZEF(J,2)
               Z = XYZEF(J,3)
               CALL REFORC( 4,NOOBVC, 3, X,Y,Z, 0D0,0D0,0D0,
     %                      LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
            ENDDO
C
C           COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS VOLUMIQUES
            COEF = CoefFo * DELTAe
            DO I=1,10
               DO J=1,10
                  S = P2P23D(I,J) * COEF
                  VE(I,1) = VE(I,1) + S * FORCE(1,J)
                  VE(I,2) = VE(I,2) + S * FORCE(2,J)
                  VE(I,3) = VE(I,3) + S * FORCE(3,J)
               ENDDO
            ENDDO
C
         ENDIF
C
      ENDIF
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES P2 SUR LE TETRAEDRE
C     INTERPOLATION DE FGAMMA AUX 6 NOEUDS DU TRIANGLE FACE
C     -------------------------------------------------------------------
      IF( IEFOCL .LE. 0 ) GOTO 20
C
C     IL EXISTE UNE SURFACE AVEC FORCE INTERNE
      DO NF=1,4
C
C        NO DE SURFACE DE LA FACE NF
         NOOB = NOOBSF(NF)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FORCE?
            IF( LTDESU( LPFORC, NOOB ) .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE NF
C              CALCUL DE LA CONTRIBUTION DE LA FACE NF A VE
C              ............................................
C              NO ELEMENTAIRE DU SOMMET 1 DE LA FACE NF
               N = NOSFTE(1,NF)
               X = XYZEF(N,1)
               Y = XYZEF(N,2)
               Z = XYZEF(N,3)
C
C              NO ELEMENTAIRE DU SOMMET 2 DE LA FACE NF
               N = NOSFTE(2,NF)
               DGL(1,1) = XYZEF(N,1) - X
               DGL(1,2) = XYZEF(N,2) - Y
               DGL(1,3) = XYZEF(N,3) - Z
C
C              NO ELEMENTAIRE DU SOMMET 3 DE LA FACE NF
               N = NOSFTE(3,NF)
               DGL(2,1) = XYZEF(N,1) - X
               DGL(2,2) = XYZEF(N,2) - Y
               DGL(2,3) = XYZEF(N,3) - Z
C
C              CALCUL DU VECTEUR NORMAL A LA FACE NF SUPPOSEE PLANE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO L=1,3
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE NF
                  N = NOSFTE(L,NF)
                  GL(1) = XYZEF(N,1)
                  GL(2) = XYZEF(N,2)
                  GL(3) = XYZEF(N,3)
                  CALL REFORC( 3, NOOB, 3,
     %                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     %                         LTDESU(LPFORC,NOOB), FORCE(1,L) )
C
C                 CALCUL DES FORCES AU MILIEU DE L'ARETE L DE LA FACE NF
                  IF( L .NE. 3 ) THEN
                     L1 = L + 1
                  ELSE
                     L1 = 1
                  ENDIF
                  N = NOSFTE(L, NF)
                  N1= NOSFTE(L1,NF)
                  GL(1) = ( XYZEF(N,1) + XYZEF(N1,1) ) * 0.5D0
                  GL(2) = ( XYZEF(N,2) + XYZEF(N1,2) ) * 0.5D0
                  GL(3) = ( XYZEF(N,3) + XYZEF(N1,3) ) * 0.5D0
                  CALL REFORC( 3, NOOB, 3,
     %                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     %                         LTDESU(LPFORC,NOOB), FORCE(1,L+3) )
               ENDDO
C
C              PRISE EN COMPTE DES COEFF Integrale( tP2 P2 ) et CoefFo
               COEF = CoefFo * DGLN
               DO I=1,6
                  NOE = NOSFTE(I,NF)
                  DO M=1,3
                     S = 0D0
                     DO J=1,6
                        S = S + P2P22D(I,J) * FORCE(1,J)
                     ENDDO
                     VE(NOE,M) = VE(NOE,M) + COEF * S
                  ENDDO
               ENDDO
C
            ENDIF
         ENDIF
C
      ENDDO
C
C     - CoGrPr Integrale  tP2 Grad P1 Pm(Se) )  dX / Rho
C     --------------------------------------------------
C     Integrale P2i dX = -1/120 i=1,4
 20   A = CoGrPr * CoefFo * DELTAe / 120D0
C     Integrale P2i dX =1/30 i=5,10
      B = A * 4D0
C
      DO K=1,3
         DO I=1,4
            S = 0D0
            DO J=1,4
               S = S + DFM1DLa(K,J) * PRESSm(NONOSO(NONOTE(J)))
            ENDDO
            VE(I,K) = VE(I,K) + A * S
         ENDDO
C
         DO I=5,10
            S = 0D0
            DO J=1,4
               S = S + DFM1DLa(K,J) * PRESSm(NONOSO(NONOTE(J)))
            ENDDO
            VE(I,K) = VE(I,K) - B * S
         ENDDO
      ENDDO
C
      RETURN
      END
