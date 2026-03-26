      SUBROUTINE F3EX2P2P2( XYZEF,  NONOTE,
     %                      NOOBSF, NUMISU, NUMASU, LTDESU,
     %                      NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                      CoGrPr, CoForc, TP2P2,
     %                      NBNOVI, PRESS0,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P2 CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 10 NOEUDS DU TETRAEDRE
C          INTERPOLATION DE FGAMMA AUX  6 NOEUDS DU TRIANGLE FACE
C          Integrale tVitP2 [ -CoGrPr GRAD P(tn) + CoForc FOmega(tn+1) ] dX
C
C ENTREES:
C --------
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C
C NOOBSF : NUMERO DES OBJETS SURFACES DES ARETES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C
C NOOBVC : NUMERO DE VOLUME DE CET ELEMENT FINI
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS VOLUMES
C
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C CoForc : COEFFICIENT DES FORCES
C TP2P2  : Integrale tP2 P2 dX sur LE TETRAEDRE UNITE
C
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES TETRAEDRES
C NTDLVI : NOMBRE DE DL EN VITESSES = 3 * (NBSOMMETS+NBMILIEUX)
C NTDLPR : NOMBRE DE DL EN PRESSION = NBSOMMETS
C NBNOVI : NOMBRE DE NOEUDS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRE DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Juin 2012
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZEF(10,3)
      DOUBLE PRECISION  PRESS0(NBNOVI)
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOOBVC,    NUMIVO, NUMAVO
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER           NBNOVI
      INTEGER           NONOTE(10), I, J, K, L, M, MN, N,
     %                  NU, NOOB, L1, N1
      DOUBLE PRECISION  VE(10,3), FORCE(3,10)
      DOUBLE PRECISION  CoGrPr, CoForc
      DOUBLE PRECISION  DELTAe, DF(3,3), TP2P2(10,10)
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
      DOUBLE PRECISION  X1, Y1, Z1, A, B, C, D, E, S
      INTRINSIC         ABS
C
C     NO DES 3 SOMMETS ET 3 MILIEUX DES ARETES DES 4 FACES DU TETRAEDRE
      INTEGER  NOSFTE(6,4)
      DATA     NOSFTE / 1,3,2, 7, 6,5,
     %                  1,4,3, 8,10,7,
     %                  1,2,4, 5, 9,8,
     %                  2,3,4, 6,10,9 /
C
C     CALCUL de  - CoGrPr Integrale tP2k     Gradk P(t) dX   ou
C                  CoGrPr Integrale tdP2/dxk       P(t) dX
C     ----------------------------------------------------
      CALL F3EX4P2P2( XYZEF,  NONOTE, CoGrPr,
     %                NBNOVI, PRESS0,  VE )
C
C     CONTRIBUTION DES EFFORTS VOLUMIQUES INTERPOLES TAYLOR-HOOD
C     ----------------------------------------------------------
      IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
C
C        VALEUR DES EFFORTS VOLUMIQUES AUX 10 NOEUDS DU TETRAEDRE
         DO J=1,10
            X1 = XYZEF(J,1)
            Y1 = XYZEF(J,2)
            Z1 = XYZEF(J,3)
            CALL REFORC( 4,NOOBVC, 3, X1,Y1,Z1, 0D0,0D0,0D0,
     %                   LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
         ENDDO
C
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
         X1 = XYZEF(1,1)
         Y1 = XYZEF(1,2)
         Z1 = XYZEF(1,3)
         DF(1,1) = XYZEF(2,1) - X1
         DF(1,2) = XYZEF(2,2) - Y1
         DF(1,3) = XYZEF(2,3) - Z1
C
         DF(2,1) = XYZEF(3,1) - X1
         DF(2,2) = XYZEF(3,2) - Y1
         DF(2,3) = XYZEF(3,3) - Z1
C
         DF(3,1) = XYZEF(4,1) - X1
         DF(3,2) = XYZEF(4,2) - Y1
         DF(3,3) = XYZEF(4,3) - Z1
C
C        LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS VOLUMIQUES
C        CoForc Integrale tP2 P2 dX {Force(t,ne)}
         D = CoForc * DELTAe
         DO I=1,10
            DO J=1,10
               S = TP2P2(I,J) * D
               VE(I,1) = VE(I,1) + S * FORCE(1,J)
               VE(I,2) = VE(I,2) + S * FORCE(2,J)
               VE(I,3) = VE(I,3) + S * FORCE(3,J)
            ENDDO
         ENDDO
C
      ENDIF
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES P2 SUR LE TETRAEDRE
C     -------------------------------------------------------------------
      A =  CoForc / 60D0
      B = -CoForc / 72D0
      C = -CoForc / 90D0
      D =  CoForc * 4D0 / 45D0
      E =  CoForc * 2D0 / 45D0
C
      DO K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FORCE?
            MN = LTDESU( LPFORC, NOOB )
            IF( MN .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE K
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A VE
C              ...........................................
C              NO ELEMENTAIRE DU SOMMET 1 DE LA FACE K
               N = NOSFTE(1,K)
               X1 = XYZEF(N,1)
               Y1 = XYZEF(N,2)
               Z1 = XYZEF(N,3)
C
C              NO ELEMENTAIRE DU SOMMET 2 DE LA FACE K
               N = NOSFTE(2,K)
               DGL(1,1) = XYZEF(N,1) - X1
               DGL(1,2) = XYZEF(N,2) - Y1
               DGL(1,3) = XYZEF(N,3) - Z1
C
C              NO ELEMENTAIRE DU SOMMET 3 DE LA FACE K
               N = NOSFTE(3,K)
               DGL(2,1) = XYZEF(N,1) - X1
               DGL(2,2) = XYZEF(N,2) - Y1
               DGL(2,3) = XYZEF(N,3) - Z1
C
C              CALCUL DU VECTEUR NORMAL A LA FACE K SUPPOSEE PLANE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO L=1,3
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE K
                  N = NOSFTE(L,K)
                  GL(1) = XYZEF(N,1)
                  GL(2) = XYZEF(N,2)
                  GL(3) = XYZEF(N,3)
                  CALL REFORC( 3, NOOB, 3,
     %                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     %                         LTDESU(LPFORC,NOOB), FORCE(1,L) )
C
C                 CALCUL DES FORCES AU MILIEU DE L'ARETE L DE LA FACE K
                  IF( L .NE. 3 ) THEN
                     L1 = L + 1
                  ELSE
                     L1 = 1
                  ENDIF
                  N = NOSFTE(L, K)
                  N1= NOSFTE(L1,K)
                  GL(1) = ( XYZEF(N,1) + XYZEF(N1,1) ) * 0.5D0
                  GL(2) = ( XYZEF(N,2) + XYZEF(N1,2) ) * 0.5D0
                  GL(3) = ( XYZEF(N,3) + XYZEF(N1,3) ) * 0.5D0
                  CALL REFORC( 3, NOOB, 3,
     %                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     %                         LTDESU(LPFORC,NOOB), FORCE(1,L+3) )
               ENDDO
C
C              PRISE EN COMPTE DES COEFF Integrale( tP2 P2 ) et CoForc
C
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
C              NUMERO ELEMENTAIRE DU NOEUD 1 DE LA FACE K
               NU = NOSFTE(1,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( A*FORCE(M,1) +B*(FORCE(M,2)+FORCE(M,3)) +C*FORCE(M,5))
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 2 DE LA FACE K
               NU = NOSFTE(2,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( A*FORCE(M,2) +B*(FORCE(M,1)+FORCE(M,3)) +C*FORCE(M,6))
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 3 DE LA FACE K
               NU = NOSFTE(3,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( A*FORCE(M,3) +B*(FORCE(M,1)+FORCE(M,2)) +C*FORCE(M,4))
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 4 DE LA FACE K
               NU = NOSFTE(4,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( C*FORCE(M,3) +D*FORCE(M,4) +E*(FORCE(M,5))+FORCE(M,6))
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 5 DE LA FACE K
               NU = NOSFTE(5,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( C*FORCE(M,1) +D*FORCE(M,5) +E*(FORCE(M,4))+FORCE(M,6))
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 6 DE LA FACE K
               NU = NOSFTE(6,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %          ( C*FORCE(M,2) +D*FORCE(M,6) +E*(FORCE(M,4))+FORCE(M,5))
               ENDDO
C
            ENDIF
         ENDIF
C
      ENDDO
C
      RETURN
      END
