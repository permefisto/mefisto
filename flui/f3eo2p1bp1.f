      SUBROUTINE F3EO2P1BP1( XYZEF,  NONOTE,
     %                       NOOBSF, NUMISU, NUMASU, LTDESU,
     %                       NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                       DELTAT, CoGrPr, VitAng,
     %                       NBNOVI, VITES0, NBSOM,  PRESS0,
     %                       VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 4 SOMMETS+BARYCENTRE DU TETRAEDRE
C          INTERPOLATION DE FGAMMA AUX 3 SOMMETS DU TRIANGLE FACE
C          dt (-CoGrPr GRAD P(tn) + FOmega(tn+1))
C          -2 Omega x W(tn) - dOmega/dt(tn+1) x R
C           - Omega x ( Omega x r )(tn+1)
C ENTREES:
C --------
C XYZEF  : 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME SOMMET DU TETRAEDRE I=1,...,4
C          NONOTE(5) = NBSOMMET + NUMERO DU TETRAEDRE
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
C DELTAT : PAS DE TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C VitAng : VITESSE ANGULAIRE DE LA ROTATION DE L'OBJET
C
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRES DES TETRAEDRES
C NTDLVI : NOMBRE DE DL EN VITESSES = 3 * (NBSOMMETS+NBTETRAEDRES)
C VITES0 : DL DES 3 COMPOSANTES DE LA VITESSE A L'INSTANT temps
C NTDLPR : NOMBRE DE DL EN PRESSION = NBSOMMETS
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRE DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Juin 2011
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZEF(4,3), TEMPS0
      DOUBLE PRECISION  VITES0(NBNOVI,3), PRESS0(NBSOM)
      DOUBLE PRECISION  DELTAT, CoGrPr, VitAng(3,2)
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOOBVC,    NUMIVO, NUMAVO
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER           NBSOM, NBNOVI
      INTEGER           NONOTE(5), I, J, K, L, M, MN, N, NU, NOOB,
     %                  NOEUDJ
      DOUBLE PRECISION  VE(5,3), FORCE(3,5), COEFPRES
      DOUBLE PRECISION  DELTA, DFM1(3,3), DF(3,3),
     %                  DFM1DLa(3,4), TP1BP1(5,5)
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
      DOUBLE PRECISION  X1, Y1, Z1, A, B, C, D, S
      INTRINSIC         ABS
C
      DOUBLE PRECISION  TP1B(5)
C     TP1B(i) = Integrale sur e chapeau Pi dx dy dz
C
C     NO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER  NOSFTE(3,4)
      DATA     NOSFTE / 1,3,2,   1,4,3,   1,2,4,  2,3,4/
C
C     MISE A ZERO GENERALE EXCEPTE LE PREMIER BLOC DIAGONAL
C     =====================================================
      DO K = 1,3
         DO I = 1,5
            VE(I,K) = 0D0
         ENDDO
      ENDDO
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
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
C     LE DETERMINANT DE DF
      DELTA = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %           + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %           + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) )
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) )
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) )
C
      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) )
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) )
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) )
C
      DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
      DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) )
      DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) )
C
C     [DFM1] [DLambda]
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
C     L'INTEGRALE  - dt CoGrPr  tVitesse  Grad Pression dx
C     ----------------------------------------------------
      DO I=1,4
         TP1B(I) = 73.D0 / 2520.D0
      ENDDO
      TP1B(5) = 16.D0 / 315.D0
C
      COEFPRES = - DELTAT * CoGrPr
      DO I=1,5
         S = TP1B(I) * COEFPRES
         DO J=1,4
C           NO GLOBAL DU NOEUD J DU TETRAEDRE
            NOEUDJ = NONOTE(J)
            DO K=1,3
               VE(I,K) = VE(I,K) + S * DFM1DLa(K,J) * PRESS0(NOEUDJ)
            ENDDO
         ENDDO
      ENDDO
C
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
C     ----------------------------------------------------
      S = DELTAT * DELTA
      A = 29836D0 / 2494800D0 * S
      B =  9046D0 / 2494800D0 * S
      C =   956D0 /  155925D0 * S
      D =  4096D0 /  155925D0 * S
C
CCC      A=  0.11959275292608627D-01
CCC      B=  0.36259419592752926D-02
CCC      C=  0.61311527978194641D-02
CCC      D=  0.26269039602372937D-01
C
      TP1BP1(1,1) = A
      TP1BP1(2,1) = B
      TP1BP1(3,1) = B
      TP1BP1(4,1) = B
      TP1BP1(5,1) = C
C
      TP1BP1(1,2) = B
      TP1BP1(2,2) = A
      TP1BP1(3,2) = B
      TP1BP1(4,2) = B
      TP1BP1(5,2) = C
C
      TP1BP1(1,3) = B
      TP1BP1(2,3) = B
      TP1BP1(3,3) = A
      TP1BP1(4,3) = B
      TP1BP1(5,3) = C
C
      TP1BP1(1,4) = B
      TP1BP1(2,4) = B
      TP1BP1(3,4) = B
      TP1BP1(4,4) = A
      TP1BP1(5,4) = C
C
      TP1BP1(1,5) = C
      TP1BP1(2,5) = C
      TP1BP1(3,5) = C
      TP1BP1(4,5) = C
      TP1BP1(5,5) = D
C
C     CONTRIBUTION DES EFFORTS VOLUMIQUES INTERPOLES BREZZI-FORTIN
C     ------------------------------------------------------------
      IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
C
C        VALEUR DES EFFORTS VOLUMIQUES AUX 4 SOMMETS DU TETRAEDRE
         DO J=1,4
            X1 = XYZEF(J,1)
            Y1 = XYZEF(J,2)
            Z1 = XYZEF(J,3)
            CALL REFORC( 4,NOOBVC, 3, X1,Y1,Z1, 0D0,0D0,0D0,
     %                   LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
         ENDDO
C
C        VALEUR DES EFFORTS VOLUMIQUES AU BARYCENTRE
         X1 = (XYZEF(1,1) + XYZEF(2,1) + XYZEF(3,1) + XYZEF(4,1))*0.25D0
         Y1 = (XYZEF(1,2) + XYZEF(2,2) + XYZEF(3,2) + XYZEF(4,2))*0.25D0
         Z1 = (XYZEF(1,3) + XYZEF(2,3) + XYZEF(3,3) + XYZEF(4,3))*0.25D0
         CALL REFORC( 4,NOOBVC, 3, X1,Y1,Z1, 0D0,0D0,0D0,
     %                LTDEVO(LPFORC,NOOBVC), FORCE(1,5) )
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS VOLUMIQUES
         DO I=1,5
            DO J=1,5
               D = TP1BP1(I,J)
               VE(I,1) = VE(I,1) + D * FORCE(1,J)
               VE(I,2) = VE(I,2) + D * FORCE(2,J)
               VE(I,3) = VE(I,3) + D * FORCE(3,J)
            ENDDO
         ENDDO
C
      ENDIF
C
C     CONTRIBUTION VITESSE ANGULAIRE INTERPOLATION AUX 5 NOEUDS DU TETRAEDRE
C     ----------------------------------------------------------------------
      IF( LTDEVO(LPVIAN,NOOBVC) .GT. 0 ) THEN
C
C        CALCULS AU TEMPS-DELTAT et TEMPS
         TEMPS0 = TEMPS
         DO J=1,5
C
C           NUMERO DU NOEUD J DU TETRAEDRE
            NOEUDJ = NONOTE(J)
C
C           VITESSE ANGULAIRE AU SOMMET J DU TETRAEDRE
            TEMPS = REAL( TEMPS0 - DELTAT )
            IF( J .LE. 4 ) THEN
C              SOMMET J
               X1 = XYZEF(J,1)
               Y1 = XYZEF(J,2)
               Z1 = XYZEF(J,3)
            ELSE
C              BARYCENTRE DU TETRAEDRE
         X1 = (XYZEF(1,1) + XYZEF(2,1) + XYZEF(3,1) + XYZEF(4,1))*0.25D0
         Y1 = (XYZEF(1,2) + XYZEF(2,2) + XYZEF(3,2) + XYZEF(4,2))*0.25D0
         Z1 = (XYZEF(1,3) + XYZEF(2,3) + XYZEF(3,3) + XYZEF(4,3))*0.25D0
            ENDIF
            DO K=1,2
               CALL REVIAN( 4, NOOBVC, X1, Y1, Z1,
     %                      LTDEVO(LPVIAN,NOOBVC), VITANG(1,K) )
               TEMPS = TEMPS0
            ENDDO
C
C           CALCUL DES FORCES AU NOEUD J DUES A LA ROTATION
C           - dOmega/dt x r    - Omega x ( Omega x r )
C           - 2 Omega x W(tn)  - g vecteur z
            FORCE(1,J) =
     %           - (VITANG(2,2)-VITANG(2,1))/DELTAT * Z1
     %           + (VITANG(3,2)-VITANG(3,1))/DELTAT * Y1
     %           + VITANG(3,2) * ( VITANG(3,2) * X1 - VITANG(1,2) * Z1 )
     %           - VITANG(2,2) * ( VITANG(1,2) * Y1 - VITANG(2,2) * X1 )
     %           - 2D0 * ( VITANG(2,2) * VITES0(NOEUDJ,3)
     %                   - VITANG(3,2) * VITES0(NOEUDJ,2) )
C
            FORCE(2,J) =
     %           - (VITANG(3,2)-VITANG(3,1))/DELTAT * X1
     %           + (VITANG(1,2)-VITANG(1,1))/DELTAT * Z1
     %           - VITANG(3,2) * ( VITANG(2,2) * Z1 - VITANG(3,2) * Y1 )
     %           + VITANG(1,2) * ( VITANG(1,2) * Y1 - VITANG(2,2) * X1 )
     %           - 2D0 * ( VITANG(3,2) * VITES0(NOEUDJ,1)
     %                   - VITANG(1,2) * VITES0(NOEUDJ,3) )
C
            FORCE(3,J) =
     %           - (VITANG(1,2)-VITANG(1,1))/DELTAT * Y1
     %           + (VITANG(2,2)-VITANG(2,1))/DELTAT * X1
     %           + VITANG(2,2) * ( VITANG(2,2) * Z1 - VITANG(3,2) * Y1 )
     %           - VITANG(1,2) * ( VITANG(3,2) * X1 - VITANG(1,2) * Z1 )
     %           - 2D0 * ( VITANG(1,2) * VITES0(NOEUDJ,2)
     %                   - VITANG(2,2) * VITES0(NOEUDJ,1) )
C
         ENDDO
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS DE ROTATION
         DO I=1,5
            DO J=1,5
               D = TP1BP1(I,J)
               VE(I,1) = VE(I,1) + D * FORCE(1,J)
               VE(I,2) = VE(I,2) + D * FORCE(2,J)
               VE(I,3) = VE(I,3) + D * FORCE(3,J)
            ENDDO
         ENDDO
C
      ENDIF
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES P1 SUR LE TETRAEDRE
C     -------------------------------------------------------------------
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
               ENDDO
C
C              PRISE EN COMPTE DES COEFF Integrale( tLa La ) et DELTAT
               DGLN = DGLN / 12D0 * DELTAT
C
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
C              NUMERO ELEMENTAIRE DU NOEUD 1 DE LA FACE K
               NU = NOSFTE(1,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %               ( FORCE(M,1) + (FORCE(M,2)+FORCE(M,3))/2D0 )
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 2 DE LA FACE K
               NU = NOSFTE(2,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %               ( FORCE(M,2) + (FORCE(M,1)+FORCE(M,3))/2D0 )
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 3 DE LA FACE K
               NU = NOSFTE(3,K)
               DO M=1,3
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(NU,M) = VE(NU,M) + DGLN *
     %               ( FORCE(M,3) + (FORCE(M,1)+FORCE(M,2))/2D0 )
               ENDDO
            ENDIF
         ENDIF
C
      ENDDO
C
      RETURN
      END
