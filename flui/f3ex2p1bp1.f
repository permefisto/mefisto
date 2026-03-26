      SUBROUTINE F3EX2P1BP1( XYZSOM, NUSOTE,
     %                       NOOBSF, NUMISU, NUMASU, LTDESU,
     %                       NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                       DELTAT, CoGrPr, NBSOM,  PRESS0,
     %                       VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 4 SOMMETS+BARYCENTRE DU TETRAEDRE
C          INTERPOLATION DE FGAMMA AUX 3 SOMMETS DU TRIANGLE FACE
C          -dt CoGrPr GRAD P(tn) + dt FOmega(tn+1)
C
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TETRAEDRISATION
C NUSOTE : NUSOTE(I) NO GLOBAL DU I-EME SOMMET DU TETRAEDRE I=1,...,4
C          NUSOTE(5) = NBSOMMET + NUMERO DU TETRAEDRE
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
C
C NTDLVI : NOMBRE DE DL EN VITESSES = 3 * (NBSOMMETS+NBTETRAEDRES)
C NTDLPR : NOMBRE DE DL EN PRESSION = NBSOMMETS
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRE DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
C
      INTEGER           NBSOM
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  PRESS0(NBSOM)
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOOBVC,    NUMIVO, NUMAVO
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER           NUSOTE(5), I, J, K, L, M, MN, N, NU, NOOB
      DOUBLE PRECISION  VE(5,3), FORCE(3,5)
      DOUBLE PRECISION  DELTAT,  CoGrPr
      DOUBLE PRECISION  DELTAe, DF(3,3), TP1BP1(5,5)
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
      DOUBLE PRECISION  X1, Y1, Z1, A, B, C, D, S
      INTRINSIC         ABS
C
C     NO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER  NOSFTE(3,4)
      DATA     NOSFTE / 1,3,2,   1,4,3,   1,2,4,  2,3,4/
C
C     CALCUL de  -dt CoGrPr Integrale tP1B GRADk P1(t) dX
C     ---------------------------------------------------
      CALL F3EX4P1BP1( XYZSOM, NUSOTE, DELTAT, CoGrPr,
     %                 NBSOM,  PRESS0, VE )
C
C     CONTRIBUTION DES EFFORTS VOLUMIQUES INTERPOLES BREZZI-FORTIN
C     ------------------------------------------------------------
      IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN

C        VALEUR DES EFFORTS VOLUMIQUES AUX 4 SOMMETS DU TETRAEDRE
         DO J=1,4
            I = NUSOTE(J)
            X1 = XYZSOM(1,I)
            Y1 = XYZSOM(2,I)
            Z1 = XYZSOM(3,I)
            CALL REFORC( 4,NOOBVC, 3, X1,Y1,Z1, 0D0,0D0,0D0,
     %                   LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
         ENDDO

C        NO DES 4 SOMMETS
         I = NUSOTE(1)
         J = NUSOTE(2)
         K = NUSOTE(3)
         L = NUSOTE(4)

C        VALEUR DES EFFORTS VOLUMIQUES AU BARYCENTRE
         X1 =(XYZSOM(1,I) +XYZSOM(1,J) +XYZSOM(1,K) +XYZSOM(1,L))*0.25D0
         Y1 =(XYZSOM(2,I) +XYZSOM(2,J) +XYZSOM(2,K) +XYZSOM(2,L))*0.25D0
         Z1 =(XYZSOM(3,I) +XYZSOM(3,J) +XYZSOM(3,K) +XYZSOM(3,L))*0.25D0
         CALL REFORC( 4,NOOBVC, 3, X1,Y1,Z1, 0D0,0D0,0D0,
     %                LTDEVO(LPFORC,NOOBVC), FORCE(1,5) )
C
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
         X1 = XYZSOM(1,I)
         Y1 = XYZSOM(2,I)
         Z1 = XYZSOM(3,I)

         DF(1,1) = XYZSOM(1,J) - X1
         DF(1,2) = XYZSOM(2,J) - Y1
         DF(1,3) = XYZSOM(3,J) - Z1
C
         DF(2,1) = XYZSOM(1,K) - X1
         DF(2,2) = XYZSOM(2,K) - Y1
         DF(2,3) = XYZSOM(3,K) - Z1
C
         DF(3,1) = XYZSOM(1,L) - X1
         DF(3,2) = XYZSOM(2,L) - Y1
         DF(3,3) = XYZSOM(3,L) - Z1
C
C        LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C        LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
         S = DELTAT * DELTAE
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
               N = NUSOTE( NOSFTE(1,K) )
               X1 = XYZSOM(1,N)
               Y1 = XYZSOM(2,N)
               Z1 = XYZSOM(3,N)
C
C              NO ELEMENTAIRE DU SOMMET 2 DE LA FACE K
               N = NUSOTE( NOSFTE(2,K) )
               DGL(1,1) = XYZSOM(1,N) - X1
               DGL(1,2) = XYZSOM(2,N) - Y1
               DGL(1,3) = XYZSOM(3,N) - Z1
C
C              NO ELEMENTAIRE DU SOMMET 3 DE LA FACE K
               N = NUSOTE( NOSFTE(3,K) )
               DGL(2,1) = XYZSOM(1,N) - X1
               DGL(2,2) = XYZSOM(2,N) - Y1
               DGL(2,3) = XYZSOM(3,N) - Z1
C
C              CALCUL DU VECTEUR NORMAL A LA FACE K SUPPOSEE PLANE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO L=1,3
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE K
                  N = NUSOTE( NOSFTE(L,K) )
                  GL(1) = XYZSOM(1,N)
                  GL(2) = XYZSOM(2,N)
                  GL(3) = XYZSOM(3,N)
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
