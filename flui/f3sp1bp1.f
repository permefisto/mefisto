      SUBROUTINE F3SP1BP1( X,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                     BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 4 SOMMETS DU TETRAEDRE
C          INTERPOLATION DE FGAMMA AUX 3 SOMMETS DU TRIANGLE
C
C ENTREES:
C --------
C X      : 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
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
C SORTIES:
C --------
C BE     : BE(19) LE SECOND MEMBRE ELEMENTAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Novembre 2008
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      INTEGER            MCN
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(4,3)
      DOUBLE PRECISION  BE(19)
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOOBVC,    NUMIVO, NUMAVO
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
C
      INTEGER           I, J, K, L, M, MN, N, NU, NOOB
      DOUBLE PRECISION  FORCE(3,5)
      DOUBLE PRECISION  ABS, DETM33, DELTA, XD, YD, ZD,
     %                  A, B, C, D
      DOUBLE PRECISION  X21, Y21, Z21, X31, Y31, Z31, X41, Y41, Z41
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
C
      DOUBLE PRECISION  TP1BP1B(5,5)
C     TP1BLa(i,j) = Integrale sur e chapeau Pi Pj dx dy dz
C
C     NO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER  NOSFTE(3,4)
      DATA     NOSFTE / 1,3,2,   1,4,3,   1,2,4,  2,3,4/
C
C     PHASE PRELIMINAIRE
C     ------------------
C     INITIALISATION DES TABLEAUX UTILISES
      DO J=1,19
        BE(J) = 0.D0
      ENDDO
C
C     1~) CONTRIBUTION DES EFFORTS VOLUMIQUES
C     ----------------------------------------
      IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
C
         X21 = X(2,1) - X(1,1)
         X31 = X(3,1) - X(1,1)
         X41 = X(4,1) - X(1,1)
C
         Y21 = X(2,2) - X(1,2)
         Y31 = X(3,2) - X(1,2)
         Y41 = X(4,2) - X(1,2)
C
         Z21 = X(2,3) - X(1,3)
         Z31 = X(3,3) - X(1,3)
         Z41 = X(4,3) - X(1,3)
C
C        CALCUL DU DETERMINANT DE LA JACOBIENNE
         DELTA = ABS( DETM33( X21, X31, X41,
     %                        Y21, Y31, Y41,
     %                        Z21, Z31, Z41  ) )
C
C        VALEUR DES EFFORTS VOLUMIQUES AUX 4 SOMMETS DU TETRAEDRE
         DO J=1,4
            XD = X(J,1)
            YD = X(J,2)
            ZD = X(J,3)
            CALL REFORC( 4,NOOBVC, 3, XD,YD,ZD, 0D0,0D0,0D0,
     %                   LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
         ENDDO
C
C        LE BARYCENTRE
         XD = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
         YD = ( X(1,2) + X(2,2) + X(3,2) + X(4,2) ) * 0.25D0
         ZD = ( X(1,3) + X(2,3) + X(3,3) + X(4,3) ) * 0.25D0
         CALL REFORC( 4,NOOBVC, 3, XD,YD,ZD, 0D0,0D0,0D0,
     %                LTDEVO(LPFORC,NOOBVC), FORCE(1,5) )
C
C        LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
C        TP1BLa(i,j) = Integrale sur e chapeau Pi Pj dx dy dz
         A = 29836D0 / 2494800D0 * DELTA
         B =  9046D0 / 2494800D0 * DELTA
         C =   956D0 /  155925D0 * DELTA
         D =  4096D0 /  155925D0 * DELTA
C
CCC      A=  0.11959275292608627D-01
CCC      B=  0.36259419592752926D-02
CCC      C=  0.61311527978194641D-02
CCC      D=  0.26269039602372937D-01
C
         TP1BP1B(1,1) = A
         TP1BP1B(2,1) = B
         TP1BP1B(3,1) = B
         TP1BP1B(4,1) = B
         TP1BP1B(5,1) = C
C
         TP1BP1B(1,2) = B
         TP1BP1B(2,2) = A
         TP1BP1B(3,2) = B
         TP1BP1B(4,2) = B
         TP1BP1B(5,2) = C
C
         TP1BP1B(1,3) = B
         TP1BP1B(2,3) = B
         TP1BP1B(3,3) = A
         TP1BP1B(4,3) = B
         TP1BP1B(5,3) = C
C
         TP1BP1B(1,4) = B
         TP1BP1B(2,4) = B
         TP1BP1B(3,4) = B
         TP1BP1B(4,4) = A
         TP1BP1B(5,4) = C
C
         TP1BP1B(1,5) = C
         TP1BP1B(2,5) = C
         TP1BP1B(3,5) = C
         TP1BP1B(4,5) = C
         TP1BP1B(5,5) = D
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS VOLUMIQUES
         DO I=1,5
            DO J=1,5
               D = TP1BP1B(I,J)
               BE(I   ) = BE(I   ) + D * FORCE(1,J)
               BE(I+ 5) = BE(I+ 5) + D * FORCE(2,J)
               BE(I+10) = BE(I+10) + D * FORCE(3,J)
            ENDDO
         ENDDO
      ENDIF
C
C     2°) CONTRIBUTION DES EFFORTS SURFACIQUES
C     ----------------------------------------
      DO 100 K=1,4
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
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C              ...........................................
C              LES 3 PREMIERS NOEUDS SONT LES SOMMETS DE LA FACE K
C              SUPPOSEE DROITE C'EST A DIRE INTERPOLEE P1
C              NO ELEMENTAIRE DU SOMMET 1 DE LA FACE K
               N = NOSFTE(1,K)
               XD = X(N,1)
               YD = X(N,2)
               ZD = X(N,3)
C
C              NO ELEMENTAIRE DU SOMMET 2 DE LA FACE K
               N = NOSFTE(2,K)
               DGL(1,1) = X(N,1) - XD
               DGL(1,2) = X(N,2) - YD
               DGL(1,3) = X(N,3) - ZD
C
C              NO ELEMENTAIRE DU SOMMET 3 DE LA FACE K
               N = NOSFTE(3,K)
               DGL(2,1) = X(N,1) - XD
               DGL(2,2) = X(N,2) - YD
               DGL(2,3) = X(N,3) - ZD
C
C              CALCUL DU VECTEUR NORMAL A LA FACE K SUPPOSEE PLANE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO L=1,3
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE K
                  N = NOSFTE(L,K)
                  GL(1) = X(N,1)
                  GL(2) = X(N,2)
                  GL(3) = X(N,3)
                  CALL REFORC( 3, NOOB, 3,
     %                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     %                         LTDESU(LPFORC,NOOB), FORCE(1,L) )
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 1 DE LA FACE K
               NU = NOSFTE(1,K)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
C              PRISE EN COMPTE DES COEFF Integrale( tLa La )
               DGLN = DGLN / 12D0
               DO M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               ( FORCE(M,1) + (FORCE(M,2)+FORCE(M,3))/2D0 )
                  N = N + 5
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 2 DE LA FACE K
               NU = NOSFTE(2,K)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               ( FORCE(M,2) + (FORCE(M,1)+FORCE(M,3))/2D0 )
                  N = N + 5
               ENDDO
C
C              NUMERO ELEMENTAIRE DU NOEUD 3 DE LA FACE K
               NU = NOSFTE(3,K)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               ( FORCE(M,3) + (FORCE(M,1)+FORCE(M,2))/2D0 )
                  N = N + 5
               ENDDO
            ENDIF
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END
