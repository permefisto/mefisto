      SUBROUTINE BGUGUTH2( NBSOM, XYZNOE, NBNOEF, NBEF, NONOEF, NONOSO,
     %                     P2DP2, Rho, DELTAT,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     NBNOVI, VITXYZ, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
C ----     - dt CoPres LAPLACIEN P = dt Div( -Fomega + Rho u. Grad u)
C          POUR LE TRIANGLE TAYLOR-HOOD
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C P2DP2  : P2DP2(i,k,j) = Integrale P2i DP2j/Dxk dx
C Rho    : DENSITE DE MASSE
C DELTAT : PAS DE TEMPS
C NOOBSF : NUMERO DE SURFACE DU FLUIDE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C VITXYZ : 2 COMPOSANTES DE LA VITESSE EN LES NBNOVI NOEUDS VITESSE
C          DU MAILLAGE
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Juin 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOVI, NBNOEF, NBEF,
     %                  NONOEF(NBEF,NBNOEF), NONOSO(NBNOVI),
     %                  NOOBSF, NUMISU, NUMASU
      DOUBLE PRECISION  DELTAT, Rho, VITXYZ(NBNOVI,2), BG(NBSOM)
C
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      DOUBLE PRECISION  DELTAe, DFM1(2,2), DFM1DLa(2,3),
     %                  S, FL, X, Y, Z, X21, Y21, X31, Y31
      INTEGER           NOSOTR(3), NS1, NS2, NS3
      EQUIVALENCE      (NOSOTR(1),NS1),(NOSOTR(2),NS2),(NOSOTR(3),NS3)
      DOUBLE PRECISION  FORCE(2,4:6)
      INTEGER           I, J, JJ, K, L, N, NSJ, NSJJ, NEF
      DOUBLE PRECISION  P2DP2(6,2,6)
C
C     INTEGRALE P2i dx SUR LE TRIANGLE UNITE  INTP2(1:3)=0D0
      DOUBLE PRECISION  INTP2(4:6)
      DATA              INTP2/    0.1666666666666666667D0,
     %   0.1666666666666666667D0, 0.1666666666666666667D0 /
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
C
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         X21 = XYZNOE(1,NS2) - XYZNOE(1,NS1)
         X31 = XYZNOE(1,NS3) - XYZNOE(1,NS1)
C
         Y21 = XYZNOE(2,NS2) - XYZNOE(2,NS1)
         Y31 = XYZNOE(2,NS3) - XYZNOE(2,NS1)
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
C        DFM1 DLambda(2,3)
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
C        ----------------------------------------------------------
         IF( LTDESU(LPFORC,NOOBSF) .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS SURFACIQUES AUX 3 MILIEUX
C           DES ARETES DU TRIANGLE (POIDS = 0 AUX SOMMETS)
            DO J=4,6
               NSJ = NONOEF(NEF,J)
               X = XYZNOE(1,NSJ)
               Y = XYZNOE(2,NSJ)
               Z = XYZNOE(3,NSJ)
               CALL REFORC( 3,NOOBSF, 2, X,Y,Z,  0D0,0D0,0D0,
     %                      LTDESU(LPFORC,NOOBSF), FORCE(1,J) )
            ENDDO
C
         ENDIF
C
         DO I=1,3
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            S = 0D0
C
            DO L=1,2
C
C              COEFFICIENTS DES EFFORTS SURFACIQUES L
               FL = 0D0
               IF( LTDESU(LPFORC,NOOBSF) .GT. 0 ) THEN
                  DO J=4,6
                     FL = FL + INTP2(J) * FORCE(L,J)
                  ENDDO
               ENDIF
C
C              - Integrale  Grad P1  Rho ( u. Grad u) dx
               DO K=1,2
                  DO J=1,6
C                    NO GLOBAL DU NOEUD J DU TRIANGLE TH
                     NSJ = NONOEF(NEF,J)
                     DO N=1,2
                        DO JJ=1,6
                           NSJJ = NONOEF(NEF,JJ)
C                          P2DP2(i,k,j) = integrale P2j dP2jj/dxn dx dy dz
                           FL = FL
     %                        - Rho * P2DP2(J,N,JJ) * VITXYZ(NSJ,K)
     %                              * DFM1(K,N) * VITXYZ(NSJJ,L)/ DELTAe
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               S = S + DFM1DLa(L,I) * FL
            ENDDO
C
C           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTR(I) ) )
            NSJ = NONOSO( NOSOTR(I) )
            BG(NSJ) = BG(NSJ) + S * DELTAT
C
         ENDDO
C
 100  CONTINUE
C
      RETURN
      END
