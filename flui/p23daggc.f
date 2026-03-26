      SUBROUTINE P23DAGGC( COEF,   NBNOVI, XYZSOM,
     %                     NBNOEF, NBEF,   NONOEF,
     %                     NBCOPG, LPLIGN, LPCOLO, PG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE MORSE GLOBALE LES
C ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) POUR UNE
C          INTERPOLATION P2 EN 3D DE FONCTIONS DE BASE LAMBDA(1:10)
C
C ENTREES:
C --------
C COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C         (4 POUR BREZZI-FORTIN & 10 POUR TAYLOR-HOOD)
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
C
C NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
C LPLIGN : POINTEUR SUR LES COEFFICIENTS DIAGONAUX      DE PG
C LPCOLO : NUMERO DES COLONNES DES COEFFICIENTS STOCKES DE PG
C
C MODIFIES:
C ---------
C PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Juin 2012
C23456---------------------------------------------------------------012
      REAL              XYZSOM(3,NBNOVI)
      INTEGER           NONOEF(NBEF,NBNOEF),
     %                  LPLIGN(0:NBNOVI),
     %                  LPCOLO(1:NBCOPG)
      DOUBLE PRECISION  COEF, PG(1:NBCOPG)
C
      DOUBLE PRECISION  AE(55), DF(3,3), DFM1(3,3), DELTA, TDFDF(3,3), S
      INTEGER           NOSOTE(10)
      EQUIVALENCE      (NS1,NOSOTE(1)),(NS2,NOSOTE(2)),(NS3,NOSOTE(3)),
     %                 (NS4,NOSOTE(4))
      INTRINSIC         ABS
C
      include"./incl/dp2dp23d.inc"
C     DOUBLE PRECISION  DP2DP23D(3,10,3,10)
C
C     MISE A ZERO DE LA MATRICE MORSE GLOBALE
      DO I = 1, NBCOPG
         PG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 10 NOEUDS DU TETRAEDRE NEF
         DO I=1,10
            NOSOTE(I) = NONOEF(NEF,I)
         ENDDO
C
C        CONSTRUCTION DE LA MATRICE DF
         S = XYZSOM(1,NS1)
         DF(1,1) = XYZSOM(1,NS2) - S
         DF(2,1) = XYZSOM(1,NS3) - S
         DF(3,1) = XYZSOM(1,NS4) - S
C
         S = XYZSOM(2,NS1)
         DF(1,2) = XYZSOM(2,NS2) - S
         DF(2,2) = XYZSOM(2,NS3) - S
         DF(3,2) = XYZSOM(2,NS4) - S
C
         S = XYZSOM(3,NS1)
         DF(1,3) = XYZSOM(3,NS2) - S
         DF(2,3) = XYZSOM(3,NS3) - S
         DF(3,3) = XYZSOM(3,NS4) - S
C
C        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3))
     %             + DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3))
     %             + DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
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
C        t[DF]-1 [DF]-1  EST UNE MATRICE SYMETRIQUE
C        CALCUL DES COEFFICIENTS DE LA MATRICE DE VISCOSITE
C        tDF-1 * DF-1 * DELTAe * CoefV
         S = COEF / DELTA
         TDFDF(1,1) = ( DFM1(1,1) * DFM1(1,1)
     %                + DFM1(2,1) * DFM1(2,1)
     %                + DFM1(3,1) * DFM1(3,1) ) * S
C
         TDFDF(2,1) = ( DFM1(1,2) * DFM1(1,1)
     %                + DFM1(2,2) * DFM1(2,1)
     %                + DFM1(3,2) * DFM1(3,1) ) * S
         TDFDF(1,2) = TDFDF(2,1)
C
         TDFDF(3,1) = ( DFM1(1,3) * DFM1(1,1)
     %                + DFM1(2,3) * DFM1(2,1)
     %                + DFM1(3,3) * DFM1(3,1) ) * S
         TDFDF(1,3) = TDFDF(3,1)
C
         TDFDF(2,2) = ( DFM1(1,2) * DFM1(1,2)
     %                + DFM1(2,2) * DFM1(2,2)
     %                + DFM1(3,2) * DFM1(3,2) ) * S
C
         TDFDF(3,2) = ( DFM1(1,3) * DFM1(1,2)
     %                + DFM1(2,3) * DFM1(2,2)
     %                + DFM1(3,3) * DFM1(3,2) ) * S
         TDFDF(2,3) = TDFDF(3,2)
C
         TDFDF(3,3) = ( DFM1(1,3) * DFM1(1,3)
     %                + DFM1(2,3) * DFM1(2,3)
     %                + DFM1(3,3) * DFM1(3,3) ) * S
C
C        AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C        LIGNE PAR LIGNE ET DE GAUCHE A DROITE JUSQU'A LA DIAGONALE
         KE = 0
         DO I = 1, 10
            DO J = 1, I
C
C              Coef Integrale tdP2i/dxk dP2j/dxl * deltae
               S = 0D0
               DO K=1,3
                  DO L=1,3
                     S = S + DP2DP23D(K,I,L,J) * TDFDF(K,L)
                  ENDDO
               ENDDO
C
C              COEFFICIENT Ae(I,J) de la MATRICE Coef DP2 DP2
               KE = KE + 1
               AE(KE) = S
C
            ENDDO
         ENDDO
C
C        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE MORSE GLOBALE
         CALL ASMEGC( 10, NOSOTE, 1, AE, AE, 1, LPLIGN, LPCOLO, PG )
C
 100  CONTINUE
C
      RETURN
      END
