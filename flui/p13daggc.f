      SUBROUTINE P13DAGGC( COEF,   NBSOM,  XYZSOM,
     %                     NBNOEF, NBELEM, NUNOEF, NONOSO,
     %                     NBCOPG, LPLIGN, LPCOLO, PG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE MORSE GLOBALE LES
C ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) POUR UNE
C          INTERPOLATION P1 EN 3D DE FONCTIONS DE BASE LAMBDA(1:4)
C
C ENTREES:
C --------
C COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C         (4 POUR BREZZI-FORTIN & 10 POUR TAYLOR-HOOD)
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
C NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM
C
C NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
C LPLIGN : POINTEUR SUR LES COEFFICIENTS DIAGONAUX      DE PG
C LPCOLO : NUMERO DES COLONNES DES COEFFICIENTS STOCKES DE PG
C
C MODIFIES:
C ---------
C PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NUNOEF(NBELEM,NBNOEF),
     %                  NONOSO(1:*),
     %                  LPLIGN(0:NBSOM),
     %                  LPCOLO(1:NBCOPG)
      DOUBLE PRECISION  COEF, PG(1:NBCOPG)
C
      DOUBLE PRECISION  AE(10), DF(3,3), DFM1(3,3), DELTA, DFM1DLa(3,4),
     %                  S, X, Y, Z, C, EPSii, EPSij
      INTEGER           NOSOTE(4)
      EQUIVALENCE      (NS1,NOSOTE(1)),(NS2,NOSOTE(2)),(NS3,NOSOTE(3)),
     %                 (NS4,NOSOTE(4))
      INTRINSIC         ABS
C
C     MISE A ZERO DE LA MATRICE MORSE GLOBALE
      DO I = 1, NBCOPG
         PG( I ) = 0D0
      ENDDO
C
      DO 100 NUELEM = 1, NBELEM
C
C        NUMERO DES 4 SOMMETS DU TETRAEDRE NUELEM
         DO I=1,4
            NOSOTE(I) = NUNOEF(NUELEM,I)
         ENDDO
C
C        CONSTRUCTION DE LA MATRICE DF
         X = XYZSOM(1,NS1)
         DF(1,1) = XYZSOM(1,NS2) - X
         DF(2,1) = XYZSOM(1,NS3) - X
         DF(3,1) = XYZSOM(1,NS4) - X
C
         Y = XYZSOM(2,NS1)
         DF(1,2) = XYZSOM(2,NS2) - Y
         DF(2,2) = XYZSOM(2,NS3) - Y
         DF(3,2) = XYZSOM(2,NS4) - Y
C
         Z = XYZSOM(3,NS1)
         DF(1,3) = XYZSOM(3,NS2) - Z
         DF(2,3) = XYZSOM(3,NS3) - Z
         DF(3,3) = XYZSOM(3,NS4) - Z
C
C        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3))
     %             + DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3))
     %             + DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )
C
         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'p13daggc: ATTENTION EF',NUELEM,
     %                ' de VOLUME*6=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*,'p13daggc: ATTENTION FE',NUELEM,
     %                ' of VOLUME*6=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF
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
C        [DFM1] [DLa]
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
C        COEF / ( 6D0 * JACOBIEN de DF )
         C = COEF / ( 6D0 * DELTA )
C
C        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 120D0
         EPSii = EPSij * 2
C
C        AE = COEF  Integrale  t([DFM1] [DLa]) ([DFM1] [DLa]) dx
         M = 0
         DO I = 1, 4
            DO J = 1, I
               S = 0D0
               DO K = 1, 3
                  S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
               ENDDO
               M = M + 1
               IF( I .EQ. J ) THEN
                  AE(M) = S * C + EPSii
               ELSE
                  AE(M) = S * C + EPSij
               ENDIF
            ENDDO
         ENDDO

         IF( NBNOEF .GT. 4 ) THEN
C           EF TAYLOR-YOUNG
C           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
C           NUNOEF SONT LES NUMEROS DES NOEUDS P2
C           NONOSO SONT LES NUMEROS DE SOMMET DES NOEUDS P2
C           MAIS LPLIGN, LPCOLO SONT DEFINIS SUR UN MAILLAGE P1
            DO K=1,4
               NOSOTE(K) = NONOSO( NOSOTE(K) )
            ENDDO
         ENDIF
C
C        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE MORSE GLOBALE
         CALL ASMEGC( 4, NOSOTE, 1, AE, AE, 1, LPLIGN, LPCOLO, PG )
C
 100  CONTINUE
C
      RETURN
      END
